#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <zstd.h>

#include "delta/Gdelta/gdelta.h"
#include "delta/Gdelta/sketch_contract.h"
#include "utils/xxhash.h"

extern "C" {
#include "delta/xdelta3/xdelta3.h"
}

namespace {

using Clock = std::chrono::steady_clock;

constexpr std::size_t kReadSize = 4U * 1024U * 1024U;
constexpr std::size_t kMinChunk = 4U * 1024U;
constexpr std::size_t kAvgChunk = 8U * 1024U;
constexpr std::size_t kMaxChunk = 32U * 1024U;
constexpr std::size_t kQueueCapacity = 256;
constexpr int kZstdLevel = 1;

constexpr std::uint64_t kFastCdcSmallMask = 0x0000d9f003530000ULL;
constexpr std::uint64_t kFastCdcLargeMask = 0x0000d90003530000ULL;

double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

class UsageError : public std::runtime_error {
public:
    explicit UsageError(const std::string& message) : std::runtime_error(message) {}
};

class VerificationError : public std::runtime_error {
public:
    explicit VerificationError(const std::string& message) : std::runtime_error(message) {}
};

enum class Scheme { SsG, SsGS, SsX, OdG, OdX };
enum class Pipeline { Sequential, Parallel };

struct Options {
    std::string input;
    Scheme scheme = Scheme::SsG;
    Pipeline pipeline = Pipeline::Parallel;
    std::size_t features = 0;
    std::size_t super_features = 0;
    bool scheme_set = false;
    bool features_set = false;
    bool super_features_set = false;
    bool json = false;
    bool verify = true;
};

const char* scheme_name(Scheme scheme) {
    switch (scheme) {
    case Scheme::SsG: return "ss-g";
    case Scheme::SsGS: return "ss-g-s";
    case Scheme::SsX: return "ss-x";
    case Scheme::OdG: return "od-g";
    case Scheme::OdX: return "od-x";
    }
    return "unknown";
}

const char* pipeline_name(Pipeline pipeline) {
    return pipeline == Pipeline::Sequential ? "sequential" : "parallel";
}

bool is_odess(Scheme scheme) {
    return scheme == Scheme::OdG || scheme == Scheme::OdX;
}

bool is_xdelta(Scheme scheme) {
    return scheme == Scheme::SsX || scheme == Scheme::OdX;
}

bool is_sketch_gdelta(Scheme scheme) {
    return scheme == Scheme::SsGS;
}

void print_usage(std::ostream& out, const char* program) {
    out << "usage: " << program
        << " --input PATH --scheme ss-g|ss-g-s|ss-x|od-g|od-x"
           " [--pipeline sequential|parallel]"
           " [--features N --super-features N] [--json] [--no-verify]\n";
}

std::size_t parse_positive(const std::string& text, const char* option) {
    if (text.empty() || !std::all_of(text.begin(), text.end(), [](char c) {
            return c >= '0' && c <= '9';
        })) {
        throw UsageError(std::string(option) + " requires a positive integer");
    }
    std::size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used, 10);
    } catch (const std::exception&) {
        throw UsageError(std::string(option) + " requires a positive integer");
    }
    if (used != text.size() || value == 0 ||
        value > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
        throw UsageError(std::string(option) + " requires a positive integer");
    }
    return static_cast<std::size_t>(value);
}

Options parse_options(int argc, char** argv) {
    Options options;
    bool input_set = false;
    bool pipeline_set = false;
    bool json_set = false;
    bool no_verify_set = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage(std::cout, argv[0]);
            std::exit(0);
        }
        if (arg == "--json") {
            if (json_set) throw UsageError("duplicate --json");
            json_set = true;
            options.json = true;
            continue;
        }
        if (arg == "--no-verify") {
            if (no_verify_set) throw UsageError("duplicate --no-verify");
            no_verify_set = true;
            options.verify = false;
            continue;
        }
        if (i + 1 >= argc) throw UsageError(arg + " requires a value");
        const std::string value(argv[++i]);
        if (arg == "--input") {
            if (input_set) throw UsageError("duplicate --input");
            input_set = true;
            options.input = value;
        } else if (arg == "--scheme") {
            if (options.scheme_set) throw UsageError("duplicate --scheme");
            options.scheme_set = true;
            if (value == "ss-g") options.scheme = Scheme::SsG;
            else if (value == "ss-g-s") options.scheme = Scheme::SsGS;
            else if (value == "ss-x") options.scheme = Scheme::SsX;
            else if (value == "od-g") options.scheme = Scheme::OdG;
            else if (value == "od-x") options.scheme = Scheme::OdX;
            else throw UsageError("invalid --scheme: " + value);
        } else if (arg == "--pipeline") {
            if (pipeline_set) throw UsageError("duplicate --pipeline");
            pipeline_set = true;
            if (value == "sequential") options.pipeline = Pipeline::Sequential;
            else if (value == "parallel") options.pipeline = Pipeline::Parallel;
            else throw UsageError("invalid --pipeline: " + value);
        } else if (arg == "--features") {
            if (options.features_set) throw UsageError("duplicate --features");
            options.features_set = true;
            options.features = parse_positive(value, "--features");
        } else if (arg == "--super-features") {
            if (options.super_features_set) throw UsageError("duplicate --super-features");
            options.super_features_set = true;
            options.super_features = parse_positive(value, "--super-features");
        } else {
            throw UsageError("unknown option: " + arg);
        }
    }

    if (!input_set || options.input.empty()) throw UsageError("--input is required");
    if (!options.scheme_set) throw UsageError("--scheme is required");
    if (is_odess(options.scheme)) {
        if (!options.features_set || !options.super_features_set) {
            throw UsageError("ODESS schemes require --features and --super-features");
        }
        if (options.features < options.super_features ||
            options.features % options.super_features != 0) {
            throw UsageError("--features must be >= and divisible by --super-features");
        }
    } else if (options.features_set || options.super_features_set) {
        throw UsageError("--features and --super-features are only valid for ODESS schemes");
    }
    return options;
}

class FileDescriptor {
public:
    explicit FileDescriptor(const std::string& path) : fd_(::open(path.c_str(), O_RDONLY | O_CLOEXEC)) {
        if (fd_ < 0) throw std::runtime_error("cannot open input '" + path + "': " + std::strerror(errno));
        struct stat st {};
        if (::fstat(fd_, &st) != 0) {
            const std::string message = "cannot stat input '" + path + "': " + std::strerror(errno);
            ::close(fd_);
            fd_ = -1;
            throw std::runtime_error(message);
        }
        if (!S_ISREG(st.st_mode)) {
            ::close(fd_);
            fd_ = -1;
            throw std::runtime_error("input must be a regular seekable file");
        }
        size_ = static_cast<std::uint64_t>(st.st_size);
    }

    ~FileDescriptor() {
        if (fd_ >= 0) ::close(fd_);
    }

    FileDescriptor(const FileDescriptor&) = delete;
    FileDescriptor& operator=(const FileDescriptor&) = delete;

    int get() const { return fd_; }
    std::uint64_t size() const { return size_; }

private:
    int fd_ = -1;
    std::uint64_t size_ = 0;
};

void pread_exact(int fd, std::uint64_t offset, std::size_t size,
                 std::vector<std::uint8_t>& data) {
    data.resize(size);
    std::size_t done = 0;
    while (done < size) {
        const ssize_t n = ::pread(fd, data.data() + done, size - done,
                                  static_cast<off_t>(offset + done));
        if (n > 0) {
            done += static_cast<std::size_t>(n);
            continue;
        }
        if (n < 0 && errno == EINTR) continue;
        if (n == 0) throw std::runtime_error("input changed or truncated during pread");
        throw std::runtime_error(std::string("pread failed: ") + std::strerror(errno));
    }
}

std::size_t fastcdc_cut(const std::uint8_t* data, std::size_t available) {
    const std::size_t limit = std::min(available, kMaxChunk);
    if (limit <= kMinChunk) return limit;

    std::uint64_t fingerprint = 0;
    const std::size_t normal = std::min(limit, kAvgChunk);
    std::size_t i = kMinChunk;
    for (; i < normal; ++i) {
        fingerprint = (fingerprint >> 8) + GEARmx[data[i]];
        if ((fingerprint & kFastCdcSmallMask) == 0) return i;
    }
    for (; i < limit; ++i) {
        fingerprint = (fingerprint >> 8) + GEARmx[data[i]];
        if ((fingerprint & kFastCdcLargeMask) == 0) return i;
    }
    return limit;
}

struct Chunk {
    std::uint64_t offset = 0;
    std::vector<std::uint8_t> data;
};

class ChunkReader {
public:
    ChunkReader(int fd, std::uint64_t expected_size) : fd_(fd), expected_size_(expected_size) {}

    template <typename Emit>
    void run(Emit emit) {
        std::vector<std::uint8_t> buffer(kReadSize + kMaxChunk);
        std::size_t carried = 0;
        std::uint64_t chunk_offset = 0;
        bool eof = false;

        while (!eof) {
            std::size_t added = 0;
            while (added < kReadSize) {
                const auto read_start = Clock::now();
                const ssize_t n = ::read(fd_, buffer.data() + carried + added, kReadSize - added);
                work_seconds_ += seconds_since(read_start);
                if (n > 0) {
                    added += static_cast<std::size_t>(n);
                    source_read_bytes_ += static_cast<std::size_t>(n);
                    continue;
                }
                if (n < 0 && errno == EINTR) continue;
                if (n < 0) throw std::runtime_error(std::string("read failed: ") + std::strerror(errno));
                eof = true;
                break;
            }

            const std::size_t valid = carried + added;
            std::size_t consumed = 0;
            while (valid - consumed >= kMaxChunk || (eof && consumed < valid)) {
                const auto start = Clock::now();
                const std::size_t length = fastcdc_cut(buffer.data() + consumed, valid - consumed);
                if (length == 0) throw std::runtime_error("FastCDC produced an empty chunk");

                Chunk chunk;
                chunk.offset = chunk_offset;
                chunk.data.resize(length);
                std::memcpy(chunk.data.data(), buffer.data() + consumed, length);
                work_seconds_ += seconds_since(start);
                consumed += length;
                chunk_offset += length;
                if (!emit(std::move(chunk))) return;
            }

            carried = valid - consumed;
            if (carried != 0) {
                const auto move_start = Clock::now();
                std::memmove(buffer.data(), buffer.data() + consumed, carried);
                work_seconds_ += seconds_since(move_start);
            }
        }

        if (source_read_bytes_ != expected_size_ || chunk_offset != expected_size_) {
            throw std::runtime_error("input size changed while processing");
        }
    }

    std::uint64_t source_read_bytes() const { return source_read_bytes_; }
    double work_seconds() const { return work_seconds_; }

private:
    int fd_;
    std::uint64_t expected_size_;
    std::uint64_t source_read_bytes_ = 0;
    double work_seconds_ = 0;
};

struct ChunkRef {
    std::uint64_t offset = 0;
    std::uint32_t size = 0;
    std::uint64_t sketch = 0;
};

std::uint64_t speed_sketch(const std::vector<std::uint8_t>& data) {
    std::uint64_t fingerprint = 0;
    std::uint64_t sketch = 0;
    for (std::uint8_t byte : data) {
        fingerprint = speedsketch_contract::gear_step(fingerprint, byte);
        const std::size_t bit = speedsketch_contract::fingerprint_class(fingerprint);
        if (speedsketch_contract::is_proxy(fingerprint)) sketch |= 1ULL << bit;
    }
    return sketch;
}

class SimilarityIndex {
public:
    SimilarityIndex(const Options& options)
        : odess_(is_odess(options.scheme)), feature_count_(options.features),
          super_feature_count_(options.super_features) {
        if (!odess_) return;
        std::mt19937 generator(922);
        std::uniform_int_distribution<std::uint32_t> distribution(
            std::numeric_limits<std::uint32_t>::min(),
            std::numeric_limits<std::uint32_t>::max());
        transpose_m_.resize(feature_count_);
        transpose_a_.resize(feature_count_);
        for (std::size_t i = 0; i < feature_count_; ++i) {
            transpose_m_[i] = ((distribution(generator) >> 1U) << 1U) + 1U;
            transpose_a_[i] = distribution(generator);
        }
        features_.resize(feature_count_);
        super_features_.resize(super_feature_count_);
        odess_tables_.resize(super_feature_count_);
        rolling_hashes_.reserve(kMaxChunk / 64U);
        group_.resize(feature_count_ / super_feature_count_);
    }

    bool find_and_insert(const std::vector<std::uint8_t>& data, const ChunkRef& current,
                         ChunkRef& base, std::uint64_t& sketch) {
        if (!odess_) {
            sketch = speed_sketch(data);
            const std::uint32_t feature = static_cast<std::uint32_t>(sketch);
            const auto found = speed_table_.find(feature);
            const bool matched = found != speed_table_.end();
            if (matched) base = found->second;
            ChunkRef stored = current;
            stored.sketch = sketch;
            if (found == speed_table_.end()) speed_table_.emplace(feature, stored);
            else found->second = stored;
            return matched;
        }

        calculate_odess(data);
        bool matched = false;
        for (std::size_t i = 0; i < super_feature_count_; ++i) {
            const auto found = odess_tables_[i].find(super_features_[i]);
            if (!matched && found != odess_tables_[i].end()) {
                base = found->second;
                matched = true;
            }
        }
        for (std::size_t i = 0; i < super_feature_count_; ++i) {
            odess_tables_[i][super_features_[i]] = current;
        }
        sketch = 0;
        return matched;
    }

private:
    void calculate_odess(const std::vector<std::uint8_t>& data) {
        std::fill(features_.begin(), features_.end(), 0);
        rolling_hashes_.clear();
        std::uint32_t hash = 0;
        for (std::uint8_t byte : data) {
            hash = static_cast<std::uint32_t>((hash >> 1U) + GEARmx[byte]);
            if ((hash & 0x2ac6U) == 0) rolling_hashes_.push_back(hash);
        }
        for (std::size_t i = 0; i < feature_count_; ++i) {
            for (std::uint32_t rolling : rolling_hashes_) {
                const std::uint32_t transformed = static_cast<std::uint32_t>(
                    static_cast<std::uint64_t>(transpose_m_[i]) * rolling + transpose_a_[i]);
                features_[i] = std::max(features_[i], transformed);
            }
        }
        for (std::size_t i = 0; i < super_feature_count_; ++i) {
            for (std::size_t j = 0; j < group_.size(); ++j) {
                group_[j] = features_[j * super_feature_count_ + i];
            }
            super_features_[i] = XXH64(group_.data(), group_.size() * sizeof(group_[0]), 0);
        }
    }

    bool odess_;
    std::size_t feature_count_;
    std::size_t super_feature_count_;
    std::unordered_map<std::uint32_t, ChunkRef> speed_table_;
    std::vector<std::unordered_map<std::uint64_t, ChunkRef>> odess_tables_;
    std::vector<std::uint32_t> transpose_m_;
    std::vector<std::uint32_t> transpose_a_;
    std::vector<std::uint32_t> features_;
    std::vector<std::uint64_t> super_features_;
    std::vector<std::uint32_t> rolling_hashes_;
    std::vector<std::uint64_t> group_;
};

struct DeduplicatedChunk {
    Chunk chunk;
    bool duplicate = false;
    std::uint64_t base_read_bytes = 0;
    double dedup_seconds = 0;
};

class Deduplicator {
public:
    explicit Deduplicator(int fd) : fd_(fd) {}

    DeduplicatedChunk process(Chunk chunk) {
        DeduplicatedChunk result;
        result.chunk = std::move(chunk);

        const auto dedup_start = Clock::now();
        const XXH64_hash_t hash = dedup_hash(result.chunk.data);
        const auto found = dedup_.find(hash);
        if (found != dedup_.end()) {
            for (const ChunkRef& candidate : found->second) {
                if (candidate.size != result.chunk.data.size()) continue;
                pread_exact(fd_, candidate.offset, candidate.size, candidate_buffer_);
                result.base_read_bytes += candidate.size;
                if (candidate_buffer_ == result.chunk.data) {
                    result.duplicate = true;
                    break;
                }
            }
        }
        if (!result.duplicate) {
            const ChunkRef current{result.chunk.offset,
                                   static_cast<std::uint32_t>(result.chunk.data.size()), 0};
            if (found == dedup_.end()) {
                dedup_.emplace(hash, std::vector<ChunkRef>{current});
            } else {
                found->second.push_back(current);
            }
        }
        result.dedup_seconds = seconds_since(dedup_start);
        return result;
    }

private:
    static XXH64_hash_t dedup_hash(const std::vector<std::uint8_t>& data) {
#ifdef SPEEDSKETCH_FORCE_XXH64_COLLISION
        (void)data;
        return 0;
#else
        return XXH64(data.data(), data.size(), 0);
#endif
    }

    int fd_;
    std::vector<std::uint8_t> candidate_buffer_;
    std::unordered_map<XXH64_hash_t, std::vector<ChunkRef>> dedup_;
};

struct ClassifiedChunk {
    Chunk chunk;
    bool duplicate = false;
    bool has_base = false;
    ChunkRef base;
    std::uint64_t sketch = 0;
    std::uint64_t base_read_bytes = 0;
    double dedup_seconds = 0;
    double sketch_seconds = 0;
};

class Sketcher {
public:
    explicit Sketcher(const Options& options) : similarity_(options) {}

    ClassifiedChunk process(DeduplicatedChunk item) {
        ClassifiedChunk result;
        result.chunk = std::move(item.chunk);
        result.duplicate = item.duplicate;
        result.base_read_bytes = item.base_read_bytes;
        result.dedup_seconds = item.dedup_seconds;
        if (result.duplicate) return result;

        const auto sketch_start = Clock::now();
        const ChunkRef current{result.chunk.offset,
                               static_cast<std::uint32_t>(result.chunk.data.size()), 0};
        result.has_base = similarity_.find_and_insert(result.chunk.data, current,
                                                      result.base, result.sketch);
        result.sketch_seconds = seconds_since(sketch_start);
        return result;
    }

private:
    SimilarityIndex similarity_;
};

enum class Choice { Duplicate, Delta, Self, Raw };

struct EncodedChunk {
    bool duplicate = false;
    bool similar = false;
    Choice choice = Choice::Raw;
    std::uint64_t raw_size = 0;
    std::uint64_t after_delta_size = 0;
    std::uint64_t final_size = 0;
    std::uint64_t source_base_read_bytes = 0;
    std::uint64_t verification_checks = 0;
    std::uint64_t verification_failures = 0;
    GdeltaStats gdelta_stats;
    double dedup_seconds = 0;
    double sketch_seconds = 0;
    double encoding_seconds = 0;
    double verification_seconds = 0;
};

struct ZstdContextDeleter {
    void operator()(ZSTD_CCtx* context) const { ZSTD_freeCCtx(context); }
};

#ifdef SPEEDSKETCH_TESTING
std::uint64_t injected_encoding_failure_offset = std::numeric_limits<std::uint64_t>::max();
#endif

class Encoder {
public:
    Encoder(int fd, const Options& options)
        : fd_(fd), options_(options), zstd_(ZSTD_createCCtx()) {
        if (!zstd_) throw std::runtime_error("cannot create ZSTD context");
    }

    EncodedChunk process(ClassifiedChunk item) {
#ifdef SPEEDSKETCH_TESTING
        if (item.chunk.offset == injected_encoding_failure_offset) {
            throw std::runtime_error("injected encoding failure");
        }
#endif
        EncodedChunk result;
        result.duplicate = item.duplicate;
        result.similar = item.has_base;
        result.raw_size = item.chunk.data.size();
        result.source_base_read_bytes = item.base_read_bytes;
        result.dedup_seconds = item.dedup_seconds;
        result.sketch_seconds = item.sketch_seconds;
        if (item.duplicate) {
            result.choice = Choice::Duplicate;
            return result;
        }

        const auto encode_start = Clock::now();
        const std::size_t self_size = zstd_compress(item.chunk.data, self_zstd_buffer_);
        std::size_t delta_choice_size = std::numeric_limits<std::size_t>::max();
        if (item.has_base) {
            pread_exact(fd_, item.base.offset, item.base.size, base_buffer_);
            result.source_base_read_bytes += item.base.size;
            if (is_xdelta(options_.scheme)) {
                xdelta_encode(item.chunk.data, base_buffer_, delta_buffer_);
                delta_choice_size = delta_buffer_.size();
            } else {
                gdelta_encode(item.chunk.data, base_buffer_, item.sketch, item.base.sketch,
                              result.gdelta_stats, delta_buffer_);
                delta_choice_size = zstd_compress(delta_buffer_, delta_zstd_buffer_);
            }
            result.after_delta_size = std::min<std::uint64_t>(result.raw_size,
                                                              delta_buffer_.size());
        } else {
            result.after_delta_size = result.raw_size;
        }
        result.encoding_seconds = seconds_since(encode_start);

        if (item.has_base && delta_choice_size < self_size &&
            delta_choice_size < item.chunk.data.size()) {
            result.choice = Choice::Delta;
            result.final_size = delta_choice_size;
        } else if (self_size < item.chunk.data.size()) {
            result.choice = Choice::Self;
            result.final_size = self_size;
        } else {
            result.choice = Choice::Raw;
            result.final_size = item.chunk.data.size();
        }

        if (options_.verify &&
            (result.choice == Choice::Delta || result.choice == Choice::Self)) {
            const auto verify_start = Clock::now();
            ++result.verification_checks;
            bool valid = false;
            if (result.choice == Choice::Self) {
                valid = verify_zstd(self_zstd_buffer_, item.chunk.data);
            } else if (is_xdelta(options_.scheme)) {
                valid = verify_xdelta(delta_buffer_, base_buffer_, item.chunk.data);
            } else {
                valid = verify_gdelta(delta_zstd_buffer_, delta_buffer_.size(), base_buffer_,
                                      item.chunk.data);
            }
            result.verification_seconds = seconds_since(verify_start);
            if (!valid) {
                ++result.verification_failures;
                throw VerificationError("selected representation verification failed at input offset " +
                                        std::to_string(item.chunk.offset));
            }
        }
        return result;
    }

private:
    std::size_t zstd_compress(const std::vector<std::uint8_t>& input,
                              std::vector<std::uint8_t>& output) {
        output.resize(ZSTD_compressBound(input.size()));
        const std::size_t size = ZSTD_compressCCtx(zstd_.get(), output.data(),
                                                   output.size(), input.data(),
                                                   input.size(), kZstdLevel);
        if (ZSTD_isError(size)) {
            throw std::runtime_error(std::string("ZSTD compression failed: ") +
                                     ZSTD_getErrorName(size));
        }
        output.resize(size);
        return size;
    }

    void gdelta_encode(const std::vector<std::uint8_t>& current,
                       const std::vector<std::uint8_t>& base,
                       std::uint64_t current_sketch,
                       std::uint64_t base_sketch,
                       GdeltaStats& stats,
                       std::vector<std::uint8_t>& delta) {
        delta.resize(std::max<std::size_t>(32, current.size() * 2 + 32));
        std::size_t size = 0;
        int status = call_gdelta(current, base, current_sketch, base_sketch,
                                 delta.data(), delta.size(), &size, &stats);
        if (status == GDELTA_OUTPUT_TOO_SMALL) {
            delta.resize(size);
            status = call_gdelta(current, base, current_sketch, base_sketch,
                                 delta.data(), delta.size(), &size, &stats);
        }
        if (status != GDELTA_OK) {
            throw std::runtime_error("Gdelta encode failed with status " + std::to_string(status));
        }
        delta.resize(size);
    }

    int call_gdelta(const std::vector<std::uint8_t>& current,
                    const std::vector<std::uint8_t>& base,
                    std::uint64_t current_sketch, std::uint64_t base_sketch,
                    std::uint8_t* output, std::size_t capacity, std::size_t* size,
                    GdeltaStats* stats) {
        if (is_sketch_gdelta(options_.scheme)) {
            return gencode_whash_into(current.data(), current.size(), base.data(), base.size(),
                                      current_sketch, base_sketch, output, capacity, size,
                                      &gdelta_workspace_, stats);
        }
        return gencode_into(current.data(), current.size(), base.data(), base.size(),
                            output, capacity, size, &gdelta_workspace_, stats);
    }

    static void xdelta_encode(const std::vector<std::uint8_t>& current,
                              const std::vector<std::uint8_t>& base,
                              std::vector<std::uint8_t>& delta) {
        std::size_t capacity = std::max<std::size_t>(1024, current.size() * 2 + 1024);
        for (int attempt = 0; attempt < 4; ++attempt) {
            delta.resize(capacity);
            usize_t size = 0;
            const int status = xd3_encode_memory(
                current.data(), static_cast<usize_t>(current.size()),
                base.data(), static_cast<usize_t>(base.size()), delta.data(), &size,
                static_cast<usize_t>(delta.size()), XD3_ADLER32);
            if (status == 0) {
                delta.resize(size);
                return;
            }
            if (status != ENOSPC) {
                throw std::runtime_error("Xdelta encode failed with status " +
                                         std::to_string(status));
            }
            capacity *= 2;
        }
        throw std::runtime_error("Xdelta output exceeded bounded retry capacity");
    }

    bool verify_zstd(const std::vector<std::uint8_t>& compressed,
                     const std::vector<std::uint8_t>& expected) {
        restored_buffer_.resize(expected.size());
        const std::size_t size = ZSTD_decompress(restored_buffer_.data(), restored_buffer_.size(),
                                                 compressed.data(), compressed.size());
        return !ZSTD_isError(size) && size == expected.size() && restored_buffer_ == expected;
    }

    bool verify_xdelta(const std::vector<std::uint8_t>& delta,
                       const std::vector<std::uint8_t>& base,
                       const std::vector<std::uint8_t>& expected) {
        restored_buffer_.resize(expected.size());
        usize_t size = 0;
        const int status = xd3_decode_memory(
            delta.data(), static_cast<usize_t>(delta.size()),
            base.data(), static_cast<usize_t>(base.size()), restored_buffer_.data(), &size,
            static_cast<usize_t>(restored_buffer_.size()), 0);
        return status == 0 && size == expected.size() && restored_buffer_ == expected;
    }

    bool verify_gdelta(const std::vector<std::uint8_t>& compressed_delta,
                       std::size_t raw_delta_size,
                       const std::vector<std::uint8_t>& base,
                       const std::vector<std::uint8_t>& expected) {
        decoded_delta_buffer_.resize(raw_delta_size);
        const std::size_t decompressed_size = ZSTD_decompress(
            decoded_delta_buffer_.data(), decoded_delta_buffer_.size(),
            compressed_delta.data(), compressed_delta.size());
        if (ZSTD_isError(decompressed_size) || decompressed_size != raw_delta_size) return false;

        restored_buffer_.resize(expected.size());
        std::size_t size = 0;
        const int status = gdecode_into(decoded_delta_buffer_.data(), decoded_delta_buffer_.size(),
                                        base.data(), base.size(), restored_buffer_.data(),
                                        restored_buffer_.size(), &size);
        return status == GDELTA_OK && size == expected.size() && restored_buffer_ == expected;
    }

    int fd_;
    const Options& options_;
    std::unique_ptr<ZSTD_CCtx, ZstdContextDeleter> zstd_;
    std::vector<std::uint8_t> base_buffer_;
    std::vector<std::uint8_t> delta_buffer_;
    std::vector<std::uint8_t> self_zstd_buffer_;
    std::vector<std::uint8_t> delta_zstd_buffer_;
    std::vector<std::uint8_t> decoded_delta_buffer_;
    std::vector<std::uint8_t> restored_buffer_;
    GdeltaWorkspace gdelta_workspace_;
};

struct Metrics {
    std::uint64_t input_bytes = 0;
    std::uint64_t chunks = 0;
    std::uint64_t unique_chunks = 0;
    std::uint64_t duplicate_chunks = 0;
    std::uint64_t similar_chunks = 0;
    std::uint64_t delta_chunks = 0;
    std::uint64_t self_compressed_chunks = 0;
    std::uint64_t raw_chunks = 0;
    std::uint64_t after_dedup_bytes = 0;
    std::uint64_t after_delta_bytes = 0;
    std::uint64_t estimated_final_payload_bytes = 0;
    std::uint64_t source_read_bytes = 0;
    std::uint64_t base_read_bytes = 0;
    std::uint64_t verification_checks = 0;
    std::uint64_t verification_failures = 0;
    std::uint64_t judgments = 0;
    std::uint64_t excluded_lookups = 0;
    std::uint64_t false_positive_lookups = 0;
    double wall_seconds = 0;
    double chunking_seconds = 0;
    double dedup_seconds = 0;
    double sketch_seconds = 0;
    double encoding_seconds = 0;
    double verification_seconds = 0;

    void add_duplicate(const DeduplicatedChunk& item) {
        ++chunks;
        ++duplicate_chunks;
        base_read_bytes += item.base_read_bytes;
        dedup_seconds += item.dedup_seconds;
    }

    void add(const EncodedChunk& item) {
        ++chunks;
        dedup_seconds += item.dedup_seconds;
        sketch_seconds += item.sketch_seconds;
        encoding_seconds += item.encoding_seconds;
        verification_seconds += item.verification_seconds;
        base_read_bytes += item.source_base_read_bytes;
        verification_checks += item.verification_checks;
        verification_failures += item.verification_failures;
        judgments += item.gdelta_stats.judgments;
        excluded_lookups += item.gdelta_stats.excluded_lookups;
        false_positive_lookups += item.gdelta_stats.false_positive_lookups;

        if (item.duplicate) {
            ++duplicate_chunks;
            return;
        }
        ++unique_chunks;
        after_dedup_bytes += item.raw_size;
        after_delta_bytes += item.after_delta_size;
        estimated_final_payload_bytes += item.final_size;
        if (item.similar) ++similar_chunks;
        switch (item.choice) {
        case Choice::Delta: ++delta_chunks; break;
        case Choice::Self: ++self_compressed_chunks; break;
        case Choice::Raw: ++raw_chunks; break;
        case Choice::Duplicate: break;
        }
    }
};

template <typename T>
class BoundedQueue {
public:
    explicit BoundedQueue(std::size_t capacity) : capacity_(capacity) {}

    bool push(T item, const std::atomic<bool>& cancelled) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [&] { return queue_.size() < capacity_ || closed_ || cancelled.load(); });
        if (closed_ || cancelled.load()) return false;
        queue_.push_back(std::move(item));
        not_empty_.notify_one();
        return true;
    }

    bool pop(T& item, const std::atomic<bool>& cancelled) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_empty_.wait(lock, [&] { return !queue_.empty() || closed_ || cancelled.load(); });
        if (cancelled.load()) return false;
        if (queue_.empty()) return false;
        item = std::move(queue_.front());
        queue_.pop_front();
        not_full_.notify_one();
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> lock(mutex_);
        closed_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    std::size_t capacity_;
    std::deque<T> queue_;
    std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    bool closed_ = false;
};

Metrics run_sequential(int fd, const Options& options, std::uint64_t input_size) {
    ChunkReader reader(fd, input_size);
    Deduplicator deduplicator(fd);
    Sketcher sketcher(options);
    Encoder encoder(fd, options);
    Metrics metrics;
    const auto wall_start = Clock::now();
    reader.run([&](Chunk chunk) {
        DeduplicatedChunk item = deduplicator.process(std::move(chunk));
        if (item.duplicate) metrics.add_duplicate(item);
        else metrics.add(encoder.process(sketcher.process(std::move(item))));
        return true;
    });
    metrics.wall_seconds = seconds_since(wall_start);
    metrics.input_bytes = input_size;
    metrics.source_read_bytes = reader.source_read_bytes();
    metrics.chunking_seconds = reader.work_seconds();
    return metrics;
}

Metrics run_parallel(int fd, const Options& options, std::uint64_t input_size) {
    ChunkReader reader(fd, input_size);
    Deduplicator deduplicator(fd);
    Sketcher sketcher(options);
    Encoder encoder(fd, options);
    Metrics metrics;
    Metrics duplicate_metrics;
    BoundedQueue<Chunk> chunks(kQueueCapacity);
    BoundedQueue<DeduplicatedChunk> deduplicated(kQueueCapacity);
    BoundedQueue<ClassifiedChunk> classified(kQueueCapacity);
    std::atomic<bool> cancelled(false);
    std::mutex error_mutex;
    std::exception_ptr error;

    const auto fail = [&](std::exception_ptr candidate) {
        bool expected = false;
        if (cancelled.compare_exchange_strong(expected, true)) {
            std::lock_guard<std::mutex> lock(error_mutex);
            error = candidate;
        }
        chunks.close();
        deduplicated.close();
        classified.close();
    };

    const auto wall_start = Clock::now();
    std::vector<std::thread> threads;
    try {
        threads.emplace_back([&] {
            try {
                reader.run([&](Chunk chunk) { return chunks.push(std::move(chunk), cancelled); });
                chunks.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });
        threads.emplace_back([&] {
            try {
                Chunk chunk;
                while (chunks.pop(chunk, cancelled)) {
                    DeduplicatedChunk item = deduplicator.process(std::move(chunk));
                    if (item.duplicate) {
                        duplicate_metrics.add_duplicate(item);
                        continue;
                    }
                    if (!deduplicated.push(std::move(item), cancelled)) break;
                }
                deduplicated.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });
        threads.emplace_back([&] {
            try {
                DeduplicatedChunk item;
                while (deduplicated.pop(item, cancelled)) {
                    if (!classified.push(sketcher.process(std::move(item)), cancelled)) break;
                }
                classified.close();
            } catch (...) {
                fail(std::current_exception());
            }
        });
        threads.emplace_back([&] {
            try {
                ClassifiedChunk item;
                while (classified.pop(item, cancelled)) {
                    metrics.add(encoder.process(std::move(item)));
                }
            } catch (...) {
                fail(std::current_exception());
            }
        });
    } catch (...) {
        fail(std::current_exception());
    }

    for (std::thread& thread : threads) {
        if (thread.joinable()) thread.join();
    }
    if (error) std::rethrow_exception(error);

    metrics.chunks += duplicate_metrics.chunks;
    metrics.duplicate_chunks += duplicate_metrics.duplicate_chunks;
    metrics.base_read_bytes += duplicate_metrics.base_read_bytes;
    metrics.dedup_seconds += duplicate_metrics.dedup_seconds;

    metrics.wall_seconds = seconds_since(wall_start);
    metrics.input_bytes = input_size;
    metrics.source_read_bytes = reader.source_read_bytes();
    metrics.chunking_seconds = reader.work_seconds();
    return metrics;
}

struct OptionalNumber {
    bool present = false;
    double value = 0;
};

OptionalNumber ratio(double numerator, double denominator) {
    if (denominator == 0) return {};
    return {true, numerator / denominator};
}

std::string json_escape(const std::string& value) {
    std::string result;
    result.reserve(value.size() + 2);
    for (unsigned char c : value) {
        switch (c) {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\b': result += "\\b"; break;
        case '\f': result += "\\f"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if (c < 0x20) {
                static const char hex[] = "0123456789abcdef";
                result += "\\u00";
                result += hex[c >> 4U];
                result += hex[c & 15U];
            } else {
                result += static_cast<char>(c);
            }
        }
    }
    return result;
}

void print_optional_json(std::ostream& out, const OptionalNumber& number) {
    if (number.present) out << number.value;
    else out << "null";
}

void print_json(const Options& options, const Metrics& metrics) {
    const double mib = 1024.0 * 1024.0;
    const double input_mib = static_cast<double>(metrics.input_bytes) / mib;
    const OptionalNumber estimated_drr = ratio(
        static_cast<double>(metrics.input_bytes),
        static_cast<double>(metrics.estimated_final_payload_bytes));
    OptionalNumber dce;
    if (metrics.after_dedup_bytes != 0) {
        dce = {true, (1.0 - static_cast<double>(metrics.after_delta_bytes) /
                           static_cast<double>(metrics.after_dedup_bytes)) * 100.0};
    }
    OptionalNumber exclusion;
    OptionalNumber false_positive;
    if (is_sketch_gdelta(options.scheme) && metrics.judgments != 0) {
        const double judgments = static_cast<double>(metrics.judgments);
        exclusion = {true, 100.0 * static_cast<double>(metrics.excluded_lookups) / judgments};
        false_positive = {
            true, 100.0 * static_cast<double>(metrics.false_positive_lookups) / judgments};
    }
    const OptionalNumber wall_rate = ratio(input_mib, metrics.wall_seconds);
    const OptionalNumber chunk_rate = ratio(input_mib, metrics.chunking_seconds);
    const OptionalNumber dedup_rate = ratio(input_mib, metrics.dedup_seconds);
    const OptionalNumber sketch_rate = ratio(input_mib, metrics.sketch_seconds);
    const OptionalNumber encoding_rate = ratio(input_mib, metrics.encoding_seconds);

    std::cout << std::setprecision(17)
              << "{\"schema_version\":1"
              << ",\"input\":\"" << json_escape(options.input) << "\""
              << ",\"scheme\":\"" << scheme_name(options.scheme) << "\""
              << ",\"pipeline\":\"" << pipeline_name(options.pipeline) << "\""
              << ",\"verify\":" << (options.verify ? "true" : "false")
              << ",\"features\":";
    if (is_odess(options.scheme)) std::cout << options.features;
    else std::cout << "null";
    std::cout << ",\"super_features\":";
    if (is_odess(options.scheme)) std::cout << options.super_features;
    else std::cout << "null";
    std::cout << ",\"input_bytes\":" << metrics.input_bytes
              << ",\"chunks\":" << metrics.chunks
              << ",\"unique_chunks\":" << metrics.unique_chunks
              << ",\"duplicate_chunks\":" << metrics.duplicate_chunks
              << ",\"similar_chunks\":" << metrics.similar_chunks
              << ",\"delta_chunks\":" << metrics.delta_chunks
              << ",\"self_compressed_chunks\":" << metrics.self_compressed_chunks
              << ",\"raw_chunks\":" << metrics.raw_chunks
              << ",\"after_dedup_bytes\":" << metrics.after_dedup_bytes
              << ",\"after_delta_bytes\":" << metrics.after_delta_bytes
              << ",\"estimated_final_payload_bytes\":" << metrics.estimated_final_payload_bytes
              << ",\"source_read_bytes\":" << metrics.source_read_bytes
              << ",\"base_read_bytes\":" << metrics.base_read_bytes
              << ",\"estimated_drr\":";
    print_optional_json(std::cout, estimated_drr);
    std::cout << ",\"dce_pct\":";
    print_optional_json(std::cout, dce);
    std::cout << ",\"exclusion_rate_pct\":";
    print_optional_json(std::cout, exclusion);
    std::cout << ",\"false_positive_rate_pct\":";
    print_optional_json(std::cout, false_positive);
    std::cout << ",\"wall_seconds\":" << metrics.wall_seconds
              << ",\"chunking_seconds\":" << metrics.chunking_seconds
              << ",\"dedup_seconds\":" << metrics.dedup_seconds
              << ",\"sketch_seconds\":" << metrics.sketch_seconds
              << ",\"encoding_seconds\":" << metrics.encoding_seconds
              << ",\"verification_seconds\":" << metrics.verification_seconds
              << ",\"wall_mib_s\":";
    print_optional_json(std::cout, wall_rate);
    std::cout << ",\"chunking_mib_s\":";
    print_optional_json(std::cout, chunk_rate);
    std::cout << ",\"dedup_mib_s\":";
    print_optional_json(std::cout, dedup_rate);
    std::cout << ",\"sketch_mib_s\":";
    print_optional_json(std::cout, sketch_rate);
    std::cout << ",\"encoding_mib_s\":";
    print_optional_json(std::cout, encoding_rate);
    std::cout << ",\"verification_checks\":" << metrics.verification_checks
              << ",\"verification_failures\":" << metrics.verification_failures
              << "}\n";
}

void print_human(const Options& options, const Metrics& metrics) {
    std::cout << std::fixed << std::setprecision(4)
              << "scheme: " << scheme_name(options.scheme) << '\n'
              << "pipeline: " << pipeline_name(options.pipeline) << '\n'
              << "verify: " << (options.verify ? "true" : "false") << '\n'
              << "input bytes: " << metrics.input_bytes << '\n'
              << "chunks: " << metrics.chunks << " (unique " << metrics.unique_chunks
              << ", duplicate " << metrics.duplicate_chunks << ")\n"
              << "similar chunks: " << metrics.similar_chunks << '\n'
              << "selected: delta " << metrics.delta_chunks << ", self "
              << metrics.self_compressed_chunks << ", raw " << metrics.raw_chunks << '\n'
              << "after dedup bytes: " << metrics.after_dedup_bytes << '\n'
              << "after delta bytes: " << metrics.after_delta_bytes << '\n'
              << "estimated final payload bytes: " << metrics.estimated_final_payload_bytes << '\n'
              << "source/base read bytes: " << metrics.source_read_bytes << '/'
              << metrics.base_read_bytes << '\n'
              << "verification: " << metrics.verification_checks << " checks, "
              << metrics.verification_failures << " failures\n"
              << "wall seconds: " << metrics.wall_seconds << '\n';
}

} // namespace

#ifndef SPEEDSKETCH_TESTING
int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const FileDescriptor input(options.input);
        Metrics metrics = options.pipeline == Pipeline::Sequential
            ? run_sequential(input.get(), options, input.size())
            : run_parallel(input.get(), options, input.size());
        if (options.json) print_json(options, metrics);
        else print_human(options, metrics);
        return 0;
    } catch (const UsageError& error) {
        std::cerr << "error: " << error.what() << '\n';
        print_usage(std::cerr, argv[0]);
        return 2;
    } catch (const VerificationError& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
#endif
