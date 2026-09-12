#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#define SPEEDSKETCH_TESTING
#define SPEEDSKETCH_FORCE_XXH64_COLLISION
#include "../speedsketch.cpp"

namespace {

void keep_embedded_cli_symbols_used() {
    (void)&parse_options;
    (void)&run_sequential;
    (void)&run_parallel;
    (void)&print_json;
    (void)&print_human;
}

class TempFile {
public:
    TempFile() {
        char path[] = "/tmp/speedsketch-core-XXXXXX";
        fd = ::mkstemp(path);
        assert(fd >= 0);
        assert(::unlink(path) == 0);
    }

    ~TempFile() { ::close(fd); }

    void append(const std::vector<std::uint8_t>& data) {
        std::size_t written = 0;
        while (written < data.size()) {
            const ssize_t size = ::write(fd, data.data() + written, data.size() - written);
            assert(size > 0);
            written += static_cast<std::size_t>(size);
        }
    }

    int fd = -1;
};

Chunk make_chunk(std::uint64_t offset, const std::vector<std::uint8_t>& data) {
    Chunk chunk;
    chunk.offset = offset;
    chunk.data = data;
    return chunk;
}

std::vector<std::uint8_t> random_bytes(std::size_t size, std::uint32_t state) {
    std::vector<std::uint8_t> data(size);
    for (std::uint8_t& byte : data) {
        state ^= state << 13U;
        state ^= state >> 17U;
        state ^= state << 5U;
        byte = static_cast<std::uint8_t>(state);
    }
    return data;
}

ClassifiedChunk make_similar_chunk(std::uint64_t offset,
                                   const std::vector<std::uint8_t>& data,
                                   std::uint64_t base_offset,
                                   const std::vector<std::uint8_t>& base) {
    ClassifiedChunk item;
    item.chunk = make_chunk(offset, data);
    item.has_base = true;
    item.base = {base_offset, static_cast<std::uint32_t>(base.size()), speed_sketch(base)};
    item.sketch = speed_sketch(data);
    return item;
}

void test_fixed_sketch_and_contract() {
    std::vector<std::uint8_t> input(256);
    for (std::size_t i = 0; i < input.size(); ++i) {
        input[i] = static_cast<std::uint8_t>(i);
    }
    assert(speed_sketch(input) == UINT64_C(0x0080000802002000));

    for (std::uint64_t fingerprint : {UINT64_C(0), UINT64_C(63), UINT64_C(64),
                                      UINT64_C(0x123456789abcdef0), UINT64_MAX}) {
        assert(speedsketch_contract::fingerprint_class(fingerprint) ==
               static_cast<std::size_t>(fingerprint & 63U));
    }

    const std::vector<std::uint8_t> current{0, 1, 2, 3, 4, 5, 6, 7};
    const std::vector<std::uint8_t> base{255, 254, 253, 252, 251, 250, 249, 248};
    std::uint64_t fingerprint = 0;
    for (std::uint8_t byte : current) {
        fingerprint = speedsketch_contract::gear_step(fingerprint, byte);
    }
    const std::size_t bit = speedsketch_contract::fingerprint_class(fingerprint);

    GdeltaWorkspace workspace;
    GdeltaStats stats;
    std::vector<std::uint8_t> delta(64);
    std::size_t delta_size = 0;
    assert(gencode_whash_into(current.data(), current.size(), base.data(), base.size(),
                              0, UINT64_C(1) << bit, delta.data(), delta.size(),
                              &delta_size, &workspace, &stats) == GDELTA_OK);
    assert(stats.judgments == 1);
    assert(stats.excluded_lookups == 1);
    assert(stats.false_positive_lookups == 0);

    const std::size_t other_bit = (bit + 1) % speedsketch_contract::kSketchBits;
    assert(gencode_whash_into(current.data(), current.size(), base.data(), base.size(),
                              0, UINT64_C(1) << other_bit, delta.data(), delta.size(),
                              &delta_size, &workspace, &stats) == GDELTA_OK);
    assert(stats.judgments == 1);
    assert(stats.excluded_lookups == 0);
    assert(stats.false_positive_lookups == 1);
}

void test_forced_xxh64_collision() {
    const std::vector<std::uint8_t> first(32U * 1024U, 0x11);
    const std::vector<std::uint8_t> different(32U * 1024U, 0x22);
    const std::vector<std::uint8_t> small(2U * 1024U, 0x33);
    TempFile file;
    file.append(first);
    file.append(different);
    file.append(small);
    file.append(small);
    file.append(first);

    Deduplicator dedup(file.fd);
    const DeduplicatedChunk one = dedup.process(make_chunk(0, first));
    const DeduplicatedChunk two = dedup.process(make_chunk(first.size(), different));
    const std::uint64_t small_offset = first.size() + different.size();
    const DeduplicatedChunk three = dedup.process(make_chunk(small_offset, small));
    const DeduplicatedChunk four = dedup.process(make_chunk(small_offset + small.size(), small));
    const DeduplicatedChunk five = dedup.process(
        make_chunk(small_offset + 2 * small.size(), first));

    assert(!one.duplicate && one.base_read_bytes == 0);
    assert(!two.duplicate && two.base_read_bytes == first.size());
    assert(!three.duplicate && three.base_read_bytes == 0);
    assert(four.duplicate && four.base_read_bytes == small.size());
    assert(five.duplicate && five.base_read_bytes == first.size());
}

void test_fastcdc_bounds() {
    static_assert(kMinChunk == 4U * 1024U, "the production FastCDC minimum is 4 KiB");
    std::vector<std::uint8_t> data(kMaxChunk, 0);
    assert(fastcdc_cut(data.data(), 1) == 1);
    assert(fastcdc_cut(data.data(), kMinChunk - 1) == kMinChunk - 1);
    assert(fastcdc_cut(data.data(), kMinChunk) == kMinChunk);
    const std::size_t near_min_cut = fastcdc_cut(data.data(), kMinChunk + 1);
    assert(near_min_cut >= kMinChunk && near_min_cut <= kMinChunk + 1);
    const std::size_t cut = fastcdc_cut(data.data(), data.size());
    assert(cut >= kMinChunk && cut <= kMaxChunk);
}

void test_chunk_reader_roundtrip() {
    std::vector<std::uint8_t> input(kReadSize + kMaxChunk + 17);
    for (std::size_t i = 0; i < input.size(); ++i) {
        input[i] = static_cast<std::uint8_t>((i * 131U + 17U) & 0xffU);
    }
    TempFile file;
    file.append(input);
    assert(::lseek(file.fd, 0, SEEK_SET) == 0);

    ChunkReader reader(file.fd, input.size());
    std::vector<std::uint8_t> restored;
    std::uint64_t expected_offset = 0;
    reader.run([&](Chunk chunk) {
        assert(chunk.offset == expected_offset);
        assert(!chunk.data.empty() && chunk.data.size() <= kMaxChunk);
        expected_offset += chunk.data.size();
        restored.insert(restored.end(), chunk.data.begin(), chunk.data.end());
        return true;
    });
    assert(reader.source_read_bytes() == input.size());
    assert(expected_offset == input.size());
    assert(restored == input);
}

void test_parallel_encoding_failure_cancels() {
    const std::vector<std::uint8_t> input(16U * 1024U * 1024U, 0x5a);
    TempFile file;
    file.append(input);
    assert(::lseek(file.fd, 0, SEEK_SET) == 0);

    Options options;
    injected_encoding_failure_offset = 0;
    bool caught = false;
    try {
        (void)run_parallel(file.fd, options, input.size());
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "injected encoding failure";
    }
    injected_encoding_failure_offset = std::numeric_limits<std::uint64_t>::max();
    assert(caught);
}

void test_duplicate_does_not_reach_encoder() {
    const std::vector<std::uint8_t> seed(kMaxChunk, 0x5a);
    assert(fastcdc_cut(seed.data(), seed.size()) == seed.size());
    std::vector<std::uint8_t> input = seed;
    input.insert(input.end(), seed.begin(), seed.end());

    for (Pipeline pipeline : {Pipeline::Sequential, Pipeline::Parallel}) {
        TempFile file;
        file.append(input);
        assert(::lseek(file.fd, 0, SEEK_SET) == 0);
        Options options;
        options.pipeline = pipeline;
        injected_encoding_failure_offset = seed.size();
        const Metrics metrics = pipeline == Pipeline::Sequential
            ? run_sequential(file.fd, options, input.size())
            : run_parallel(file.fd, options, input.size());
        assert(metrics.chunks == 2);
        assert(metrics.unique_chunks == 1);
        assert(metrics.duplicate_chunks == 1);
    }
    injected_encoding_failure_offset = std::numeric_limits<std::uint64_t>::max();
}

void test_encoder_scratch_reuse() {
    const std::vector<std::uint8_t> base_big = random_bytes(32U * 1024U, 7);
    const std::vector<std::uint8_t> base_small = random_bytes(2U * 1024U, 11);
    std::vector<std::uint8_t> current_big = base_big;
    std::vector<std::uint8_t> current_small = base_small;
    current_big[current_big.size() / 2] ^= 1;
    current_small[current_small.size() / 2] ^= 1;

    for (Scheme scheme : {Scheme::SsG, Scheme::SsX}) {
        TempFile file;
        file.append(base_big);
        file.append(base_small);
        Options options;
        options.scheme = scheme;
        Encoder encoder(file.fd, options);
        for (const ClassifiedChunk& source : {
                 make_similar_chunk(100000, current_big, 0, base_big),
                 make_similar_chunk(200000, current_small, base_big.size(), base_small),
                 make_similar_chunk(300000, current_big, 0, base_big)}) {
            const EncodedChunk result = encoder.process(source);
            assert(result.choice == Choice::Delta);
            assert(result.verification_checks == 1);
            assert(result.verification_failures == 0);
        }
    }

    TempFile file;
    Options options;
    Encoder encoder(file.fd, options);
    const std::vector<std::uint8_t> self_big(32U * 1024U, 0x5a);
    const std::vector<std::uint8_t> self_small(2U * 1024U, 0x5a);
    for (const std::vector<std::uint8_t>* data : {&self_big, &self_small, &self_big}) {
        ClassifiedChunk item;
        item.chunk = make_chunk(0, *data);
        const EncodedChunk result = encoder.process(std::move(item));
        assert(result.choice == Choice::Self);
        assert(result.verification_checks == 1);
        assert(result.verification_failures == 0);
    }
}

}  // namespace

int main() {
    keep_embedded_cli_symbols_used();
    test_fixed_sketch_and_contract();
    test_forced_xxh64_collision();
    test_fastcdc_bounds();
    test_chunk_reader_roundtrip();
    test_parallel_encoding_failure_cancels();
    test_duplicate_does_not_reach_encoder();
    test_encoder_scratch_reuse();
    std::cout << "core tests passed\n";
}
