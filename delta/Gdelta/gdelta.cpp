#include "gdelta.h"
#include "gear_matrix.h"
#include "sketch_contract.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace {

constexpr std::size_t kWordSize = WordSize;
constexpr std::size_t kLegacyCapacity = ChunkSize;
constexpr std::size_t kLegacyInitialCapacity = INIT_BUFFER_SIZE;
constexpr std::uint8_t kVarintPayloadMask = 0x7f;
constexpr std::uint8_t kHeadLengthMask = 0x3f;
// The original on-disk format shifts the continuation by FLAGLEN (one bit),
// even though the first byte physically reserves six bits for length.
constexpr unsigned kLegacyHeadValueBits = 1;

struct DeltaUnit {
    bool copy = false;
    std::uint64_t length = 0;
    std::uint64_t offset = 0;
};

struct Reader {
    const std::uint8_t* data;
    std::size_t size;
    std::size_t cursor;
};

struct DeltaLayout {
    std::size_t instruction_offset;
    std::size_t instruction_size;
    std::size_t literal_offset;
    std::size_t output_size;
};

bool valid_buffer(const std::uint8_t* buffer, std::size_t size) {
    return buffer != nullptr || size == 0;
}

bool checked_add(std::size_t left, std::size_t right, std::size_t* result) {
    if (right > std::numeric_limits<std::size_t>::max() - left) {
        return false;
    }
    *result = left + right;
    return true;
}

bool equal_word(const std::uint8_t* left, const std::uint8_t* right) {
    std::uint64_t left_word;
    std::uint64_t right_word;
    std::memcpy(&left_word, left, sizeof(left_word));
    std::memcpy(&right_word, right, sizeof(right_word));
    return left_word == right_word;
}

std::uint64_t gear_fingerprint(const std::uint8_t* data) {
    std::uint64_t fingerprint = 0;
    for (std::size_t i = 0; i < kWordSize; ++i) {
        fingerprint = speedsketch_contract::gear_step(fingerprint, data[i]);
    }
    return fingerprint;
}

std::uint64_t roll_fingerprint(std::uint64_t fingerprint, std::uint8_t next) {
    return speedsketch_contract::gear_step(fingerprint, next);
}

bool sketch_bit_matches(std::uint64_t fingerprint,
                        std::uint64_t new_hash,
                        std::uint64_t base_hash) {
    return speedsketch_contract::sketch_bits_match(fingerprint, new_hash, base_hash);
}

void write_varint(std::vector<std::uint8_t>& output, std::uint64_t value) {
    do {
        std::uint8_t byte = static_cast<std::uint8_t>((value & kVarintPayloadMask) << 1);
        value >>= 7;
        if (value != 0) {
            byte |= 1;
        }
        output.push_back(byte);
    } while (value != 0);
}

std::size_t write_varint(std::uint8_t output[10], std::uint64_t value) {
    std::size_t size = 0;
    do {
        std::uint8_t byte = static_cast<std::uint8_t>((value & kVarintPayloadMask) << 1);
        value >>= 7;
        if (value != 0) {
            byte |= 1;
        }
        output[size++] = byte;
    } while (value != 0);
    return size;
}

void write_unit(std::vector<std::uint8_t>& instructions, const DeltaUnit& unit) {
    const bool more = unit.length > 1;
    const std::uint8_t head =
        static_cast<std::uint8_t>((unit.copy ? 1 : 0) |
                                  (more ? 2 : 0) |
                                  ((unit.length & 1) << 2));
    instructions.push_back(head);
    if (more) {
        write_varint(instructions, unit.length >> kLegacyHeadValueBits);
    }
    if (unit.copy) {
        write_varint(instructions, unit.offset);
    }
}

bool read_byte(Reader& input, std::uint8_t* value) {
    if (input.cursor >= input.size) {
        return false;
    }
    *value = input.data[input.cursor++];
    return true;
}

bool read_varint(Reader& input, std::uint64_t* value) {
    std::uint64_t result = 0;
    unsigned shift = 0;
    for (unsigned part = 0; part < 10; ++part, shift += 7) {
        std::uint8_t byte;
        if (!read_byte(input, &byte)) {
            return false;
        }
        const std::uint64_t payload = byte >> 1;
        if (payload > (std::numeric_limits<std::uint64_t>::max() >> shift)) {
            return false;
        }
        result |= payload << shift;
        if ((byte & 1) == 0) {
            *value = result;
            return true;
        }
    }
    return false;
}

bool read_unit(Reader& instructions, DeltaUnit* unit) {
    std::uint8_t head;
    if (!read_byte(instructions, &head)) {
        return false;
    }

    unit->copy = (head & 1) != 0;
    unit->length = (head >> 2) & kHeadLengthMask;
    unit->offset = 0;
    if ((head & 2) != 0) {
        std::uint64_t remaining;
        if (!read_varint(instructions, &remaining) ||
            remaining > (std::numeric_limits<std::uint64_t>::max() >>
                         kLegacyHeadValueBits)) {
            return false;
        }
        unit->length |= remaining << kLegacyHeadValueBits;
    }
    if (unit->copy && !read_varint(instructions, &unit->offset)) {
        return false;
    }
    return unit->length != 0;
}

void emit_copy(GdeltaWorkspace& workspace, std::size_t offset, std::size_t length) {
    if (length == 0) {
        return;
    }
    write_unit(workspace.instructions,
               {true, static_cast<std::uint64_t>(length),
                static_cast<std::uint64_t>(offset)});
}

void emit_literal(GdeltaWorkspace& workspace,
                  const std::uint8_t* input,
                  std::size_t offset,
                  std::size_t length) {
    if (length == 0) {
        return;
    }
    write_unit(workspace.instructions,
               {false, static_cast<std::uint64_t>(length), 0});
    workspace.data.insert(workspace.data.end(), input + offset, input + offset + length);
}

int build_hash_table(const std::uint8_t* base,
                     std::size_t begin,
                     std::size_t end,
                     GdeltaWorkspace& workspace) {
    const std::size_t length = end - begin;
    if (length < kWordSize) {
        workspace.hash_table.clear();
        return GDELTA_OK;
    }

    std::size_t desired;
    if (!checked_add(length, 10, &desired)) {
        return GDELTA_INVALID_ARGUMENT;
    }
    std::size_t hash_size = 1;
    while (hash_size < desired) {
        if (hash_size > std::numeric_limits<std::size_t>::max() / 2) {
            return GDELTA_INVALID_ARGUMENT;
        }
        hash_size *= 2;
    }
    workspace.hash_table.assign(hash_size, 0);

    const std::size_t mask = hash_size - 1;
    std::size_t position = begin;
    std::uint64_t fingerprint = gear_fingerprint(base + position);
    while (position + kWordSize < end) {
        const std::size_t index = static_cast<std::size_t>(fingerprint) & mask;
        // offset+1 reserves zero as the empty sentinel and makes base offset zero usable.
        workspace.hash_table[index] = static_cast<std::uint32_t>(position + 1);
        ++position;
        fingerprint = roll_fingerprint(fingerprint, base[position + kWordSize - 1]);
    }
    return GDELTA_OK;
}

int inspect_delta(const std::uint8_t* delta,
                  std::size_t delta_size,
                  std::size_t base_size,
                  DeltaLayout* layout) {
    Reader header{delta, delta_size, 0};
    std::uint64_t instruction_size_u64;
    if (!read_varint(header, &instruction_size_u64) ||
        instruction_size_u64 > std::numeric_limits<std::size_t>::max()) {
        return GDELTA_INVALID_DELTA;
    }
    const std::size_t instruction_size = static_cast<std::size_t>(instruction_size_u64);
    if (instruction_size > delta_size - header.cursor) {
        return GDELTA_INVALID_DELTA;
    }

    const std::size_t literal_offset = header.cursor + instruction_size;
    Reader instructions{delta + header.cursor, instruction_size, 0};
    Reader literals{delta + literal_offset, delta_size - literal_offset, 0};
    std::size_t output_size = 0;
    while (instructions.cursor < instructions.size) {
        DeltaUnit unit;
        if (!read_unit(instructions, &unit) ||
            unit.length > std::numeric_limits<std::size_t>::max()) {
            return GDELTA_INVALID_DELTA;
        }
        const std::size_t length = static_cast<std::size_t>(unit.length);
        if (unit.copy) {
            if (unit.offset > base_size || length > base_size - unit.offset) {
                return GDELTA_INVALID_DELTA;
            }
        } else {
            if (length > literals.size - literals.cursor) {
                return GDELTA_INVALID_DELTA;
            }
            literals.cursor += length;
        }
        if (!checked_add(output_size, length, &output_size)) {
            return GDELTA_INVALID_DELTA;
        }
    }
    if (literals.cursor != literals.size) {
        return GDELTA_INVALID_DELTA;
    }

    *layout = {header.cursor, instruction_size, literal_offset, output_size};
    return GDELTA_OK;
}

int encode_impl(const std::uint8_t* new_buf,
                std::size_t new_size,
                const std::uint8_t* base_buf,
                std::size_t base_size,
                bool use_sketch,
                std::uint64_t new_hash,
                std::uint64_t base_hash,
                std::uint8_t* delta_buf,
                std::size_t delta_capacity,
                std::size_t* delta_size,
                GdeltaWorkspace* workspace,
                GdeltaStats* stats) {
    if (stats != nullptr) {
        *stats = {};
    }
    if (delta_size == nullptr || workspace == nullptr ||
        !valid_buffer(new_buf, new_size) || !valid_buffer(base_buf, base_size) ||
        (delta_buf == nullptr && delta_capacity != 0) ||
        new_size > std::numeric_limits<std::uint32_t>::max() ||
        base_size > std::numeric_limits<std::uint32_t>::max()) {
        return GDELTA_INVALID_ARGUMENT;
    }
    *delta_size = 0;

    try {
        workspace->data.clear();
        workspace->instructions.clear();
        if (workspace->data.capacity() < new_size) {
            workspace->data.reserve(new_size);
        }
        const std::size_t instruction_reserve = new_size + 32;
        if (workspace->instructions.capacity() < instruction_reserve) {
            workspace->instructions.reserve(instruction_reserve);
        }

        const std::size_t common_size = std::min(new_size, base_size);
        std::size_t prefix = 0;
        while (prefix + kWordSize <= common_size &&
               equal_word(base_buf + prefix, new_buf + prefix)) {
            prefix += kWordSize;
        }
        while (prefix < common_size && base_buf[prefix] == new_buf[prefix]) {
            ++prefix;
        }
        if (prefix <= 16) {
            prefix = 0;
        }

        std::size_t suffix = 0;
        const std::size_t suffix_limit = common_size - prefix;
        while (suffix + kWordSize <= suffix_limit &&
               equal_word(base_buf + base_size - suffix - kWordSize,
                          new_buf + new_size - suffix - kWordSize)) {
            suffix += kWordSize;
        }
        while (suffix < suffix_limit &&
               base_buf[base_size - suffix - 1] == new_buf[new_size - suffix - 1]) {
            ++suffix;
        }
        if (suffix <= 16) {
            suffix = 0;
        }

        emit_copy(*workspace, 0, prefix);

        const std::size_t base_end = base_size - suffix;
        const std::size_t new_end = new_size - suffix;
        int status = build_hash_table(base_buf, prefix, base_end, *workspace);
        if (status != GDELTA_OK) {
            return status;
        }

        std::size_t position = prefix;
        std::size_t literal_begin = position;
        if (!workspace->hash_table.empty()) {
            const std::size_t hash_mask = workspace->hash_table.size() - 1;
            bool have_fingerprint = false;
            std::uint64_t fingerprint = 0;
            while (position + kWordSize <= new_end) {
                if (!have_fingerprint) {
                    fingerprint = gear_fingerprint(new_buf + position);
                    have_fingerprint = true;
                }

                bool lookup = true;
                if (use_sketch) {
                    if (stats != nullptr) {
                        ++stats->judgments;
                    }
                    lookup = sketch_bit_matches(fingerprint, new_hash, base_hash);
                    if (!lookup && stats != nullptr) {
                        ++stats->excluded_lookups;
                    }
                }

                bool matched = false;
                std::size_t base_offset = 0;
                if (lookup) {
                    const std::uint32_t entry =
                        workspace->hash_table[static_cast<std::size_t>(fingerprint) & hash_mask];
                    if (entry != 0) {
                        base_offset = static_cast<std::size_t>(entry - 1);
                        matched = base_offset + kWordSize <= base_end &&
                                  std::memcmp(new_buf + position,
                                              base_buf + base_offset,
                                              kWordSize) == 0;
                    }
                    if (use_sketch && !matched && stats != nullptr) {
                        ++stats->false_positive_lookups;
                    }
                }

                if (!matched) {
                    ++position;
                    if (position + kWordSize <= new_end) {
                        fingerprint = roll_fingerprint(
                            fingerprint, new_buf[position + kWordSize - 1]);
                    }
                    continue;
                }

                const std::size_t literal_length = position - literal_begin;
                std::size_t backward = 0;
                while (backward < base_offset && backward < literal_length &&
                       base_buf[base_offset - backward - 1] ==
                           new_buf[position - backward - 1]) {
                    ++backward;
                }
                emit_literal(*workspace, new_buf, literal_begin,
                             literal_length - backward);
                std::size_t match_length = kWordSize;
                while (base_offset + match_length + kWordSize <= base_end &&
                       position + match_length + kWordSize <= new_end &&
                       equal_word(base_buf + base_offset + match_length,
                                  new_buf + position + match_length)) {
                    match_length += kWordSize;
                }
                while (base_offset + match_length < base_end &&
                       position + match_length < new_end &&
                       base_buf[base_offset + match_length] ==
                           new_buf[position + match_length]) {
                    ++match_length;
                }
                emit_copy(*workspace, base_offset - backward,
                          match_length + backward);
                position += match_length;
                literal_begin = position;
                have_fingerprint = false;
            }
        }
        emit_literal(*workspace, new_buf, literal_begin, new_end - literal_begin);
        emit_copy(*workspace, base_size - suffix, suffix);

        std::uint8_t header[10];
        const std::size_t header_size =
            write_varint(header, workspace->instructions.size());
        std::size_t required;
        if (!checked_add(header_size, workspace->instructions.size(), &required) ||
            !checked_add(required, workspace->data.size(), &required)) {
            return GDELTA_INVALID_ARGUMENT;
        }
        *delta_size = required;
        if (required > delta_capacity || (required != 0 && delta_buf == nullptr)) {
            return GDELTA_OUTPUT_TOO_SMALL;
        }

        std::size_t cursor = 0;
        std::memcpy(delta_buf + cursor, header, header_size);
        cursor += header_size;
        if (!workspace->instructions.empty()) {
            std::memcpy(delta_buf + cursor, workspace->instructions.data(),
                        workspace->instructions.size());
            cursor += workspace->instructions.size();
        }
        if (!workspace->data.empty()) {
            std::memcpy(delta_buf + cursor, workspace->data.data(), workspace->data.size());
        }
        return GDELTA_OK;
    } catch (const std::bad_alloc&) {
        return GDELTA_ALLOCATION_FAILED;
    } catch (const std::length_error&) {
        return GDELTA_ALLOCATION_FAILED;
    }
}

template <typename Encode>
int legacy_encode(Encode encode,
                  std::uint32_t new_size,
                  std::uint8_t** delta_buf,
                  std::uint32_t* delta_size) {
    if (delta_buf == nullptr || delta_size == nullptr) {
        return GDELTA_INVALID_ARGUMENT;
    }
    *delta_size = 0;
    bool allocated = false;
    std::size_t capacity = kLegacyCapacity;
    if (*delta_buf == nullptr) {
        const std::size_t conservative = static_cast<std::size_t>(new_size) * 2 + 32;
        capacity = std::max(kLegacyInitialCapacity, conservative);
        *delta_buf = static_cast<std::uint8_t*>(std::malloc(capacity));
        if (*delta_buf == nullptr) {
            return GDELTA_ALLOCATION_FAILED;
        }
        allocated = true;
    }

    std::size_t actual = 0;
    int status = encode(*delta_buf, capacity, &actual);
    if (status == GDELTA_OUTPUT_TOO_SMALL && allocated) {
        void* resized = std::realloc(*delta_buf, actual);
        if (resized == nullptr) {
            std::free(*delta_buf);
            *delta_buf = nullptr;
            return GDELTA_ALLOCATION_FAILED;
        }
        *delta_buf = static_cast<std::uint8_t*>(resized);
        capacity = actual;
        status = encode(*delta_buf, capacity, &actual);
    }
    if (status != GDELTA_OK || actual > std::numeric_limits<std::uint32_t>::max() ||
        actual > static_cast<std::size_t>(INT_MAX)) {
        if (allocated) {
            std::free(*delta_buf);
            *delta_buf = nullptr;
        }
        return status == GDELTA_OK ? GDELTA_INVALID_ARGUMENT : status;
    }
    *delta_size = static_cast<std::uint32_t>(actual);
    return static_cast<int>(actual);
}

}  // namespace

int gencode_into(const std::uint8_t* new_buf, std::size_t new_size,
                 const std::uint8_t* base_buf, std::size_t base_size,
                 std::uint8_t* delta_buf, std::size_t delta_capacity,
                 std::size_t* delta_size, GdeltaWorkspace* workspace,
                 GdeltaStats* stats) {
    return encode_impl(new_buf, new_size, base_buf, base_size, false, 0, 0,
                       delta_buf, delta_capacity, delta_size, workspace, stats);
}

int gencode_whash_into(const std::uint8_t* new_buf, std::size_t new_size,
                       const std::uint8_t* base_buf, std::size_t base_size,
                       std::uint64_t new_hash, std::uint64_t base_hash,
                       std::uint8_t* delta_buf, std::size_t delta_capacity,
                       std::size_t* delta_size, GdeltaWorkspace* workspace,
                       GdeltaStats* stats) {
    return encode_impl(new_buf, new_size, base_buf, base_size, true,
                       new_hash, base_hash, delta_buf, delta_capacity,
                       delta_size, workspace, stats);
}

int gdecode_into(const std::uint8_t* delta_buf, std::size_t delta_size,
                 const std::uint8_t* base_buf, std::size_t base_size,
                 std::uint8_t* out_buf, std::size_t out_capacity,
                 std::size_t* out_size) {
    if (out_size == nullptr || !valid_buffer(delta_buf, delta_size) ||
        !valid_buffer(base_buf, base_size) ||
        (out_buf == nullptr && out_capacity != 0) || delta_size == 0) {
        return GDELTA_INVALID_ARGUMENT;
    }
    *out_size = 0;

    DeltaLayout layout;
    const int status = inspect_delta(delta_buf, delta_size, base_size, &layout);
    if (status != GDELTA_OK) {
        return status;
    }
    *out_size = layout.output_size;
    if (layout.output_size > out_capacity ||
        (layout.output_size != 0 && out_buf == nullptr)) {
        return GDELTA_OUTPUT_TOO_SMALL;
    }

    Reader instructions{delta_buf + layout.instruction_offset,
                        layout.instruction_size, 0};
    Reader literals{delta_buf + layout.literal_offset,
                    delta_size - layout.literal_offset, 0};
    std::size_t output_cursor = 0;
    while (instructions.cursor < instructions.size) {
        DeltaUnit unit;
        if (!read_unit(instructions, &unit)) {
            *out_size = 0;
            return GDELTA_INVALID_DELTA;
        }
        const std::size_t length = static_cast<std::size_t>(unit.length);
        if (unit.copy) {
            std::memcpy(out_buf + output_cursor,
                        base_buf + static_cast<std::size_t>(unit.offset), length);
        } else {
            std::memcpy(out_buf + output_cursor, literals.data + literals.cursor, length);
            literals.cursor += length;
        }
        output_cursor += length;
    }
    *out_size = output_cursor;
    return GDELTA_OK;
}

int gencode(const std::uint8_t* new_buf, std::uint32_t new_size,
            const std::uint8_t* base_buf, std::uint32_t base_size,
            std::uint8_t** delta_buf, std::uint32_t* delta_size) {
    thread_local GdeltaWorkspace workspace;
    return legacy_encode(
        [&](std::uint8_t* output, std::size_t capacity, std::size_t* actual) {
            return gencode_into(new_buf, new_size, base_buf, base_size,
                                output, capacity, actual, &workspace);
        },
        new_size, delta_buf, delta_size);
}

int gencodeWHash(const std::uint8_t* new_buf, std::uint32_t new_size,
                 const std::uint8_t* base_buf, std::uint32_t base_size,
                 std::uint8_t** delta_buf, std::uint32_t* delta_size,
                 std::uint64_t new_hash, std::uint64_t base_hash) {
    thread_local GdeltaWorkspace workspace;
    return legacy_encode(
        [&](std::uint8_t* output, std::size_t capacity, std::size_t* actual) {
            return gencode_whash_into(new_buf, new_size, base_buf, base_size,
                                      new_hash, base_hash, output, capacity,
                                      actual, &workspace);
        },
        new_size, delta_buf, delta_size);
}

int gdecode(const std::uint8_t* delta_buf, std::uint32_t delta_size,
            const std::uint8_t* base_buf, std::uint32_t base_size,
            std::uint8_t** out_buf, std::uint32_t* out_size) {
    if (out_buf == nullptr || out_size == nullptr) {
        return GDELTA_INVALID_ARGUMENT;
    }
    *out_size = 0;

    bool allocated = false;
    std::size_t capacity = kLegacyCapacity;
    if (*out_buf == nullptr) {
        std::size_t required = 0;
        const int query = gdecode_into(delta_buf, delta_size, base_buf, base_size,
                                       nullptr, 0, &required);
        if (query != GDELTA_OUTPUT_TOO_SMALL && query != GDELTA_OK) {
            return query;
        }
        capacity = std::max<std::size_t>(required, 1);
        *out_buf = static_cast<std::uint8_t*>(std::malloc(capacity));
        if (*out_buf == nullptr) {
            return GDELTA_ALLOCATION_FAILED;
        }
        allocated = true;
    }

    std::size_t actual = 0;
    const int status = gdecode_into(delta_buf, delta_size, base_buf, base_size,
                                    *out_buf, capacity, &actual);
    if (status != GDELTA_OK || actual > std::numeric_limits<std::uint32_t>::max() ||
        actual > static_cast<std::size_t>(INT_MAX)) {
        if (allocated) {
            std::free(*out_buf);
            *out_buf = nullptr;
        }
        return status == GDELTA_OK ? GDELTA_INVALID_ARGUMENT : status;
    }
    *out_size = static_cast<std::uint32_t>(actual);
    return static_cast<int>(actual);
}
