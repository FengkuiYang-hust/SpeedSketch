#ifndef GDELTA_GDELTA_H
#define GDELTA_GDELTA_H

#include <cstddef>
#include <cstdint>
#include <vector>

// Kept for the legacy gdeltaLshift.cpp translation unit.
#define ChunkSize (300 * 1024)
#define INIT_BUFFER_SIZE (128 * 1024)
#define FPTYPE uint64_t
#define WordSize 8
#define BaseSampleRate 2
#define PRINT_PERF 0
#define DEBUG_UNITS 0

const int hashLength = 64;

enum GdeltaStatus {
    GDELTA_OK = 0,
    GDELTA_INVALID_ARGUMENT = -1,
    GDELTA_OUTPUT_TOO_SMALL = -2,
    GDELTA_INVALID_DELTA = -3,
    GDELTA_ALLOCATION_FAILED = -4,
};

struct GdeltaStats {
    std::size_t judgments = 0;
    std::size_t excluded_lookups = 0;
    std::size_t false_positive_lookups = 0;
};

// Keep one workspace per calling thread. Its allocations are retained between
// calls and grow only when a larger chunk requires more capacity.
struct GdeltaWorkspace {
    std::vector<std::uint8_t> data;
    std::vector<std::uint8_t> instructions;
    std::vector<std::uint32_t> hash_table;
};

// On success, *_size is the produced length. If the caller buffer is too
// small, GDELTA_OUTPUT_TOO_SMALL is returned and *_size is the required length.
// Any supplied stats object is reset at the start of each encode call.
int gencode_into(const std::uint8_t* new_buf, std::size_t new_size,
                 const std::uint8_t* base_buf, std::size_t base_size,
                 std::uint8_t* delta_buf, std::size_t delta_capacity,
                 std::size_t* delta_size, GdeltaWorkspace* workspace,
                 GdeltaStats* stats = nullptr);

int gencode_whash_into(const std::uint8_t* new_buf, std::size_t new_size,
                       const std::uint8_t* base_buf, std::size_t base_size,
                       std::uint64_t new_hash, std::uint64_t base_hash,
                       std::uint8_t* delta_buf, std::size_t delta_capacity,
                       std::size_t* delta_size, GdeltaWorkspace* workspace,
                       GdeltaStats* stats = nullptr);

int gdecode_into(const std::uint8_t* delta_buf, std::size_t delta_size,
                 const std::uint8_t* base_buf, std::size_t base_size,
                 std::uint8_t* out_buf, std::size_t out_capacity,
                 std::size_t* out_size);

// Legacy allocation-capable wrappers. New code should use the bounded APIs
// above; successful legacy calls continue to return the produced byte count.
int gencode(const std::uint8_t* new_buf, std::uint32_t new_size,
            const std::uint8_t* base_buf, std::uint32_t base_size,
            std::uint8_t** delta_buf, std::uint32_t* delta_size);

int gencodeWHash(const std::uint8_t* new_buf, std::uint32_t new_size,
                 const std::uint8_t* base_buf, std::uint32_t base_size,
                 std::uint8_t** delta_buf, std::uint32_t* delta_size,
                 std::uint64_t new_hash, std::uint64_t base_hash);

int gdecode(const std::uint8_t* delta_buf, std::uint32_t delta_size,
            const std::uint8_t* base_buf, std::uint32_t base_size,
            std::uint8_t** out_buf, std::uint32_t* out_size);

#endif
