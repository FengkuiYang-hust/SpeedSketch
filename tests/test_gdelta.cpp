#include "gdelta.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

namespace {

using Bytes = std::vector<std::uint8_t>;

Bytes random_bytes(std::size_t size, std::uint32_t seed) {
    std::mt19937 generator(seed);
    Bytes bytes(size);
    for (auto& byte : bytes) {
        byte = static_cast<std::uint8_t>(generator());
    }
    return bytes;
}

Bytes encode(const Bytes& input,
             const Bytes& base,
             bool sketch,
             std::uint64_t input_hash,
             std::uint64_t base_hash,
             GdeltaWorkspace& workspace,
             GdeltaStats* stats = nullptr) {
    Bytes delta(input.size() * 2 + 32);
    std::size_t delta_size = 0;
    const int status = sketch
        ? gencode_whash_into(input.data(), input.size(), base.data(), base.size(),
                             input_hash, base_hash, delta.data(), delta.size(),
                             &delta_size, &workspace, stats)
        : gencode_into(input.data(), input.size(), base.data(), base.size(),
                       delta.data(), delta.size(), &delta_size, &workspace, stats);
    assert(status == GDELTA_OK);
    delta.resize(delta_size);
    return delta;
}

void assert_roundtrip(const Bytes& input,
                      const Bytes& base,
                      bool sketch,
                      std::uint64_t input_hash,
                      std::uint64_t base_hash) {
    GdeltaWorkspace workspace;
    GdeltaStats stats{99, 99, 99};
    const Bytes delta = encode(input, base, sketch, input_hash, base_hash,
                               workspace, &stats);

    Bytes output(input.size());
    std::size_t output_size = 0;
    assert(gdecode_into(delta.data(), delta.size(), base.data(), base.size(),
                        output.data(), output.size(), &output_size) == GDELTA_OK);
    assert(output_size == input.size());
    assert(output == input);

    if (!sketch) {
        assert(stats.judgments == 0);
        assert(stats.excluded_lookups == 0);
        assert(stats.false_positive_lookups == 0);
    }

    std::size_t required = 0;
    assert(gencode_into(input.data(), input.size(), base.data(), base.size(),
                        nullptr, 0, &required, &workspace) ==
           GDELTA_OUTPUT_TOO_SMALL);
    assert(required != 0);

    if (!input.empty()) {
        output_size = 0;
        assert(gdecode_into(delta.data(), delta.size(), base.data(), base.size(),
                            output.data(), input.size() - 1, &output_size) ==
               GDELTA_OUTPUT_TOO_SMALL);
        assert(output_size == input.size());
    }
}

void test_roundtrips() {
    Bytes base = random_bytes(32 * 1024, 7);
    Bytes identical = base;
    Bytes single_change = base;
    single_change[single_change.size() / 2] ^= 0x5a;
    Bytes random = random_bytes(base.size(), 19);

    for (bool sketch : {false, true}) {
        assert_roundtrip(identical, base, sketch, 0x1234, 0x1234);
        assert_roundtrip(single_change, base, sketch, 0x1234, 0x1235);
        assert_roundtrip(random, base, sketch, 0, ~std::uint64_t{0});
    }
}

void test_sketch_stats_and_reuse() {
    const Bytes base = random_bytes(16 * 1024, 23);
    const Bytes unrelated = random_bytes(base.size(), 29);
    GdeltaWorkspace workspace;
    GdeltaStats stats;

    encode(unrelated, base, true, 0, ~std::uint64_t{0}, workspace, &stats);
    assert(stats.judgments > 0);
    assert(stats.excluded_lookups == stats.judgments);
    assert(stats.false_positive_lookups == 0);
    assert(std::any_of(workspace.hash_table.begin(), workspace.hash_table.end(),
                       [](std::uint32_t entry) { return entry != 0; }));

    const auto data_capacity = workspace.data.capacity();
    const auto instruction_capacity = workspace.instructions.capacity();
    const auto hash_capacity = workspace.hash_table.capacity();
    encode(unrelated, base, true, 0, 0, workspace, &stats);
    assert(stats.judgments > 0);
    assert(stats.excluded_lookups == 0);
    assert(stats.false_positive_lookups > 0);
    assert(workspace.data.capacity() == data_capacity);
    assert(workspace.instructions.capacity() == instruction_capacity);
    assert(workspace.hash_table.capacity() == hash_capacity);
}

void test_offset_zero_sentinel() {
    Bytes base = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Bytes input = {42, 1, 2, 3, 4, 5, 6, 7, 8};
    GdeltaWorkspace workspace;
    const Bytes delta = encode(input, base, false, 0, 0, workspace);
    assert(std::find(workspace.hash_table.begin(), workspace.hash_table.end(), 1) !=
           workspace.hash_table.end());
    assert(delta.size() < input.size());

    Bytes output(input.size());
    std::size_t output_size = 0;
    assert(gdecode_into(delta.data(), delta.size(), base.data(), base.size(),
                        output.data(), output.size(), &output_size) == GDELTA_OK);
    assert(output == input);
}

void test_stats_status_and_retry() {
    const Bytes base = random_bytes(4096, 37);
    const Bytes input = random_bytes(4096, 41);
    GdeltaWorkspace workspace;
    GdeltaStats stats{99, 99, 99};
    std::size_t required = 0;
    assert(gencode_whash_into(input.data(), input.size(), base.data(), base.size(),
                              0, ~std::uint64_t{0}, nullptr, 0, &required,
                              &workspace, &stats) == GDELTA_OUTPUT_TOO_SMALL);
    assert(required > 0);
    assert(stats.judgments > 0);
    assert(stats.excluded_lookups == stats.judgments);
    assert(stats.false_positive_lookups == 0);
    const GdeltaStats query_stats = stats;

    Bytes delta(required);
    std::size_t actual = 0;
    stats = {99, 99, 99};
    assert(gencode_whash_into(input.data(), input.size(), base.data(), base.size(),
                              0, ~std::uint64_t{0}, delta.data(), delta.size(), &actual,
                              &workspace, &stats) == GDELTA_OK);
    assert(actual == required);
    assert(stats.judgments == query_stats.judgments);
    assert(stats.excluded_lookups == query_stats.excluded_lookups);
    assert(stats.false_positive_lookups == query_stats.false_positive_lookups);

    stats = {99, 99, 99};
    assert(gencode_whash_into(input.data(), input.size(), base.data(), base.size(),
                              0, 0, delta.data(), delta.size(), nullptr,
                              &workspace, &stats) == GDELTA_INVALID_ARGUMENT);
    assert(stats.judgments == 0);
    assert(stats.excluded_lookups == 0);
    assert(stats.false_positive_lookups == 0);
}

void test_truncated_and_corrupt_delta() {
    const Bytes base = random_bytes(4096, 31);
    Bytes input = base;
    input[100] ^= 1;
    GdeltaWorkspace workspace;
    const Bytes delta = encode(input, base, false, 0, 0, workspace);
    Bytes output(input.size());

    for (std::size_t size = 1; size < delta.size(); ++size) {
        std::size_t output_size = 0;
        assert(gdecode_into(delta.data(), size, base.data(), base.size(),
                            output.data(), output.size(), &output_size) < 0);
    }

    // Instruction stream of two bytes: copy one byte from offset base.size().
    const Bytes invalid_copy = {4, 5, 18};
    std::size_t output_size = 0;
    assert(gdecode_into(invalid_copy.data(), invalid_copy.size(),
                        base.data(), 9, output.data(), output.size(),
                        &output_size) == GDELTA_INVALID_DELTA);
}

}  // namespace

int main() {
    test_roundtrips();
    test_sketch_stats_and_reuse();
    test_offset_zero_sentinel();
    test_stats_status_and_retry();
    test_truncated_and_corrupt_delta();
    std::cout << "gdelta tests passed\n";
}
