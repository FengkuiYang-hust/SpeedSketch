#ifndef SPEEDSKETCH_SKETCH_CONTRACT_H
#define SPEEDSKETCH_SKETCH_CONTRACT_H

#include <cstddef>
#include <cstdint>

// Defined once by gear_matrix.h in the Gdelta translation unit.
extern std::uint64_t GEARmx[256];

namespace speedsketch_contract {

constexpr std::size_t kSketchBits = 64;

constexpr std::uint64_t kSketchMasks[kSketchBits] = {
    12448694272ULL, 15669919744ULL, 51774488576ULL, 37782290432ULL,
    60934848512ULL, 11257511936ULL, 14193524736ULL, 28487712768ULL,
    66706210816ULL, 31977373696ULL, 15300820992ULL, 56992202752ULL,
    17095983104ULL,   822083584ULL, 66706210816ULL, 27665629184ULL,
    40282095616ULL, 51959037952ULL, 58200162304ULL, 46640660480ULL,
    49023025152ULL, 63149441024ULL, 64307068928ULL, 25618808832ULL,
    57797509120ULL, 26105348096ULL, 28890365952ULL, 54626615296ULL,
    62679678976ULL, 63786975232ULL, 31574720512ULL, 55448698880ULL,
    23051894784ULL,  9110028288ULL, 58032390144ULL, 64256737280ULL,
    57310969856ULL, 55348035584ULL, 10636754944ULL, 16374562816ULL,
    49056579584ULL, 14344519680ULL, 11962155008ULL, 25283264512ULL,
    54660169728ULL, 67796729856ULL, 50482642944ULL, 48133832704ULL,
    29142024192ULL, 45113933824ULL, 34024194048ULL, 40114323456ULL,
    19746783232ULL, 19528679424ULL, 25182601216ULL, 15586033664ULL,
    10619977728ULL, 53687091200ULL, 23521656832ULL, 65162706944ULL,
    43637538816ULL, 61706600448ULL, 39409680384ULL, 21407727616ULL,
};

inline std::uint64_t gear_step(std::uint64_t fingerprint, std::uint8_t byte) {
    return (fingerprint >> 8) + GEARmx[byte];
}

inline std::size_t fingerprint_class(std::uint64_t fingerprint) {
    return static_cast<std::size_t>(fingerprint & (kSketchBits - 1));
}

inline bool is_proxy(std::uint64_t fingerprint) {
    return (fingerprint & kSketchMasks[fingerprint_class(fingerprint)]) == 0;
}

inline bool sketch_bits_match(std::uint64_t fingerprint,
                              std::uint64_t left,
                              std::uint64_t right) {
    return ((left ^ right) & (std::uint64_t{1} << fingerprint_class(fingerprint))) == 0;
}

}  // namespace speedsketch_contract

#endif
