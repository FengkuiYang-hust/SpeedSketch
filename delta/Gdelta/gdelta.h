#ifndef GDELTA_GDELTA_H
#define GDELTA_GDELTA_H
using namespace std;
#include <iostream>
#include <cstdint>
#include <bitset>

/*****Parameter*****/
#define ChunkSize (300 * 1024)
#define INIT_BUFFER_SIZE (128 * 1024)
#define FPTYPE uint64_t
//#define FPTYPE uint32_t
#define WordSize 8
// #define SkipStep 2
// #define SkipOn  // 如果有2^(skipstep)位的不匹配，则快速跳过后面的bytes
#define BaseSampleRate 2
// #define ReverseMatch
/*****Parameter*****/

#define PRINT_PERF 0
#define DEBUG_UNITS 0

const int hashLength = 64;

int gencode(uint8_t *newBuf, uint32_t newSize, uint8_t *baseBuf,
            uint32_t baseSize, uint8_t **deltaBuf, uint32_t *deltaSize);

int gencodeWHash(uint8_t *newBuf, uint32_t newSize, uint8_t *baseBuf,
            uint32_t baseSize, uint8_t **deltaBuf, uint32_t *deltaSize, uint64_t newHash, uint64_t baseHash);


int gdecode(uint8_t *deltaBuf, uint32_t deltaSize, uint8_t *baseBuf,
            uint32_t baseSize, uint8_t **outBuf, uint32_t *outSize);



#endif // GDELTA_GDELTA_H
