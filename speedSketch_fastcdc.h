/**
 * This is an implementation of fastCDC
 * The origin paper is Wen Xia, Yukun Zhou, Hong Jiang, Dan Feng, Yu Hua, Yuchong Hu, Yucheng Zhang, Qing Liu, "FastCDC: a Fast and Efficient Content-Defined Chunking Approach for Data Deduplication", in Proceedings of USENIX Annual Technical Conference (USENIX ATC'16), Denver, CO, USA, June 22–24, 2016, pages: 101-114
 * and Wen Xia, Xiangyu Zou, Yukun Zhou, Hong Jiang, Chuanyi Liu, Dan Feng, Yu Hua, Yuchong Hu, Yucheng Zhang, "The Design of Fast Content-Defined Chunking for Data Deduplication based Storage Systems", IEEE Transactions on Parallel and Distributed Systems (TPDS), 2020
 *
 */ 

#include <openssl/md5.h>
#include <openssl/sha.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <zlib.h>
#include "uthash.h"
#include "Gdelta/gear_matrix.h"

// typedef std::bitset<num_masks> MYHASH;
typedef uint64_t MYHASH;
typedef uint64_t SFP_type;
typedef uint64_t OHASH_FPTYPE;
const MYHASH subOhashMask = (MYHASH(0xffffffff));
uint64_t mask_mask = 0x0000000fff000000ULL;	// 18
// uint64_t mask_mask = 0x00000f0fff0f0000ULL;

constexpr int num_masks = 64;
const uint32_t mask = 0x2fa6;

int movebitlength;

uint64_t Masks[64] = {
    12448694272,       15669919744,       51774488576,       37782290432,       60934848512,       11257511936,       14193524736,       28487712768,
    66706210816,       31977373696,       15300820992,       56992202752,       17095983104,         822083584,       66706210816,       27665629184,
    40282095616,       51959037952,       58200162304,       46640660480,       49023025152,       63149441024,       64307068928,       25618808832,
    57797509120,       26105348096,       28890365952,       54626615296,       62679678976,       63786975232,       31574720512,       55448698880,
    23051894784,        9110028288,       58032390144,       64256737280,       57310969856,       55348035584,       10636754944,       16374562816,
    49056579584,       14344519680,       11962155008,       25283264512,       54660169728,       67796729856,       50482642944,       48133832704,
    29142024192,       45113933824,       34024194048,       40114323456,       19746783232,       19528679424,       25182601216,       15586033664,
    10619977728,       53687091200,       23521656832,       65162706944,       43637538816,       61706600448,       39409680384,       21407727616,
};

void initMasks(){
	uint64_t tempMasks[num_masks];
	for (int i = 0; i< num_masks; ++i){
		tempMasks[i] = GEARmx[4*i+3] & mask_mask;
	}	// 使用GEARmx进行掩码得到的mask，最后的结果也会与GEAEmx 有关系

// static const uint64_t tempMasks[64] = {
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,
// 	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask,	dis(gen) &  mask_mask
// };

// static const uint64_t tempMasks[num_masks] = {
	// 23051894784,        9110028288,       58032390144,       64256737280,       57310969856,       55348035584,       10636754944,       16374562816,
	// 49056579584,       14344519680,       11962155008,       25283264512,       54660169728,       67796729856,       50482642944,       48133832704,
	// 29142024192,       45113933824,       34024194048,       40114323456,       19746783232,       19528679424,       25182601216,       15586033664,
	// 10619977728,       53687091200,       23521656832,       65162706944,       43637538816,       61706600448,       39409680384,       21407727616,
	// 2617245696,       16424894464,       31239176192,       25048383488,       60632858624,       29947330560,       57780731904,       15065939968,
	// 48888807424,        5452595200,       54408511488,       48922361856,       37379637248,       64793608192,       45986349056,       40114323456,
	// 26944208896,        9948889088,       50616860672,       58200162304,        4076863488,       47479521280,       13019119616,       51472498688,
	// 12482248704,       43637538816,       48804921344,       50230984704,       65531805696,       30349983744,       26675773440,       14864613376
// };		// 结合两个DCE很高的mask

	std::memcpy(Masks, tempMasks, sizeof(Masks));
}


#define SymbolCount 256
#define SeedLength 64
#define CacheSize 1024 * 1024 * 1024

#define ORIGIN_CDC 1
#define ROLLING_2Bytes 2
#define NORMALIZED_CDC 3
#define NORMALIZED_2Bytes 4

// Rolling2Bytes Mask
uint32_t FING_GEAR_08KB_ls = 0xd9300353 << 1;
uint32_t FING_GEAR_02KB_ls = 0xd9000353 << 1;
uint32_t FING_GEAR_32KB_ls = 0xd9f00353 << 1;

uint64_t LEARv2[256];

uint64_t FING_GEAR_08KB_ls_64 = 0x0000d93003530000 << 1;
uint64_t FING_GEAR_02KB_ls_64 = 0x0000d90003530000 << 1;
uint64_t FING_GEAR_32KB_ls_64 = 0x0000d9f003530000 << 1;
uint64_t FING_GEAR_08KB_64 = 0x0000d93003530000;

uint64_t FING_GEAR_02KB_64 = 0x0000d90003530000;
uint64_t FING_GEAR_32KB_64 = 0x0000d9f003530000;

// global variants
struct timeval tmStart, tmEnd;
struct chunk_info *users = NULL;

float totalTm = 0;
int chunk_dist[30];
uint32_t g_global_matrix[SymbolCount];
uint32_t g_global_matrix_left[SymbolCount];
uint32_t expectCS;
uint32_t Mask_15;
uint32_t Mask_11;
uint64_t Mask_11_64, Mask_15_64;

uint32_t MinSize = 8192 / 4;
uint32_t MinSize_divide_by_2;
uint32_t MaxSize = 8192 * 4;
int sameCount = 0;
int tmpCount = 0;
int smalChkCnt = 0;  //记录小于8KB的分块

// init function
void fastCDC_init(void);

int (*chunking) (unsigned char*p, int n);

// origin fastcdc function
int cdc_origin_64(unsigned char *p, int n);

// fastcdc with once rolling 2 bytes 
int rolling_data_2byes_64(unsigned char *p, int n);

// normalized fastcdc
int normalized_chunking_64(unsigned char *p, int n);

// normalized fastcdc with once rolling 2 bytes
int normalized_chunking_2byes_64(unsigned char *p, int n);

// fix size chunking 
int fix_size_chunking(unsigned char *p, int n);

int normalized_chunking_64_with_features(unsigned char *p, int n, uint64_t &bitFeature);

// functions
void fastCDC_init(void) {
    // initMasks();

    movebitlength = sizeof(MYHASH) * 8 / WordSize;
    if (sizeof(OHASH_FPTYPE) * 8 % WordSize != 0)
        movebitlength++;

    unsigned char md5_digest[16];
    uint8_t seed[SeedLength];
    for (int i = 0; i < SymbolCount; i++) {

        for (int j = 0; j < SeedLength; j++) {
            seed[j] = i;
        }

        g_global_matrix[i] = 0;
        MD5(seed, SeedLength, md5_digest);
        memcpy(&(g_global_matrix[i]), md5_digest, 4);
        g_global_matrix_left[i] = g_global_matrix[i] << 1;
    }

    // 64 bit init
    for (int i = 0; i < SymbolCount; i++) {
        LEARv2[i] = GEARmx[i] << 1;
    }

    MinSize = 8192 / 4;
    MaxSize = 8192 * 4;    // 32768;
    Mask_15 = 0xf9070353;  //  15个1
    Mask_11 = 0xd9000353;  //  11个1
    Mask_11_64 = 0x0000d90003530000;
    Mask_15_64 = 0x0000f90703530000;
    MinSize_divide_by_2 = MinSize / 2;
}

int normalized_chunking_64(unsigned char *p, int n) {
    uint64_t fingerprint = 0, digest;
    MinSize = 6 * 1024;
    int i = MinSize, Mid = 8 * 1024;

    // the minimal subChunk Size.
    if (n <= MinSize)  
        return n;

    if (n > MaxSize)
        n = MaxSize;
    else if (n < Mid)
        Mid = n;

    while (i < Mid) {
        fingerprint = (fingerprint >> movebitlength) + (GEARmx[p[i]]);
        // fingerprint = (fingerprint << 1) + (GEARmx[p[i]]);

        if ((!(fingerprint & FING_GEAR_32KB_64))) {
            return i;
        }

        i++;
    }

    while (i < n) {
        fingerprint = (fingerprint >> movebitlength) + (GEARmx[p[i]]);
        // fingerprint = (fingerprint << 1) + (GEARmx[p[i]]);

        if ((!(fingerprint & FING_GEAR_02KB_64))) {
            return i;
        }

        i++;
    }

    return n;
}

int normalized_chunking_2byes_64(unsigned char *p, int n) {
    uint64_t fingerprint = 0, digest;
    MinSize = 6 * 1024;
    int i = MinSize / 2, Mid = 8 * 1024;

    // the minimal subChunk Size.
    if (n <= MinSize) 
        return n;

    if (n > MaxSize)
        n = MaxSize;
    else if (n < Mid)
        Mid = n;

    while (i < Mid / 2) {
        int a = i * 2;
        fingerprint = (fingerprint << 2) + (LEARv2[p[a]]);

        if ((!(fingerprint & FING_GEAR_32KB_ls_64))) {
            return a;
        }

        fingerprint += GEARmx [p[a + 1]];  

        if ((!(fingerprint & FING_GEAR_32KB_64))) {
            return a + 1;
        }

        i++;
    }

    while (i < n / 2) {
        int a = i * 2;
        fingerprint = (fingerprint << 2) + (LEARv2[p[a]]);

        if ((!(fingerprint & FING_GEAR_02KB_ls_64))) {
            return a;
        }

        fingerprint += GEARmx[p[a + 1]];

        if ((!(fingerprint & FING_GEAR_02KB_64))) {
            return a + 1;
        }

        i++;
    }

    return n;
}

int rolling_data_2byes_64(unsigned char *p, int n) {
    uint64_t fingerprint = 0, digest;
    int i = MinSize_divide_by_2;

    // the minimal subChunk Size.
    if (n <= MinSize) 
        return n;

    if (n > MaxSize) n = MaxSize;

    while (i < n / 2) {
        int a = i * 2;
        fingerprint = (fingerprint << 2) + (LEARv2[p[a]]);

        if ((!(fingerprint & FING_GEAR_08KB_ls_64))) {
            return a;
        }

        fingerprint += GEARmx[p[a + 1]];

        if ((!(fingerprint & FING_GEAR_08KB_64))) {
            return a + 1;
        }

        i++;
    }

    return n;
}

int cdc_origin_64(unsigned char *p, int n) {
    uint64_t fingerprint = 0, digest;
    int i = MinSize;
    // return n;
    // the minimal subChunk Size.
    if (n <= MinSize)  
        return n;

    if (n > MaxSize) n = MaxSize;

    while (i < n) {
        fingerprint = (fingerprint << 1) + (GEARmx[p[i]]);
        if ((!(fingerprint & FING_GEAR_08KB_64))) {
            return i;
        }
        i++;
    }

    return n;
}

int fix_size_chunking(unsigned char *p, int n) {
    if (n <= 8192)  
        return n;

    if (n > 8192) n = 8192;

    return n;
}


int normalized_chunking_64_with_features(unsigned char *p, int n, uint64_t &bitFeature) {
    uint64_t fingerprint = 0, digest;
    MinSize = 6 * 1024;
    uint64_t local_bitFeature = 0ULL;
    int i = 0, Min = 6 * 1024, Mid = 8 * 1024;
    uint32_t bit_mask = sizeof(uint64_t) * 8 - 1;
    uint32_t highbit = 0;

    // the minimal subChunk Size.
    if (n <= Min)  
        Min = n;
    else if (n > MaxSize)
        n = MaxSize;
    else if (n < Mid)
        Mid = n;

    while (i < Min) {
        fingerprint = (fingerprint >> movebitlength) + (GEARmx[p[i]]);

        // highbit = (fingerprint >> 58) & bit_mask;	// 左移使用低位，涉及更多的byte 58 = sizeof(uint64_t) * 8 - 6;
		highbit = (fingerprint >> 4) & bit_mask;		// 右移使用高位，涉及更多的byte
		if ((!(fingerprint & (Masks[highbit])))){
			local_bitFeature |= (1ULL << highbit);
			// bitMaskHash[highbit] = true;
            // highbitArray[highbit] = 1;
		}

        i++;
    }
    if (i < MinSize){
        bitFeature = local_bitFeature;
        return i;
    }

    while (i < Mid) {
        fingerprint = (fingerprint >> movebitlength) + (GEARmx[p[i]]);

        if ((!(fingerprint & FING_GEAR_32KB_64))) {
            bitFeature = local_bitFeature;
            return i;
        }

        // highbit = (fingerprint >> 58) & bit_mask;	// 左移使用低位，涉及更多的byte 58 = sizeof(uint64_t) * 8 - 6;
		highbit = (fingerprint >> 4) & bit_mask;		// 右移使用高位，涉及更多的byte
		if ((!(fingerprint & (Masks[highbit])))){
			local_bitFeature |= (1ULL << highbit);
			// bitMaskHash[highbit] = true;
            // highbitArray[highbit] = 1;
		}

        i++;
    }

    while (i < n) {
        fingerprint = (fingerprint >> movebitlength) + (GEARmx[p[i]]);

        if ((!(fingerprint & FING_GEAR_02KB_64))) {
            bitFeature = local_bitFeature;
            return i;
        }

        // highbit = (fingerprint >> 58) & bit_mask;	// 左移使用低位，涉及更多的byte 58 = sizeof(uint64_t) * 8 - 6;
		highbit = (fingerprint >> 4) & bit_mask;		// 右移使用高位，涉及更多的byte
		if ((!(fingerprint & Masks[highbit]))){
			local_bitFeature |= (1ULL << highbit);
			// bitMaskHash[highbit] = true;
            // highbitArray[highbit] = 1;
		}

        i++;
    }
    
    bitFeature = local_bitFeature;
    return n;
}