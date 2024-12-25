
#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <iostream>
// #include <cstring>
#include <iomanip>
#include <random>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <pthread.h>
#include <sys/time.h>
#include <bitset>
#include <zstd.h> 
#include "tools.h"
#include "xxhash.h"
#include "lz4.c"
#include "Gdelta/gear_matrix.h"
#include "Gdelta/gdeltaLshift.cpp"
extern "C" {
	#include "xdelta3/xdelta3.c"
	#include "speedSketch_fastcdc.h"
}

#define ENCODEWHASH 0
#define XDELTA 1
#define UNIFIED 0   // A value of 0 enables the following definition to take effect

#define OHASH_FEATURES 0
#define ODESS_FEATURES 0
#define FINESSE_FEATURES 0   // TODO: do nothing now...
#define NTRANSFORM_FEATURES 1

#define USING_ZSTD 0
#define USING_LZ4 1

std::random_device rd;
std::mt19937 gen(rd()), gen1;
std::uniform_int_distribution<uint64_t> full_uint64_t;
std::uniform_int_distribution<uint32_t> full_uint32_t;
uint64_t Twopower = 1ULL << 32;
std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);

const int INF = INT32_MAX;

using namespace std;

#ifndef BLOCK_SIZE
#define BLOCK_SIZE (8 * 1024)
#endif
#define COMPRESSION

// constexpr int FEATURE_NUM = 1;
// constexpr int SFP_NUM = 1;
constexpr int GRAD_FEATURE_NUM = 1;

constexpr uint32_t LocalMinSize = 8192 / 4;
constexpr uint32_t LocalMaxSize = 8192 * 4;

// char block[LocalMaxSize];
char aux_block[LocalMaxSize];
char compressed[2 * LocalMaxSize];
char delta_compressed[2 * LocalMaxSize];

typedef struct Node {
	uint32_t id;									/* current chunk id */
    SFP_type *sfp;                  /* sketch */
    MYHASH gradFeatures[GRAD_FEATURE_NUM];         /* sketch */
    uint64_t block_offset;
    uint32_t block_size;
} node;


uint32_t CalculateSimThreshold(MYHASH *srcGradFeatures, MYHASH *dstGradFeatures)
{
	uint32_t similarity = 0;
	MYHASH and_result = srcGradFeatures[0] & dstGradFeatures[0];	// total 1-bit shared bits
	// uint32_t and_result = srcGradFeatures[0] ^ dstGradFeatures[0];	// hamming distance

	// while (and_result > 0){
	// 	++similarity;
	// 	and_result &= (and_result - 1);
	// }
	similarity = __builtin_popcountll(and_result);	// 

	// similarity = and_result.count();	// std::bitset count the number of 1-bit

	// return (num_masks - similarity);
	return (similarity);
}


class OhashUnique {
public:
	OhashUnique(char* name, char* tag, int feature_num, int sfp_num);
	~OhashUnique();
	long long getN() { return N; }
    uint32_t getSimilaritySize() {  return similarity[0].size();   }
    void show();
    int fastcdc_chunking(char* filename);
    int fastcdc_chunking_cmp(std::vector<uint64_t> &chunks_offset, char* filename);

	uint64_t fileSize = 0;

private:
    void calOHash(unsigned char* ptr, uint32_t block_size, uint64_t &bitMaskHash, SFP_type* sfp);
    void calOdessFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
    void calFinesseFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
    void calNTransformFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
	int count1InHash(MYHASH hash);

    char fileName[100];
    char* tag_;
	int batch_num_;
    int FEATURE_NUM;
    int SFP_NUM;

	uint64_t N;
    uint32_t W = 32;
    const uint32_t A = 37, MOD = 1000000007;
	const uint64_t Apower = 1;
	const uint64_t Twopower = (1ULL << 32);

    uint64_t dedupTime = 0, dedupCnt = 0, dedupReduce = 0;
    uint64_t deltaTime = 0, deltaCnt = 0, deltaReduce = 0, deltaSize = 0;
    uint64_t selfTime = 0,  selfCnt = 0,  selfReduce = 0, selfSize = 0;
    uint64_t noCompCnt = 0, noCompSize = 0;

    uint64_t total = 0, compCnt = 0;
	uint64_t chunking_time = 0, produce_feature_time = 0;
	uint64_t readBlockTime = 0, readBlockTime4Delta = 0;
	uint64_t searchTime = 0;
	uint64_t temp_count = 0;
	FILE *seq_file = NULL;
	FILE *aux_file = NULL;

    std::unordered_map<XXH64_hash_t, uint32_t> dedup;
	std::unordered_map<SFP_type, node*> *similarity;
    uint32_t *features;

	uint32_t* TRANSPOSE_M;
	uint32_t* TRANSPOSE_A;

    uint32_t count1[num_masks];
	uint32_t maskHightBit[num_masks];
};

OhashUnique::OhashUnique(char* name, char* tag, int feature_num, int sfp_num) {
	tag_ = tag;
    FEATURE_NUM = feature_num;
    SFP_NUM = sfp_num;
    // std::cout << "----------------- file preprocessing -----------------" << std::endl;

	// std::cout << "the mask_mask is :"  << std::setw(16) << std::setfill('0') << std::hex << mask_mask << ", and the mask of using is :" << std::endl;
	// std::cout << "the subOhashMask is : " << std::setfill('0') << subOhashMask << std::endl;
	// // recover decimal 
    // std::cout << std::dec << std::setw(0) << std::setfill(' ');
	// for(int i = 0; i < num_masks; i+=8){
	// 	for(int j = 0; j < std::min(8, (num_masks - i)); ++j){
	// 		std::cout << Masks[i + j] << ","  << std::setw(3) << "\t";
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl;

    seq_file = fopen(name, "r+");
    aux_file = fopen(name, "r+");
    if ((seq_file == NULL) || (aux_file == NULL)) {
        perror("Fail to open file");
        exit(-1);
    }
    fseeko64(aux_file, 0, SEEK_END);   // seek to the end of the file
    fileSize = ftello64(aux_file);
	
	char* file_n = strrchr(name, '/');
    if (file_n != nullptr) {
        file_n++;  
	}

	TRANSPOSE_M = new uint32_t[FEATURE_NUM];
	TRANSPOSE_A = new uint32_t[FEATURE_NUM];
    similarity = new std::unordered_map<SFP_type, node*>[SFP_NUM];
    features = new uint32_t[FEATURE_NUM];

    for (int i = 0; i < FEATURE_NUM; ++i) {
        TRANSPOSE_M[i] = ((full_uint32_t(gen1) >> 1) << 1) + 1;
        TRANSPOSE_A[i] = full_uint32_t(gen1);
    }

    std::memset(maskHightBit, 0, sizeof(maskHightBit));
	std::memset(count1, 0, sizeof(count1));

    std::array<char, 100> timeString;
    GetNowTime(timeString);
	std::string time_str(timeString.data(), 10);
}

void OhashUnique::show() {
    std::cout << "---the results of statistics---" << std::endl;

    std::cout << "Total Block Number : " << std::setw(12) << N << std::endl;
    std::cout << "DedupNum :           " << std::setw(12) << dedupCnt      << " \tDeltaCompressNum : " << deltaCnt << 
            "\tlosslessCompressNum : " << selfCnt << std::endl;
    std::cout << "DedupReduce :        " << std::setw(12) << dedupReduce   << " \tDeltaReduce :      " << deltaReduce << 
            "\tlosslessReduce :      " << selfReduce << std::endl;
    std::cout << "DedupSize :          " << std::setw(12) << dedupReduce   << " \tDeltaSize :        " << deltaSize << 
            "\tlosslessSize :        " << selfSize << std::endl;
    std::cout << "AfterCompress Size :" << std::setw(12) << total 		 << " \tTotal Size :       " << fileSize << std::endl;
    std::cout << "noCompNum : \t" << noCompCnt << "\tnoCompSize : \t" << noCompSize << std::endl;

    double persetage = (double)total / fileSize;
    std::cout << "Compress Ratio :  " << std::setw(12) << persetage <<" (" << 1 / persetage <<")" << std::endl;
    std::cout << "Delta Compression Efficiency (DCE) : " << std::setw(12) << 100.0 * deltaReduce / deltaSize << " %" << std::endl;
    std::cout << "Read Block Time :  " << readBlockTime / 1000000  << "s "
                        << readBlockTime % 1000000 / 1000 << "ms " 
                        << readBlockTime % 1000 << "us" << std::endl;
    std::cout << "Delta Block Time : " << readBlockTime4Delta / 1000000  << "s "
                        << readBlockTime4Delta % 1000000 / 1000 << "ms " 
                        << readBlockTime4Delta % 1000 << "us" << std::endl;
    std::cout << "Search Time :      " << searchTime / 1000000  << "s "
                        << searchTime % 1000000 / 1000 << "ms " 
                        << searchTime % 1000 << "us" << std::endl;
    std::cout << "Dedup Time :       " << dedupTime / 1000000  << "s "
                        << dedupTime % 1000000 / 1000 << "ms " 
                        << dedupTime % 1000 << "us" << std::endl;
    std::cout << "Delta Time :       " << deltaTime / 1000000  << "s "
                        << deltaTime % 1000000 / 1000 << "ms " 
                        << deltaTime % 1000 << "us" << std::endl;
    std::cout << "Self Time :        " << selfTime / 1000000  << "s "
                        << selfTime % 1000000 / 1000 << "ms " 
                        << selfTime % 1000 << "us" << std::endl;
    std::cout << "Chunking time :            " << chunking_time / 1000000 << "s " 
                        << chunking_time % 1000000 / 1000 << "ms " 
                        << chunking_time % 1000 << "us" << std::endl;
    std::cout << "Produce Features time :    " << produce_feature_time / 1000000 << "s " 
                        << produce_feature_time % 1000000 / 1000 << "ms " 
                        << produce_feature_time % 1000 << "us" << std::endl;
}

OhashUnique::~OhashUnique() {
	// std::cout << "the count of maskBit is selected :" << std::endl;
	// for(int i=0; i< num_masks; ++i){
	// 	std::cout << i+1 << ": " <<maskHightBit[i] << ",\t";
	// }
	// std::cout << std::endl;
	// std::cout << "the numbers of count \"1\" is contained in a bitMaskHash with total unique blocks : " << std::endl;
	uint32_t total_blocks = 0;
	for(int i=0; i< num_masks; i+=8){
        for(int j = 0; j < std::min(8, (num_masks - i)); ++j){
			// std::cout << i+j+1 << ": " <<count1[i+j] << "," << std::setw(5) << "\t"  ;
            total_blocks += count1[i];
		}
		// std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << "total_blocks : " << total_blocks << std::endl;
	std::cout << "temp count is : " << temp_count << std::endl;
	if(seq_file) delete seq_file;
	if(aux_file) delete aux_file;
	if(TRANSPOSE_M) delete[] TRANSPOSE_M;
	if(TRANSPOSE_A) delete[] TRANSPOSE_A;
    std::cout << "delete nodes in similarity " << std::endl;
    for(int i = 0; i < SFP_NUM; ++i){
        for (auto &sim : similarity[i]){
            if(sim.second){
                free(sim.second->sfp);
                delete sim.second;
                sim.second = NULL;
            }
        }
    }
    std::cout << "delete similarity " << std::endl;
    if(similarity) delete[] similarity;
    std::cout << "delete features " << std::endl;
    if(features) delete[] features;
}

void OhashUnique::calOHash(unsigned char *ptr, uint32_t block_size, uint64_t &bitMaskHash, SFP_type* sfp) {
	struct timeval start_time, end_time;
	MYHASH local_bitMaskHash = 0;
	OHASH_FPTYPE fingerPrint = 0;
	uint64_t bit_mask = (num_masks - 1);
	uint64_t highbit = 0;
    uint32_t sfpMoveBit = sizeof(uint32_t) * 8 / SFP_NUM;
	// uint32_t i = 0;
	// for (; i < WordSize - 1; ++i){
	// 	fingerPrint = (fingerPrint >> movebitlength) + GEARmx[ptr[i]];
	// }
	for (uint64_t i = 0; i < block_size; ++i){
        // Lshift or Rshift？Lshift->high bit fail; Rshift low bit fail
		fingerPrint = (fingerPrint >> movebitlength) + GEARmx[ptr[i]];	

		// highbit = (fingerPrint >> (num_masks - 6)) & bit_mask;	// Lshift will affecte more bytes
		highbit = (fingerPrint >> 4) & bit_mask;
		if (unlikely(!(fingerPrint & Masks[highbit]))){
			++(maskHightBit[highbit]);
			local_bitMaskHash |= (1ULL << highbit);
			// bitMaskHash[highbit] = true;
		}
	}
    for(int i = 0; i < SFP_NUM; ++i){
		// feature[i] = ((gradFeatures[0]) & (subOhashMask)).to_ullong();
		sfp[i] = (local_bitMaskHash >> (i * sfpMoveBit)) & (subOhashMask);
		// feature[i] = (gradFeatures[0] >> 32) & (subOhashMask);
	}

    bitMaskHash = local_bitMaskHash;
	return;
}

void OhashUnique::calOdessFeature(unsigned char* ptr, int block_size, SFP_type* sfp){
	for (int i = 0; i< FEATURE_NUM; ++i){
		features[i] = 0;
	}
    // Characterizing
    vector<uint32_t> rolling_hashes;

	uint32_t hash = 0;
	for (uint32_t i = 0; i < block_size; ++i) {
		hash = (hash >> 1) + GEARmx[ptr[i]]; 
		// Sampling
		if ((hash & mask) == 0){
	        rolling_hashes.push_back(hash);
		}
    }

    // Selecting
    int hashes_size = rolling_hashes.size();
    for (uint32_t i = 0; i < FEATURE_NUM; ++i) {
		for (uint32_t j = 0; j < hashes_size; j++) {
			uint32_t temp = ((uint64_t) TRANSPOSE_M[i] * (uint64_t) rolling_hashes[j] + (uint64_t) TRANSPOSE_A[i]) % Twopower;
			if (features[i] < temp) {
				features[i] = temp;
			}
		}
    }

    for (int i = 0; i < SFP_NUM; ++i) {
		uint64_t temp[FEATURE_NUM / SFP_NUM];
		for (int j = 0; j < FEATURE_NUM / SFP_NUM; ++j) {
			temp[j] = features[j * SFP_NUM + i];
		}
		sfp[i] = XXH64(temp, sizeof(uint64_t) * FEATURE_NUM / SFP_NUM, 0);
	}
}

void OhashUnique::calFinesseFeature(unsigned char* ptr, int block_size, SFP_type* sfp){
    
}
void OhashUnique::calNTransformFeature(unsigned char* ptr, int block_size, SFP_type* sfp){
    for (int i = 0; i < FEATURE_NUM; ++i){
        features[i] = 0;
    }

	int64_t fp = 0;

	for (int m = 0; m < W; m++) {
		fp *= (uint64_t) A;
		fp += (uint64_t) (unsigned char) ptr[m];
		fp %= (uint64_t) MOD;

	}

	for (int m = 0; m < block_size - W + 1; m++) {

		for (int i = 0; i < FEATURE_NUM; i++) {
			uint32_t temp = ((uint64_t) TRANSPOSE_M[i] * (uint64_t) fp + (uint64_t) TRANSPOSE_A[i]) % Twopower;
			if (features[i] < temp) {
				features[i] = temp;
			}
		}

		fp -= ((uint64_t) ptr[m] * (uint64_t) Apower) % (uint64_t) MOD;
		while (fp < 0) fp += MOD;

		if (m != block_size - W) {
			fp *= (uint64_t) A;
            fp += (uint64_t) (unsigned char) ptr[m + W];
            fp %= (uint64_t) MOD;
		}
	}

    for (int i = 0; i < SFP_NUM; ++i) {
		uint64_t temp[FEATURE_NUM / SFP_NUM];
		for (int j = 0; j < FEATURE_NUM / SFP_NUM; ++j) {
			temp[j] = features[j * SFP_NUM + i];
		}
		sfp[i] = XXH64(temp, sizeof(uint64_t) * FEATURE_NUM / SFP_NUM, 0);
	}
}

int OhashUnique::count1InHash(MYHASH hash){
	int count = 0;
	// while(hash){
	// 	++count;
	// 	hash &= (hash - 1);
	// }
	count = __builtin_popcountll((unsigned long long)hash); // need gcc or clang > C11
	// count = hash.count();
	++(count1[count]);
	return count;
}


int OhashUnique::fastcdc_chunking(char* filename){
	struct timeval start_time, end_time;
    std::printf("using fastcdc chunking, the function is normalized_chunking_64 with feature\n");
    size_t readStatus = 0;
    uint32_t offset = 0, chunkLength = 0, readFlag = 0;
	uint64_t global_offset = 0;
    
    uint32_t chunk_num = 0, end = CacheSize - 1;
    unsigned char *fileCache = (unsigned char *)malloc(CacheSize);
    memset(fileCache, 0, sizeof(CacheSize));

    uint32_t comp_self = INF, dcomp_lsh = INF;
    uint64_t tempself = 0, tempdelta = 0;

    fastCDC_init();
    gettimeofday(&tmStart, NULL);
    gettimeofday(&start_time, NULL);
    readStatus = fread(fileCache, 1, CacheSize, seq_file);
    gettimeofday(&end_time, NULL);
    readBlockTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                     (end_time.tv_usec - start_time.tv_usec));

	for (;;) {
        chunk_num += 1;
        comp_self = INF, dcomp_lsh = INF;
        MYHASH bitFeatures;
        SFP_type *sfp = (SFP_type *)malloc(sizeof(SFP_type) * SFP_NUM);

        gettimeofday(&start_time, NULL);
#if UNIFIED
		// ----------------Unified Blocking and Fingerprint Calculation
        chunkLength = normalized_chunking_64_with_features(fileCache + offset, CacheSize - offset + 1, bitFeatures); 

        // calOdessFeature(fileCache + offset, chunkLength, sfp);
        feature[0] = (bitFeatures & subOhashMask);
        gettimeofday(&end_time, NULL);
        chunking_time += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                    (end_time.tv_usec - start_time.tv_usec));  
        produce_feature_time += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                    (end_time.tv_usec - start_time.tv_usec));  

        count1InHash(bitFeatures);
        // ----------------Unified Blocking and Fingerprint Calculation
#else
        // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
        chunkLength = normalized_chunking_64(fileCache + offset, CacheSize - offset + 1); 
        gettimeofday(&end_time, NULL);
        chunking_time += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                    (end_time.tv_usec - start_time.tv_usec));  
        // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
#endif
        // --------------------post-deduplication delta compression------------------------
        gettimeofday(&start_time, NULL);
        XXH64_hash_t h = XXH64(fileCache + offset, chunkLength, 0);
		if (!(dedup.count(h))) { // non-deduplication
			dedup[h]=chunk_num - 1;
            gettimeofday(&end_time, NULL);
            dedupTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                        (end_time.tv_usec - start_time.tv_usec));
#if UNIFIED
            ;
#else
            // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
            gettimeofday(&start_time, NULL);
    #if OHASH_FEATURES
            calOHash(fileCache + offset, chunkLength, bitFeatures, sfp);
    #endif
    #if ODESS_FEATURES
            calOdessFeature(fileCache + offset, chunkLength, sfp);
    #endif
    #if NTRANSFORM_FEATURES
            calNTransformFeature(fileCache + offset, chunkLength, sfp);
    #endif
    #if FINESSE_FEATURES
            calFinesseFeature(fileCache + offset, chunkLength, sfp);
    #endif
            gettimeofday(&end_time, NULL);
            produce_feature_time += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                    (end_time.tv_usec - start_time.tv_usec));  

            count1InHash(bitFeatures);
            // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
#endif
            // TODO losslesscompression
            gettimeofday(&start_time, NULL);

#if USING_ZSTD
            comp_self = ZSTD_compress(compressed, 2 * LocalMaxSize, (uint8_t *)(fileCache + offset), chunkLength, 3);
#elif USING_LZ4
            comp_self = LZ4_compress_default((char *)(fileCache + offset), compressed, chunkLength, 2 * LocalMaxSize);
#endif
            comp_self = (comp_self < 0) ? INF : comp_self;
            gettimeofday(&end_time, NULL);
            tempself = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                     (end_time.tv_usec - start_time.tv_usec));
            selfTime += tempself;

            bool find_similarity = false;
            for(int i = 0; i < SFP_NUM; ++i){
                if(similarity[i].count(sfp[i])){
                    node* temp_node = similarity[i][sfp[i]];
                    // TODO delta compression 
                    // requirement: baseblock offset, newblock
                    // gettimeofday(&start_time, NULL);
                    fseeko64(aux_file, temp_node->block_offset, SEEK_SET);
                    readStatus = fread(aux_block, 1, temp_node->block_size, aux_file);
                    // gettimeofday(&end_time, NULL);
                    // readBlockTime4Delta = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                    //      (end_time.tv_usec - start_time.tv_usec));

                    gettimeofday(&start_time, NULL);
                    uint32_t localDeltaSize = 0;
                    ++temp_count;

    #if ENCODEWHASH
                    dcomp_lsh = gencodeWHash((uint8_t *)(fileCache + offset), chunkLength, 
                                        (uint8_t *)aux_block, temp_node->block_size, 
                                        (uint8_t **)&delta_compressed, &localDeltaSize, bitFeatures, temp_node->gradFeatures[0]);
    #elif XDELTA

                    dcomp_lsh = xdelta3_compress((char *)(fileCache + offset), chunkLength, 
                                    (char *)aux_block, temp_node->block_size, 
                                    delta_compressed, 1);
    #else
                    dcomp_lsh = gencode((uint8_t *)(fileCache + offset), chunkLength, 
                                        (uint8_t *)aux_block, temp_node->block_size, 
                                        (uint8_t **)&delta_compressed, &localDeltaSize);
    #endif
                    gettimeofday(&end_time, NULL);
                    tempdelta = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                                        (end_time.tv_usec - start_time.tv_usec));
                    deltaTime += tempdelta;

                    // replace the block in similarity
                    temp_node->id = chunk_num - 1;
                    memcpy(temp_node->sfp, sfp, sizeof(sfp));
                    temp_node->gradFeatures[0] = bitFeatures;
                    temp_node->block_offset = global_offset;
                    temp_node->block_size = chunkLength;     
                    for (int s = 0; s < SFP_NUM; ++s){
                        similarity[s][sfp[s]] = temp_node;
                    }
                    free(sfp);
                    find_similarity = true;
                    break;
                }
            }
            if(!find_similarity){
                node* sfpNode = new node;
                sfpNode->id = chunk_num - 1;
                sfpNode->sfp = sfp;
                // memcpy(sfpNode->sfp, sfp, sizeof(sfp));
                sfpNode->gradFeatures[0] = bitFeatures;
                sfpNode->block_offset = global_offset;
                sfpNode->block_size = chunkLength;
                for (int s = 0; s < SFP_NUM; ++s){
                    similarity[s][sfp[s]] = sfpNode;
                }
            }
//             if(similarity.count(feature[0])){  // resemblance
//                 // TODO delta compression 
//                 // requirement: baseblock offset, newblock
//                 // gettimeofday(&start_time, NULL);
//                 fseeko64(aux_file, similarity[feature[0]]->block_offset, SEEK_SET);
// 		        readStatus = fread(aux_block, 1, similarity[feature[0]]->block_size, aux_file);
//                 // gettimeofday(&end_time, NULL);
//                 // readBlockTime4Delta = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
//                 //      (end_time.tv_usec - start_time.tv_usec));

//                 gettimeofday(&start_time, NULL);
//                 uint32_t localDeltaSize = 0;
//                 ++temp_count;
// #if ENCODEWHASH
//                 dcomp_lsh = gencodeWHash((uint8_t *)(fileCache + offset), chunkLength, 
//                                     (uint8_t *)aux_block, similarity[feature[0]]->block_size, 
//                                     (uint8_t **)&delta_compressed, &localDeltaSize, bitFeatures, similarity[feature[0]]->gradFeatures[0]);
// #else
//                 dcomp_lsh = gencode((uint8_t *)(fileCache + offset), chunkLength, 
//                                     (uint8_t *)aux_block, similarity[feature[0]]->block_size, 
//                                     (uint8_t **)&delta_compressed, &localDeltaSize);
// 			    // dcomp_lsh = xdelta3_compress((char *)(fileCache + offset), chunkLength, 
//                 //                     (char *)aux_block, similarity[feature[0]]->block_size, 
//                 //                     delta_compressed, 1);
// #endif
//                 gettimeofday(&end_time, NULL);
//                 tempdelta = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
//                                     (end_time.tv_usec - start_time.tv_usec));
//                 deltaTime += tempdelta;

//                 // replace the block in similarity
//                 similarity[feature[0]]->id = chunk_num - 1;
//                 similarity[feature[0]]->sfp[0] = feature[0];
//                 similarity[feature[0]]->gradFeatures[0] = bitFeatures;
//                 similarity[feature[0]]->block_offset = global_offset;
//                 similarity[feature[0]]->block_size = chunkLength;
// 		    }else{                      // non-resemblance
//                 node* sfpNode = new node{chunk_num - 1, feature[0], bitFeatures, global_offset, chunkLength};
//                 // sfpNode->id = chunk_num - 1;
//                 // sfpNode->sfp[0] = feature;
//                 // sfpNode->gradFeatures[0] = bitFeatures;
//                 // sfpNode->block_offset = global_offset;
//                 // sfpNode->block_size = chunkLength;
//                 similarity[feature[0]] = sfpNode;
//             }
            // TODO compare the size of lossless and delta
            if (std::min(comp_self, chunkLength) > dcomp_lsh) { // delta compress
                total += dcomp_lsh;
                deltaReduce += (chunkLength - dcomp_lsh);
                deltaSize += chunkLength;
                deltaCnt += 1;

            }else if(comp_self < chunkLength){    // self compress
                total += comp_self;
                selfReduce += (chunkLength - comp_self);
                selfSize += chunkLength;
                selfCnt += 1;

            }else{ // no compress
                total += chunkLength;
                noCompSize += (chunkLength);
                noCompCnt += 1;

            }

		}else{                  // deduplication
            gettimeofday(&end_time, NULL);
            dedupTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                        (end_time.tv_usec - start_time.tv_usec));
            dedupCnt += 1;
			dedupReduce += chunkLength;
            // dedup[h].push_back(i);
            ;
        }

        // --------------------post-deduplication delta compression------------------------
        
        if(chunkLength <= 1024 * 4) ++smalChkCnt;
        offset += chunkLength;
        global_offset += chunkLength;

        if (CacheSize - offset < MaxSize) {
            gettimeofday(&start_time, NULL);
            memcpy(fileCache, fileCache + offset + 1, CacheSize - offset);
            readStatus = fread(fileCache + CacheSize - offset, 1, offset, seq_file);
            end = CacheSize - 1 - offset + readStatus;
            gettimeofday(&end_time, NULL);
            readBlockTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                     (end_time.tv_usec - start_time.tv_usec));

            if (readStatus < offset + 1) {
                // all the files are read
                readFlag = 1;
            }

            offset = 0;
        }

        if (offset >= end && readFlag == 1) 
            break;
        if (++compCnt % 5000 == 0) {
			// fprintf(stderr, "%d/%lld\t\t%d/%d\r", i, N, lsh.getRebuildMonitor1(), lsh.getRebuildMonitor2());
			std::cerr <<"processingSize : " << global_offset / (1024 * 1024) << " / " << fileSize / (1024 * 1024) << " MiB, ChunkNum : " << chunk_num << "\r" << std::flush;
		}
    }
    std::cout << std::endl;

    gettimeofday(&tmEnd, NULL);
    totalTm += ((tmEnd.tv_sec - tmStart.tv_sec) * 1000000 + 
                        (tmEnd.tv_usec - tmStart.tv_usec));
    N = chunk_num;

    std::printf("totalTime is %f s\n", totalTm / 1000000);
    std::printf("chunknum is %d\n", chunk_num);
    std::printf("small chunknum is %d\n", smalChkCnt);

    // clear the items
    std::free(fileCache);
    // fileCache = NULL;
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
		cerr << "usage: "<< argv[0] << " [input_file] [feature_num] [sfp_num] [tag] \n";
        cerr << "example: "<< argv[0] << " /data/data1 64 1000 tag\n";
		exit(0);
	}
    
    std::array<char, 100> timeString;
    GetNowTime(timeString);
	cout << timeString.data() << endl;
    std::cout << std::endl;
    char* filename = argv[1];
    int feature_num = atoi(argv[2]);
    int sfp_num = atoi(argv[3]);
    char* tag = argv[4];

    std::cout << "----------------- execute graph method-----------------" << std::endl;
    std::cout << "the executable program:  \t" << argv[0] << std::endl;
    std::cout << "process file or dir is:  \t" << filename << std::endl;
    std::cout << "FEATURE_NUM: \t" << feature_num << std::endl;
    std::cout << "SFP_NUM: \t" << sfp_num << std::endl;
    std::cout << "BLOCK_SIZE: \t" << "using fastcdc..." << std::endl;
    std::cout << "UNIFIED: \t" << UNIFIED << std::endl;
    std::cout << "-OHASH_FEATURES: \t" << OHASH_FEATURES << std::endl;
    std::cout << "-ODESS_FEATURES: \t" << ODESS_FEATURES << std::endl;
    std::cout << "-FINESSE_FEATURES: \t" << FINESSE_FEATURES << std::endl;
    std::cout << "-NTRANSFORM_FEATURES: \t" << NTRANSFORM_FEATURES << std::endl;
    std::cout << "ENCODEWHASH: \t" << ENCODEWHASH << std::endl;
    std::cout << "XDELTA: \t" << XDELTA << std::endl;
    std::cout << std::endl;
    std::cout << "USING_ZSTD: \t" << USING_ZSTD << std::endl;
    std::cout << "USING_LZ4: \t" << USING_LZ4 << std::endl;
    std::cout << "Runing tag: \t" << tag << std::endl;
    

	OhashUnique* ohash = new OhashUnique(filename, tag, feature_num, sfp_num);
	
    std::cout << "-------start chunking, calculate bithash and sfp, data reduction------" << std::endl;
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);

    ohash->fastcdc_chunking(filename);

    std::cout << "-------complete chunking, calculate bithash and sfp, data reduction------" << std::endl;
    std::cout << "graph->similarity size : " << ohash->getSimilaritySize() << std::endl;
	gettimeofday(&end_time, NULL);
	long long calTime = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);


    std::cout << "-------show all results------" << std::endl;                              
                           
    std::cout << "Total Time :  " << calTime / 1000000  << "s "
                        << calTime % 1000000 / 1000 << "ms " 
                        << calTime % 1000 << "us" << std::endl;
	ohash->show();
	
    delete ohash;

    GetNowTime(timeString);
	cout << timeString.data() << endl;
}
