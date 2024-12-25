
#include <stdint.h>
#include <stdlib.h>
#include <thread>
#include <condition_variable>
#include <queue>
#include <vector>
#include <atomic>
#include <iostream>
#include <fstream>
// #include <cstring>
#include <iomanip>
#include <random>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <pthread.h>
#include <sched.h>
#include <sys/time.h>
#include <bitset>
#include <zstd.h> 
#include <lz4.h>
// #include "lz4.c"
#include "tools.h"
#include "xxhash.h"
#include "Gdelta/gear_matrix.h"
#include "Gdelta/gdelta.cpp"
extern "C" {
	#include "xdelta3/xdelta3.c"
	#include "speedSketch_fastcdc.h"
}
// #define FEATURE_WITH_DATA

#define GDELTA 0
#define GDELTAWHASH 1
#define XDELTA 0
#define UNIFIED 0   // A value of 0 enables the following definition to take effect

#define OHASH_FEATURES 1
#define ODESS_FEATURES 0
#define FINESSE_FEATURES 0   // TODO: do nothing now...
#define NTRANSFORM_FEATURES 0

#define USING_ZSTD 1
#define USING_LZ4 0

std::random_device rd;
std::mt19937 gen(rd()), gen1(922);
std::uniform_int_distribution<uint64_t> full_uint64_t;
std::uniform_int_distribution<uint32_t> full_uint32_t;
// const uint64_t Twopower = 1ULL << 32;
std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);

const int INF = INT32_MAX;

using namespace std;
using std::atomic_bool;

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
// char aux_block[LocalMaxSize];

typedef struct Node {
	// uint32_t id;								   /* current chunk id */
    uint32_t ref = 0;
    // SFP_type *sfp = nullptr;                          /* similar fp */
    MYHASH gradFeatures;         /* sketch */
    uint64_t block_offset;
    uint32_t block_size;
    // uint8_t *block_data = nullptr;      /* more effcient with data*/
} Node;


uint32_t CalculateSimThreshold(MYHASH srcGradFeatures, MYHASH dstGradFeatures)
{
	uint32_t sim = 0;
	// MYHASH and_result = srcGradFeatures[0] & dstGradFeatures[0];	// the number of shared 1-bit
	MYHASH and_result = ~(srcGradFeatures ^ dstGradFeatures);	// anti-hamming distance

	// while (and_result > 0){
	// 	++sim;
	// 	and_result &= (and_result - 1);
	// }
	sim = __builtin_popcountll(and_result);

	// similarity = and_result.count();	// std::bitset to count the number of 1-bits

	// return (num_masks - sim);
	return (sim);
}


class OhashUnique {
public:
	OhashUnique(char* name, char* tag, int feature_num, int sfp_num);
	~OhashUnique();
	long long getN() { return N; }
    uint32_t getSimilaritySize() {  return similarity[0].size();   }
    void show();
    int start_processing();

	uint64_t fileSize = 0;

private:
    void calOHash(unsigned char* ptr, uint32_t block_size, uint64_t &bitMaskHash, SFP_type* sfp);
    void calOdessFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
    void calFinesseFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
    void calNTransformFeature(unsigned char* ptr, int block_size, SFP_type* sfp);
	int count1InHash(MYHASH hash);

    void doChunkingT();
    void doFeatureCalculateionT();
    void doCompressionT();

    char fileName[100];
    char* tag_;
	int batch_num_;
    int FEATURE_NUM;
    int SFP_NUM;
    uint32_t sfpMoveBit = 0;

	uint64_t N;
    uint32_t W = 32;
    const uint32_t A = 37, MOD = 1000000007;
	uint64_t Apower = 1;
	const uint64_t Twopower = (1ULL << 32);

    std::queue<Node*> blocks; // chunks queue
    std::queue<std::pair<Node*, Node*>> similarBlocks; // sketch queue
    std::mutex mutex_blocks;
    std::mutex mutex_similarBlocks;
    std::condition_variable cond_blocks;
    std::condition_variable cond_similarBlocks;
    std::atomic<bool> finished_features;
    std::atomic<bool> finished_chunking;
    
    uint64_t dedupTime = 0, dedupCnt = 0, dedupReduce = 0;
    uint64_t deltaTime = 0, deltaCnt = 0, deltaReduce = 0, deltaSize = 0;
    uint64_t selfTime = 0,  selfCnt = 0,  selfReduce = 0, selfSize = 0;
    uint64_t noCompCnt = 0, noCompSize = 0;

    uint64_t total = 0, compCnt = 0;
	uint64_t chunkingTime = 0, produceFeatureTime = 0, featureIndexingTime = 0, findSimilarChunkTime = 0, insertFeatureTime = 0;
	uint64_t readBlockTime = 0, readBlockTime4Delta = 0;
	uint64_t searchTime = 0;
	uint64_t compressionCount = 0;
    uint64_t findSimilar = 0;
    uint64_t uniqueCount = 0;

    std::ifstream *chunk_file;
    std::ifstream *feature_file;
    std::ifstream *comp_file;

	// FILE *seq_file = NULL;
	// FILE *aux_file = NULL;

    std::unordered_map<XXH64_hash_t, uint32_t> dedup;
	std::unordered_map<SFP_type, Node*> *similarity;
    uint32_t *features;
	uint64_t bit_mask = (num_masks - 1);

	uint32_t* TRANSPOSE_M;
	uint32_t* TRANSPOSE_A;

    uint32_t count1[num_masks] = {};
	uint32_t maskHightBit[num_masks] = {};
    std::pair<uint64_t, uint64_t> similarDCE[num_masks+1] = {};
};


void OhashUnique::doChunkingT(){
    // 获取当前线程 ID
    pthread_t thread = pthread_self();
    
    // 创建一个 CPU 集合
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset); // clear cpu set
    CPU_SET(0, &cpuset); // performed with CPU 0

    if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) != 0) {
        std::cerr << "Error, Thread doChunkingT() setting thread affinity" << std::endl;
        return;
    }

    int cpu = sched_getcpu();
    std::cout << "Thread doChunkingT() is running on CPU: " << cpu << std::endl;


    int testCount = 0, testCount2 = 0;
    struct timeval start_time, end_time;
    std::printf("using fastcdc chunking, the function is normalized_chunking_64 with feature\n");
    size_t readStatus = 0;
    uint32_t offset = 0, chunkLength = 0, readFlag = 0;
	uint64_t global_offset = 0;
    
    uint32_t chunk_num = 0, end = CacheSize - 1;
    unsigned char *fileCache = (unsigned char *)malloc(CacheSize);
    memset(fileCache, 0, sizeof(CacheSize));

    fastCDC_init();
    gettimeofday(&tmStart, NULL);
    gettimeofday(&start_time, NULL);
    chunk_file->read((char *)fileCache,CacheSize);

    gettimeofday(&end_time, NULL);
    readBlockTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                     (end_time.tv_usec - start_time.tv_usec));

	for (;;) {
        // std::cout << "chunking : " << ++testCount << std::endl;
        ++chunk_num;
        gettimeofday(&start_time, NULL);
        // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
        chunkLength = normalized_chunking_64(fileCache + offset, CacheSize - offset + 1); 
        gettimeofday(&end_time, NULL);
        chunkingTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                    (end_time.tv_usec - start_time.tv_usec)); 
        // ----------------original fastcdc calcultate the bithash and feature out of chunking processing
        gettimeofday(&start_time, NULL);
        XXH64_hash_t h = XXH64(fileCache + offset, chunkLength, 0);
		if (!(dedup.count(h))) { // non-deduplication
            ++uniqueCount;
            Node *node = new Node;
            // node->id = chunk_num - 1;
            node->ref = SFP_NUM; // new chunk must update sfp and map
            node->block_offset = global_offset;
            node->block_size = chunkLength;
            // node->sfp = (SFP_type *)malloc(sizeof(SFP_type) * SFP_NUM);
            // node->block_data = (uint8_t *)malloc(chunkLength);
            // memcpy(node->block_data, fileCache + offset, chunkLength);

			dedup[h]= chunk_num - 1;
            gettimeofday(&end_time, NULL);
            dedupTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                        (end_time.tv_usec - start_time.tv_usec));

            {
                std::lock_guard<std::mutex> lock(mutex_blocks);
                blocks.push(node);
            // std::cout << "push blocks in chunking : " << ++testCount2 << std::endl;

            }
            cond_blocks.notify_one();

		}else{                  // deduplication
            gettimeofday(&end_time, NULL);
            dedupTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                        (end_time.tv_usec - start_time.tv_usec));
            dedupCnt += 1;
			dedupReduce += chunkLength;
            // dedup[h].push_back(i);
        }

        if(chunkLength <= 1024 * 4) ++smalChkCnt;
        offset += chunkLength;
        global_offset += chunkLength;

        if (CacheSize - offset < MaxSize) {
            gettimeofday(&start_time, NULL);
            memcpy(fileCache, fileCache + offset + 1, CacheSize - offset);
            chunk_file->read((char *)fileCache + CacheSize - offset, offset);
            readStatus = chunk_file->gcount();
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
        // std::cout << "chunking : " << testCount << std::endl;


        if (offset >= end && readFlag == 1) {
            finished_chunking.store(true);
            cond_blocks.notify_all();
            break;
        }
        if (++compCnt % 100000 == 0) {
			// fprintf(stderr, "%d/%lld\t\t%d/%d\r", i, N, lsh.getRebuildMonitor1(), lsh.getRebuildMonitor2());
			std::cerr <<"processingSize : " << global_offset / (1024 * 1024) << " / " << fileSize / (1024 * 1024) << " MiB, ChunkNum : " << chunk_num << "\r" << std::flush;
		}
    }
    std::cout << std::endl;

    gettimeofday(&tmEnd, NULL);
    totalTm += ((tmEnd.tv_sec - tmStart.tv_sec) * 1000000 + 
                        (tmEnd.tv_usec - tmStart.tv_usec));
    N = chunk_num;
    std::cout << "Read Block Time :  " << readBlockTime / 1000000  << "s "
                    << readBlockTime % 1000000 / 1000 << "ms " 
                    << readBlockTime % 1000 << "us" << std::endl;

    std::printf("totalTime of chunking is %f s\n", totalTm / 1000000);
    std::printf("chunknum is %d\n", chunk_num);
    std::printf("small chunknum is %d\n", smalChkCnt);

    // clear the items
    std::free(fileCache);
    // fileCache = NULL;
}

void OhashUnique::doFeatureCalculateionT(){
    // 获取当前线程 ID
    pthread_t thread = pthread_self();
    
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(1, &cpuset); // binding to CPU 1

    if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) != 0) {
        std::cerr << "Error, Thread doFeatureCalculateionT() setting thread affinity" << std::endl;
        return;
    }

    int cpu = sched_getcpu();
    std::cout << "Thread doFeatureCalculateionT() is running on CPU: " << cpu << std::endl;

    int testCount = 0, testCount2 = 0;
    struct timeval start_time, end_time;
    // size_t readStatus = 0;

    SFP_type *sfp = new SFP_type[SFP_NUM];
    uint8_t *data_pointer = (uint8_t *)malloc(LocalMaxSize);

    while (true) {
        Node* block;
        {
            std::unique_lock<std::mutex> lock(mutex_blocks);
            cond_blocks.wait(lock, [this] { return !blocks.empty() || finished_chunking.load(); });
            if (blocks.empty() && finished_chunking.load()) {
                finished_features.store(true);
                cond_similarBlocks.notify_all();
                break;
            }
            block = std::move(blocks.front());
            blocks.pop();
        }

        feature_file->seekg(block->block_offset);
        feature_file->read((char *)data_pointer, block->block_size);

        gettimeofday(&start_time, NULL);
#if OHASH_FEATURES
        calOHash(data_pointer, block->block_size, block->gradFeatures, sfp);
#elif ODESS_FEATURES
        calOdessFeature(data_pointer, block->block_size, sfp);
#elif NTRANSFORM_FEATURES
        calNTransformFeature(data_pointer, block->block_size, sfp);
#elif FINESSE_FEATURES
        calFinesseFeature(data_pointer, block->block_size, sfp);
#endif
        gettimeofday(&end_time, NULL);
        produceFeatureTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec)); 
    	++(count1[count1InHash(block->gradFeatures)]);

        Node *current_node = new Node(*block);  // a new node for compression
        // current_node->sfp = (SFP_type *)malloc(sizeof(SFP_type) * SFP_NUM);
        // memcpy(current_node->sfp, sfp, sizeof(SFP_type) * SFP_NUM); // alloc memory in chunking
        // free(data_pointer);
        // block->block_data = nullptr;
// #ifdef FEATURE_WITH_DATA
//         block->block_data = data_pointer;
//         current_node->block_data = (uint8_t *)malloc(current_node->block_size);
//         memcpy(current_node->block_data, data_pointer, current_node->block_size);
// #endif
        Node *aux_node = nullptr;  // a new aux_node for compression

        bool find_similarity = false;
        gettimeofday(&start_time, NULL);

        for(int i = 0; i < SFP_NUM; ++i){
            
            if(similarity[i].count(sfp[i])){
                // std::cout << "push block in featureCalculateion : " << ++testCount2 << std::endl;
                if(!aux_node)    ++findSimilar;
                gettimeofday(&end_time, NULL);
                findSimilarChunkTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec));
                gettimeofday(&start_time, NULL);
                Node *temp_node = similarity[i][sfp[i]];

                if(!aux_node)    aux_node = new Node(*temp_node);
                gettimeofday(&end_time, NULL);
                readBlockTime4Delta += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec));
                gettimeofday(&start_time, NULL);

                // replace the block in similarity
                similarity[i][sfp[i]] = block; // update block
                if ((temp_node->ref - 1)) { // ref > 1, it will be used by other chunk
                    temp_node->ref -= 1;    // reduce ref
                    continue;
                }

                delete temp_node;
                temp_node = nullptr;

                gettimeofday(&end_time, NULL);
                insertFeatureTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec));
                gettimeofday(&start_time, NULL);
                find_similarity = true;
                // std::cout << "push block in featureCalculateion : " << testCount2 << std::endl;

                // break;
            }else{
                similarity[i][sfp[i]] = block; // update block
            }
            
        }
        gettimeofday(&end_time, NULL);
        featureIndexingTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec)); 
        gettimeofday(&start_time, NULL);

        // std::cout << "processing block in featureCalculateion : " << testCount << std::endl;

        // if(!find_similarity){
        //     for (int s = 0; s < SFP_NUM; ++s){
        //         similarity[s][sfp[s]] = block;
        //     }
        // }
        gettimeofday(&end_time, NULL);
        insertFeatureTime += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
                (end_time.tv_usec - start_time.tv_usec)); 

        {
            std::lock_guard<std::mutex> lock(mutex_similarBlocks);
            similarBlocks.push(std::make_pair(current_node, aux_node));
            // std::cout << "push block in featureCalculateion : " << testCount2 << std::endl;

        }
        cond_similarBlocks.notify_one(); // notify compression thread
    }
    delete[] sfp;
    free(data_pointer);
}

void OhashUnique::doCompressionT(){
    pthread_t thread = pthread_self();
    
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset); // 
    CPU_SET(2, &cpuset); // binding to CPU 2

    if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) != 0) {
        std::cerr << "Error, Thread doCompression() setting thread affinity" << std::endl;
        return;
    }

    int cpu = sched_getcpu();
    std::cout << "Thread doCompression() is running on CPU: " << cpu << std::endl;

    int testCount = 0, testCount2 = 0;
    size_t readStatus;
    int doCompressionCount = 0;
    char srcData[LocalMaxSize], auxData[LocalMaxSize];
    char compressed[2 * LocalMaxSize], delta_compressed[2 * LocalMaxSize];
    uint8_t *delta_ptr = (uint8_t *)delta_compressed; 
    uint8_t **ptr_to_delta_ptr = &delta_ptr;
    struct timeval start_time, end_time;
    std::pair<Node*, Node*> compBlocks;
    Node *src_node = nullptr, *aux_node = nullptr;
    
    while(true){
        {
            std::unique_lock<std::mutex> lock(mutex_similarBlocks);
            cond_similarBlocks.wait(lock, [this] { return !similarBlocks.empty() || finished_features.load(); });
            if (similarBlocks.empty() && finished_features.load()) {
                return;
            }
            compBlocks = std::move(similarBlocks.front());
            similarBlocks.pop();
        }
        ++compressionCount;
        src_node = compBlocks.first;
        aux_node = compBlocks.second;
        
        gettimeofday(&start_time, NULL);
        if (aux_node){
            comp_file->seekg(aux_node->block_offset);
            comp_file->read(auxData, aux_node->block_size);
        }
        comp_file->seekg(src_node->block_offset);
        comp_file->read(srcData, src_node->block_size);

        gettimeofday(&end_time, NULL);
        readBlockTime4Delta += ((end_time.tv_sec - start_time.tv_sec) * 1000000 +
        (end_time.tv_usec - start_time.tv_usec));
        gettimeofday(&start_time, NULL);

        // std::cout << "processing block in Compression : " << ++testCount << std::endl;

        int comp_self = INF, comp_delta = INF;
        uint64_t  tempself = 0, tempdelta = 0;

#if USING_ZSTD
        comp_self = ZSTD_compress(compressed, 2 * LocalMaxSize, (uint8_t *)(srcData), src_node->block_size, 1);
#elif USING_LZ4
        comp_self = LZ4_compress_default((char *)(srcData), compressed, src_node->block_size, 2 * LocalMaxSize);
#endif
        comp_self = (comp_self < 0) ? INF : comp_self;
		gettimeofday(&end_time, NULL);

        gettimeofday(&end_time, NULL);
        tempself = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                    (end_time.tv_usec - start_time.tv_usec));
        selfTime += tempself;

        gettimeofday(&start_time, NULL);
        uint32_t localDeltaSize = 0;

        if (aux_node != nullptr){

#if GDELTA
        comp_delta = gencode((uint8_t *)(srcData), src_node->block_size, 
                                (uint8_t *)(auxData), aux_node->block_size, 
                                ptr_to_delta_ptr, &localDeltaSize);
    #if USING_ZSTD
        comp_delta = ZSTD_compress(compressed, 2 * LocalMaxSize, (uint8_t *)(delta_compressed), localDeltaSize, 1);
    #elif USING_LZ4
        comp_delta = LZ4_compress_default((char *)(delta_compressed), compressed, localDeltaSize, 2 * LocalMaxSize);
    #endif
#elif GDELTAWHASH
            comp_delta = gencodeWHash((uint8_t *)(srcData),    src_node->block_size, 
                    (uint8_t *)(auxData), aux_node->block_size, 
                    ptr_to_delta_ptr, &localDeltaSize, src_node->gradFeatures, aux_node->gradFeatures);
        #if USING_ZSTD
            comp_delta = ZSTD_compress(compressed, 2 * LocalMaxSize, (uint8_t *)(delta_compressed), localDeltaSize, 1);
        #elif USING_LZ4
            comp_delta = LZ4_compress_default((char *)(delta_compressed), compressed, localDeltaSize, 2 * LocalMaxSize);
        #endif
            // int similarity = CalculateSimThreshold(src_node->gradFeatures, aux_node->gradFeatures);
            // similarDCE[similarity].first += comp_delta;
            // similarDCE[similarity].second += src_node->block_size;
#elif XDELTA

            comp_delta = xdelta3_compress((char *)(srcData), src_node->block_size, 
                            (char *)(auxData), aux_node->block_size, 
                            delta_compressed, 1);
#else
        std::cout << " unknown compression algorithm" << std::endl;

#endif
            gettimeofday(&end_time, NULL);
            tempdelta = ((end_time.tv_sec - start_time.tv_sec) * 1000000 + 
                                (end_time.tv_usec - start_time.tv_usec));
            deltaTime += tempdelta;
        }

        if (std::min(comp_self, int(src_node->block_size)) > comp_delta) { // delta compress
            total += comp_delta;
            deltaReduce += (src_node->block_size - comp_delta);
            deltaSize += src_node->block_size;
            deltaCnt += 1;
#if OHASH_FEATURES
            int similarity = CalculateSimThreshold(src_node->gradFeatures, aux_node->gradFeatures);
            similarDCE[similarity].first += (src_node->block_size - comp_delta);
            similarDCE[similarity].second += src_node->block_size;
#endif
        }else if(comp_self < src_node->block_size){    // self compress
            total += comp_self;
            selfReduce += (src_node->block_size - comp_self);
            selfSize += src_node->block_size;
            selfCnt += 1;

        }else{ // no compress
            total += src_node->block_size;
            noCompSize += src_node->block_size;
            noCompCnt += 1;
        }

        if (aux_node != nullptr){
            // free(aux_node->sfp);
            // free(aux_node->block_data);
            delete aux_node;
            aux_node = nullptr;
        }

        // free(src_node->sfp);
        // free(src_node->block_data);
        delete src_node;
        src_node = nullptr;
        // std::cout << "processing block in Compression : " << testCount << std::endl;
        if(((++doCompressionCount) % 50000 == 0) && finished_chunking.load()){
            std::cerr <<"compressionProgress : " << doCompressionCount  << " / " << uniqueCount << "\r" << std::flush;
        }

    }
    std::cout << std::endl;
}

OhashUnique::OhashUnique(char* name, char* tag, int feature_num, int sfp_num) {
	tag_ = tag;
    FEATURE_NUM = feature_num;
    SFP_NUM = sfp_num;
    full_uint32_t = std::uniform_int_distribution<uint32_t>(std::numeric_limits<uint32_t>::min(), std::numeric_limits<uint32_t>::max());

    chunk_file = new std::ifstream(name, std::ios::binary);
    feature_file = new std::ifstream(name, std::ios::binary);
    comp_file = new std::ifstream(name, std::ios::binary);
        // move to the end of file
    feature_file->seekg(0, std::ios::end);
    fileSize = feature_file->tellg();

    if ((!chunk_file->is_open()) || (!feature_file->is_open()) || (!comp_file->is_open())) {
        std::cerr << "Failed to open file " << name << " for reading. Error code: " << errno << std::endl;
        throw std::runtime_error("File open error");
    }
	
	char* file_n = strrchr(name, '/');
    if (file_n != nullptr) {
        file_n++; // Remove directory separators
	}

	TRANSPOSE_M = new uint32_t[FEATURE_NUM];
	TRANSPOSE_A = new uint32_t[FEATURE_NUM];
    similarity = new std::unordered_map<SFP_type, Node*>[SFP_NUM];
    features = new uint32_t[FEATURE_NUM];

    for (int i = 0; i < FEATURE_NUM; ++i) {
        TRANSPOSE_M[i] = ((full_uint32_t(gen1) >> 1) << 1) + 1;
        TRANSPOSE_A[i] = full_uint32_t(gen1);
    }
    for (int i = 0; i < W - 1; ++i) {
        Apower *= A;
        Apower %= MOD;
    }
    if(SFP_NUM > 1){
        sfpMoveBit = sizeof(uint32_t) * 8 / (SFP_NUM - 1);
    }
    std::cout << "sfpMoveBit : " << sfpMoveBit << std::endl;

    std::array<char, 100> timeString;
    GetNowTime(timeString);
	std::string time_str(timeString.data(), 10);
    finished_features.store(false);
    finished_chunking.store(false);
}

void OhashUnique::show() {
    std::cout << "---the results of statistics---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Total Block Number : " << std::setw(12) << N << std::endl;
    std::cout << "DedupNum :           " << std::setw(12) << dedupCnt      << " \tDeltaCompressNum : " << deltaCnt << 
            "\tlosslessCompressNum : " << selfCnt << std::endl;
    std::cout << "DedupReduce :        " << std::setw(12) << dedupReduce   << " \tDeltaReduce :      " << deltaReduce << 
            "\tlosslessReduce :      " << selfReduce << std::endl;
    std::cout << "DedupSize :          " << std::setw(12) << dedupReduce   << " \tDeltaSize :        " << deltaSize << 
            "\tlosslessSize :        " << selfSize << std::endl;
    std::cout << "AfterCompress Size :" << std::setw(12) << total 		 << " \tTotal Size :       " << fileSize << std::endl;
    std::cout << "noCompNum : \t" << noCompCnt << "\tnoCompSize : \t" << noCompSize << std::endl;
    std::cout << "Unique Chunk Nums : " << uniqueCount << std::endl;

    double persetage = (double)total / fileSize;
    std::cout << "Compress Ratio :  " << std::setw(12) << persetage <<" (" << 1 / persetage <<")" << std::endl;
    std::cout << "Delta Compression Efficiency (DCE) : " << std::setw(12) << 100.0 * deltaReduce / deltaSize << " %" << std::endl;


    // std::cout << "Search Time :      " << searchTime / 1000000  << "s "
    //                     << searchTime % 1000000 / 1000 << "ms " 
    //                     << searchTime % 1000 << "us" << std::endl;
    std::cout << "Dedup Time :       " << dedupTime / 1000000  << "s "
                        << dedupTime % 1000000 / 1000 << "ms " 
                        << dedupTime % 1000 << "us" << std::endl;
    std::cout << "Delta Time :       " << deltaTime / 1000000  << "s "
                        << deltaTime % 1000000 / 1000 << "ms " 
                        << deltaTime % 1000 << "us" << std::endl;
    std::cout << "Self Time :        " << selfTime / 1000000  << "s "
                        << selfTime % 1000000 / 1000 << "ms " 
                        << selfTime % 1000 << "us" << std::endl;
    std::cout << "Chunking Time :            " << chunkingTime / 1000000  << "s "
                        << chunkingTime % 1000000 / 1000 << "ms " 
                        << chunkingTime % 1000 << "us" << std::endl;
    std::cout << "Producing Features time :    " << produceFeatureTime / 1000000 << "s " 
                        << produceFeatureTime % 1000000 / 1000 << "ms " 
                        << produceFeatureTime % 1000 << "us" << std::endl;
    std::cout << "Feature Indexing time :    " << featureIndexingTime / 1000000 << "s " 
                        << featureIndexingTime % 1000000 / 1000 << "ms " 
                        << featureIndexingTime % 1000 << "us" << std::endl;
    std::cout << " -Finding Similar Chunks time :    " << findSimilarChunkTime / 1000000 << "s " 
                        << findSimilarChunkTime % 1000000 / 1000 << "ms " 
                        << findSimilarChunkTime % 1000 << "us" << std::endl;
    std::cout << " -Inserting Chunks time :    " << insertFeatureTime / 1000000 << "s " 
                        << insertFeatureTime % 1000000 / 1000 << "ms " 
                        << insertFeatureTime % 1000 << "us" << std::endl;
    std::cout << "Reading Delta Block Time : " << readBlockTime4Delta / 1000000  << "s "
                        << readBlockTime4Delta % 1000000 / 1000 << "ms " 
                        << readBlockTime4Delta % 1000 << "us" << std::endl;
}

OhashUnique::~OhashUnique() {
	std::cout << "compressionCount is : " << compressionCount << std::endl;
    std::cout << "findSimilar : " << findSimilar << std::endl;
#if OHASH_FEATURES
    for (int i = 0; i< num_masks+1; i+=8){
        for(int j = 0; j < std::min(8, (num_masks + 1 - i)); ++j){
            std::cout << 100.0 * similarDCE[i+j].first / similarDCE[i+j].second << ",\t" ;
        }
        std::cout << std::endl;
    }
#endif
    std::cout << std::endl;
    if(chunk_file->is_open()) chunk_file->close();
    if(feature_file->is_open()) feature_file->close();
    if(comp_file->is_open()) comp_file->close();

	if(TRANSPOSE_M) delete[] TRANSPOSE_M;
	if(TRANSPOSE_A) delete[] TRANSPOSE_A;
    std::cout << "delete nodes in similarity " << std::endl;
    uint64_t totalMemory = 0;
    totalMemory += (dedup.size() * (sizeof(XXH64_hash_t) + sizeof(uint32_t)));
    int testCount = 0;
    for(int i = 0; i < SFP_NUM; ++i){
        std::cout << "the size of similarity " << i+1 << ": " << similarity[i].size() << std::endl;

#if OHASH_FEATURES
            totalMemory += ((sizeof(SFP_type) / 2 + sizeof(uint32_t)) * SFP_NUM) * similarity[i].size();
#else
            totalMemory += ((sizeof(Node) + sizeof(uint32_t)) * SFP_NUM) * similarity[i].size();
#endif
        // std::cout << "delete the node in " << i << "th similar table" << std::endl;
        testCount = 0;
        // std::cout << "the size of similairy is " << similarity[i].size() << std::endl;
        for (auto &sim : similarity[i]){
            // std::cout << "delete the node in " << i << "th similar table" << "\ttestCount is: " << ++testCount << std::endl;
            if ((sim.second->ref - 1)) {
                sim.second->ref -= 1;
                sim.second = nullptr;
                continue;
            }
            delete sim.second;
            sim.second = nullptr;
        }
    }
    std::cout << "total usage of memory is : " << totalMemory / (1024 * 1024) << "MB "
                                                << (totalMemory % (1024 * 1024)) / (1024) << "KB " 
                                                << (totalMemory % 1024) << "B" << std::endl;
    delete[] similarity;
    std::cout << "delete similarity " << std::endl;
    std::cout << "delete features " << std::endl;
    if(features) delete[] features;
}

void OhashUnique::calOHash(unsigned char *ptr, uint32_t block_size, uint64_t &bitMaskHash, SFP_type* sfp) {
	// struct timeval start_time, end_time;
	MYHASH local_bitMaskHash = 0;
	OHASH_FPTYPE fingerPrint = 0;
	uint64_t highbit = 0;
	// uint32_t i = 0;
	// for (; i < WordSize - 1; ++i){
	// 	fingerPrint = (fingerPrint >> movebitlength) + GEARmx[ptr[i]];
	// }
	for (uint64_t i = 0; i < block_size; ++i){
		fingerPrint = (fingerPrint >> movebitlength) + GEARmx[ptr[i]];

		// highbit = (fingerPrint >> (num_masks - 6)) & bit_mask;
		highbit = (fingerPrint) & bit_mask;
		if (unlikely(!(fingerPrint & Masks[highbit]))){
			++(maskHightBit[highbit]);
			local_bitMaskHash |= (1ULL << highbit);
			// bitMaskHash[highbit] = true;
		}
	}
    // std::cout << "bitMaskHash : " << std::bitset<64>(local_bitMaskHash) << "\t";
    for(int i = 0; i < SFP_NUM; ++i){
		// feature[i] = ((gradFeatures[0]) & (subOhashMask)).to_ullong();
		sfp[i] = (local_bitMaskHash >> (i * sfpMoveBit)) & (subOhashMask);
		// feature[i] = (gradFeatures[0] >> 32) & (subOhashMask);
        // std::cout<< "sfp[" << i << "] = " << std::bitset<32>(sfp[i]) << "\t";
	}
    // std::cout<< std::endl;

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
    for (int i = 0; i < SFP_NUM; ++i){
        sfp[i] = 0;
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
	count = __builtin_popcountll((unsigned long long)hash);
	// count = hash.count();
	return count;
}

int OhashUnique::start_processing(){

    std::thread block_thread(&OhashUnique::doChunkingT,this);
    std::thread feature_thread(&OhashUnique::doFeatureCalculateionT,this);
    std::thread compress_thread(&OhashUnique::doCompressionT,this);

    // 等待线程完成
    block_thread.join();
    feature_thread.join();
    compress_thread.join();

    // doChunkingT();
    // doFeatureCalculateionT();
    // doCompressionT();

    return 0;
}


int main(int argc, char* argv[]) {
    if (argc != 5) {
		cerr << "usage: "<< argv[0] << " [input_file] [feature_num] [sfp_num] [tag] \n";
        cerr << "example: "<< argv[0] << " /data/data1 12 3 tag\n";
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
    std::cout << std::endl;
    std::cout << "UNIFIED: \t" << UNIFIED << std::endl;
    std::cout << "-OHASH_FEATURES: \t" << OHASH_FEATURES << std::endl;
    std::cout << "-ODESS_FEATURES: \t" << ODESS_FEATURES << std::endl;
    std::cout << "-FINESSE_FEATURES: \t" << FINESSE_FEATURES << std::endl;
    std::cout << "-NTRANSFORM_FEATURES: \t" << NTRANSFORM_FEATURES << std::endl;
    std::cout << std::endl;
    std::cout << "GDELTAWHASH: \t" << GDELTAWHASH << std::endl;
    std::cout << "XDELTA: \t" << XDELTA << std::endl;
    std::cout << "GDELTA: \t" << GDELTA << std::endl;
    std::cout << std::endl;
    std::cout << "USING_ZSTD: \t" << USING_ZSTD << std::endl;
    std::cout << "USING_LZ4: \t" << USING_LZ4 << std::endl;

    std::cout << "Runing tag: \t" << tag << std::endl;

	OhashUnique* ohash = new OhashUnique(filename, tag, feature_num, sfp_num);
	
    std::cout << "-------start chunking, calculate bithash and sfp, data reduction------" << std::endl;
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);

    ohash->start_processing();

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
