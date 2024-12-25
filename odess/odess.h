#include <tuple>
#include <limits>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <array>
#include "tools.h"
#include "xxhash.h"

#include "Gdelta/gear_matrix.h"
#define G GEARmx
// Odess
using namespace std;

class Odess {
	private:
	std::mt19937 gen1, gen2;
	std::uniform_int_distribution<uint32_t> full_uint32_t;

	int W;
	int SF_NUM, FEATURE_NUM;
	uint32_t mask = 0x2ac6;	// 0x2a6(32): 0010 1010 0110;  0x2ac6(128):0010 1010 0110 1100; 0x2a6c5(512):0010 1010 0110 1100 0101

	uint32_t* TRANSPOSE_M;
	uint32_t* TRANSPOSE_A;

	const uint32_t A = 37, MOD = 1000000007;
	uint64_t Apower = 1;
	uint64_t Twopower = 1ULL << 32;

	uint32_t* feature;
	uint64_t* superfeature;

	std::map<uint64_t, std::vector<int>>* sfTable;

	void calFeatures(unsigned char* ptr, int blcok_size);
    uint32_t GearHash(const unsigned char* data, uint32_t len, uint32_t hash);
	void calSuperFeatures();

	public:
    uint64_t produce_time = 0, index_time = 0, insert_time = 0, calSF_time = 0;
	uint64_t memoryOverhead = 0;
	Odess(int _W, int _SF_NUM, int _FEATURE_NUM) {
		gen1 = std::mt19937(922);
		gen2 = std::mt19937(314159);
		full_uint32_t = std::uniform_int_distribution<uint32_t>(std::numeric_limits<uint32_t>::min(), std::numeric_limits<uint32_t>::max());

		W = _W;
		SF_NUM = _SF_NUM;
		FEATURE_NUM = _FEATURE_NUM;

		feature = new uint32_t[FEATURE_NUM];
		superfeature = new uint64_t[SF_NUM];

		sfTable = new std::map<uint64_t, std::vector<int>>[SF_NUM];
		// 
// 		TRANSPOSE_M = new uint32_t[30]{
//     53306, 15204, 54395, 65535, 12320, 30849, 21870, 38763, 64878, 41396,
//     27992, 22161, 33929, 39401, 61717, 42322, 51550, 11966, 51985, 359, 76,
//     46662, 49675, 30359, 21866, 32644, 25889, 64697, 13204, 47875
// };
// 		TRANSPOSE_A = new uint32_t[30]{
// 	24595, 12584, 44078, 2835, 57400, 11786, 15474, 49521, 44583, 43399, 46690,
//     12756, 63638, 11352, 15407, 22270, 58848, 42951, 18874, 58239, 23522,
//     48324, 29813, 55683, 42855, 7206, 61887, 505, 58277, 42830
// 		};
		//
		TRANSPOSE_M = new uint32_t[FEATURE_NUM];
		TRANSPOSE_A = new uint32_t[FEATURE_NUM];

		for (int i = 0; i < FEATURE_NUM; ++i) {
			TRANSPOSE_M[i] = ((full_uint32_t(gen1) >> 1) << 1) + 1;
			TRANSPOSE_A[i] = full_uint32_t(gen1);
		}
		for (int i = 0; i < W - 1; ++i) {
			Apower *= A;
			Apower %= MOD;
		}
	}
	~Odess() {
		delete[] TRANSPOSE_M;
		delete[] TRANSPOSE_A;
		delete[] feature;
		delete[] superfeature;
		delete[] sfTable;
	}
	int request(unsigned char* ptr, int blcok_size);
	void insert(int label);
	void calMemory();
	uint64_t getFeatureTime(){ return produce_time + calSF_time;};
	uint64_t getIndexingTime(){ return index_time + insert_time;};
	uint64_t getMemoryOverhead()	{return memoryOverhead;}
};

// Gear rolling hash
uint32_t Odess::GearHash(const unsigned char* data, uint32_t len, uint32_t hash) {
    // uint32_t hash_l = 0;
    for (uint32_t i = 0; i < len; ++i) {
        // hash_l = hash_l * 33 + data[i];  // readPaper
        hash = (hash << 1) + G[data[i]];    // Ddelta 
        // hash_l += G[data[i]];          // readPaper
        // hash_l <<= 1;
    }
    return hash;
}

// Odess feature generation
void Odess::calFeatures(unsigned char* ptr, int block_size){
    for (int i = 0; i < FEATURE_NUM; ++i) {
		// feature[i] = UINT32_MAX;
		feature[i] = 0;
	}
    // Characterizing
    vector<uint32_t> rolling_hashes;

    // for (uint32_t i = 0; i <= block_size - W; ++i) {
	// 	uint32_t hash = 0;
    // 	uint32_t rolling_hash = GearHash(ptr + i, W, hash);
	// 	// Sampling
	// 	if ((rolling_hash & mask) == 0){
	// 		rolling_hashes.push_back(rolling_hash);
	// 	}
    // }

	uint32_t hash = 0;
	for (uint32_t i = 0; i < block_size; ++i) {
		hash = (hash >> 1) + G[ptr[i]]; 
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
			if (feature[i] < temp) {
				feature[i] = temp;
			}
		}
    }
}

void Odess::calSuperFeatures(){
	for (int i = 0; i < SF_NUM; ++i) superfeature[i] = 0;

	for (int i = 0; i < SF_NUM; ++i) {
		uint64_t temp[FEATURE_NUM / SF_NUM];
		for (int j = 0; j < FEATURE_NUM / SF_NUM; ++j) {
			temp[j] = feature[j * SF_NUM + i];
		}
		superfeature[i] = XXH64(temp, sizeof(uint64_t) * FEATURE_NUM / SF_NUM, 0);
	}
}

int Odess::request(unsigned char* ptr, int block_size) {
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);

	calFeatures(ptr, block_size);

	gettimeofday(&end_time, NULL);
	produce_time += (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);	
	gettimeofday(&start_time, NULL);

	calSuperFeatures();

	gettimeofday(&end_time, NULL);
	calSF_time += (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
	gettimeofday(&start_time, NULL);

	// uint32_t r = full_uint32_t(gen2) % SF_NUM;
	std::vector<std::vector<int>> results;
	for (int i = 0; i < SF_NUM; ++i) {
		// int index = (r + i) % SF_NUM;
		int index = i;
		if (sfTable[index].count(superfeature[index])) {
			gettimeofday(&end_time, NULL);
			index_time += (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
			return sfTable[index][superfeature[index]].back();
		}
	}
	gettimeofday(&end_time, NULL);
	index_time += (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
	return -1;
}

// insert "prev calculated" sf: label
void Odess::insert(int label) {
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);
	for (int i = 0; i < SF_NUM; ++i) {
		sfTable[i][superfeature[i]].push_back(label);
	}
	gettimeofday(&end_time, NULL);
	long long temp = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
	insert_time += temp;
}

void Odess::calMemory() {
	for (int i = 0; i < SF_NUM; i++) {
		for (auto it = sfTable[i].begin(); it != sfTable[i].end(); it++) {
			memoryOverhead += sizeof(uint64_t) + sizeof(int) * it->second.size();
		}
	}
	int memory_MB = memoryOverhead / (1024 * 1024);
	int memory_KB = (memoryOverhead % (1024 * 1024)) / 1024;
	int memory_B = memoryOverhead % 1024;
	cout << "Memory Overhead : " << memory_MB << "MB " << memory_KB << "KB "<< memory_B << "B" << endl;
	// cout << "Calculate FingerPrint Time : " << produce_time  << " us"<< endl;
	// cout << "Search Reference Block Time : " << index_time << " us" << endl;
	// cout << "Insert Time : " << insert_time << " us" << endl;

	cout << "Calculate FingerPrint Time : " << produce_time / 1000000  << "s " 
											<< produce_time % 1000000 / 1000 << "ms " 
											<< produce_time % 1000 << "us" << endl;
	cout << "Calculate SuperFeature Time : " << calSF_time / 1000000  << "s " 
											<< calSF_time % 1000000 / 1000 << "ms " 
											<< calSF_time % 1000 << "us" << endl;
	cout << "Search Reference Block Time : " << index_time / 1000000  << "s " 
											<< index_time % 1000000 / 1000 << "ms " 
											<< index_time % 1000 << "us" << endl;
	cout << "Insert Time : " 				<< insert_time / 1000000  << "s " 
											<< insert_time % 1000000 / 1000 << "ms " 
											<< insert_time % 1000 << "us" << endl;
}

