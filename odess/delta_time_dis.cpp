#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <bitset>
#include <map>
#include <cmath>
#include <algorithm>

#include "edelta/src/edelta.cc"
#include "Gdelta/gdelta.cpp"
#include "odess.h"
// #include <zdlib.h>
// #include <zdconf.h>
#include "zdelta/zdlib.h"
#include "zdelta/zdconf.h"
#include <zstd.h> 
#include <lz4.h>
extern "C" {
	#include "xdelta3/xdelta3.c"
}
#define INF 987654321
using namespace std;


typedef struct COMPRESS_PARAM{
    uint64_t total = 0;
    
    uint64_t deltaTime = 0;
	uint64_t deltaSize = 0;
    uint64_t deltaReduce = 0;
    uint64_t deltaCnt = 0;
    
    uint64_t selfTime = 0;
	uint64_t selfSize = 0;
    uint64_t selfReduce = 0;
    uint64_t selfCnt = 0;

}COMPRESS_PARAM;

void statisticResult(int comp_self, int delta_size, uint64_t self_time, uint64_t delta_time, COMPRESS_PARAM& compress_param){
    if (min(comp_self, BLOCK_SIZE) > delta_size) { // delta compress
        compress_param.total += delta_size;
        compress_param.deltaTime += delta_time;
        compress_param.deltaSize += (BLOCK_SIZE);
        compress_param.deltaReduce += (BLOCK_SIZE - delta_size);
        compress_param.deltaCnt += 1;
    }
    else {
        if (comp_self < BLOCK_SIZE) { // self compress
            compress_param.total += comp_self;
            compress_param.selfTime += self_time;
            compress_param.selfSize += (BLOCK_SIZE);
            compress_param.selfReduce += (BLOCK_SIZE - comp_self);
            compress_param.selfCnt += 1;
        }
        else { // no compress
            compress_param.total += BLOCK_SIZE;
        }
    }
}

void printResults(long long N, COMPRESS_PARAM& compress_param){
    printf("Final size: %llu (%.2lf%%)\n", compress_param.total, (double)compress_param.total * 100 / N / BLOCK_SIZE);

	cout << "Delta Time : " << compress_param.deltaTime / 1000000  << "s "
                            << compress_param.deltaTime % 1000000 / 1000 << "ms " 
                            << compress_param.deltaTime % 1000 << "us";
	cout << "\tAverage Delta Time : " << (double)compress_param.deltaTime / compress_param.deltaCnt << std::endl;
    cout << "Self Time : " << compress_param.selfTime / 1000000  << "s "
                            << compress_param.selfTime % 1000000 / 1000 << "ms " 
                            << compress_param.selfTime % 1000 << "us" << endl;	

	cout << "Delta Size : " << compress_param.deltaSize << "\tSelf Size : " << compress_param.selfSize << endl;
	cout << "Delta Reduce : " << compress_param.deltaReduce << "\tSelf Reduce: " << compress_param.selfReduce << endl;
	cout << "Delta Count : " << compress_param.deltaCnt << "\t Self Count : " << compress_param.selfCnt << endl;
}

int main(int argc, char* argv[]) {
	if (argc != 5) {
		cerr << "usage: "<< argv[0] << " [input_file] [window_size] [SF_NUM] [FEATURE_NUM]\n";
		cerr << "this method reads block from hardware one by one\n";
		exit(0);
	}
    cout << "---------------------delta time distribution-------------------------------" << endl;
	cout  << "this method reads block from hardware one by one\n";

	std::array<char, 100> timeString;
	GetNowTime(timeString);
	cout << timeString.data() << endl;
	struct timeval total_time_s, total_time_e;
	struct timeval request_time_s, request_time_e;

	int W = atoi(argv[2]);
	int SF_NUM = atoi(argv[3]);
	int FEATURE_NUM = atoi(argv[4]);
	char* filename = argv[1];
	cout << "input file: \t" << filename << endl;
	cout << "windows size: \t" << W << endl;
	cout << "SF Number: \t" << SF_NUM << endl;
	cout << "Feature Number: \t" << FEATURE_NUM << endl;

	char compressed[2 * BLOCK_SIZE];
	char delta_compressed[2 * BLOCK_SIZE];
	uint8_t *delta_ptr = (uint8_t *)delta_compressed;
    uint8_t **ptr_to_delta_ptr = &delta_ptr;
	
	File_OP file_t(argv[1]);
	File_OP file_s(argv[1]);
	long long N = file_t.getBlockNum();
    cout << "the number of the blocks of the file is: " << N << endl;
	
	char* file_n = strrchr(filename, '/');
    if (file_n != nullptr) {
        file_n++;   
	}

    cout << "using CDF..." << endl;
    char time_filename[timeString.size() + strlen(file_n) + strlen("_all_delta_time_dis.csv")];
    strncpy(time_filename, timeString.data(), 11);
    time_filename[11] = '\0';
    strcat(time_filename, file_n);
    strcat(time_filename, "_all_delta_time_dis.csv");
    // cout << "time_file : " << time_filename << "size: "  << strlen(time_filename) << endl;

    std::ofstream time_file(time_filename, std::ios::out);
    cout << "the file of time is saved to: " << time_filename <<endl;
    time_file << "method," << "block_reduction_ratio," << "time(us)" << std::endl;

	map<XXH64_hash_t, int> dedup;
	Odess lsh(W, SF_NUM, FEATURE_NUM); // parameter

	unsigned long long total = 0;
	uint64_t readBlockTime = 0, dedupTime = 0;
	uint64_t deltaTime = 0,  selfTime = 0;
	uint64_t deltaReduce = 0, selfReduce = 0;
	uint64_t deltaCnt = 0, selfCnt = 0, dedupCnt = 0, insertCnt = 0;
	uint64_t dedupSize = 0;
	uint64_t requestTime = 0;

	gettimeofday(&total_time_s, NULL);

	// char block[BLOCK_SIZE +1];
	// block[BLOCK_SIZE] = '\0';
	std::array<char, BLOCK_SIZE> block;
	std::array<char, BLOCK_SIZE> block_ref;


    COMPRESS_PARAM xdelta_param;
    COMPRESS_PARAM zdelta_param;
    COMPRESS_PARAM edelta_param;
    COMPRESS_PARAM gdelta_param;

	for (int i = 0; i < N / 3; ++i) {
		struct timeval start_time, end_time;

		// RECIPE r;

		gettimeofday(&start_time, NULL);
		memcpy(block.data(), file_t.nextBlock(), BLOCK_SIZE);
		gettimeofday(&end_time, NULL);
		readBlockTime += (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
		
        long long tempdedup = 0;
		gettimeofday(&start_time, NULL);
		XXH64_hash_t h = XXH64(block.data(), BLOCK_SIZE, 0);

		if (dedup.count(h)) { // deduplication
			dedupCnt += 1;
			dedupSize += BLOCK_SIZE; //
			// set_ref(r, dedup[h]);
			// set_flag(r, 0b10);
			// f.recipe_insert(r);
			continue;
		}

		dedup[h] = i;
		gettimeofday(&end_time, NULL);
        tempdedup = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
		dedupTime += tempdedup;
        time_file << "dedup," << 0.0 << "," << tempdedup << std::endl;
		
		long long temp_self = 0;
		gettimeofday(&start_time, NULL);
		// int comp_self = LZ4_compress_default(block.data(), compressed, BLOCK_SIZE, 2 * BLOCK_SIZE);
		int comp_self = ZSTD_compress(compressed, 2 * BLOCK_SIZE, (uint8_t *)(block.data()), BLOCK_SIZE, 1);
		gettimeofday(&end_time, NULL);
		temp_self = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);

		int dcomp_xdelta = INF, dcomp_lsh_ref;
        int dcomp_zdelta = INF;
        int dcomp_edelta = INF;
		int dcomp_gdelta = INF;
		uint32_t gdeltaDeltaSize = INF;

        time_file << "lossless," << std::fixed << std::setprecision(4) << (double)comp_self / BLOCK_SIZE << "," << temp_self << std::endl;
        
		gettimeofday(&request_time_s, NULL);
		dcomp_lsh_ref = lsh.request((unsigned char*)block.data(), BLOCK_SIZE);

		long long temp_xdelta = 0, temp_zdelta = 0, temp_edelta = 0, temp_gdelta = 0;
		if (dcomp_lsh_ref != -1) {
            memcpy(block_ref.data(), file_s.getBlockFromSeek(dcomp_lsh_ref), BLOCK_SIZE);

			gettimeofday(&start_time, NULL);
			dcomp_xdelta = xdelta3_compress(block.data(), BLOCK_SIZE, block_ref.data(), BLOCK_SIZE, delta_compressed, 1);
			gettimeofday(&end_time, NULL);
			temp_xdelta = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
            time_file << "xdelta," << std::fixed << std::setprecision(4) << (double)dcomp_xdelta / BLOCK_SIZE << "," << temp_xdelta << std::endl;
			
            gettimeofday(&start_time, NULL);
            zd_compress((Bytef*)block_ref.data(), BLOCK_SIZE, (Bytef*)block.data(), BLOCK_SIZE, (Bytef*)delta_compressed, (uLongf*)&dcomp_zdelta);
			gettimeofday(&end_time, NULL);
            temp_zdelta = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
            time_file << "zdelta," << std::fixed << std::setprecision(4) << (double)dcomp_zdelta / BLOCK_SIZE << "," << temp_zdelta << std::endl;

            gettimeofday(&start_time, NULL);
            EDeltaEncode((uint8_t*)block.data(), BLOCK_SIZE, (uint8_t*)block_ref.data(), BLOCK_SIZE, (uint8_t*)delta_compressed, (uint32_t*)&dcomp_edelta);
            gettimeofday(&end_time, NULL);
            temp_edelta = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
            time_file << "edelta," << std::fixed << std::setprecision(4) << (double)dcomp_edelta / BLOCK_SIZE << "," << temp_edelta << std::endl;

            gettimeofday(&start_time, NULL);
            dcomp_gdelta = gencode((uint8_t*)block.data(), BLOCK_SIZE, (uint8_t*)block_ref.data(), BLOCK_SIZE, 
                                        ptr_to_delta_ptr, &gdeltaDeltaSize);
            gettimeofday(&end_time, NULL);
            temp_gdelta = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
            time_file << "gdelta," << std::fixed << std::setprecision(4) << (double)dcomp_gdelta / BLOCK_SIZE << "," << temp_gdelta << std::endl;

		}

        statisticResult(comp_self, dcomp_xdelta, temp_self, temp_xdelta, xdelta_param);
        statisticResult(comp_self, dcomp_zdelta, temp_self, temp_zdelta, zdelta_param);
        statisticResult(comp_self, dcomp_edelta, temp_self, temp_edelta, edelta_param);
		statisticResult(comp_self, dcomp_gdelta, temp_self, temp_gdelta, gdelta_param);

		// if (min(comp_self, BLOCK_SIZE) > dcomp_xdelta) { // delta compress
		// 	total += dcomp_xdelta;
		// 	deltaTime += temp_xdelta;
		// 	deltaReduce += (BLOCK_SIZE - dcomp_xdelta);
		// 	deltaCnt += 1;
		// }
		// else {
		// 	if (comp_self < BLOCK_SIZE) { // self compress
		// 		total += comp_self;
		// 		selfTime += temp_self;
		// 		selfReduce += (BLOCK_SIZE - comp_self);
		// 		selfCnt += 1;
		// 	}
		// 	else { // no compress
		// 		total += BLOCK_SIZE;
		// 	}
		// }


		lsh.insert(i);
		insertCnt++;
		gettimeofday(&request_time_e, NULL);
		requestTime += (request_time_e.tv_sec - request_time_s.tv_sec) * 1000000 + 
						(request_time_e.tv_usec - request_time_s.tv_usec);
								// f.recipe_insert(r);

		if (i % 1000 == 0) {
			fprintf(stderr, "%d/%lld\r", i, N);
		}
	}
	// f.recipe_write();

	cout << "N-Transform" << endl;
	gettimeofday(&total_time_e, NULL);
	long long total_time = (total_time_e.tv_sec - total_time_s.tv_sec) * 1000000 + (total_time_e.tv_usec - total_time_s.tv_usec);
	printf("Trace: %s\n", argv[1]);
	printf("method: Odess, W = %d, SF = %d, feature = %d\n", W, SF_NUM, FEATURE_NUM);
    cout << "------------------------ results -------------------------" << endl;

	cout << "Total time(us): " << total_time <<"us" << endl;
	cout << "Total time: " << (total_time_e.tv_sec - total_time_s.tv_sec) << "s "
							<< (total_time_e.tv_usec - total_time_s.tv_usec) / 1000 << "ms "
							<< (total_time_e.tv_usec - total_time_s.tv_usec) % 1000 << "us " << endl;
	cout << "The Number of Blocks : " << N << "\tRequest Time : " << requestTime / 1000000  << "s "
                                                                << requestTime % 1000000 / 1000 << "ms " 
                                                                << requestTime % 1000 << "us" << endl;
    cout << "Read Blocks Time : " << readBlockTime / 1000000  << "s "
                                << readBlockTime % 1000000 / 1000 << "ms " 
                                << readBlockTime % 1000 << "us" << endl;
    cout << "Dedup Time : " << dedupTime / 1000000  << "s "
                            << dedupTime % 1000000 / 1000 << "ms " 
                            << dedupTime % 1000 << "us" << endl;
	cout << "Dedup Size : " << dedupSize << "\tDedup Reduce : " << dedupSize << endl;
	cout << "Dedup Count : " << dedupCnt << "\t Insert SF Count : " << insertCnt << endl;

    cout << "---------------------the result of xdelta----------------------" << endl;
    printResults(N, xdelta_param);

    cout << "---------------------the result of zdelta----------------------" << endl;
    printResults(N, zdelta_param);
    
    cout << "---------------------the result of edelta----------------------" << endl;
    printResults(N, edelta_param);

	cout << "---------------------the result of gdelta----------------------" << endl;
    printResults(N, gdelta_param);

    cout << "---------------------the overhead of memory----------------------" << endl;

	lsh.calMemory();
	std::cout << "Average Characterzing Time : " << (double)lsh.getFeatureTime()/insertCnt << std::endl;
	std::cout << "Average Indexing Time : " << (double)lsh.getIndexingTime()/insertCnt << std::endl;
    cout << "-----------------------------end----------------------------" << endl;

}
