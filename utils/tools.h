#pragma once
#ifndef __DEFINE_H__
#define __DEFINE_H__
#include <stdint.h>
#include <string.h>

#ifndef BLOCK_SIZE
	#define BLOCK_SIZE (8 * 1024)
#endif

#if (defined(__GNUC__) && (__GNUC__ >= 3)) || (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 800)) || defined(__clang__)
#  define expect(expr,value)    (__builtin_expect ((expr),(value)) )
#else
#  define expect(expr,value)    (expr)
#endif

#ifndef likely
#define likely(expr)     expect((expr) != 0, 1)
#endif
#ifndef unlikely
#define unlikely(expr)   expect((expr) != 0, 0)
#endif

#include <array>
#include <sys/time.h>
#include <ctime>
void GetNowTime(std::array<char, 100>& timeString){
	std::time_t now = std::time(nullptr);
    std::tm* timeinfo = std::localtime(&now);
    std::strftime(timeString.data(), timeString.size(), "%Y-%m-%d_%H:%M:%S", timeinfo);
}

#include <iostream>
class File_OP {
    public:
    File_OP(char* file){
        read_file_size(file);
        block = new char[BLOCK_SIZE];

        blocks = nullptr;
        temp_blocks = nullptr;
        batchNums = -1;
        seekForBlock = 0;
        seekFile(0);
    }
    File_OP(char* file, int b_nums){
        batchNums = b_nums;
        read_file_size(file);
        
        blocks = new char *[batchNums];
        temp_blocks = new char[batchNums * BLOCK_SIZE];
        for(int i = 0; i < batchNums +1; i++){
            blocks[i] = temp_blocks + i * BLOCK_SIZE;
        }
        block = nullptr;
        seekForBlock = 0;
        seekFile(0);
    }

    void read_file_size(char* file);
    long long getBlockNum() { return N; }
    long long get_file_size() { return fileSize;   };
    char* getFirstBlock()   {   return getBlockFromSeek(0); };
    char** getFirstBlocks() {   return getBlocksFromSeek(0); };
    char* getBlockFromSeek(int i);
    char** getBlocksFromSeek(int i);
    char* nextBlock();
    char** nextBlocks();
    void seekFile(long long s);
    long long getSeekNum() {    return seekForBlock;   }

    ~File_OP(){
        if (block){
            delete[] block;
            block = nullptr;
        }        
        if(temp_blocks){
            delete[] temp_blocks;
            temp_blocks =nullptr;
        }
        if (blocks){
            delete[] blocks;
            blocks = nullptr;
        }

        if(fp){
            fclose(fp);
            fp = nullptr;
        }
    }


    private:
    char fileName[100];
    long long N;
    long long fileSize;
    int batchNums;
    char* block;
    char** blocks;
    char* temp_blocks;
    FILE* fp = nullptr;
    long long seekForBlock;
};

void File_OP::seekFile(long long s){
    seekForBlock = s;
    uint64_t offset = (uint64_t) seekForBlock * BLOCK_SIZE;
    if(fp)  fseeko64(fp, offset, SEEK_SET);
    else {
        fp = fopen(fileName, "rb"); 
        fseeko64(fp, offset, SEEK_SET);
        }  
}

void File_OP::read_file_size(char* file){
	// printf("Trace: %s\n", file);
	sprintf(fileName, "%s", file);
    if(!fp) fp = fopen(fileName, "rb");
	fseeko64(fp, 0L, SEEK_END);  
	fileSize = ftell(fp);
	// cout << "file size : " << size << endl;
	N = fileSize / (BLOCK_SIZE);
	// cout << "Total Block : " << N << endl;
}

char* File_OP::getBlockFromSeek(int i) { 
	// FILE * fp = fopen(fileName, "rb");    
    seekFile(i);
	return nextBlock();
}

char** File_OP::getBlocksFromSeek(int i) { 
	// FILE * fp = fopen(fileName, "rb");
    seekFile(i);
	return nextBlocks();
}


char* File_OP::nextBlock(){
    if (seekForBlock >= N){
        std::cout << "read file finish, the last Block is the " << N << "th Block." << std::endl;
        return block;
    }

	if(fp) int now = fread(block, 1, BLOCK_SIZE, fp);
    else{
        fp = fopen(fileName, "rb");
	    fseeko64(fp, 0L, (uint64_t) seekForBlock * BLOCK_SIZE);
        int now = fread(block, 1, BLOCK_SIZE, fp);
    }
    seekForBlock += 1;

    return block;
}

char** File_OP::nextBlocks(){
    long long read_nums = batchNums;
    if (seekForBlock >= N){
        std::cout << "read file finish, EOF" << std::endl;
        return blocks;
    }
    if (seekForBlock >= N-batchNums){
        read_nums = N - seekForBlock;
        std::cout << "The number of the last batch is: " << read_nums << std::endl;
    }
    // uint64_t offset = (uint64_t) seekForBlocks * BLOCK_SIZE;
    // fseeko64(fp, offset, SEEK_SET);
    
	if(fp) int now = fread(blocks[0], BLOCK_SIZE, read_nums, fp);
    else{
        fp = fopen(fileName, "rb");
	    fseeko64(fp, 0L, (uint64_t) seekForBlock * BLOCK_SIZE);
        int now = fread(blocks[0], BLOCK_SIZE, read_nums, fp);
    }
    // for(int j = 0; j < read_nums; j++){
    //     int byteRead = fread(blocks[j], BLOCK_SIZE, 1, fp);
    // }
    seekForBlock += read_nums;
    // cout << "the blocks of this read is: " << read_nums <<endl;
    return blocks;  
}


#endif // __DEFINE_H__