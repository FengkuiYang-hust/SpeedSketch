
CXX = g++
CLANGXX=clang++
CXX_FLAG = -pthread -std=c++14 -mcx16 -Wno-invalid-offsetof -latomic -lcrypto -lz -lzstd -L./lib/ -I./delta/ -I./utils/
OPT_FLAG = -O3

tag=_test

SpeedSketchSeqSRC = ./speedSketch_seq.cpp

SpeedSketchThreadSRC = ./speedSketch_threads.cpp

speedSketch_seq:
	$(CXX) ${SpeedSketchSeqSRC} -o ./speedSketch_seq${tag} $(CXX_FLAG) $(OPT_FLAG) $(MALLOC_FLAG)

speedSketch_threads:
	$(CXX) ${SpeedSketchThreadSRC} -o ./speedSketch_threads${tag} $(CXX_FLAG) $(OPT_FLAG) $(MALLOC_FLAG)


clean:
	rm -rf speedSketch_*${tag}


	
