
CXX = g++
CLANGXX=clang++
CXX_FLAG = -pthread -std=c++17 -mcx16 -Wno-invalid-offsetof -latomic -lcrypto -lz -lzstd -lzd -L../lib/ -I../delta/ -I../utils/
OPT_FLAG = -O2


DeltaTimeDis = ./delta_time_dis.cpp ../delta/edelta/util/*.cc

all_deltaTimeDis:
	$(CXX) $(DeltaTimeDis) -o ./all_deltaTimeDis$(tag) $(CXX_FLAG) $(OPT_FLAG) $(MALLOC_FLAG)

clean:
	rm -rf odess${tag}


	
