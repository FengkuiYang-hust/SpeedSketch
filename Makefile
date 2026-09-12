CXX ?= g++
CC ?= gcc

TARGET ?= speedsketch
TEST_TARGET ?= test_gdelta
CORE_TEST_TARGET ?= test_core
BUILD_DIR ?= .build
OPT_FLAGS ?= -O3
EXTRA_FLAGS ?=

CPPFLAGS := -I. -Iutils -Idelta/Gdelta -Idelta/xdelta3 -D_FILE_OFFSET_BITS=64
CXXFLAGS := -std=c++14 -pthread -Wall -Wextra -Wpedantic $(OPT_FLAGS) $(EXTRA_FLAGS)
CFLAGS := -std=c99 $(OPT_FLAGS) $(EXTRA_FLAGS)
LDLIBS := -pthread -lzstd

APP_OBJECTS := $(BUILD_DIR)/speedsketch.o $(BUILD_DIR)/gdelta.o $(BUILD_DIR)/xdelta3.o
TEST_OBJECTS := $(BUILD_DIR)/test_gdelta.o $(BUILD_DIR)/gdelta.o
CORE_TEST_OBJECTS := $(BUILD_DIR)/test_core.o $(BUILD_DIR)/gdelta.o $(BUILD_DIR)/xdelta3.o
DEPS := $(APP_OBJECTS:.o=.d) $(TEST_OBJECTS:.o=.d) $(CORE_TEST_OBJECTS:.o=.d)

.PHONY: all check asan tsan clean

all: $(TARGET)

$(TARGET): $(APP_OBJECTS)
	$(CXX) $(EXTRA_FLAGS) $^ $(LDLIBS) -o $@

$(TEST_TARGET): $(TEST_OBJECTS)
	$(CXX) $(EXTRA_FLAGS) $^ -pthread -o $@

$(CORE_TEST_TARGET): $(CORE_TEST_OBJECTS)
	$(CXX) $(EXTRA_FLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/speedsketch.o: speedsketch.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/gdelta.o: delta/Gdelta/gdelta.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/xdelta3.o: delta/xdelta3/xdelta3.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -w -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_gdelta.o: tests/test_gdelta.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_core.o: tests/test_core.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

check: $(TARGET) $(TEST_TARGET) $(CORE_TEST_TARGET)
	./$(TEST_TARGET)
	timeout 30 ./$(CORE_TEST_TARGET)
	python3 tests/test_cli.py --binary ./$(TARGET)

asan:
	ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 $(MAKE) \
		TARGET=speedsketch_asan TEST_TARGET=test_gdelta_asan CORE_TEST_TARGET=test_core_asan \
		BUILD_DIR=.build/asan \
		OPT_FLAGS='-O1 -g' EXTRA_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' check

tsan:
	$(MAKE) \
		TARGET=speedsketch_tsan TEST_TARGET=test_gdelta_tsan CORE_TEST_TARGET=test_core_tsan \
		BUILD_DIR=.build/tsan \
		OPT_FLAGS='-O1 -g' EXTRA_FLAGS='-fsanitize=thread -fno-omit-frame-pointer' \
		speedsketch_tsan test_gdelta_tsan test_core_tsan
	TSAN_OPTIONS=halt_on_error=1 setarch $(shell uname -m) -R ./test_gdelta_tsan
	TSAN_OPTIONS=halt_on_error=1 timeout 30 setarch $(shell uname -m) -R ./test_core_tsan
	TSAN_OPTIONS=halt_on_error=1 setarch $(shell uname -m) -R \
		python3 tests/test_cli.py --binary ./speedsketch_tsan

clean:
	$(RM) speedsketch speedsketch_asan speedsketch_tsan
	$(RM) test_gdelta test_gdelta_asan test_gdelta_tsan
	$(RM) test_core test_core_asan test_core_tsan
	$(RM) -r .build

-include $(DEPS)
