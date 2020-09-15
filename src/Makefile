WARN_FLAGS = \
	-std=c++11 -O3 -g  \
	 -I../builds/sdsl-lite/include \
	  -I../src/BWTIL/extern/bitvector/include \
		 -L../builds/sdsl-lite/lib -lsdsl -ldivsufsort -ldivsufsort64 -lz

BIN = RLZI.x
OBJECTS = RLZI.o
BUILD_DIR = ../bin
# OUTF = $(BUILD_DIR)/RLZI.x
#BIN2 = test2
#OBJECTS2 = test2.o

#all : $(BIN) $(BIN2)
all: $(BUILD_DIR)/$(BIN)

$(BUILD_DIR)/$(BIN): $(OBJECTS)
	$(MKDIR_P) $(dir $@)
	$(CXX) -o $@ $(OBJECTS) $(WARN_FLAGS)

#$(BIN2): $(OBJECTS2)
	#$(CXX) -o $@ $(OBJECTS2) $(WARN_FLAGS)

RLZI.o: RLZI.cpp optparse.h
	$(CXX) $(WARN_FLAGS) -c $< -o $@

# Check includes in multiple files.
# check: test.o
# 	echo '#include "optparse.h"' > tmp_second.cc
# 	$(CXX) -o tmp_test test.o tmp_second.cc
# 	rm tmp_second.cc
# 	./tmp_test
# 	rm tmp_test

clean:
	rm -f *.o $(BUILD_DIR)/*

.PHONY: clean

MKDIR_P = mkdir -p
