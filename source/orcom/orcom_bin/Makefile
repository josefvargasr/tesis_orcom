.PHONY: orcom_bin

all: orcom_bin

ifndef CXX
	#CXX = g++ -pg -flto -fwhole-program
	CXX = nvcc -pg 
endif
CXX = nvcc -pg 
CXX_FLAGS += -O3 -DNDEBUG 
CXX_FLAGS += -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
CXX_FLAGS += -std=c++11

CXX_OBJS = BinOperator.o \
	BinFile.o \
	DnaPacker.o \
	FastqStream.o \
	FileStream.o

CXX_OBJS2 = BinModule.o \
	BinOperator.o \
	BinFile.o \
	DnaPacker.o \
	DnaCategorizer.o \
	DnaParser.o \
	FastqStream.o \
	FileStream.o


CXX_LIBS += -lz
CXX_LIBS += -lpthread

.cpp.o:
	$(CXX) $(CXX_FLAGS) -c $< -o $@

orcom_bin: $(CXX_OBJS)
	$(CXX) $(CXX_FLAGS) -c DnaParser.cu -o DnaParser.o
	$(CXX) $(CXX_FLAGS) -c BinModule.cu -o BinModule.o
	$(CXX) $(CXX_FLAGS) -c DnaCategorizer.cu -o DnaCategorizer.o
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp $(CXX_OBJS2) $(CXX_LIBS)
	#strip $@
	#./orcom_bin e -isample_1.fastq -osample_1.bin
	#gprof orcom_bin gmon.out > analysis.txt
	
clean:
	-rm -f $(CXX_OBJS)
	-rm -f orcom_bin
	#-rm -f sample_1.bin.bmeta
	#-rm -f sample_1.bin.bdna
	#-rm -f gmon.out
	#-rm -f analysis.txt
