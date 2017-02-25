.PHONY: clean perf_microBRE_3 valgrind_microBRE_3 gprof_sample_1 gprof_microBRE_3 

all: fastqFile_small

perf_microBRE_3:
	sudo perf record ./orcom_bin e -iMicroBRE_3.fastq -oMicroBRE_3.bin
	sudo perf report

valgrind_microBRE_3:
	valgrind --tool=callgrind ./orcom_bin e -iMicroBRE_3.fastq -oMicroBRE_3.bin
	kcachegrind callgrind.out.*

gprof_sample_1: 
	-rm -f sample_1.bin.bmeta
	-rm -f sample_1.bin.bdna
	-rm -f gmon.out
	-rm -f analysis_sample_1.txt
	./orcom_bin e -isample_1.fastq -osample_1.bin
	gprof orcom_bin gmon.out > analysis_sample_1.txt
	nano analysis_sample_1.txt

gprof_microBRE_3: 
	-rm -f MicroBRE_3.bin.bmeta
	-rm -f MicroBRE_3.bin.bdna
	-rm -f gmon.out
	-rm -f analysis_sample_1.txt
	./orcom_bin e -iMicroBRE_3.fastq -oMicroBRE_3.bin
	gprof orcom_bin gmon.out > analysis_MicroBRE_3.txt
	nano analysis_MicroBRE_3.txt

#103 MB
m: 	
	./orcom_bin e -iMicroBRE_3.fastq -oMicroBRE_3.bin -t1
	#gprof orcom_bin gmon.out > analysis_MicroBRE_3.txt
	#nano analysis_MicroBRE_3.txt

#1.6 GB
n:	
	./orcom_bin e -iERR251006.fastq -oERR251006.bin -t1
	#gprof orcom_bin gmon.out > analysis_ERR251006.txt
	#nano analysis_ERR251006.txt

#8.0 GB
o:	
	./orcom_bin e -iERR160123_1.fastq -oERR160123_1.bin -t1
	#gprof orcom_bin gmon.out > analysis_ERR251006.txt
	#nano analysis_ERR251006.txt

#17 GB
p:	
	./orcom_bin e -iSRR608815.fastq -oSRR608815.bin -t1
	#gprof orcom_bin gmon.out > analysis_ERR251006.txt
	#nano analysis_ERR251006.txt

fastqFile_small: 
	./orcom_bin e -iFastqFile_small.fastq -oFastqFile_small.bin -t1
	#gprof orcom_bin gmon.out > analysis_FastqFile_small.txt
	#nano analysis_FastqFile_small.txt

clean:
	-rm -f sample_1.bin.bmeta
	-rm -f sample_1.bin.bdna
	-rm -f MicroBRE_3.bin.bmeta
	-rm -f MicroBRE_3.bin.bdna
	-rm -f gmon.out
	-rm -f analysis_MicroBRE_3.txt
	-rm -f analysis_MicroBRE_3.txt
