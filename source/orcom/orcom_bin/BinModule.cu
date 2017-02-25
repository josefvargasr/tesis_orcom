/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#include "Globals.h"

#include <vector>
#include <iostream>
#include <sys/time.h>

#include "BinModule.h"
#include "FastqStream.h"
#include "DnaParser.h"
#include "DnaPacker.h"
#include "DnaCategorizer.h"
#include "BinFile.h"
#include "BinOperator.h"
#include "DnaBlockData.h"
#include "Exception.h"
#include "Thread.h"

// CUDA
//
#ifdef CUDA
	//extern char* h_reads;
	extern uint32* arr_minim; 
#endif

extern uint32 countReads;
//clock_t start2, end1, start_1, end_2;
//extern clock_t end_1, start_2;

struct timeval start2, end1, start_1, end_2;
extern struct timeval end_1, start_2;
int diff_1, diff_2, diff_loop;

void BinModule::Fastq2Bin(const std::vector<std::string> &inFastqFiles_, const std::string &outBinFile_,
						  uint32 threadNum_,  bool compressedInput_, bool verboseMode_)
{
	// TODO: try/catch to free resources
	//
	IFastqStreamReader* fastqFile = NULL;
	if (compressedInput_)
		fastqFile = new MultiFastqFileReaderGz(inFastqFiles_);		// Abre el archivo Fastq Comprimido
	else
		fastqFile = new MultiFastqFileReader(inFastqFiles_);		// Abre el archivo Fastq sin Comprimir


	BinFileWriter binFile;
	binFile.StartCompress(outBinFile_, config);		// Prepara las configuraciones Iniciales del proceso de compresion

	const uint32 minimizersCount = config.minimizer.TotalMinimizersCount();
	if (threadNum_ > 1)
	{
		FastqChunkPool* fastqPool = NULL;
		FastqChunkQueue* fastqQueue = NULL;
		BinaryPartsPool* binPool = NULL;
		BinaryPartsQueue* binQueue = NULL;

		FastqChunkReader* fastqReader = NULL;
		BinChunkWriter* binWriter = NULL;

		const uint32 partNum = threadNum_ * 4;
		fastqPool = new FastqChunkPool(partNum, config.fastqBlockSize);
		fastqQueue = new FastqChunkQueue(partNum, 1);

		binPool = new BinaryPartsPool(partNum, minimizersCount);
		binQueue = new BinaryPartsQueue(partNum, threadNum_);

		fastqReader = new FastqChunkReader(fastqFile, fastqQueue, fastqPool);
		binWriter = new BinChunkWriter(&binFile, binQueue, binPool);

		// launch stuff
		//
		mt::thread readerThread(mt::ref(*fastqReader));

		std::vector<IOperator*> operators;
		operators.resize(threadNum_);

#ifdef USE_BOOST_THREAD
		boost::thread_group opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new BinEncoder(config.minimizer, config.catParams,
										  fastqQueue, fastqPool,
										  binQueue, binPool);
			opThreadGroup.create_thread(mt::ref(*operators[i]));
		}

		(*binWriter)();

		readerThread.join();
		opThreadGroup.join_all();


#else
		std::vector<mt::thread> opThreadGroup;

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			operators[i] = new BinEncoder(config.minimizer, config.catParams,
										  fastqQueue, fastqPool, binQueue, binPool);
			opThreadGroup.push_back(mt::thread(mt::ref(*operators[i])));
		}

		(*binWriter)();

		readerThread.join();

		for (mt::thread& t : opThreadGroup)
		{
			t.join();
		}

#endif

		for (uint32 i = 0; i < threadNum_; ++i)
		{
			delete operators[i];
		}

		TFREE(binWriter);
		TFREE(fastqReader);

		TFREE(binQueue);
		TFREE(binPool);
		TFREE(fastqQueue);
		TFREE(fastqPool);
	}
	else
	{
		DnaParser parser;
		DnaCategorizer categorizer(config.minimizer, config.catParams);
		DnaPacker packer(config.minimizer);

		DataChunk fastqChunk(config.fastqBlockSize);
		std::vector<DnaRecord> records;										// Vector de objetos DnaRecords, que guarda los Reads leidos
		records.resize(1 << 10);
		//h_reads = (char*)malloc(10 * sizeof(char));
		DnaBinBlock dnaBins(minimizersCount);
		BinaryBinBlock binBins;
		DataChunk dnaBuffer;

#ifdef CUDA
		printf("\n----- Start  CUDA -----\n");
#else
		printf("\n----- Start ORCOM -----\n");
#endif				
		countReads=0;
		diff_1=0;
		diff_2=0;
		
		struct timeval start;
		gettimeofday(&start, NULL);
		gettimeofday(&end1, NULL);
		
		while (fastqFile->ReadNextChunk(&fastqChunk))						// Lee Partes del archivo Fastq
		{
			gettimeofday(&start_1, NULL);
			uint64 recordsCount = 0;
			parser.ParseFrom(fastqChunk, dnaBuffer, records, recordsCount);	// Toma las partes leidas y las pone en un vector de objetos, cada objeto contiene un Read
			
			ASSERT(recordsCount > 0);
			categorizer.Categorize(records, recordsCount, dnaBins);			// Busca los Minimizers y los almacena en Bins

			packer.PackToBins(dnaBins, binBins);							// Empaqueta los Bins

			binFile.WriteNextBlock(&binBins);								// Escribe el paquete de Bins
		
			gettimeofday(&end_2, NULL);
			diff_1 += (end_1.tv_sec - start_1.tv_sec);
			diff_2 += (end_2.tv_sec - start_2.tv_sec);
		}
		gettimeofday(&start2, NULL);
		struct timeval end;
		gettimeofday(&end, NULL);
		diff_loop = double(end.tv_sec - start.tv_sec);
		

		printf("Number of Minimizers processed: %d\n",countReads);
		printf("Number of Minimizers in file: %d\n",countReads/2);
#ifdef CUDA
		printf("\n----- End  CUDA -----\n");
#else
		printf("\n----- End ORCOM -----\n");
#endif
	}

	binFile.FinishCompress();

	if (verboseMode_)
	{
		std::vector<uint64> recordCounts;
		binFile.GetBinStats(recordCounts);

		std::cout << "Signatures count: " << recordCounts.size() << std::endl;
		std::cout << "Records distribution in bins by signature:\n";
		for (uint32 i = 0; i < recordCounts.size(); ++i)
		{
			if (recordCounts[i] > 0)
				std::cout << i << " : " << recordCounts[i] << '\n';
		}
		std::cout << std::endl;
	}

	delete fastqFile;
}


void BinModule::Bin2Dna(const std::string &inBinFile_, const std::string &outDnaFile_)
{
	// TODO: try/catch to free resources
	//
	BinFileReader binFile;

	binFile.StartDecompress(inBinFile_, config);
	uint32 minimizersCount = config.minimizer.TotalMinimizersCount();

	DnaFileWriter dnaFile(outDnaFile_);
	DataChunk fastqChunk(config.fastqBlockSize >> 1);			// WARNING! --- here can be a BUG
	DnaPacker packer(config.minimizer);
	DnaParser parser;

	DnaBinBlock dnaBins(minimizersCount);
	BinaryBinBlock binBins;
	DataChunk dnaBuffer;

	while (binFile.ReadNextBlock(&binBins))
	{
		packer.UnpackFromBins(binBins, dnaBins, dnaBuffer);
		parser.ParseTo(dnaBins, fastqChunk);

		dnaFile.WriteNextChunk(&fastqChunk);
	}

	dnaFile.Close();
	binFile.FinishDecompress();
}
