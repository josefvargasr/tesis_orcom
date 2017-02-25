/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#include "Globals.h"

#include <string.h>
#include <algorithm>
#include <sys/time.h>

#include "DnaCategorizer.h"
#include "DnaBlockData.h"

// CUDA
// -----------------------------------
uint32 countReads = 0;

double total_f1 = 0, total_f2 = 0, total_f3 = 0, diff_f1 = 0, diff_f2 = 0, diff_f3 = 0;
uint32 countCategorize = 0;

struct timeval start_2, end_1, start_bin, start_bin2, end_bin, end_bin2;
int diff_bin, diff_bin2, diff_f33, diff__f1 = 0, total_f33 = 0;

#ifdef CUDA
	__device__ char d_symbolIdxTable[128];
	__device__ uint32 d_nBinValue;
	__device__ uint32 d_maxLongMinimValue;
	__device__ uint32 d_maxShortMinimValue;			
	__device__ uint8 d_params_signatureLen;
	__device__ uint8 d_params_skipZoneLen; 
	uint64 h_dnaSize;	
	char* h_reads;
	char* d_reads;char* d_readsRC;
	uint32* h_posReads;
	uint32* d_posReads;
	uint16* h_lenReads;
	uint16* d_lenReads;
	uint32* h_arr_minim;
	uint32* d_arr_minim;
	#define numStream	8	
	#define STREAMS		1	
#endif
// -----------------------------------

DnaCategorizer::DnaCategorizer(const MinimizerParameters& params_, const CategorizerParameters& catParams_)
	:	params(params_)
	,	catParams(catParams_)
	,	maxShortMinimValue(1 << (2 * params.signatureSuffixLen))
	,	maxLongMinimValue(1 << (2 * params.signatureLen))
	,	nBinValue(maxLongMinimValue)
{
	ASSERT(params.signatureSuffixLen <= params.signatureLen);

	std::fill(symbolIdxTable, symbolIdxTable + 128, -1);		// Llena el vector symbolIdxTable con '-1'
	for (uint32 i = 0; i < 5; ++i)
		symbolIdxTable[(int32)params.dnaSymbolOrder[i]] = i;	// En la posicion de cada letra (ACGTN) poner su respectivo Id
																// A = symbolIdxTable[65] = 0; C = symbolIdxTable[67] = 1; G = symbolIdxTable[71] = 2;
																// T = symbolIdxTable[84] = 3; N = symbolIdxTable[78] = 4;
	freqTable.resize(maxShortMinimValue, 0);
}


void DnaCategorizer::Categorize(std::vector<DnaRecord>& records_, uint64 recordsCount_, DnaBinBlock& bin_)
{
	
	ASSERT(recordsCount_ > 0);
	ASSERT(recordsCount_ <= records_.size());

	// clear bins
	//
	for (uint32 i = 0; i < bin_.stdBins.Size(); ++i)
		bin_.stdBins[i].Clear();
	bin_.nBin.Clear();

	std::fill(freqTable.begin(), freqTable.end(), 0);

	// process records
	//
	countCategorize++;
	
	clock_t start_f1 = clock();

	gettimeofday(&end_1, NULL);
	DistributeToBins(records_, recordsCount_, bin_.stdBins, bin_.nBin);		// Lee cada Read del vector records_, encuentra su minimizador y guarda el Read en el Bin correspondiente

	gettimeofday(&start_2, NULL);

	clock_t end_f1 = clock();	
	diff_f1 = 0;
	diff_f1 = double(end_f1 - start_f1)/CLOCKS_PER_SEC;
	total_f1 += diff_f1;

	// sort the bins
	//

	FindMinimizerPositions(bin_.stdBins);
	
	
	DnaRecordComparator comparator(params.signatureLen - params.signatureSuffixLen);


	for (uint32 i = 0; i < bin_.stdBins.Size(); ++i)
	{
		DnaBin& db = bin_.stdBins[i];
		if (db.Size() > 0)
			std::sort(db.Begin(), db.End(), comparator);
	}
		
	
	printf("F1 Time: %fs, Total: %fs ---- DistributeToBins\n", diff_f1, total_f1);
	printf("F2 Time: %fs, Total: %fs ---- FindMinimizer\n", diff_f2, total_f2);
	printf("F3 Time: %fs, Total: %fs ---- PrepReads, countCategorize: %d\n", diff_f3, total_f3, countCategorize);	
	printf("F33 Time: %ds, Total: %ds ---- GPU\n", total_f33, diff_f33);

}

#ifdef CUDA

#define HANDLE_ERROR(err) handleError(err, __FILE__, __LINE__)

static void handleError(cudaError_t err, const char *file, int line){
	if(err != cudaSuccess){
		printf("Cuda error at %s:%d: %s\n", file, line, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}

void getCudaMemInfo(){
	size_t free_byte, total_byte;

	cudaMemGetInfo(&free_byte, &total_byte);

	double free = (double) free_byte;
	double total = (double) total_byte;
	double used = total - free;

	printf("Cuda Memory. Used: %fMB, Free: %fMB, Total: %fMB\n", used/1024/1024, free/1024/1024, total/1024/1024);
}

#define N_THREADS	512//256//512//1024
#define N_BLOCKS	32767//16384//32767//65535

__device__ bool d_IsMinimizerValid(uint32_t minim_, uint32_t mLen_)
{
        if(minim_ == 0)
                minim_ = 0;

        const uint32_t excludeLen = 3;
        const uint32_t excludeMask = 0x3F;
        const uint32_t symbolExcludeTable[] = {0x00, 0x15, 0x2A, 0x3F};

        minim_ &= (1 << (2 * mLen_)) - 1;
        bool hasInvalidSeq = false;
        for (uint32_t i = 0; !hasInvalidSeq && i <= (mLen_ - excludeLen); ++i){
                uint32_t x = minim_ & excludeMask;

                for(uint32_t j = 0; j < 4; ++j){
                        hasInvalidSeq |= (x == symbolExcludeTable[j]);
                }

                minim_ >>= 2;
        }

        return !hasInvalidSeq;
}

__device__ uint32 d_ComputeMinimizer(const char* dna_, uint32_t mLen_)
{
        uint32_t r = 0;

        for(uint32_t i = 0; i < mLen_; ++i){

                if(dna_[i] == 'N')
                        return d_nBinValue;

                r <<= 2;
                r |= d_symbolIdxTable[(uint32_t)dna_[i]];
        }

        return r;
}

__device__ uint32 d_FindMinimizer(char * d_dna_, uint16 d_len_)
{
	uint32 minimizer = d_maxLongMinimValue;
	#if EXP_USE_RC_ADV
		const int32 ibeg = 0; //rec_.reverse ? d_params_skipZoneLen : 0;
		const int32 iend = d_len_ - d_params_signatureLen + 1 - (rec_.reverse ? 0 : d_params_skipZoneLen);
	#else	
		const int32 ibeg = 0;
		const int32 iend = d_len_ - d_params_signatureLen + 1 - d_params_skipZoneLen;
	#endif
	for (int32 i = ibeg; i < iend; ++i){
		uint32 m = d_ComputeMinimizer(d_dna_ + i, d_params_signatureLen);		
		if (m < minimizer && d_IsMinimizerValid(m, d_params_signatureLen))		
			minimizer = m;
	}
	
	if (minimizer >= d_maxLongMinimValue)
		return d_nBinValue;

	return minimizer & (d_maxShortMinimValue - 1);
}

__device__ char* d_ComputeRC(char* d_reads_, char* d_readsRC_, uint16 d_len_)
{
	const char rcCodes[24] = {-1,'T',-1,'G',-1,-1,-1,'C',
				  -1,-1,-1,-1,-1,-1,'N',-1,
				  -1,-1,-1,-1,'A',-1,-1,-1,
				};
		
	for(uint32 i = 0; i < d_len_; ++i){
		d_readsRC_[d_len_-1-i] = rcCodes[(int32)d_reads_[i] - 64];
	}
	
	return d_readsRC_;
}

#ifndef STREAMS
__global__ void d_DistributeToBins(uint32 n_reads_, uint32* d_arr_minim_, uint16* d_lenReads_, uint32* d_posReads_, char* d_reads, char* d_readsRC)
#else
__global__ void d_DistributeToBins(uint32 n_reads_, uint32* d_arr_minim_, uint16* d_lenReads_, uint32* d_posReads_, char* d_reads, char* d_readsRC, uint32 offset, uint32 sizeOffset2)
#endif
{
#ifdef STREAMS
	uint32 tid = offset +  threadIdx.x + (blockDim.x * blockIdx.x);
	uint32 tid2 = threadIdx.x + (blockDim.x * blockIdx.x);
	if(tid < n_reads_ && tid2 < sizeOffset2 )
#else
	uint32 tid = threadIdx.x + (blockDim.x * blockIdx.x);
	if(tid < n_reads_ )
#endif
		
	{
		d_arr_minim_[tid] = d_FindMinimizer(&d_reads[d_posReads_[tid]], d_lenReads_[tid]);
		d_arr_minim_[tid + n_reads_] = d_FindMinimizer(d_ComputeRC(&d_reads[d_posReads_[tid]], &d_readsRC[d_posReads_[tid]], d_lenReads_[tid]), d_lenReads_[tid]);
	}
}

#endif

// todo: split into rev and non-rev
void DnaCategorizer::DistributeToBins(std::vector<DnaRecord>& records_, uint64 recordsCount_, DnaBinCollection& bins_, DnaBin& nBin_)
{	
	struct timeval start__f1;
	gettimeofday(&start__f1, NULL);

	char revBuffer[1024];		// TODO: make size constant depending on the record max len
	DnaRecord rcRec;		//Objeto que guardar el complemento reverso del Read
	rcRec.dna = revBuffer;
	rcRec.reverse = true;		//Indica que el objeto guarda el Read de forma inversa
	
#ifdef CUDA
	uint32 n_reads = 0;
	uint16 h_len_ = records_[0].len;
	if(countReads == 0){
		printf("Longitud del primer Read: %d\n", h_len_);
		// Alojar variables globales en GPU
		uint32* dd_maxLongMinimValue;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_maxLongMinimValue, d_maxLongMinimValue));
		HANDLE_ERROR(cudaMemcpy(dd_maxLongMinimValue, &maxLongMinimValue, sizeof(uint32), cudaMemcpyHostToDevice));
		uint32* dd_maxShortMinimValue;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_maxShortMinimValue, d_maxShortMinimValue));
		HANDLE_ERROR(cudaMemcpy(dd_maxShortMinimValue, &maxShortMinimValue, sizeof(uint32), cudaMemcpyHostToDevice));
		uint8* dd_params_signatureLen;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_params_signatureLen, d_params_signatureLen));
		HANDLE_ERROR(cudaMemcpy(dd_params_signatureLen, &params.signatureLen, sizeof(uint8), cudaMemcpyHostToDevice));
		uint8* dd_params_skipZoneLen;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_params_skipZoneLen, d_params_skipZoneLen));
		HANDLE_ERROR(cudaMemcpy(dd_params_skipZoneLen, &params.skipZoneLen, sizeof(uint8), cudaMemcpyHostToDevice));
		char* dd_symbolIdxTable;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_symbolIdxTable, d_symbolIdxTable));
		HANDLE_ERROR(cudaMemcpy(dd_symbolIdxTable, symbolIdxTable, 128 * sizeof(char), cudaMemcpyHostToDevice));
		uint32* dd_nBinValue;
		HANDLE_ERROR(cudaGetSymbolAddress((void**)&dd_nBinValue, d_nBinValue));
		HANDLE_ERROR(cudaMemcpy(dd_nBinValue, &nBinValue, sizeof(uint32), cudaMemcpyHostToDevice));
		
		HANDLE_ERROR(cudaMalloc((void**)&d_readsRC, h_dnaSize * sizeof(char)));
		HANDLE_ERROR(cudaMalloc((void**)&d_reads, h_dnaSize * sizeof(char)));
		HANDLE_ERROR(cudaMalloc((void**)&d_lenReads, recordsCount_ * sizeof(uint16)));
		h_arr_minim = (uint32*)malloc(2 * recordsCount_ * sizeof(uint32));
		h_posReads = (uint32*)malloc(recordsCount_ * sizeof(uint32));
		HANDLE_ERROR(cudaMalloc((void**)&d_posReads, recordsCount_ * sizeof(uint32)));
		HANDLE_ERROR(cudaMalloc((void**)&d_arr_minim, 2 * recordsCount_ * sizeof(uint32)));
	}	
#endif
	
	struct timeval end__f1;
	gettimeofday(&end__f1, NULL);
	diff__f1 = (end__f1.tv_sec - start__f1.tv_sec);

	if (params.tryReverseCompliment)					// Si se usa el Read en forma reversa
	{
	#ifdef CUDA
		clock_t start_f3 = clock();

		uint32 pos = 0;
		for(uint32 i = 0; i < recordsCount_; ++i){
			h_posReads[i] = pos;
			pos += records_[i].len;
			h_lenReads[i] = records_[i].len;
		}

		n_reads = recordsCount_;
		countReads += n_reads * 2;
		clock_t end_f3 = clock();	
		diff_f3 = 0;
		diff_f3 = double(end_f3 - start_f3)/CLOCKS_PER_SEC;
		total_f3 += diff_f3;

		struct timeval start_f33;
		gettimeofday(&start_f33, NULL);

	#ifndef STREAMS

		// Copiar Datos a la GPU
		HANDLE_ERROR(cudaMemcpy(d_reads, h_reads, h_dnaSize * sizeof(char), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(d_lenReads, h_lenReads, recordsCount_ * sizeof(uint16), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(d_posReads, h_posReads, recordsCount_ * sizeof(uint32), cudaMemcpyHostToDevice));

		// Lanzar Kernel
		d_DistributeToBins<<<((n_reads) + N_THREADS-1)/N_THREADS, N_THREADS>>>(n_reads, d_arr_minim, d_lenReads, d_posReads, d_reads, d_readsRC);
		
		HANDLE_ERROR(cudaGetLastError());
		HANDLE_ERROR(cudaDeviceSynchronize());
		//getCudaMemInfo();
	
		// Copiar Datos a la CPU
		HANDLE_ERROR(cudaMemcpy(h_arr_minim, d_arr_minim, 2 * recordsCount_ * sizeof(uint32), cudaMemcpyDeviceToHost));

	#else

		const int n = (n_reads) + (N_THREADS-1);
		const int nStreams = numStream;

		// Crear Streams
		cudaStream_t stream[numStream];
		for (uint32 i = 0; i < nStreams; ++i)
			cudaStreamCreate(&stream[i]);
		
		// Copiar Datos a la GPU y Lanzar Kernel por cada Stream
		for (uint32 i = 0; i < nStreams; ++i)
		{	
			uint32 offset1 = i * (h_dnaSize/nStreams);
			uint32 offset2 = i * (recordsCount_/nStreams);
			uint32 sizeOffset1;
			uint32 sizeOffset2;
			uint32 bSize = (((n)/N_THREADS)/nStreams) + 1;

			if(i == nStreams - 1){
				sizeOffset1 = h_dnaSize/nStreams;
				sizeOffset1 = sizeOffset1 + (h_dnaSize - (sizeOffset1 * nStreams));
				sizeOffset2 = recordsCount_/nStreams;
				sizeOffset2 = sizeOffset2 + (recordsCount_ - (sizeOffset2 * nStreams));
			}else{
				sizeOffset1 = h_dnaSize/nStreams;
				sizeOffset2 = recordsCount_/nStreams;
			}

			HANDLE_ERROR(cudaMemcpyAsync(&d_reads[offset1], &h_reads[offset1], sizeOffset1 * sizeof(char), cudaMemcpyHostToDevice, stream[i]));
			HANDLE_ERROR(cudaMemcpyAsync(&d_lenReads[offset2], &h_lenReads[offset2], sizeOffset2 * sizeof(uint16), cudaMemcpyHostToDevice, stream[i]));
			HANDLE_ERROR(cudaMemcpyAsync(&d_posReads[offset2], &h_posReads[offset2], sizeOffset2 * sizeof(uint32), cudaMemcpyHostToDevice, stream[i]));
			d_DistributeToBins<<< bSize, N_THREADS, 0, stream[i]>>>(n_reads, d_arr_minim, d_lenReads, d_posReads, d_reads, d_readsRC, offset2, sizeOffset2);
		}	
		
		for (uint32 i = 0; i < nStreams; ++i)
			HANDLE_ERROR(cudaStreamSynchronize(stream[i]));
		
		// Copiar Datos a la CPU
		for (uint32 i = 0; i < nStreams; ++i)	
		{
			uint32 sizeOffset2;
			uint32 offset2 = i * 2 * (recordsCount_/nStreams);

			if(i == nStreams - 1){
				sizeOffset2 = recordsCount_/nStreams;
				sizeOffset2 = sizeOffset2 + (recordsCount_ - (sizeOffset2 * nStreams));
			}else{
				sizeOffset2 = recordsCount_/nStreams;
			}	
			HANDLE_ERROR(cudaMemcpyAsync(&h_arr_minim[offset2], &d_arr_minim[offset2], 2 * sizeOffset2 * sizeof(uint32), cudaMemcpyDeviceToHost, stream[i]));
		}

		
		// Destruir Streams
		for (uint32 i = 0; i < nStreams; ++i)
			cudaStreamDestroy(stream[i]);
		
		
	#endif
		struct timeval end_f33;
		gettimeofday(&end_f33, NULL);

		total_f33 = (end_f33.tv_sec - start_f33.tv_sec);
		diff_f33 += total_f33;

		gettimeofday(&start_bin, NULL);
				
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			DnaRecord& rec = records_[i];				//Objeto que guardar el Read de forma directa
			rec.reverse = false;
	
			rec.ComputeRC(rcRec);
			
			const uint32 minimizerFwd = h_arr_minim[(i)];			
			const uint32 minimizerRev = h_arr_minim[(i) + (recordsCount_)];
			uint32 minimizer = 0;
			bool reverse = false;
		
			if (minimizerFwd <= minimizerRev)			
			{
				minimizer = minimizerFwd;
			}
			else
			{
				minimizer = minimizerRev;
				reverse = true;
			}

			// store record to bin
			//
			if (minimizer != nBinValue)				// !TODO --- find here minimizer pos
			{
				if (reverse)
				{
					rec.reverse = true;
					std::copy(rcRec.dna, rcRec.dna + rec.len, rec.dna);
				}

				bins_[minimizer].Insert(rec);			// Guarda el read en el Bin corresspondiente dependiedo del minimizador
			}
			else
			{
				rec.reverse = false;
				rec.minimizerPos = 0;
				nBin_.Insert(rec);				// Guarda el read en el Bin corresspondiente si no tiene minimizador
			}
		}
		
		gettimeofday(&end_bin, NULL);
  		int diffa = (end_bin.tv_sec - start_bin.tv_sec);
		diff_bin += diffa;
	
	#else
		
		
		for (uint32 i = 0; i < recordsCount_; ++i){
			DnaRecord& rec = records_[i];				//Objeto que guardar el Read de forma directa
			rec.reverse = false;

			ASSERT(rec.len > 0);

			rec.ComputeRC(rcRec);					//Se encarga de buscar el complemento inverso del Read

			// find and select minimizers
			//
			const uint32 minimizerFwd = FindMinimizer(rec);		// Busca el minimizador del read en forma directa
			const uint32 minimizerRev = FindMinimizer(rcRec);	// Busca el minimizador del complemento inverso del Read

			uint32 minimizer = 0;
			bool reverse = false;

			if (minimizerFwd <= minimizerRev)			// Verifica si elige minimizerFwd o minimizerRev
			{
				minimizer = minimizerFwd;
			}
			else
			{
				minimizer = minimizerRev;
				reverse = true;
			}

			// store record to bin
			//
			if (minimizer != nBinValue)				// !TODO --- find here minimizer pos
			{
				if (reverse)
				{
					rec.reverse = true;
					std::copy(rcRec.dna, rcRec.dna + rec.len, rec.dna);
				}

				bins_[minimizer].Insert(rec);			// Guarda el read en el Bin corresspondiente dependiedo del minimizador
			}
			else
			{
				rec.reverse = false;
				rec.minimizerPos = 0;
				nBin_.Insert(rec);				// Guarda el read en el Bin corresspondiente si no tiene minimizador
			}
		}
	#endif

	}
	else
	{
		printf("NO tryReverseCompliment\n");
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			DnaRecord& r = records_[i];
			ASSERT(r.len > 0);
			ASSERT(!r.reverse);
			uint32 minimizer = FindMinimizer(r);			// Encuentra el minimizador

			if (minimizer != nBinValue)				// !TODO --- find here minimizer pos
			{
				bins_[minimizer].Insert(r);			// Guarda el read en el Bin corresspondiente dependiedo del minimizador
			}
			else
			{
				r.minimizerPos = 0;
				nBin_.Insert(r);				// Guarda el read en el Bin corresspondiente si no tiene minimizador
			}
		}
	}

	// re-balance bins
	//
	gettimeofday(&start_bin2, NULL);

	for (uint32 i = 0; i < bins_.Size(); ++i)
	{
		DnaBin& db = bins_[i];
		if (db.Size() == 0 || db.Size() >= catParams.minBlockBinSize)
			continue;

		for (uint32 j = 0; j < db.Size(); ++j)
		{
			DnaRecord& r = db[j];
			r.minimizerPos = 0;

			// un-reverse the record
			if (r.reverse)
			{
				r.ComputeRC(rcRec);
				r.reverse = false;
				std::copy(rcRec.dna, rcRec.dna + r.len, r.dna);
			}

			std::map<uint32, uint16> mins = FindMinimizers(r);
			std::map<uint32, uint16>::iterator mit = mins.begin();

			for ( ; mit != mins.end(); ++mit)
			{
				const uint32 m = mit->first;
				if (bins_[m].Size() >= catParams.minBlockBinSize)// && bins_[m].Size() < maxBinSize)
				{
					r.minimizerPos = mit->second;
					bins_[m].Insert(r);
					break;
				}
			}

			// only one minimizer or we did not find appropriate bin
			if (mit == mins.end())
			{
				nBin_.Insert(r);
			}
		}
		db.Clear();
	}

	gettimeofday(&end_bin2, NULL);
	int diffb = (end_bin2.tv_sec - start_bin2.tv_sec);
	diff_bin2 += diffb;
		
}


void DnaCategorizer::FindMinimizerPositions(DnaBinCollection& bins_)
{
	
	for (uint32 i = 0; i < params.TotalMinimizersCount(); ++i)
	{
		if (bins_[i].Size() == 0)
			continue;

		char minString[64] = {0};
		params.GenerateMinimizer(i, minString);

		for (uint32 j = 0; j < bins_[i].Size(); ++j)
		{
			DnaRecord& r = bins_[i][j];
	#if EXP_USE_RC_ADV
			const char* beg = r.dna + (r.reverse ? params.skipZoneLen : 0);
			const char* end = r.dna + r.len - (r.reverse ? params.skipZoneLen : 0);
			const char* mi = std::search(beg, end, minString, minString + params.signatureSuffixLen);

			ASSERT(mi != r.dna + r.len);

			r.minimizerPos = mi - r.dna;
			ASSERT((!r.reverse && r.minimizerPos < r.len - params.skipZoneLen) ||
					(r.minimizerPos >= params.skipZoneLen) );
	#else
			const char* mi = std::search(r.dna, r.dna + r.len, minString, minString + params.signatureSuffixLen);

			ASSERT(mi != r.dna + r.len);

			r.minimizerPos = mi - r.dna;
			ASSERT(r.minimizerPos < r.len - params.skipZoneLen);
	#endif

		}
	}
}

uint32 DnaCategorizer::FindMinimizer(DnaRecord &rec_)
{
	
	uint32 minimizer = maxLongMinimValue;
	ASSERT(rec_.len >= params.signatureLen - params.skipZoneLen + 1);
	
	// Calculo de la posicion donde se comienza a iterar y numero de iteraciones que se realizarn en el Read para encontrar el minimizer
	#if EXP_USE_RC_ADV
		const int32 ibeg = rec_.reverse ? params.skipZoneLen : 0;
		const int32 iend = rec_.len - params.signatureLen + 1 - (rec_.reverse ? 0 : params.skipZoneLen);
	#else
		const int32 ibeg = 0;
		const int32 iend = rec_.len - params.signatureLen + 1 - params.skipZoneLen;
	#endif

	for (int32 i = ibeg; i < iend; ++i){
		uint32 m = ComputeMinimizer(rec_.dna + i, params.signatureLen);		
		if (m < minimizer && IsMinimizerValid(m, params.signatureLen))		
			minimizer = m;
	}

	countReads++;

	//printf("%d, Minimizer: %d\n", countReads, minimizer);
	
	if (minimizer >= maxLongMinimValue)
		return nBinValue;

	return minimizer & (maxShortMinimValue - 1);
}

std::map<uint32, uint16> DnaCategorizer::FindMinimizers(DnaRecord &rec_)
{
	ASSERT(rec_.len >= params.signatureLen - params.skipZoneLen + 1);

	// find all
	std::map<uint32, uint16> signatures;
	for (int32 i = 0; i < rec_.len - params.signatureLen + 1 - params.skipZoneLen; ++i)
	{
		uint32 m = ComputeMinimizer(rec_.dna + i, params.signatureLen);

		if (IsMinimizerValid(m, params.signatureLen) && m < maxLongMinimValue)
		{
			m &= (maxShortMinimValue - 1);
			if (signatures.count(m) == 0)
				signatures[m] = i + (params.signatureLen - params.signatureSuffixLen);
		}
	}

	return signatures;
}


uint32 DnaCategorizer::ComputeMinimizer(const char* dna_, uint32 mLen_)
{
	uint32 r = 0;						// Minimizador Retornado en binario

	for (uint32 i = 0; i < mLen_; ++i)			// Recorre toda la cadena de caracteres
	{
		if (dna_[i] == 'N')				// Si alguna letra es 'N' retorna un valor predeterminado
			return nBinValue;

		ASSERT(dna_[i] >= 'A' && dna_[i] <= 'T');	// Si el caracter es correcto continua 
		r <<= 2;					// Desplaza 'r' dos bits a la izquierda
		r |= symbolIdxTable[(uint32)dna_[i]];		// Hace r = r OR symbolIdxTable[ACGT]
	}

	return r;						// Variable que contiene el minimizador hallado en binario
}


bool DnaCategorizer::IsMinimizerValid(uint32 minim_, uint32 mLen_)
{
	if (minim_ == 0)
		minim_ = 0;

	const uint32 excludeLen = 3;
	const uint32 excludeMask = 0x3F;
	const uint32 symbolExcludeTable[] = {0x00, 0x15, 0x2A, 0x3F};	//0, 21, 42, 63

	minim_ &= (1 << (2*mLen_)) - 1;
	bool hasInvalidSeq = false;

	for (uint32 i = 0; !hasInvalidSeq && i <= (mLen_ - excludeLen); ++i)
	{
		uint32 x = minim_ & excludeMask;

		for (uint32 j = 0; j < 4; ++j)
			hasInvalidSeq |= (x == symbolExcludeTable[j]);

		minim_ >>= 2;
	}

	return !hasInvalidSeq;
}
