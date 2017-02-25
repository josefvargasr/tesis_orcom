/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#include "Globals.h"

#include <string.h>
#include <algorithm>

#include "DnaCategorizer.h"
#include "DnaBlockData.h"


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
	DistributeToBins(records_, recordsCount_, bin_.stdBins, bin_.nBin);		// Lee cada Read del vector records_, encuentra su minimizador y guarda el Read en el Bin correspondiente

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
}


// todo: split into rev and non-rev
void DnaCategorizer::DistributeToBins(std::vector<DnaRecord>& records_, uint64 recordsCount_, DnaBinCollection& bins_, DnaBin& nBin_)
{
	char revBuffer[1024];	// TODO: make size constant depending on the record max len
	DnaRecord rcRec;		//Objeto que guardar el complemento reverso del Read
	rcRec.dna = revBuffer;
	rcRec.reverse = true;		//Indica que el objeto guarda el Read de forma inversa

	if (params.tryReverseCompliment)							// Si se usa el Read en forma reversa
	{
		for (uint32 i = 0; i < recordsCount_; ++i)				
		{
			DnaRecord& rec = records_[i];						//Objeto que guardar el Read de forma directa
			rec.reverse = false;

			ASSERT(rec.len > 0);

			rec.ComputeRC(rcRec);								//Se encarga de buscar el complemento inverso del Read

			// find and select minimizers
			//
			const uint32 minimizerFwd = FindMinimizer(rec);		// Busca el minimizador del read en forma directa
			const uint32 minimizerRev = FindMinimizer(rcRec);	// Busca el minimizador del complemento inverso del Read
			uint32 minimizer = 0;
			bool reverse = false;

			if (minimizerFwd <= minimizerRev)					// Verifica si elige minimizerFwd o minimizerRev
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
			if (minimizer != nBinValue)								// !TODO --- find here minimizer pos
			{
				if (reverse)
				{
					rec.reverse = true;
					std::copy(rcRec.dna, rcRec.dna + rec.len, rec.dna);
				}

				bins_[minimizer].Insert(rec);					// Guarda el read en el Bin corresspondiente dependiedo del minimizador
			}
			else
			{
				rec.reverse = false;
				rec.minimizerPos = 0;
				nBin_.Insert(rec);								// Guarda el read en el Bin corresspondiente si no tiene minimizador
			}
		}
	}
	else
	{
		for (uint32 i = 0; i < recordsCount_; ++i)
		{
			DnaRecord& r = records_[i];
			ASSERT(r.len > 0);
			ASSERT(!r.reverse);
			uint32 minimizer = FindMinimizer(r);				// Encuentra el minimizador

			if (minimizer != nBinValue)										// !TODO --- find here minimizer pos
			{
				bins_[minimizer].Insert(r);						// Guarda el read en el Bin corresspondiente dependiedo del minimizador
			}
			else
			{
				r.minimizerPos = 0;
				nBin_.Insert(r);								// Guarda el read en el Bin corresspondiente si no tiene minimizador
			}
		}
	}

	// re-balance bins
	//
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
	
	for (int32 i = ibeg; i < iend; ++i)										// Itera hasta encontrar el minimizador en el read
	{
		uint32 m = ComputeMinimizer(rec_.dna + i, params.signatureLen);		// Encuentra el minimizador en el Read en cada iteracin y lo convierte en binario con sus respectivos Id

		if (m < minimizer && IsMinimizerValid(m, params.signatureLen))		// verifica si el minimizador hallado es valido
			minimizer = m;
	}

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
	uint32 r = 0;		// Minimizador Retornado en binario

	for (uint32 i = 0; i < mLen_; ++i)				// Recorre toda la cadena de caracteres
	{
		if (dna_[i] == 'N')							// Si alguna letra es 'N' retorna un valor predeterminado
			return nBinValue;

		ASSERT(dna_[i] >= 'A' && dna_[i] <= 'T');	// Si el caracter es correcto continua 
		r <<= 2;									// Desplaza 'r' dos bits a la izquierda
		r |= symbolIdxTable[(uint32)dna_[i]];		// Hace r = r OR symbolIdxTable[ACGT]
	}

	return r;										// Variable que contiene el minimizador hallado en binario
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
