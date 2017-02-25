/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#ifndef H_DNAPACKER
#define H_DNAPACKER

#include "Globals.h"
#include "Params.h"
#include "BinBlockData.h"


class DnaPacker
{
public:
	DnaPacker(const MinimizerParameters& params_);
	/*
		Funci�n:
			Inicializa variables de la clase DnaPacker, se enumera cada letra (A, C, G, T, N)
		Entradas: 
			params_: Objeto de la clase MinimizerParameters
	*/
	
	void PackToBins(const DnaBinBlock& dnaBins_, BinaryBinBlock& binBins_);
	/*
		Funci�n:
			Guarda los paquetes de Bins para ser escritos en el archivo posteriormente
		Entradas: 
			dnaBins_: Bins sin empaquetar
			binBins_: Paquete de Bins que ser� escrito en el archivo
	*/
	void UnpackFromBins(const BinaryBinBlock& binBins_, DnaBinBlock& dnaBins_, DataChunk& dnaChunk_);
	void UnpackFromBin(const BinaryBinBlock& binBin_, DnaBin& dnaBin_, uint32 minimizerId_,
					   DataChunk& dnaChunk_, bool append_ = false);

private:
	struct BinPackSettings
	{
		bool hasConstLen;
		uint32 minLen;
		uint32 maxLen;
		uint32 suffixLen;
		uint32 bitsPerLen;

		BinPackSettings()
			:	hasConstLen(false)
			,	minLen((uint32)-1)
			,	maxLen(0)
			,	suffixLen(0)
			,	bitsPerLen(0)
		{}
	};

	static const uint32 LenBits = 8;

	const MinimizerParameters params;

	char dnaToIdx[128];
	char idxToDna[8];


	void PackToBin(const DnaBin& dnaBin_, BitMemoryWriter& metaWriter_, BitMemoryWriter& dnaWriter_,
				   BinaryBinDescriptor& binDesc_, bool cutMinimizer_ = false);

	void UnpackFromBin(const BinaryBinDescriptor& desc_, DnaBin& dnaBin_, DataChunk& dnaChunk_,
					   BitMemoryReader& metaReader_, BitMemoryReader& dnaReader_, uint32 minimizerId_);

	void StoreNextRecord(const DnaRecord& rec_, BitMemoryWriter& binWriter_,
						 BitMemoryWriter& dnaWriter_, const BinPackSettings& settings_);

	bool ReadNextRecord(BitMemoryReader& binReader_, BitMemoryReader& dnaReader_,
						DnaRecord& rec_, const BinPackSettings& settings_);

};


#endif // H_DNAPACKER
