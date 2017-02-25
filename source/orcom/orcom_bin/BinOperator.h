/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#ifndef H_BINOPERATOR
#define H_BINOPERATOR

#include "Globals.h"
#include "DataPool.h"
#include "DataQueue.h"
#include "BinBlockData.h"
#include "Params.h"

typedef BinaryBinBlock BinaryPart;

// operators for multi threaded processing
//
typedef TDataPool<BinaryPart> BinaryPartsPool;
typedef TDataPool<DataChunk> FastqChunkPool;

typedef TDataQueue<BinaryPart> BinaryPartsQueue;
typedef TDataQueue<DataChunk> FastqChunkQueue;


// TODO: readers and writers can be fully templatized
//
class FastqChunkReader : public IOperator
{
public:
	FastqChunkReader(IFastqStreamReader* partsStream_, FastqChunkQueue* partsQueue_, FastqChunkPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}
	/*
		Función:
			Inicializar variables de la clase FastqChunkReader, las cuales se usan para leer el archivo Fastq
		Entradas: 
			partsStream_: Objeto de la clase IFastqStreamReader del archivo Fastq a leer
			partsQueue_: Objeto de la clase FastqChunkQueue 
			partsPool_: Objeto de la clase FastqChunkPool 
	*/

	void Run();

private:
	IFastqStreamReader* partsStream;
	FastqChunkQueue* partsQueue;
	FastqChunkPool* partsPool;
};


class BinChunkWriter : public IOperator
{
public:
	BinChunkWriter(BinFileWriter* partsStream_, BinaryPartsQueue* partsQueue_, BinaryPartsPool* partsPool_)
		:	partsStream(partsStream_)
		,	partsQueue(partsQueue_)
		,	partsPool(partsPool_)
	{}
	/*
		Función:
			Inicializar variables de la clase BinChunkWriter, las cuales se usan para escribir el archivo Bin
		Entradas: 
			partsStream_: Objeto de la clase BinFileWriter del archivo Bin a escribir
			partsQueue_: Objeto de la clase BinaryPartsQueue 
			partsPool_: Objeto de la clase BinaryPartsPool 
	*/

	void Run();

private:
	BinFileWriter* partsStream;
	BinaryPartsQueue* partsQueue;
	BinaryPartsPool* partsPool;
};


class BinEncoder : public IOperator
{
public:
	BinEncoder(const MinimizerParameters& params_,
			   const CategorizerParameters& catParams_,
			   FastqChunkQueue* fqPartsQueue_, FastqChunkPool* fqPartsPool_,
			   BinaryPartsQueue* binPartsQueue_, BinaryPartsPool* binPartsPool_)
		:	params(params_)
		,	catParams(catParams_)
		,	fqPartsQueue(fqPartsQueue_)
		,	fqPartsPool(fqPartsPool_)
		,	binPartsQueue(binPartsQueue_)
		,	binPartsPool(binPartsPool_)
	{}

protected:
	const MinimizerParameters params;
	const CategorizerParameters catParams;

	FastqChunkQueue* fqPartsQueue;
	FastqChunkPool* fqPartsPool;
	BinaryPartsQueue* binPartsQueue;
	BinaryPartsPool* binPartsPool;

	void Run();
};


#endif // H_BINOPERATOR
