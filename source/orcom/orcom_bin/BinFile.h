/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#ifndef H_BINFILE
#define H_BINFILE

#include "Globals.h"

#include <vector>

#include "FileStream.h"
#include "Params.h"


struct BinFileHeader
{
	static const uint32 ReservedBytes = 8;
	static const uint32 HeaderSize = 4*8 + ReservedBytes;

	uint64 footerOffset;
	uint64 recordsCount;
	uint64 blockCount;
	uint64 footerSize;
	uchar reserved[ReservedBytes];
};


struct BinFileFooter
{
	static const uint32 ParametersSize = sizeof(BinModuleConfig);	// 16 B

	BinModuleConfig params;
	std::vector<uint32> metaBinSizes;
	std::vector<uint32> dnaBinSizes;

	std::vector<uint32> recordsCounts;
	std::vector<uint32> rawDnaSizes;

	void Clear()
	{
		metaBinSizes.clear();
		dnaBinSizes.clear();
		recordsCounts.clear();
		rawDnaSizes.clear();
	}
};


class BinFileWriter
{
public:
	BinFileWriter();
	/*
		Función:
			Funcion que inicialza parametros de la clase BinFileWriter.
	*/
	~BinFileWriter();

	void StartCompress(const std::string& filename_, const BinModuleConfig& params_);
	/*
		Función:
			Se encarga de crear los archivos Bin de salida (.bmeta y .bdna) e inicializar parametros para la compresión
		Entradas: 
			filename_: Nombre del archivo Bin de salida
			params_: Parametros de compresion
	*/

	void WriteNextBlock(const BinaryBinBlock* block_);
	/*
		Función:
			Función que escribe los paquetes de Bins en el archivo final
		Entradas: 
			block_: Paquete de Bins que se escribirá en el archivo final
	*/
	void FinishCompress();
	/*
		Función:
			Funcion que se encarga de terminar el proceso de compresion, cerrando los archivos y liberando los
			espacios de memoria necesarios
	*/

	void GetBinStats(std::vector<uint64>& recordsCounts_);
	/*
		Función:
			Se encarga de recolectar información extra acerca del proceso de compresión
		Entradas: 
			recordsCounts_: Vector donde se guardará la información recolectada
	*/

protected:
	IDataStreamWriter* metaStream;
	IDataStreamWriter* dnaStream;

	BinFileHeader fileHeader;
	BinFileFooter fileFooter;

	uint64 currentBlockId;
	uint64 minimizersCount;

	void WriteFileHeader();
	void WriteFileFooter();
};


class BinFileReader
{
public:
	BinFileReader();
	~BinFileReader();

	void StartDecompress(const std::string& fileName_, BinModuleConfig& params_);

	bool ReadNextBlock(BinaryBinBlock* block_);
	void FinishDecompress();

	uint64 BlockCount() const
	{
		return fileHeader.blockCount;
	}

protected:
	IDataStreamReader* metaStream;
	IDataStreamReader* dnaStream;
	BinFileHeader fileHeader;
	BinFileFooter fileFooter;

	uint64 currentBlockId;
	uint64 minimizersCount;
	uint64 currentNonEmptyBlockId;

	void ReadFileHeader();
	void ReadFileFooter();
};

#endif // H_BINFILE
