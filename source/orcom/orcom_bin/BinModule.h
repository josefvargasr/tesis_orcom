/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#ifndef H_BINMODUL
#define H_BINMODULE

#include "Globals.h"

#include <string>
#include <vector>

#include "Params.h"


class BinModule
{
public:
	void Fastq2Bin(const std::vector<std::string>& inFastqFiles_, const std::string& outBinFile_,
				   uint32 threadNum_ = 1, bool compressedInput_ = false, bool verboseMode_ = false);
	/*
		Función:
			Función principal del proceso de compresion, permite leer archivos fastq y escribir 
			archivos bin con sus respctivos datos.
		Entradas: 
			inFastqFiles_: Nombre del archivo fastq a leer
			outBinFile_: Nombre del archivo bin a escribir
			threadNum_: Numero de threads a utilizar
			compressedInput_: Define si los datos a leer están comprimidos
			verboseMode_: Define si se da información detallada del proceso
		Salida:
			int: Define si se ejecuta correctamente la función
	*/
	void Bin2Dna(const std::string& inBinFile_, const std::string& outDnaFile_);

	void SetModuleConfig(const BinModuleConfig& config_)
	{
		config = config_;
	}

private:
	BinModuleConfig config;
};


#endif // H_BINMODULE
