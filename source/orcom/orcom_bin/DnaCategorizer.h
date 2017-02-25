/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#ifndef H_DNACATEGORIZER
#define H_DNACATEGORIZER

#include "Globals.h"

#include <vector>
#include <map>

#include "DnaRecord.h"
#include "Collections.h"
#include "Params.h"


class DnaCategorizer
{
public:
	DnaCategorizer(const MinimizerParameters& params_, const CategorizerParameters& catParams_);
	/*
		Funci�n:
			Inicializa variables de la clase DnaCategorizer
		Entradas: 
			params_: Objeto de la clase MinimizerParameters
			catParams_: Objeto de la clase CategorizerParameters
	*/

	void Categorize(std::vector<DnaRecord>& records_, uint64 recordsCount_, DnaBinBlock& bin_);
	/*
		Funci�n:
			Funci�n encargada de hacer el procesamiento de los datos, limpia los Bins, procesa los Reads en busca de minimizadores y los organiza en los Bins 
		Entradas: 
			records_: Vector de objetos DnaRecord, el cual contienen una lista de Reads obtenidos del archivo Fastq
			recordsCount_: Cantidad de Reads a procesar
			bin_: Objeto DnaBinBlock en el cual se almacenar�n los Reads de forma ordenada
	*/

private:
	const MinimizerParameters& params;
	const CategorizerParameters catParams;

	const uint32 maxShortMinimValue;
	const uint32 maxLongMinimValue;
	const uint32 nBinValue;

	char symbolIdxTable[128];		// Almacena el Id de las letras (ACGTN = 01234) en su respectiva posici�n segun su numero en ASCII

	std::vector<uint64> freqTable;

	void DistributeToBins(std::vector<DnaRecord>& records_, uint64 recordsCount_, DnaBinCollection& bins_, DnaBin& nBin_);
	/*
		Funci�n:
			Lee cada Read del vector records_, encuentra su minimizador y guarda el Read en el Bin correspondiente
		Entradas: 
			records_: Vector de objetos DnaRecord, el cual contienen una lista de Reads obtenidos del archivo Fastq
			recordsCount_: Cantidad de Reads a procesar
			bins_: Objeto DnaBinCollection en el cual se almacenar�n los Reads de forma ordenada
			nBin_: 
	*/
	
	void FindMinimizerPositions(DnaBinCollection& bins_);

	uint32 FindMinimizer(DnaRecord& rec_);
	/*
		Funci�n:
			Busca el minimizador y su validez 
		Entradas: 
			rec_: Objeto que contiene el Read obtenido del archivo Fastq
		Salida:
			uint32: minimizer & (maxShortMinimValue - 1);
	*/
	
	
	std::map<uint32, uint16> FindMinimizers(DnaRecord &rec_);
	uint32 ComputeMinimizer(const char* dna_, uint32 mLen_);
	/*
		Funci�n:
			Encuentra el minimizador en el Read y lo convierte en binario con sus respectivos Id
		Entradas: 
			dna_: Read leido, cadena de caracteres que contiene informaci�n genomica
			mLen_: Longitud del minimizador, default=8
		Salida:
			uint32: Minimizador
	*/

	bool IsMinimizerValid(uint32 minim_, uint32 mLen_);
	/*
		Funci�n:
			Fuinci�n que define si un minimizador es valido
		Entradas: 
			minim_: Minimizador hallado en forma binaria
			mLen_: Longitud del mimimizador, default=8
		Salida:
			bool: Variable que define si el minimizador es valido
	*/
};


#endif // H_DNACATEGORIZER
