/*
  This file is a part of ORCOM software distributed under GNU GPL 2 licence.
  Homepage:	http://sun.aei.polsl.pl/orcom
  Github:	http://github.com/lrog/orcom

  Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski
*/

#include "Globals.h"

#include <iostream>
#include <string.h>
#include <sys/time.h>

#include "main.h"
#include "BinModule.h"
#include "Utils.h"
#include "Thread.h"


uint32 InputArguments::AvailableCoresNumber = mt::thread::hardware_concurrency();
uint32 InputArguments::DefaultThreadNumber = MIN(8, InputArguments::AvailableCoresNumber);

//extern clock_t end1, start2;
//clock_t start, end, start1, end2;
extern int diff_1, diff_2, diff_loop, diff_bin, diff_bin2, diff_f33, diff__f1;

struct timeval  start, end, start1, end2;
extern struct timeval end1, start2;

int main(int argc_, const char* argv_[])
{	
	gettimeofday(&start, NULL);
	gettimeofday(&start1, NULL);

	diff_bin=0;
	diff_bin2=0;
	diff_f33=0;
	if (argc_ < 1 + 3 || (argv_[1][0] != 'e' && argv_[1][0] != 'd'))
	{
		usage();
		return -1;
	}

	

	InputArguments args;
	if (!parse_arguments(argc_, argv_, args))
		return -1;
	
	int ret = 0;
	if (args.mode == InputArguments::EncodeMode){
		int ret = fastq2bin(args);

		gettimeofday(&end, NULL);
		gettimeofday(&end2, NULL);

		

		printf(" ====== main Time: %ds ======\n", end.tv_sec - start.tv_sec);
		printf(" ====== main Time: %dm%ds ======\n\n", (end.tv_sec - start.tv_sec)/60, (end.tv_sec - start.tv_sec)%60);

		int diff1 = (end1.tv_sec - start1.tv_sec);
		int diff2 = (end2.tv_sec - start2.tv_sec);

		//printf("diff1: %ds, diff2: %ds, diff_1: %ds, diff_2: %ds, diff_bin: %fs, diff_bin2: %ds, diff_f33: %ds, diff__f1: %ds\n", diff1, diff2, diff_1, diff_2, diff_bin, diff_bin2, diff_f33, diff__f1);

		printf(" ====== Total Time by GPU: %ds \n", diff_f33);
		printf(" ====== Total Time Without GPU: %ds \n", diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1);
		printf(" ====== Total Time With GPU: %ds \n\n", diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1 + diff_f33);

		printf(" ====== Total Time by GPU: %dm%ds \n", (diff_f33)/60, (diff_f33)%60);
		printf(" ====== Total Time Without GPU: %dm%ds  \n", (diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1)/60, (diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1)%60);
		printf(" ====== Total Time With GPU: %dm%ds  \n\n", (diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1 + diff_f33)/60, (diff1 + diff2 + diff_1 + diff_2 + diff_bin + diff_bin2 + diff__f1 + diff_f33)%60);

		return ret;
	}
	return bin2dna(args);

}


void usage()
{
	std::cerr << "Overlapping Reads COmpression with Minimizers\n";
	std::cerr << "orcom_bin - DNA records binning tool\n";
	std::cerr << "Version: " << APP_VERSION << '\n';
	std::cerr << "Authors: Sebastian Deorowicz, Szymon Grabowski and Lucas Roguski\n\n";

	std::cerr << "usage:\n\torcom_bin <e|d> [options]\n";
	std::cerr << "options:\n";

	std::cerr << "\t-i<file>\t: input file" << '\n';
	std::cerr << "\t-f\"<f1> <f2> ... <fn>\": input file list" << '\n';
	std::cerr << "\t-g\t\t: input compressed in .gz format\n";
	std::cerr << "\t-o<f>\t\t: output files prefix" << '\n';

	std::cerr << "\t-p<n>\t\t: signature length, default: " << MinimizerParameters::DefaultSignatureLen << '\n';
	std::cerr << "\t-s<n>\t\t: skip-zone length, default: " << MinimizerParameters::DefaultskipZoneLen << '\n';
	std::cerr << "\t-b<n>\t\t: FASTQ input buffer size (in MB), default: " << (BinModuleConfig::DefaultFastqBlockSize >> 20) << '\n';
	std::cerr << "\t-t<n>\t\t: worker threads number, default: " << InputArguments::DefaultThreadNumber << '\n';
	std::cerr << "\t-v\t\t: verbose mode, default: false\n";

#if (DEV_TWEAK_MODE)
	std::cerr << "\t-l<n>\t\t: signature suffix len, default: " << MinimizerParameters::DefaultSignatureSuffixLen << '\n';
	std::cerr << "\t-m<n>\t\t: mimimum block bin size, default: " << CategorizerParameters::DefaultMinimumPartialBinSize << '\n';
	std::cerr << "\t-r\t\t: disable reverse-compliment strategy";
#endif

}


int fastq2bin(const InputArguments& args_)
{
	try
	{
		BinModule module;

		module.SetModuleConfig(args_.config);
		module.Fastq2Bin(args_.inputFiles, args_.outputFile, args_.threadsNum, args_.compressedInput, args_.verboseMode);
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	return 0;
}


int bin2dna(const InputArguments& args_)
{
	try
	{
		BinModule module;
		module.Bin2Dna(args_.inputFiles[0], args_.outputFile);
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	return 0;
}


bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_)
{
	outArgs_.mode = (argv_[1][0] == 'e') ? InputArguments::EncodeMode : InputArguments::DecodeMode;

	MinimizerParameters& pars = outArgs_.config.minimizer;

	// parse params
	//
	for (int i = 2; i < argc_; ++i)
	{
		const char* param = argv_[i];
		if (param[0] != '-')
			continue;

		int pval = -1;
		int len = strlen(param);
		if (len > 2 && len < 10)
			pval = to_num((const uchar*)param + 2, len - 2);
		const char* str = param + 2;
		const uint32 slen = len - 2;

		switch (param[1])
		{
			case 'i':	outArgs_.inputFiles.push_back(std::string(str, str + slen));	break;
			case 'o':	outArgs_.outputFile.assign(str, str + slen);					break;
			case 'p':	pars.signatureLen = pval;										break;
			case 's':	pars.skipZoneLen = pval;										break;
			case 'b':	outArgs_.config.fastqBlockSize = (uint64)pval << 20;			break;
			case 'g':	outArgs_.compressedInput = true;								break;
			case 't':	outArgs_.threadsNum = pval;										break;
			case 'v':	outArgs_.verboseMode = true;									break;
			case 'f':
			{
				int beg = 2;
				for (int i = 2; i < len-1; ++i)
				{
					if (param[i] == ' ' || param[i] == '\n')
					{
						outArgs_.inputFiles.push_back(std::string(param + beg, i - beg));
						beg = i + 1;
					}
				}
				outArgs_.inputFiles.push_back(std::string(param + beg, len - beg));

				break;
			}

#if (DEV_TWEAK_MODE)
			case 'l':	pars.signatureSuffixLen = pval;									break;
			case 'r':	pars.tryReverseCompliment = false;								break;
			case 'm':	outArgs_.config.catParams.minBlockBinSize = pval;				break;
#endif
		}
	}

	if (pars.signatureSuffixLen != pars.signatureLen && !DEV_TWEAK_MODE)
		pars.signatureSuffixLen = pars.signatureLen;


	// check params
	//
	if (outArgs_.inputFiles.size() == 0)
	{
		std::cerr << "Error: no input file(s) specified\n";
		return false;
	}

	if (outArgs_.outputFile.length() == 0)
	{
		std::cerr << "Error: no output file specified\n";
		return false;
	}

	if (outArgs_.threadsNum == 0 || outArgs_.threadsNum > 64)
	{
		std::cerr << "Error: invalid number of threads specified\n";
		return false;
	}

	return true;
}
