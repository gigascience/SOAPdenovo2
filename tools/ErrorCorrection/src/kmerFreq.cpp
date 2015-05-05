/*
 * kmerFreq.cpp
 *
 *  Created on: 2011-8-10
 *      Author: Aqua
 */

#include"kmerFreq.h"

using namespace std;
using namespace aqua;
using namespace kmerfreq;

int kmerfreq::conKmerSize(17), kmerfreq::spaKmerSize(16), kmerfreq::spaSize(9);
string prefix("output");

static void usage(char ** argv)
{
	cerr
	<<argv[0]<<" "<<argv[1]<<" [PARAMETERS] <List of Read Files>"<<'\n'
	<<"    -k <int>     Consecutive Kmer Size, Default: "<<conKmerSize<<'\n'
	<<"    -K <int>     Space Kmer Solid Part Size, Default: "<<spaKmerSize<<'\n'
	<<"    -S <int>     Space Kmer Space Part Size, Default: "<<spaSize<<'\n'
	<<"    -p <string>  Output prefix, Default: "<<prefix<<'\n'
	<<"    -v           Increase Verbosity, 3 times max"<<'\n'
	<<"    -h           This help"<<'\n';

	exit(EXIT_FAILURE);
}

int Kmerfreq(int argc, char ** argv)
{
	char strBuf[1024] = "";
	int status;
	if(argc < 3)
		usage(argv);
	while((status = getopt(argc-1, argv+1, "vk:K:S:p:h")) !=  -1)
	{
		switch(status)
		{
			case 'k': conKmerSize = atoi(optarg); break;
			case 'K': spaKmerSize = atoi(optarg); break;
			case 'S': spaSize = atoi(optarg); break;
			case 'p': prefix = optarg; break;
			case 'v': ModifyVerbosity(); break;
			case 'h': usage(argv); break;
			default: usage(argv); break;
		}
	}
	string readFNList;
	if(optind < (argc - 1))
		readFNList = ((argv+1)[optind++]);
	else
		usage(argv);

	mvnv(2, "User defined parameters:");
	mvnv(2, "     Consecutive kmer size: %d", conKmerSize);
	mvnv(2, "Space kmer solid-part size: %d", spaKmerSize);
	mvnv(2, "Space kmer space-part size: %d", spaSize);
	mvnv(2, "             Output prefix: %s", prefix.c_str());
	mvnv(2, "                 Verbosity: %d", (int)log2(verbosity));
	mvnv(2, "Kmer Scene:\n%s", PrintKmerScene(conKmerSize, spaKmerSize, spaSize).c_str());

	status = 0;
	if(spaKmerSize % 2 != 0)
	{
		mvnv(1, "Space kmer solid part should be even number");
		++status;
	}
	if( (spaKmerSize + spaSize) % 2 == 0 )
	{
		mvnv(1, "Space kmer length is not even number");
		++status;
	}
	if( ((spaKmerSize / 2) + spaSize) > conKmerSize )
	{
		mvnv(1, "Consecutive kmer cannot cover the space of space kmer");
		++status;
	}
	if(status)
	{
		if(spaKmerSize % 2 != 0)
			--spaKmerSize;
		spaSize = conKmerSize - (spaKmerSize / 2);
		if( (spaKmerSize + spaSize) % 2 == 0 )
			--spaSize;

		mvnv(2, "New parameters:");
		mvnv(2, "     Consecutive kmer size: %d", conKmerSize);
		mvnv(2, "Space kmer solid-part size: %d", spaKmerSize);
		mvnv(2, "Space kmer space-part size: %d", spaSize);
		mvnv(2, "Kmer Scene:\n%s", PrintKmerScene(conKmerSize, spaKmerSize, spaSize).c_str());
	}

	uint64_t * freqArray = new uint64_t[256];
	bzero((void *)freqArray, sizeof(uint64_t) * 256);

	uint8_t * conTable, * spaTable;
	uint64_t kmerCount(0), kmerEffectiveCount(0), kmerTheoryCount(0);


	/*File Pre-opening Area*/
	gzFile conO = gzopen(RegisterFile(prefix + ".freq.gz").c_str(), "wb");		//Consecutive kmerfreq Table
	if(conO == Z_NULL)
		perrdie("Error opening %s", (prefix + ".freq.gz").c_str());

	gzFile spaO = gzopen(RegisterFile(prefix + ".sfreq.gz").c_str(), "wb");		//Space kmerfreq Table
	if(spaO == Z_NULL)
		perrdie("Error opening %s", (prefix + ".sfreq.gz").c_str());

	fstream conStatO;
	OpenFile(conStatO, RegisterFile(prefix + ".freqstat").c_str(), "wc-");		//Consecutive kmerfreq Statistics
	fstream spaStatO;
	OpenFile(spaStatO, RegisterFile(prefix + ".sfreqstat").c_str(), "wc-");		//Space kmerfreq Statistics
	/*File Pre-opening Area End*/

	{
		RecordOneTime();
		mvnv(1, "Consecutive kmer table construction in progress...");
		//Consecutive kmer frequency construction
		conTable = MakeKmerfreqFromFilelist(readFNList, conKmerSize, 0, kmerTheoryCount,\
											kmerCount, kmerEffectiveCount, freqArray);
		PrintDuration();
	}

	{
		RecordOneTime();
		mvnv(1, "Writing consecutive kmer table...")
		uint64_t rtBufSize;
		for(uint64_t i(0); i < kmerTheoryCount; i += IOBUFSIZE)
		{
			rtBufSize = ((kmerTheoryCount - i) > IOBUFSIZE) ? IOBUFSIZE : (kmerTheoryCount - i);
			gzwrite(conO, (void const *)(conTable + i), rtBufSize);
		}
		gzclose(conO);
		PrintDuration();
	}
	munmap(conTable, kmerTheoryCount);

	/*Statistics printing abbreviation*/
#ifndef PRTSTAT
#define PRTSTAT(x) \
	strBuf[0] = '\0';\
	sprintf(strBuf, "Consecutive_Kmer_Size\t%d\n", conKmerSize);\
	sprintf(strBuf, "%sSpace_Kmer_Solid\t%d\n", strBuf, spaKmerSize);\
	sprintf(strBuf, "%sSpace_Kmer_Space\t%d\n", strBuf, spaSize);\
	sprintf(strBuf, "%sKmer_Max\t%lu\n", strBuf, kmerTheoryCount);\
	sprintf(strBuf, "%sKmer_Count\t%lu\n", strBuf, kmerCount);\
	sprintf(strBuf, "%sKmer_Effective\t%lu\n", strBuf, kmerEffectiveCount);\
	x<<strBuf<<'\n';\
	x<<"KmerFreq\tKmerType\tKmerCount\tAccumulate%\n";\
	{\
		uint64_t accumKmerCount(0);\
		double accumKmerPercent = 0.;\
		\
		for(uint i(1); i < 255; ++i)\
		{\
			accumKmerCount += freqArray[i] * i;\
			accumKmerPercent = (double)accumKmerCount / kmerCount;\
			strBuf[0] = '\0';\
			sprintf(strBuf, "%d\t%lu\t%lu\t%f\n", i, freqArray[i], freqArray[i] * i, accumKmerPercent);\
			x << strBuf;\
		}\
		strBuf[0] = '\0';\
		sprintf(strBuf, "254+\t%lu\t%lu\t%f\n", freqArray[255], kmerCount-accumKmerCount, 1.);\
		x << strBuf;\
	}\
	x<<flush;
#endif
	/*Statistics printing abbreviation done*/

	PRTSTAT(conStatO);

	bzero((void *)freqArray, sizeof(uint64_t) * 256);
	kmerTheoryCount = kmerCount = kmerEffectiveCount = 0;

	{
		RecordOneTime();
		mvnv(1, "Space kmer table construction in progress...");
		//Space kmer frequency construction
		spaTable = MakeKmerfreqFromFilelist(readFNList, spaKmerSize, spaSize, kmerTheoryCount,\
											kmerCount, kmerEffectiveCount, freqArray);
		PrintDuration();
	}

	{
		RecordOneTime();
		mvnv(1, "Writing space kmer table...")
		uint64_t rtBufSize;
		for(uint64_t i(0); i < kmerTheoryCount; i += IOBUFSIZE)
		{
			rtBufSize = ((kmerTheoryCount - i) > IOBUFSIZE) ? IOBUFSIZE : (kmerTheoryCount - i);
			gzwrite(spaO, (void const *)(spaTable + i), rtBufSize);
		}
		gzclose(spaO);
		PrintDuration();
	}
	munmap(spaTable, kmerTheoryCount);

	PRTSTAT(spaStatO);

	mvnv(1, "Kmer freqeuency construction done");

	return status = 0;
}
