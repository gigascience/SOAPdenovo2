/*
 * basicOperation.cpp
 *
 *  Created on: 2011-8-10
 *      Author: Aqua
 */

#include"seqKmer.h"

using namespace std;
using namespace aqua;

const int highFreqDefault = 255;
const int sleepWait = 10;

const char bases[5] = {'A', 'C', 'T', 'G', 'N'};

const char cBases[5] = {'T', 'G', 'A', 'C', 'N'};

const uint8_t bitPos[8] = {128, 64, 32, 16, 8, 4, 2, 1};

const uint64_t bitLeft[4] =
{
	0x0000000000000000llu,
	0x4000000000000000llu,
	0x8000000000000000llu,
	0xc000000000000000llu
};

uint8_t* MapKmerfreq2Mem(const string & kmerFreqFN, const uint kmerSize, uint64_t & kmerAmount,\
						const int lowFreqCutoff, const int highFreqCutoff, uint8_t *& invalidKmerTable)
{
	int status;
	kmerAmount =  ((0x2llu << (kmerSize * 2 - 1)));
	string bitKmerFreqLowFN;
	string bitKmerFreqHighFN;
	uint8_t * freq(NULL);
	invalidKmerTable = NULL;
	{
		stringstream ssl, ssh;
		ssl << kmerFreqFN << ".bin_low_" << kmerSize << "_" << lowFreqCutoff;
		ssh << kmerFreqFN << ".bin_high_" << kmerSize << "_" << highFreqCutoff;
		bitKmerFreqLowFN = ssl.str();
		bitKmerFreqHighFN = ssh.str();
	}
	if(!FileExist(bitKmerFreqLowFN.c_str(), false) || ((highFreqCutoff > 0 && highFreqCutoff < highFreqDefault) && !FileExist(bitKmerFreqHighFN.c_str(),false)))
	{
		mvnv(1, "Transforming kmer frequency table %s to authenticated kmer table %s", kmerFreqFN.c_str(), bitKmerFreqLowFN.c_str());
		int Ol = open(RegisterFile(bitKmerFreqLowFN).c_str(), O_CREAT | O_WRONLY, 0755);
		if(Ol == -1)
			perrdie("Error opening %s", bitKmerFreqLowFN.c_str());
		int Oh = open(RegisterFile(bitKmerFreqHighFN).c_str(), O_CREAT | O_WRONLY, 0755);
		if(Oh == -1)
			perrdie("Error opening %s", bitKmerFreqHighFN.c_str());
		RecordOneTime();
		uint64_t theoryTotal(0), count(0), effectiveTotal(0); double lowFreqPercent(0.);
		freq = ReadKmerfreq2Mem(kmerFreqFN.c_str(), kmerSize, theoryTotal, count, effectiveTotal, lowFreqCutoff, highFreqCutoff, invalidKmerTable, lowFreqPercent);
		mvnv(1, " Table theory maximum: %lu", theoryTotal);
		mvnv(1, "     Total kmer count: %lu", count);
		mvnv(1, "      Total kmer hits: %lu", effectiveTotal);
		mvnv(1, "Low Freq Kmer Percent: %f", lowFreqPercent*100);
		mvnv(1, "Writing to disk...");
		{
			const size_t IOBUFSIZE(1000000);
			uint64_t rtBufSize;
			const size_t kmerTheoryCount = sizeof(uint8_t) * (kmerAmount/8+1);
			for(uint64_t i(0); i < kmerTheoryCount; i += IOBUFSIZE)
			{
				rtBufSize = ((kmerTheoryCount - i) > IOBUFSIZE) ? IOBUFSIZE : (kmerTheoryCount - i);
				status = write(Ol, (const void *)(freq + i), rtBufSize);
			}
			close(Ol);
			if(highFreqCutoff > 0 && highFreqCutoff < highFreqDefault)
			{
				for(uint64_t i(0); i < kmerTheoryCount; i += IOBUFSIZE)
				{
					rtBufSize = ((kmerTheoryCount - i) > IOBUFSIZE) ? IOBUFSIZE : (kmerTheoryCount - i);
					status = write(Oh, (const void *)(invalidKmerTable + i), rtBufSize);
				}
			}
			close(Oh);
		}
		PrintDuration();
		munmap(freq, sizeof(uint8_t) * (kmerAmount/8+1));
		if(highFreqCutoff > 0 && highFreqCutoff < highFreqDefault)
			munmap(invalidKmerTable, sizeof(uint8_t) * (kmerAmount/8+1));
		freq = invalidKmerTable = NULL;
	}

	int fdes;
	struct stat sb;
	fdes = open(bitKmerFreqLowFN.c_str(), O_RDONLY);
	if(fdes < 0)
		die("Error opening: %s", bitKmerFreqLowFN.c_str());
	if(fstat (fdes, &sb) == -1)
		die("Error stating: %s", bitKmerFreqLowFN.c_str());
	if(!S_ISREG (sb.st_mode))
		die("Non regular file: %s", bitKmerFreqLowFN.c_str());
	if(sb.st_size != (off_t)(sizeof(uint8_t) * (kmerAmount/8+1)))
	{
		off_t previousSize = sb.st_size;
		sleep(sleepWait);
		if(fstat (fdes, &sb) == -1)
			die("Error stating: %s", bitKmerFreqLowFN.c_str());
		if(sb.st_size - previousSize == 0)
		{
			die("Error loading memory mapping file, please delete %s and %s and run again.", bitKmerFreqLowFN.c_str(), bitKmerFreqHighFN.c_str());
		}
		else
		{
			mvnv(1, "Waiting another process to finish memory mapping file construction.");
			do
			{
				sleep(sleepWait);
				if(fstat (fdes, &sb) == -1)
					die("Error stating: %s", bitKmerFreqLowFN.c_str());
			} while(sb.st_size != (off_t)(sizeof(uint8_t) * (kmerAmount/8+1)));
		}
	}
	freq = (uint8_t *)mmap(0, sizeof(uint8_t) * (kmerAmount/8+1), PROT_READ, MAP_SHARED, fdes, 0);
	if((freq == NULL) || (freq == MAP_FAILED))
		die("Failed to map %s to memory", bitKmerFreqLowFN.c_str());
	if(close(fdes) == -1)
		die("Failed closing file: %s", bitKmerFreqLowFN.c_str());

	if(highFreqCutoff > 0 && highFreqCutoff < highFreqDefault)
	{
		int fdes2;
		struct stat sb2;
		fdes2 = open(bitKmerFreqHighFN.c_str(), O_RDONLY);
		if(fdes2 < 0)
			die("Error opening: %s", bitKmerFreqHighFN.c_str());
		if(fstat (fdes2, &sb2) == -1)
			die("Error stating: %s", bitKmerFreqHighFN.c_str());
		if(!S_ISREG (sb2.st_mode))
			die("Non regular file: %s", bitKmerFreqHighFN.c_str());
		if(sb2.st_size != (off_t)(sizeof(uint8_t) * (kmerAmount/8+1)))
		{
			off_t previousSize = sb2.st_size;
			sleep(sleepWait);
			if(fstat (fdes2, &sb2) == -1)
				die("Error stating: %s", bitKmerFreqHighFN.c_str());
			if(sb2.st_size - previousSize == 0)
			{
				die("Error loading memory mapping file, please delete %s and %s and run again.", bitKmerFreqLowFN.c_str(), bitKmerFreqHighFN.c_str());
			}
			else
			{
				mvnv(1, "Waiting another process to finish memory mapping file construction.");
				do
				{
					sleep(sleepWait);
					if(fstat (fdes2, &sb2) == -1)
						die("Error stating: %s", bitKmerFreqHighFN.c_str());
				} while(sb2.st_size != (off_t)(sizeof(uint8_t) * (kmerAmount/8+1)));
			}
			}
		invalidKmerTable = (uint8_t *)mmap(0, sizeof(uint8_t) * (kmerAmount/8+1), PROT_READ, MAP_SHARED, fdes2, 0);
		if((invalidKmerTable == NULL) || (invalidKmerTable == MAP_FAILED))
			die("Failed to map %s to memory", bitKmerFreqHighFN.c_str());
		if(close(fdes2) == -1)
			die("Failed closing file: %s", bitKmerFreqHighFN.c_str());
	}

	return freq;
}

uint8_t* ReadKmerfreq2Mem(const string & kmerFreqFN, const uint kmerSize, uint64_t & kmerAmount, \
						  uint64_t & kmerCount, uint64_t & kmerEffectiveCount, \
						  int lowFreqCutoff, int highFreqCutoff, uint8_t *& invalidKmerTable, double & lowFreqPercent)
{
	//kmerAmount =  ((0x2llu << (kmerSize * 2 - 1)) - 1);
	kmerAmount =  ((0x2llu << (kmerSize * 2 - 1)));

	mvnv(3, "freqTable size: %lu", kmerAmount);
	uint8_t * freq = (uint8_t*)mmap(NULL, sizeof(uint8_t) * (kmerAmount/8+1), PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE, -1, 0);
	if(freq == MAP_FAILED)
		perrdie("\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:");
	//uint8_t * freq = new uint8_t[kmerAmount / 8 + 1];
	//bzero((void *)freq, sizeof(uint8_t) * (kmerAmount / 8 + 1));
	madvise((void *)freq, (kmerAmount/8+1), MADV_RANDOM);

	if(highFreqCutoff > 0 && highFreqCutoff < highFreqDefault)
	{
		invalidKmerTable = (uint8_t*)mmap(NULL, sizeof(uint8_t) * (kmerAmount/8+1), PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE, -1, 0);
		if(invalidKmerTable == MAP_FAILED)
			perrdie("\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:");
		//bzero((void *)invalidKmerTable, sizeof(uint8_t) * (kmerAmount / 8 + 1));
		madvise((void *)invalidKmerTable, (kmerAmount/8+1), MADV_RANDOM);
	}
	else
		invalidKmerTable = NULL;

	//Read frequencies from binary gz file
	uint64_t bufSize = 1e6;
	uint8_t *buf = new uint8_t[bufSize];
	gzFile inFile = gzopen(kmerFreqFN.c_str(), "rb");

	uint64_t lowFreqNum(0);
	uint64_t countNum(0);
	uint64_t i(0);
	uint64_t idx(0), rc_idx(0);
	do
	{
		countNum = gzread(inFile, (void *)buf, bufSize);
		for (uint64_t j(0); j < countNum; ++j)
		{
			idx = i+j;

			if (buf[j] > 0)
			{
				kmerCount += buf[j];
				++kmerEffectiveCount;

				if(buf[j] > lowFreqCutoff)
				{
					freq[idx/8] |= bitPos[idx%8];
					rc_idx = getRevCompKbit(idx, kmerSize);
					freq[rc_idx/8] |= bitPos[rc_idx%8];
					if(invalidKmerTable != NULL)
					{
						if(buf[j] > highFreqCutoff)
						{
							invalidKmerTable[idx/8] |= bitPos[idx%8];
							invalidKmerTable[rc_idx/8] |= bitPos[rc_idx%8];
						}
					}
				}
				else
				{
					++lowFreqNum;
				}
			}
		}
		i += countNum;
	}
	while (countNum == bufSize);
	mvnv(2, "Imported array length: %lu/%lu", i, kmerAmount);

	gzclose(inFile);

	lowFreqPercent = (double)lowFreqNum / (double)kmerEffectiveCount;

	return freq;
}

uint8_t * MakeKmerfreqFromFilelist(const string & seqListFN, const uint kmerSize, const uint spaceSize, \
									 uint64_t & kmerAmount, uint64_t & kmerCount, uint64_t & kmerEffectiveCount, \
									 uint64_t * freqAry)
{
	mvnv(1, "Start");
	
	uint64_t kmerMasker((0x2llu << (kmerSize * 2 - 1)) - 1);
	uint64_t headMasker, tailMasker;
	uint headOffset;
	if(spaceSize)
	{
		headOffset = (kmerSize / 2);
		tailMasker = ((0x2llu << ((kmerSize/2) * 2 - 1)) - 1 - 3);
		headMasker = (((0x2llu << ((kmerSize - (kmerSize/2)) * 2 - 1)) - 1 - 3) << (headOffset * 2));
	}
	kmerAmount =  ((0x2llu << (kmerSize * 2 - 1)));
	kmerCount = 0;
	kmerEffectiveCount = 0;

	mvnv(3, "freqTable size: %lu", kmerAmount);
	uint8_t * freq = (uint8_t*)mmap(NULL, sizeof(uint8_t) * kmerAmount, PROT_READ|PROT_WRITE|PROT_EXEC, MAP_ANON|MAP_PRIVATE|MAP_NORESERVE, -1, 0);
	if(freq == MAP_FAILED)
		perrdie("\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:");
	//uint8_t * freq = new uint8_t[kmerAmount];
	bzero((void *)freq, sizeof(uint8_t) * kmerAmount);
	madvise((void *)freq, kmerAmount, MADV_RANDOM);
	
	mvnv(1, "Memory assigned");
  
	VSTR seqFNs;
	ImportToVector(seqListFN.c_str(), seqFNs, '\n', true);

	int spaceKmerSize, headPartSize, endPartSize, endPartStart;
	string kmerRaw, kmerHead, kmerEnd;
	uint64_t kbit, rckbit;
	spaceKmerSize = kmerSize + spaceSize;
	endPartSize = kmerSize / 2;
	headPartSize = kmerSize - endPartSize;
	endPartStart = headPartSize + spaceSize;

	for (uint i(0); i < seqFNs.size(); ++i)
	{
		mvnv(1, "Importing file %d: %s", i, seqFNs[i].c_str());

		igzstream I;
		OpenFileGZ(I, seqFNs[i].c_str(), "rc" );

		string read, readHead;
		while ( getline( I, readHead, '\n' ) )
		{
			if (readHead[0] == '>')
			{
				getline( I, read );
			}
			else if(readHead[0] == '@')
			{
				getline(I, read);
				I.ignore(100000, '\n');
				I.ignore(100000, '\n');
			}
			else
			{
				continue;
			}

			if (!PrepRead(read) ) { continue; }

			if(read.size() < kmerSize + spaceSize)
				continue;

			const uint kLimit(read.size() - kmerSize - spaceSize);
			for (uint k(0); k <= kLimit; ++k)
			{
				if(spaceSize)
				{
					if(k)
					{
						kbit = GetNextKbitSpace(kbit, headMasker, tailMasker, headOffset, read[k+headPartSize-1], read[k+kmerSize+spaceSize-1]);
					}
					else
					{
						kmerRaw = read.substr(k, spaceKmerSize);
						kmerHead = kmerRaw.substr(0, headPartSize);
						kmerEnd =  kmerRaw.substr(endPartStart,endPartSize);
						kbit = Seq2Bit(kmerHead + kmerEnd);
					}
				}
				else
				{
					if(k)
						kbit = GetNextKbit(kbit, kmerMasker, read[k+kmerSize-1]);
					else
						kbit = Seq2Bit(read.substr(k, kmerSize));
					mvnv(4, "kbit: %lu", kbit);
				}

				if (freq[kbit] < 255)
					++freq[kbit];
				
				++kmerCount;
			}
		}
	}

	mvnv(1, "Kmer imported");

	int addFreq;
	for (kbit = 0; kbit < kmerAmount; ++kbit)
	{
		if(!freq[kbit])
			continue;
		++kmerEffectiveCount;
		addFreq = freq[kbit];
		rckbit = getRevCompKbit(kbit, kmerSize);

		if(kbit < rckbit && freq[rckbit] > 0)
		{
			addFreq += freq[rckbit];
			if(addFreq > 255)
				addFreq = 255;
			freq[kbit] = addFreq;
			freq[rckbit] = 0;
		}

		++freqAry[addFreq];
	}

	mvnv(1, "Done");

	return freq;
}

string PrintKmerScene(const int ck, const int sksolid, const int skspace)
{
	string scene;
	scene += string(32, '+'); scene += "\n";
	scene += "Consecutive Kmer\n";
	scene += (string(32-ck, ' ') + string(ck, '|') + "\n");
	scene += (string(32-sksolid-skspace, ' ') + string((sksolid+1)/2, '|') + string(skspace, '-') + string(sksolid/2, '|') + "\n");
	scene += "Space Kmer\n";
	scene += string(32, '+'); scene += "\n";

	return scene;
}
