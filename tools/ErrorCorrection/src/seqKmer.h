/*
 * basicOperation.h
 *
 *  Created on: 2011-8-10
 *      Author: Aqua
 */

#pragma once
#ifndef BASICOPERATION_H_AQUA_
#define BASICOPERATION_H_AQUA_

#include <string>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "general.h"

#define base2int(base)        (char)(((base)&0x06)>>1)
#define int2base(seq)         "ACTG"[seq]              //0123
#define int2compbase(seq)     "TGAC"[seq]              //2301
#define int_comp(seq)         (char)(seq^0x02)         //(char)((0x4E>>((seq)<<1))&0x03)
#define base2compbase(base)   (char)int2compbase((int)base2int(base))
//0123<-->2301
/*
ASCII:
A  00    1000001
C  01    1000011
T  10    1010100
G  11    1000111
			 ||
   0x06  0000110
*/

extern const int highFreqDefault;
extern const char bases[5];
extern const char compBases[5];
extern const uint8_t bitPos[8];
extern const uint64_t bitLeft[4];

static bool validBaseArray[] =
{
	/*      0, 1, 2, 3, 4,  5, 6, 7, 8, 9*/
	/*00*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,
	/*10*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,
	/*20*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,
	/*30*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,
	/*40*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,

	/*50*/	0, 0, 0, 0, 0,	0, 0, 0, 0, 0,
	/*60*/	0, 0, 0, 0, 0,	1, 0, 1, 0, 0,
	/*70*/	0, 1, 0, 0, 0,	0, 0, 0, 1, 0,
	/*80*/	0, 0, 0, 0, 1,	0, 0, 0, 0, 0,
	/*90*/	0, 0, 0, 0, 0,	0, 0, 1, 0, 1,

	/*100*/	0, 0, 0, 1, 0,	0, 0, 0, 0, 0,
	/*110*/	1, 0, 0, 0, 0,	0, 1, 0, 0, 0,
	/*120*/	0, 0, 0, 0, 0,	0, 0, 0
};

/* Common Functions*/
static inline int PrepRead (string & read, const unsigned int trimHead = 0, const unsigned int trimTail = 0, const uint nLimit = 0)
{
	unsigned int n(0);
	const unsigned int size(read.size());
	if(trimHead||trimTail)
		read = read.substr(trimHead, read.size() - trimTail - trimHead);
	for(unsigned int i(0); i < size; ++i)
	{
		if(!validBaseArray[(int)read[i]])
			read[i] = 'N';

		read[i] = toupper(read[i]);

		if(read[i] == 'N')
			if(++n > nLimit)
				return 0;
	}

	return (int)read.size();
}

static inline void RevCompl ( const string & seq, string & nseq )
{
	nseq.resize(seq.size());

	for ( int i(seq.size() - 1), j(0); i >= 0; --i, ++j )
		nseq[j] = base2compbase(seq[i]);
}

static inline uint64_t GetNextKbit ( const uint64_t currentKbit, const uint64_t mask, const char nextBase )
{
	return ( ( ( currentKbit << 2 ) & mask ) | base2int(nextBase) );
}

static inline uint64_t GetNextKbitSpace ( const uint64_t currentKbit, const uint64_t headMask, const uint64_t tailMask,\
									const uint headOffset, const char headNextBase, const char tailNextBase)
{
//	uint64_t head = (((currentKbit << 2) & headMask) | (base2int(headNextBase) << (headOffset * 2)));
//	uint64_t tail = (((currentKbit << 2) & tailMask) | (base2int(tailNextBase)));
//	fprintf(stderr, "headMask/tailMask/headOffset/head/tail: %llX/%llX/%llX/%llX/%llX\n", headMask, tailMask, headOffset, head, tail);

	return ((((currentKbit << 2) & headMask) | (((uint64_t)base2int(headNextBase)) << (headOffset * 2))) | (((currentKbit << 2) & tailMask) | (base2int(tailNextBase))));
}

static inline uint64_t Seq2Bit ( const string & seq )
{
	uint64_t bit(0);
	for ( size_t i(0); i < seq.size(); ++i )
		bit = (( bit << 2 ) | base2int(seq[i]));

	return bit;
}

static inline void Bit2Seq ( const uint64_t bit, const int kmerSize, string & seq )
{
	seq.resize(kmerSize);

	for ( int i(0); i < kmerSize; ++i )
		seq[i] = int2base(( bit >> ( kmerSize - 1 - i ) * 2 ) & 0x3);
}

static inline uint64_t getRevCompKbit ( uint64_t kbit, const uint ksize )
{
	kbit ^= 0xAAAAAAAAAAAAAAAALLU;
	kbit = ( ( kbit & 0x3333333333333333LLU ) <<  2 ) | ( ( kbit & 0xCCCCCCCCCCCCCCCCLLU ) >>  2 );
	kbit = ( ( kbit & 0x0F0F0F0F0F0F0F0FLLU ) <<  4 ) | ( ( kbit & 0xF0F0F0F0F0F0F0F0LLU ) >>  4 );
	kbit = ( ( kbit & 0x00FF00FF00FF00FFLLU ) <<  8 ) | ( ( kbit & 0xFF00FF00FF00FF00LLU ) >>  8 );
	kbit = ( ( kbit & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( kbit & 0xFFFF0000FFFF0000LLU ) >> 16 );
	kbit = ( ( kbit & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( kbit & 0xFFFFFFFF00000000LLU ) >> 32 );
	return (kbit >> ( 64 - ( ksize << 1 ) ));
}

static inline int getFreq ( const uint8_t * freq, uint64_t idx )
{
	return ( ( freq[idx / 8] >> ( 7 - idx % 8 ) ) & 0x1u);
}

string PrintKmerScene(const int ck, const int sksolid, const int skspace);

/* Common Functions Done*/

//Kmer Construction
uint8_t * MakeKmerfreqFromFilelist(const string & seqListFN, const uint kmerSize, const uint spaceSize, \
								   uint64_t & kmerAmount, uint64_t & kmerCount, uint64_t & kmerEffectiveCount, \
								   uint64_t * freqAry);

//Read Kmer frequency from binary gzipped file into memory, store frequency value in 1 bit in memory, range 0-1
uint8_t* MapKmerfreq2Mem(const string & kmerFreqFN, const uint kmerSize, uint64_t & kmerAmount,\
						const int lowFreqCutoff, const int highFreqCutoff, uint8_t *& invalidKmerTable);

uint8_t* ReadKmerfreq2Mem(const string & kmerFreqFN, const uint kmerSize, uint64_t & kmerAmount, \
						  uint64_t & kmerCount, uint64_t & kmerEffectiveCount, \
						  const int lowFreqCutoff, const int highFreqCutoff, uint8_t *& invalidKmerTable, double & lowFreqPercent);

#endif /* BASICOPERATION_H_AQUA_ */
