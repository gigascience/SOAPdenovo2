/*
 * correct.h
 *
 *  Created on: 2011-8-13
 *      Author: Aqua
 */

//#define _CORR_DEBUG
#ifdef _CORR_DEBUG
#define mvnvd(x, ...) {fprintf(stderr, ##__VA_ARGS__); fprintf(stderr, "\n");}
#else
#define mvnvd(...) NULL
#endif

#pragma once
#ifndef CORRECT_H_AQUA_
#define CORRECT_H_AQUA_

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <sstream>
#include "kmerFreq.h"
#include "connect.h"
#include "seqKmer.h"
#include "general.h"

const int FQ(1), FA(2), QUALMAX(38);

struct FreqContReg_s {
	FreqContReg_s() : start(0), end(0), status(0) {}
	FreqContReg_s(const int startT, const int endT, const int statusT) :\
				start(startT), end(endT), status(statusT) {}
	void modify(const int, const int, const int);

	int start;
	int end;
	int status;
};

struct TreeNode_s {
	TreeNode_s();
	TreeNode_s(const uint, const uint, const uint, const uint);
	void modify(const uint, const uint, const uint, const uint);

	uint32_t pointer:26, base:2, change:3, same:1;
};

struct CorrGlobalVar_s {
	explicit CorrGlobalVar_s() : conKmerSize(kmerfreq::conKmerSize), spaKmerSize(kmerfreq::spaKmerSize), spaSize(kmerfreq::spaSize),\
			conLowFreqCutoff(10), spaLowFreqCutoff(10), conHighFreqCutoff(highFreqDefault), spaHighFreqCutoff(highFreqDefault), minHighFreqLen(10), maxChangePerRead(2), trimAllow(0), maxNode(0x1000000),\
			minOutputReadLen(50), corrByOvl(0), joinRead(0), toDelete(0), commentOn(0), minInterval(0), conn_qual_base(connect::qual_base),\
			conn_qual_offset(connect::qual_offset), conn_MATCH_CUTOFF(connect::MATCH_CUTOFF), conn_DIF(connect::DIF),\
			conn_Min(connect::Min), conn_Max(connect::Max), conn_Mean(connect::Mean), conKmerTable(0), spaKmerTable(0), conInvalidKmerTable(0), spaInvalidKmerTable(0), \
			conKmerMasker(0), spaKmerHeadMasker(0), spaKmerTailMasker(0), spaHeadOffset(0) {}
	void print(const char *);

	int conKmerSize;
	int spaKmerSize;
	int spaSize;
	int conLowFreqCutoff;
	int spaLowFreqCutoff;
	int conHighFreqCutoff;
	int spaHighFreqCutoff;
	int minHighFreqLen;
	int maxChangePerRead;
	int trimAllow;
	uint maxNode;
	int minOutputReadLen;
	int corrByOvl;
	int joinRead;
	int toDelete;
	int commentOn;
	int minInterval;
	int & conn_qual_base;
	int & conn_qual_offset;
	double & conn_MATCH_CUTOFF;
	double & conn_DIF;
	int & conn_Min;
	int & conn_Max;
	int & conn_Mean;

	uint8_t * conKmerTable;
	uint8_t * spaKmerTable;
	uint8_t * conInvalidKmerTable;
	uint8_t * spaInvalidKmerTable;

	string conKmerTableFN;
	string spaKmerTableFN;
	string readsListFN;

	vector<string> readsFNs;

	uint64_t conKmerMasker;
	uint64_t spaKmerHeadMasker;
	uint64_t spaKmerTailMasker;
	uint spaHeadOffset;
};

struct Read_s {
	explicit Read_s() : fastCorrectCount(0), bbCorrectCount(0),\
						shouldDelete(0), trimLeft(0), trimRight(0) {}

	int fastCorrectCount;
	int bbCorrectCount;
	int shouldDelete;
	int trimLeft;
	int trimRight;

	string head;
	string read;
	string qualHead;
	string qual;

	vector<int> fastCorrectRecord;
	vector<char> fastCorrectCharRec;
	vector<int> bbCorrectRecord;
	vector<char> bbCorrectCharRec;
	string comment;

	string altAllele;
};

int Correct(int, char **);

#endif /* CORRECT_H_AQUA_ */
