/*
 * connect.h
 *
 *  Created on: 2011-8-16
 *      Author: Aqua
 */

#ifndef CONNECT_H_AQUA_
#define CONNECT_H_AQUA_

#include <string>
#include <fstream>
#include <map>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <getopt.h>
#include <algorithm>
#include "gzstream.h"

using namespace std;

namespace connect
{
extern double MATCH_CUTOFF, DIF, B_CUTOFF;
extern int Min, Max, Mean, QualityOutput, Nfilter;
extern int tHead, tTail;
extern int qual_base, qual_offset;
}

struct ConnectArgs_s
{
	ConnectArgs_s(): total_pairs(0), connect_pairs(0), unconnect_pairs(0), low_quality_pairs(0), low_match_rate_pairs(0), confuse_pairs(0),\
					 Max_mismatch_num(0), Mismatch_array(NULL), mismatch_num_distr(NULL) {}

	long int total_pairs, connect_pairs, unconnect_pairs,\
			low_quality_pairs, low_match_rate_pairs, confuse_pairs;
	string read_a, read_b, read_head_a, read_head_b,\
			a_id, b_id, a_s, b_s, a_q, b_q, rread_b;

	int Max_mismatch_num;
	int * Mismatch_array;
	map<int, uint64_t> InsertSize_distr;
	uint64_t * mismatch_num_distr;
};

void Set1Read4Connect(ConnectArgs_s &, const string &, const string &, const string &, const string &,\
		const string &, const string &, const string &, const string &);

int alignConnect (ConnectArgs_s &, const string &, const string &, const int = 0);

void ConnectInit(ConnectArgs_s &);

void ConnectClean(ConnectArgs_s &);

void OutputConnectStat(ConnectArgs_s &, fstream &);

int Analyze4Correction (ConnectArgs_s &, int, ogzstream &);

int Connect (int, char **);

#endif /* CONNECT_H_AQUA_ */
