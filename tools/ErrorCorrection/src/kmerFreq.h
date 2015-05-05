/*
 * kmerFreq.h
 *
 *  Created on: 2011-8-10
 *      Author: Aqua
 */

#pragma once
#ifndef KMERFREQ_H_AQUA_
#define KMERFREQ_H_AQUA_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <zlib.h>
#include <getopt.h>
#include "seqKmer.h"
#include "gzstream.h"

#define IOBUFSIZE 0x100000

namespace kmerfreq
{
extern int conKmerSize, spaKmerSize, spaSize;
}

int Kmerfreq(int, char **);

#endif /* KMERFREQ_H_AQUA_ */
