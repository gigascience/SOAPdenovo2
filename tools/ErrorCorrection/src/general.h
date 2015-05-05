/*
 *      Filename: general.h
 *
 *
 *      Description:
 *         Basic functions
 *
 *  	Created on: Feb 8, 2010
 *      Author: Ruibang Luo, BGI
 *
 *     	History:
 *         1.
 */

#pragma once
#ifndef GENERAL_H_AQUA_
#define GENERAL_H_AQUA_

#ifndef __cplusplus
#error Please use C++ complier!
#endif

#include<fstream>
#include"gzstream.h"
#include<string>
#include<vector>
#include<boost/lexical_cast.hpp>
#include<boost/utility.hpp>
#include<boost/static_assert.hpp>
#include<stdio.h>
#include<unistd.h>
#include<ctime>
#include<deque>

//#define _ENABLE_BOOST_REGEX
#ifdef _ENABLE_BOOST_REGEX
#include<boost/regex.hpp>
#endif

using namespace std;
using namespace boost;

//Useful Variables*************************************************************
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define FN_SIZE 2048
//*****************************************************************************

//Types************************************************************************
typedef vector<string> VSTR;
typedef vector<VSTR>   VVSTR;

typedef unsigned int       uint;
typedef unsigned char      uchar;
typedef unsigned short     ushort;
typedef unsigned long      ulong;
typedef unsigned long long ullong;

typedef unsigned int       ubyte4;
typedef unsigned char      ubyte1;
typedef unsigned short     ubyte2;
typedef unsigned long long ubyte8;

template <ubyte4 i, ubyte1 j, ubyte2 k, ubyte8 l> void CheckSize() {
  BOOST_STATIC_ASSERT(sizeof(i) == 4);
  BOOST_STATIC_ASSERT(sizeof(j) == 1);
  BOOST_STATIC_ASSERT(sizeof(k) == 2);
  BOOST_STATIC_ASSERT(sizeof(l) == 8);
}

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned int       DWORD;

typedef unsigned char      u8_t;
typedef unsigned short     u16_t;
typedef unsigned int       u32_t;
typedef unsigned long long u64_t;

typedef char * chptr;

//*****************************************************************************

//Debugging********************************************************************
namespace aqua
{
	//Verbose system
	//Verbosity is separated into 4 levels: 0, 1, 2, 3
	#define VERBOSITY_BOTTOM 0
	#define VERBOSITY_TOP 4
	extern int verbosity;
	int ModifyVerbosity(const int = 0);
	#define verboseBufSize 16384
	extern char verboseStr[verboseBufSize];
}
#define ModVerboseStrAndVerbose(level, ...) \
			{\
				if(aqua::verbosity >> level)\
				{\
					snprintf(aqua::verboseStr, verboseBufSize, ##__VA_ARGS__);\
					cerr<<'['<<__FUNCTION__<<"]: "<<aqua::verboseStr<<endl;\
				}\
			}
#define mvnv(level, ...) ModVerboseStrAndVerbose(level, ##__VA_ARGS__)
#define die(...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			cerr<<"Program terminated."<<endl;\
			exit(EXIT_FAILURE);\
		}
#define sigdie(sig, ...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			cerr<<"Program terminated."<<endl;\
			exit(sig);\
		}
#define perrdie(...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			perror("");\
			cerr<<"Program terminated."<<endl;\
			exit(EXIT_FAILURE);\
		}
#define mk \
{\
	fprintf(stderr, "DBG Marker @ %s:%d\n", __FUNCTION__, __LINE__);\
}
/* Please define _EXAM_ASSERT_TEST_ to enable ASSERTMENT */
#ifdef _VERY_GOSSIP_

	#define GOSSIP(message)\
			{\
				cerr<<(message)<<endl;\
			}

#else

	#define GOSSIP(message)

#endif /* _VERY_GOSSIP_ */

#ifdef _EXAM_ASSERT_TEST_

	#define DEBUG_INFO (string("[") + __FILE__ + "]" + "[" + __FUNCTION__ + "]" + "[" + lexical_cast<string>(__LINE__) + "]")

	void exam_assert(const char *, const char *, uint);
	void exam_assert(const char *, const char *, uint, const char *);
	#define ASSERT(condition, ...) \
			{\
				if(!(condition))\
					exam_assert(__FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);\
			}

	#define DBG(message)\
			{\
				cerr<<(message)<<endl;\
			}

	#define NR(message)\
			{\
				cerr<<DEBUG_INFO<<" Will never reach here. "<<(message)<<endl;\
			}

#else

	#define ASSERT(condition, ...)
	#define DEBUG_INFO
	#define DBG(message)
	#define NR

#endif /* _EXAM_ASSERT_TEST_ */
//IO operators and manipulators************************************************
namespace aqua
{
ostream & line(ostream &);
ostream & endt(ostream &);

template<typename T>
void Cin(string warning, T & data)
{
	cerr<<warning<<": ";
	cin>>data;
}

}
//*****************************************************************************

//*****************************************************************************

	//File_operation***************************************************************
namespace aqua
{
const char * _JudgeStdio(const char *, const string &);
_Ios_Openmode _DefineOpenMode(const string &);
int FileExist(const char *, const bool);
char * DirExist(const char *);
string GetCWD();
/* File open operators*/
/*
+: append
b: binary
r: input
w: output
-: truncate
c: check for existence, exit if error.
*/
int OpenFile(fstream &, const char *, const string);
int OpenFileGZ(gzstreambase &, const char *, const string);
ubyte8 FileSizeC(FILE *);
ubyte8 FileSizeCPP(fstream &);
}
//*****************************************************************************

//Regex************************************************************************
#ifdef _ENABLE_BOOST_REGEX
namespace aqua
{
class RegexBasic : boost::noncopyable
{
public:
	explicit RegexBasic(const string & reg_str) : reg(reg_str), success(false) {}

protected:
	boost::regex reg;
	bool success;
	boost::smatch matched;
};

class RegexSearch : public RegexBasic
{
public:
	explicit RegexSearch(const string & target_str, const string reg_str) : RegexBasic(reg_str)
	{
		success = boost::regex_search(target_str, matched, reg);
	}

	bool is_success();

	smatch::const_reference operator [](uint);

	smatch::size_type size();
};
}
#endif

namespace aqua
{
int StrToVector(VSTR &, string &, const string &);
int ImportToVector(const char *, VSTR &, const char, const bool);
int ImportToVectorGZ(const char *, VSTR &, const char, const bool);
int ImportToVVSTR(char *, VVSTR &, const char, const string &, const bool);
#ifdef _ENABLE_BOOST_REGEX
int SplitByRegex(VSTR &, string &, const char *);
#endif
int SplitBySS(VSTR &, string &, const char);
}
//*****************************************************************************

namespace aqua
{
string SuitFilenameToShell(string);
}

namespace aqua
{
char * GetCurrentTimeString();
void PrintCurrentTimeString();
void RecordOneTime();
void PrintDuration();
}

namespace aqua
{
string RegisterFile(const string &);
void checkFileAndDeleteZeroSize();
string GetFilenamePart(const string &);
}

#endif /* GENERAL_H_ */
