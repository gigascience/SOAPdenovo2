/*
 *      Filename: general.cpp
 *
 *      Description:
 *         Basic functions
 *
 *  	Created on: Feb 8, 2010
 *      Author: Ruibang Luo
 *
 */

#define _EXAM_ASSERT_TEST_
#include<iostream>
#include"general.h"
#include<cstring>
#include<sstream>

#ifdef _ENABLE_BOOST_REGEX
#include<boost/regex.hpp>
#endif

#define VEC_INIT_SIZE 0xFFU

int aqua::verbosity(0x1 << 1);
static char cwd[FN_SIZE];
static time_t rawTime;
static deque<time_t> timeQueue;
static vector<string> exitSizeCheckList;

int aqua::FileExist(const char * filename, const bool aggresive)
{
	ifstream I(filename, ios::in | ios::binary);
	if(!I)
	{
		if(aggresive)
		{
			cerr<<"File missing: "<<filename<<endl;
			exit(EXIT_FAILURE);
		}
		return 0;
	}
	I.close();
	return 1;
}

const char * DirExistORCWD(const char * prefix)
{
	char sysCMD[FN_SIZE];
	snprintf(sysCMD, FN_SIZE - 1, "[ -d %s/ ]", prefix);
	if(!system(sysCMD))
	{
		return prefix;
	}
	else
	{
		return NULL;
	}
}

string aqua::GetCWD()
{
	char* status;
	status = getcwd(cwd, FN_SIZE - 1);
	return string(cwd);
}

/* File open operators*/
/*
+: append
b: binary
r: input
w: output
-: truncate
c: check for existence, exit if error.
*/

const char * stdinPath("/dev/stdin");
const char * stdoutPath("/dev/stdout");

const char * aqua::_JudgeStdio(const char * file_name, const string & str)
{
	if(strlen(file_name) == 1 && file_name[0] == '-')
	{
		if(str.find("r") != string::npos)
		{
			//cerr<<"Input redirect to STDIN"<<endl;
			return stdinPath;
		}
		else if(str.find("w") != string::npos)
		{
			//cerr<<"Output redirect to STDOUT"<<endl;
			return stdoutPath;
		}
		else
			return file_name;
	}

	return file_name;
}

ostream & aqua::line(ostream & O)
{
  return O<<"\n--------------------------------------------------------------------------------\n";
}

ostream & aqua::endt(ostream & O)
{
  return O<<"\t";
}

inline _Ios_Openmode aqua::_DefineOpenMode(const string & str)
{
	_Ios_Openmode open_mode(ios::trunc ^ ios::trunc);

	for(uint i(0); i < str.size(); ++i)
	{
		switch(str[i])
		{
			case '+': open_mode |= ios::app; break;
			case 'b': open_mode |= ios::binary; break;
			case 'r': open_mode |= ios::in; break;
			case 'w': open_mode |= ios::out; break;
			case '-': open_mode |= ios::trunc; break;
			//default: break;
		}
	}
	ASSERT( (!((open_mode & ios::app) && (open_mode & ios::trunc))), "Append and Truncate at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::out))), "Input and Output at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::trunc))), "Input and Truncate at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::app))), "Input and Append at the same time.");

	return open_mode;
}

int aqua::OpenFile(fstream & FH, const char * file_name, const string str)
{
	FH.open(_JudgeStdio(file_name, str), _DefineOpenMode(str));
	return FileExist(file_name, (str.find('c') != string::npos));
}

int aqua::OpenFileGZ(gzstreambase & FH, const char * file_name, const string str)
{
	FH.open(_JudgeStdio(file_name, str), _DefineOpenMode(str));
	return FileExist(file_name, (str.find('c') != string::npos));
}

ubyte8 aqua::FileSizeC(FILE * fh)
{
	ubyte8 cur(ftell(fh));
	fseek(fh, 0, SEEK_END);
	ubyte8 size(ftell(fh));
	fseek(fh, cur, SEEK_SET);

	return size;
}

ubyte8 aqua::FileSizeCPP(fstream & fh)
{
	ubyte8 cur(fh.tellg());
	fh.seekg(0, ios_base::end);
	ubyte8 size(fh.tellg());
	fh.seekg(cur, ios_base::beg);

	return size;
}

int aqua::StrToVector(VSTR & vec_str, string & whole_str, const string & reg_str)
{
	ASSERT((reg_str.size() != 0), "Split Regex equal to 0!");
	SplitBySS(vec_str, whole_str, reg_str.c_str()[0]);

	/*if(vec_str.size() && vec_str[vec_str.size() - 1].empty())
		vec_str.resize(vec_str.size() - 1);*/

	return vec_str.size();
}

int aqua::ImportToVector(const char * filename, VSTR & vec_str, const char line_cut = '\n', const bool ambiguous = true)
{
	fstream I;
	if(!OpenFile(I, filename, (ambiguous ? "rc" : "r")))
		return 0;

	uint pos(0);
	vec_str.resize(VEC_INIT_SIZE);
	while(true)
	{
		getline(I, vec_str[pos], line_cut);
		if(!I)
		{
			vec_str.resize(pos);
			break;
		}
		if(!vec_str[pos].empty() && vec_str[pos][0] == '#') continue;
		++pos;
		if(pos == vec_str.size())
			vec_str.resize(vec_str.size() << 1);
	}

	return vec_str.size();
}

int aqua::ImportToVectorGZ(const char * filename, VSTR & vec_str, const char line_cut = '\n', const bool ambiguous = true)
{
	igzstream I;
	if(!OpenFileGZ(I, filename, (ambiguous ? "rc" : "r")))
		return 0;

	uint pos(0);
	vec_str.resize(VEC_INIT_SIZE);
	while(true)
	{
		getline(I, vec_str[pos], line_cut);
		if(!I)
		{
			vec_str.resize(pos);
			break;
		}
		++pos;
		if(pos == vec_str.size())
			vec_str.resize(vec_str.size() << 1);
	}

	return vec_str.size();
}

int aqua::ImportToVVSTR(char * filename, VVSTR & vec_vec_str, const char line_cut, const string & reg_str, const bool ambiguous)
{
	fstream I;
	if(!OpenFile(I, filename, (ambiguous ? "rc" : "r")))
		return 0;

	uint pos(0);
	vec_vec_str.resize(VEC_INIT_SIZE);
	string tmp;
	while(true)
	{
		getline(I, tmp, line_cut);
		if(!I)
		{
			vec_vec_str.resize(pos);
			break;
		}
		StrToVector(vec_vec_str[pos], tmp, reg_str);
		++pos;
		if(pos == vec_vec_str.size())
			vec_vec_str.resize(vec_vec_str.size() << 1);
	}

	return vec_vec_str.size();
}

int aqua::SplitBySS(VSTR & vec, string & str, const char reg)
{

	stringstream ss(str);
	vec.resize(VEC_INIT_SIZE);
	register uint pos(0);
	for(;;)
	{
		getline(ss, vec[pos], reg);
		if(!ss)
		{
			vec.resize(pos);
			break;
		}
		++pos;
		if(pos == vec.size())
			vec.resize(vec.size() << 1);
	}

	return vec.size();
}

#ifdef _ENABLE_BOOST_REGEX
smatch::const_reference aqua::RegexSearch::operator [](uint i)
{
	return matched[i];
}

smatch::size_type aqua::RegexSearch::size()
{
	return matched.size();
}

bool aqua::RegexSearch::is_success()
{
	return success;
}

int aqua::SplitByRegex(VSTR & vec, string & str, const char * reg)
{
	regex e(reg);
	regex_split(back_inserter(vec), str, e);
	return vec.size();
}
#endif

void exam_assert(const char * file_name, const char * func_name, uint line_no)
{
	cerr<<endl<<DEBUG_INFO<<":Assert failure"<<endl;
	abort();
}

void exam_assert(const char * file_name, const char * func_name, uint line_no, const char * reason)
{
	cerr<<endl<<DEBUG_INFO<<":Assert failure\nReason: "<<endl<<reason<<endl;
	abort();
}

int aqua::ModifyVerbosity(const int tmp)
{
	if(tmp)
	{
		if(tmp < VERBOSITY_BOTTOM || tmp > VERBOSITY_TOP)
		{
			return -1;
		}
		aqua::verbosity = (0x1 << tmp);
	}
	else
	{
		aqua::verbosity <<= 1;
	}

	return aqua::verbosity;
}

char aqua::verboseStr[verboseBufSize];

string aqua::SuitFilenameToShell(string name)
{
	for(register uint i(0); i < name.size();)
	{
		i = name.find_first_of("<>*&$\\;|\'\"[]=:", i);
		name.insert(i, "\\");
		i += 2;
	}

	return name;
}

char * aqua::GetCurrentTimeString()
{
	time(&rawTime);
	return ctime(&rawTime);
}

void aqua::PrintCurrentTimeString()
{
	fprintf(stderr, "Now: %s", aqua::GetCurrentTimeString());
}

void aqua::RecordOneTime()
{
	timeQueue.push_back(time(0));
}

void aqua::PrintDuration()
{
	if(timeQueue.size())
	{
		fprintf(stderr, "Duration: %lu sec(s)\n", (unsigned long)time(0) - (unsigned long)timeQueue.back());
		timeQueue.pop_back();
	}
}

string aqua::RegisterFile(const string & tmp)
{
	exitSizeCheckList.push_back(tmp);
	return tmp;
}

void aqua::checkFileAndDeleteZeroSize()
{
	for(uint i(0); i < exitSizeCheckList.size(); ++i)
	{
		FILE * fh = fopen(exitSizeCheckList[i].c_str(), "r");
		if(fh != NULL)
		{
			if(aqua::FileSizeC(fh) == 0)
			{
				fclose(fh);
				if(remove(exitSizeCheckList[i].c_str()) != 0)
					mvnv(3, "Error Cleaning Zero Size File: %s", exitSizeCheckList[i].c_str());
			}
		}
	}
}

string aqua::GetFilenamePart(const string & a_file)
{
	return ((a_file.find("/")==string::npos)?(a_file):(a_file.substr(a_file.rfind("/")+1, string::npos)));
}
