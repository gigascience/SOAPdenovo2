/*
 * ec.cpp
 *
 *  Created on: 2011-8-10
 *      Author: Aqua
 */

#include<iostream>
#include<cstring>
#include<cstdio>
#include<csignal>
#include<cstdlib>
#include<cstdarg>
#include<sys/mman.h>
#include"kmerFreq.h"
#include"general.h"
#include"connect.h"
#include"correct.h"

using namespace std;
using namespace aqua;

const char * AUTHOR  = "Aqua Lo";
const char * VERSION = "0.04";
const char * CONTACT = "rbluo@cs.hku.hk";

static void usage(char ** argv)
{
	cerr
	<< "Version:       " << VERSION << "\n"
	<< "Compile Date:  " << __DATE__ << " " << __TIME__ << "\n"
	<< "Author:        " << AUTHOR << "\n"
	<< "Contact:       " << CONTACT << "\n"
	<< "\n"
	<< argv[0] << " <MODE> [PARAMETERS]\n"
	<<"    connect       Connect PE reads"<<'\n'
	<<"    kmerfreq      Build Kmer Frequency Table"<<'\n'
	<<"    correct       Reads Error Correction"<<'\n';

	exit(EXIT_FAILURE);
}

struct sigaction sigact;

static void signalHandler(int sig)
{
	if (sig == SIGINT)
	{
		mvnv(0, "Caught signal for Ctrl+C");
		mvnv(0, "Cleaning and terminating ...");
		exit(0);
	}
}

void initSignals(void)
{
	struct sigaction sigact;
	sigact.sa_handler = signalHandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, (struct sigaction *)NULL);
}

void cleanup(void)
{
	sigemptyset(&sigact.sa_mask);
	munlockall();
	checkFileAndDeleteZeroSize();
	PrintCurrentTimeString();
}

int main(int argc, char ** argv)
{
	RecordOneTime();

	atexit(cleanup);
	initSignals();

	mlockall(MCL_FUTURE);
	{
		uint8_t * tmp = (uint8_t*)mmap(NULL, sizeof(uint8_t) * 131072, PROT_READ|PROT_WRITE|PROT_EXEC, MAP_ANON|MAP_PRIVATE|MAP_NORESERVE, -1, 0);
		if(tmp == MAP_FAILED)
			perrdie("\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:");
		munmap(tmp, sizeof(uint8_t) * 131072);
	}

	if(argc < 2)
		usage(argv);

	if(strcmp(argv[1], "connect") == 0)
	{
		Connect(argc, argv);
	}
	else if(strcmp(argv[1], "kmerfreq") == 0)
	{
		Kmerfreq(argc, argv);
	}
	else if(strcmp(argv[1], "correct") == 0)
	{
		Correct(argc, argv);
	}
	else
	{
		die("Irregonized command: %s", argv[1])
	}

	PrintDuration();
	return EXIT_SUCCESS;
}
