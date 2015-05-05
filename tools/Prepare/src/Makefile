# Generated automatically from Makefile.in by configure.
SHELL = /bin/sh

exec_prefix = .
bindir = $(exec_prefix)/bin
libdir = 
mandir = 

CC = gcc
CCOPT = -O3 -minline-all-stringops -funsafe-loop-optimizations -fprefetch-loop-arrays -funroll-loops -fomit-frame-pointer -fno-bounds-check -Wall -W #-pg -g -D_DEBUG -O0 -m32 -march=i386 -mtune=i386
#CCOPT = -O0 -g -D_DEBUG -Wall
LIBS = -lm -lpthread
INCDIRS = -Iinc/
CFLAGS = ${CCOPT} ${INCDIRS}

all: finalFusion
SRCS1 = searchPath.c scaffold.c check.c seq.c bundle.c potential.c\
	loadGraph.c mem_manager.c attachPEinfo.c newhash.c\
	output_scaffold.c orderContig.c connect.c hashFunction.c\
	readseq1by1.c fib.c fibHeap.c stack.c kmer.c prepare.c
OBJS1 = searchPath.o scaffold.o check.o seq.o bundle.o potential.o\
	loadGraph.o mem_manager.o attachPEinfo.o newhash.o\
	output_scaffold.o orderContig.o connect.o hashFunction.o\
	readseq1by1.o fib.o fibHeap.o stack.o kmer.o prepare.o

SRCS2 = prlHashCtg.c prlRead2Ctg.c map.c localAsm.c\
	lib.c darray.c prlReadFillGap.c read2scaf.c
OBJS2 = prlHashCtg.o prlRead2Ctg.o map.o localAsm.o\
	lib.o darray.o prlReadFillGap.o read2scaf.o


SRCS3 = main.c 
OBJS3 = main.o

.c.o :
	$(CC) $(CFLAGS) -c $< 

finalFusion:  $(OBJS1) $(OBJS2) $(OBJS3)
	$(CC) $(CCOPT) -o $@ $^ $(LIBS)

install: all
	/bin/cp finalFusion $(bindir) 
		
# removes whatever can be built with make except xypost
clean:
	/bin/rm -f ${exec_prefix}/*.o

