/*
 * correct.cpp
 *
 *  Created on: 2011-8-13
 *      Author: Aqua
 */

#include "correct.h"

using namespace std;
using namespace aqua;

void FreqContReg_s::modify(const int startT, const int endT, const int statusT)
{
	this->start = startT;
	this->end = endT;
	this->status = statusT;
}

TreeNode_s::TreeNode_s()
{
	this->pointer = 0;
	this->base = 0;
	this->change = 0;
	this->same = 0;
}

TreeNode_s::TreeNode_s(const uint pointer_t, const uint base_t, const uint change_t, const uint same_t)
{
	this->pointer = pointer_t;
	this->base = base_t;
	this->change = change_t;
	this->same = same_t;
}

void TreeNode_s::modify(const uint pointer_t, const uint base_t, const uint change_t, const uint same_t)
{
	this->pointer = pointer_t;
	this->base = base_t;
	this->change = change_t;
	this->same = same_t;
}

void CorrGlobalVar_s::print(const char * cstr)
{
	mvnv(0, "%s", cstr);
	mvnv(0, "                 Consecutive kmer size: %d", this->conKmerSize);
	mvnv(0, "            Space kmer solid-part size: %d", this->spaKmerSize);
	mvnv(0, "            Space kmer space-part size: %d", this->spaSize);
	mvnv(0, "Consecutive Kmer Frequency lower limit: %d", this->conLowFreqCutoff);
	mvnv(0, "Consecutive Kmer Frequency upper limit: %d", this->conHighFreqCutoff);
	mvnv(0, "      Space Kmer Frequency lower limit: %d", this->spaLowFreqCutoff);
	mvnv(0, "      Space Kmer Frequency upper limit: %d", this->spaHighFreqCutoff);
	mvnv(0, " Min consecutive high freq kmers/block: %d", this->minHighFreqLen);
	mvnv(0, "     Maximum base alternation per read: %d", this->maxChangePerRead);
	mvnv(0, "Maximum trimmed bases allowed per side: %d", this->trimAllow);
	mvnv(0, "      Min read length after correction: %d", this->minOutputReadLen);
	mvnv(0, "Min base interval between 2 correction: %d", this->minInterval);
	mvnv(0, "         Delete corretion failed reads: %s", this->toDelete?"y":"n");
	mvnv(0, "        Connect reads after correction: %s", this->joinRead?"y":"n");
	mvnv(0, "                 Correct using overlap: %s", this->corrByOvl?"y":"n");
	if(this->joinRead || this->corrByOvl)
	{
		mvnv(0, "                      Min align length: %d", this->conn_Min);
		mvnv(0, "                      Max align length: %d", this->conn_Max);
		mvnv(0, "                     Mean align length: %d", this->conn_Mean);
		mvnv(0, "     Identity cutoff of aligned region: %.2f", this->conn_MATCH_CUTOFF);
		mvnv(0, "  Min diff best and second best result: %.2f", this->conn_DIF);
	}
	mvnv(0, "                             Verbosity: %d", (int)log2(verbosity));
	mvnv(0, "Kmer Scene:\n%s", PrintKmerScene(this->conKmerSize, this->spaKmerSize, this->spaSize).c_str());
}

static void GetMaxHighfreqBlock(const vector<FreqContReg_s> & cont, FreqContReg_s & max)
{
	max.start = max.end = -1;
	max.status = 1;
	if(cont.empty())
	{
		mvnv(0, "vector<FreqContReg_s> cont.size() == empty");
		return;
	}
	int highStart(0), highEnd(0);
	uint i(0);
	while(i < cont.size())
	{
		while(!cont[i].status && i < cont.size())
			++i;
		if(i < cont.size())
			highStart = cont[i].start;
		while(cont[i].status && i < cont.size())
		{
			highEnd = cont[i].end;
			++i;
		}
		if((highEnd - highStart) > (max.end - max.start))
		{
			max.start = highStart;
			max.end = highEnd;
		}
	}
}

static void GetConKmerfreqBlock(const CorrGlobalVar_s & args, Read_s & read, vector<FreqContReg_s> & conFreqContReg)
{
	const int totalKmer(read.read.size() - args.conKmerSize + 1);

	uint64_t kbit;
	int kbitFreq, lowStart, lowEnd, highStart, highEnd;
	int i(0);
	int tr(0);
	while(i < totalKmer)
	{
		lowStart = i;
		while(i < totalKmer)
		{
#define GETKBIT\
			if(tr)\
				tr = 0;\
			else if(i)\
				kbit = GetNextKbit(kbit, args.conKmerMasker, read.read[i+args.conKmerSize-1]);\
			else\
				kbit = Seq2Bit(read.read.substr(i, args.conKmerSize));\
			kbitFreq = getFreq(args.conKmerTable, kbit);
			//
			GETKBIT;

			if(kbitFreq)
				break;
			else
				++i;
		}
		lowEnd = i;
		tr = 1;

		if(lowEnd > lowStart)
			conFreqContReg.push_back(FreqContReg_s(lowStart+1, lowEnd, 0));

		highStart = i;
		while(i < totalKmer)
		{
			GETKBIT;
#undef GETKBIT

			if(kbitFreq)
				++i;
			else
				break;
		}
		highEnd = i;
		tr = 1;
		if(highEnd > highStart)
			conFreqContReg.push_back(FreqContReg_s(highStart+1, highEnd, 1));
	}

	if(args.commentOn)
	{
		stringstream ss;
		ss<<"M";
		for(uint i(0); i < conFreqContReg.size(); ++i)
		{
			if(conFreqContReg[i].status)
				ss<<"H";
			else
				ss<<"L";
			ss<<conFreqContReg[i].start<<"E"<<conFreqContReg[i].end;
		}
		read.comment += ss.str();
	}
}

static void GetAllKmerfreqBlock(const CorrGlobalVar_s & args, Read_s & read, vector<FreqContReg_s> & freqContReg)
{
	int totalKmer(read.read.size() - args.spaKmerSize - args.spaSize + 1);
	int spaHeadSize, spaTailSize;
	spaTailSize = args.spaKmerSize / 2;
	spaHeadSize = args.spaKmerSize - spaTailSize;

	string conKmerSeq, spaKmerSeq, spaSeq, spaHeadSeq, spaTailSeq;
	uint64_t conKbit, spaKbit;
	int conKbitFreq, spaKbitFreq, lowStart, lowEnd, highStart, highEnd;
	int i(0);
	int tr(0);
	while(i < totalKmer)
	{
		lowStart = i;
		while(i < totalKmer)
		{
#define GETKBIT\
			if(tr)\
				tr = 0;\
			else if(i)\
			{\
				conKbit = GetNextKbit(conKbit, args.conKmerMasker, read.read[i+args.conKmerSize-1]);\
				spaKbit = GetNextKbitSpace(spaKbit, args.spaKmerHeadMasker, args.spaKmerTailMasker, args.spaHeadOffset, read.read[i+spaHeadSize-1], read.read[i+args.spaKmerSize+args.spaSize-1]);\
			}\
			else\
			{\
				conKmerSeq = read.read.substr(i, args.conKmerSize);\
				conKbit = Seq2Bit(conKmerSeq);\
				\
				spaSeq = read.read.substr(i, args.spaKmerSize + args.spaSize);\
				spaHeadSeq = spaSeq.substr(0, spaHeadSize);\
				spaTailSeq = spaSeq.substr(spaHeadSize+args.spaSize, spaTailSize);\
				spaKbit = Seq2Bit(spaHeadSeq+spaTailSeq);\
			}\
			conKbitFreq = getFreq(args.conKmerTable, conKbit);\
			spaKbitFreq = getFreq(args.spaKmerTable, spaKbit);
			//
			GETKBIT;
			if(!conKbitFreq || !spaKbitFreq)
				++i;
			else
				break;
		}
		lowEnd = i;
		tr = 1;

		if(lowEnd > lowStart)
			freqContReg.push_back(FreqContReg_s(lowStart+1, lowEnd, 0));

		highStart = i;
		while(i < totalKmer)
		{
			GETKBIT;
#undef GETKBIT
			if(conKbitFreq && spaKbitFreq)
				i++;
			else
				break;
		}
		highEnd = i;
		tr = 1;
		if(highEnd > highStart)
			freqContReg.push_back(FreqContReg_s(highStart+1, highEnd, 1));
	}

	if(args.commentOn)
	{
		stringstream ss;
		ss<<"A";
		for(uint i(0); i < freqContReg.size(); ++i)
		{
			if(freqContReg[i].status)
				ss<<"H";
			else
				ss<<"L";
			ss<<freqContReg[i].start<<"E"<<freqContReg[i].end;
		}
		read.comment += ss.str();
	}
}

static int ConCorrect1BaseFast(const CorrGlobalVar_s & args, Read_s & read, FreqContReg_s & freqContReg)
{
	int highFreqKmerCount(0), possibleBaseCount(0);
	char errorBase(read.read[freqContReg.end-1]);
	char correctBase(errorBase);
	uint64_t kbit; int kbitFreq; int kbitInvalid(0);

	for(uint i(0); i < 4; ++i)
	{
		if(errorBase != bases[i])
		{
			read.read[freqContReg.end-1] = bases[i];
			highFreqKmerCount = 0;
			for(int j(freqContReg.start-1), k(0); j < freqContReg.end; ++j, ++k)
			{
				if(k)
					kbit = GetNextKbit(kbit, args.conKmerMasker, read.read[j+args.conKmerSize-1]);
				else
					kbit = Seq2Bit(read.read.substr(j, args.conKmerSize));
				kbitFreq = getFreq(args.conKmerTable, kbit);
				if(kbitFreq)
				{
					if(args.conInvalidKmerTable != NULL)
						kbitInvalid = getFreq(args.conInvalidKmerTable, kbit);
					if(!kbitInvalid)
						++highFreqKmerCount;
				}
				else
					break;
			}
			if(highFreqKmerCount == (freqContReg.end - freqContReg.start + 1))
			{
				if(read.altAllele.empty() || (read.altAllele[freqContReg.end-1] == 'N') || (bases[i] == read.altAllele[freqContReg.end-1]))
				{
					correctBase = bases[i];
					++possibleBaseCount;
				}
			}
			if(args.commentOn)
			{
				stringstream ss;
				ss<<"C"<<(freqContReg.end-1)<<"B"<<bases[i]<<"F"<<highFreqKmerCount;
				read.comment += ss.str();
			}
		}
	}

	if(possibleBaseCount > 1)
	{
		read.read[freqContReg.end-1] = errorBase;
		return 0;
	}
	else
	{
		if(correctBase != errorBase)
		{
			read.read[freqContReg.end-1] = correctBase;
			return 1;
		}
		else
		{
			read.read[freqContReg.end-1] = errorBase;
			return 0;
		}
	}
	return 0;
}


/*
static int AllCorrect1BaseFast(const CorrGlobalVar_s & args, Read_s & read, FreqContReg_s & freqContReg)
{
	int highFreqKmerCount(0), possibleBaseCount(0);
	char errorBase(read.read[freqContReg.end-1]);
	char correctBase(errorBase);
	int conKbitFreq, spaKbitFreq;

	int spaHeadSize, spaTailSize, spaKmerCutStartPoint;
	spaTailSize = args.spaKmerSize / 2;
	spaHeadSize = args.spaKmerSize - spaTailSize;
	string spaRawSeq, spaHeadSeq, spaTailSeq;

	for(uint i(0); i < 4; ++i)
	{
		if(errorBase != bases[i])
		{
			read.read[freqContReg.end-1] = bases[i];
			highFreqKmerCount = 0;
			for(int j(freqContReg.start-1); j < freqContReg.end; ++j)
			{
				conKbitFreq = getFreq(args.conKmerTable, Seq2Bit(read.read.substr(j+args.spaSize, args.conKmerSize)));

				if(j >= (freqContReg.start - 1 + spaTailSize))
					spaKmerCutStartPoint = j + args.spaSize;
				else
					spaKmerCutStartPoint = j;
				spaRawSeq = read.read.substr(spaKmerCutStartPoint, (args.spaKmerSize + args.spaSize));
				spaHeadSeq = spaRawSeq.substr(0, spaHeadSize);
				spaTailSeq = spaRawSeq.substr(spaHeadSize+args.spaSize, spaTailSize);
				spaKbitFreq = getFreq(args.spaKmerTable, Seq2Bit(spaHeadSeq + spaTailSeq));

				if(conKbitFreq && spaKbitFreq)
					++highFreqKmerCount;
				else
					break;
			}
			if(highFreqKmerCount == (freqContReg.end - freqContReg.start + 1 - args.spaSize))
			{
				if(read.altAllele.empty() || (read.altAllele[freqContReg.end-1] == 'N') || (bases[i] == read.altAllele[freqContReg.end-1]))
				{
					correctBase = bases[i];
					++possibleBaseCount;
				}
			}

			if(args.commentOn)
			{
				stringstream ss;
				ss<<"S"<<(freqContReg.end-1)<<"B"<<bases[i]<<"F"<<highFreqKmerCount;
				read.comment += ss.str();
			}
		}
	}

	if(possibleBaseCount > 1)
	{
		read.read[freqContReg.end-1] = errorBase;
		return 0;
	}
	else
	{
		if(correctBase != errorBase)
		{
			read.read[freqContReg.end-1] = correctBase;
			return 1;
		}
		else
		{
			read.read[freqContReg.end-1] = errorBase;
			return 0;
		}
	}
	return 0;
}
*/


static uint GetHighFreqBlockFromCont(const CorrGlobalVar_s & args, const vector<FreqContReg_s> & cont, vector<FreqContReg_s> & high)
{
	for(uint i(0); i < cont.size(); ++i)
		if(cont[i].status)
			if((cont[i].end - cont[i].start + 1) >= args.minHighFreqLen)
			{
				high.push_back(FreqContReg_s(cont[i].start, cont[i].end, (int)(i+1)));
			}
	return high.size();
}


static inline uint64_t getKmerLeftward(const int kmerSize, const vector<TreeNode_s> & nodeVec, const uint curBase, uint64_t startPointKbit, uint pos)
{
	uint64_t kbit(curBase & 0x3);
	int i(1);
	while((pos > 0) && (i < kmerSize))
	{
		kbit = ((kbit<<2) | nodeVec[pos].base);
		pos = nodeVec[pos].pointer;
		++i;
	}
	if(i < kmerSize)
		kbit = ((kbit << ((kmerSize-i)*2)) | (startPointKbit >> ((i-1)*2)));
	return kbit;
}


static inline uint64_t getKmerRightward(const int kmerSize, const vector<TreeNode_s> & nodeVec, const uint curBase, uint64_t startPointKbit, uint pos)
{
	uint64_t kbit(bitLeft[curBase & 0x3]);
	int i(1);
	while((pos > 0) && (i < kmerSize))
	{
		kbit = ((kbit>>2) | (bitLeft[nodeVec[pos].base]));
		pos = nodeVec[pos].pointer;
		++i;
	}
	while(i < kmerSize)
	{
		kbit = ((kbit>>2) | (bitLeft[startPointKbit & 0x3]));
		startPointKbit >>= 2;
		++i;
	}
	return (kbit >> (64 - kmerSize * 2));
}


static int CorrectBasesRightwardBB(const CorrGlobalVar_s & args, Read_s & read, const int chkStart, const int chkEnd, const int canTrim, int lengthNeedTrim)
{
#ifdef _CORR_DEBUG
	if(chkStart < 0 || chkEnd < 0)
	{
		mvnv(4, "chkStart < 0 || chkEnd < 0");
		return 0;
	}
#endif
	int numCorrected(0);
	uint64_t conStartPointKbit(Seq2Bit(read.read.substr(chkStart-args.conKmerSize, args.conKmerSize-1)));
	uint64_t spaStartPointKbit(Seq2Bit(read.read.substr(chkStart-args.spaKmerSize-args.spaSize, args.spaKmerSize+args.spaSize-1)));

	vector<TreeNode_s> nodeVec;
	vector<uint> nodeCurIdx;
	nodeVec.push_back(TreeNode_s());
	nodeCurIdx.push_back(0);
	uint nodeVecPos(0);
	uint64_t conKbit, spaKbit;
	int conKbitFreq, spaKbitFreq;
	int conKbitInvalid(0), spaKbitInvalid(0);
	uint isEqualBase, changedCount;

	int cycleCur(chkStart);
	uint maxAllowChange(args.maxChangePerRead - read.fastCorrectCount - read.bbCorrectCount);
	if(maxAllowChange > 3)
		maxAllowChange = 3;
	while(cycleCur <= chkEnd)
	{
		vector<uint> nodeTmpIdx;

		for(uint i(0); i < nodeCurIdx.size(); ++i)
		{
			uint & parentPos(nodeCurIdx[i]);
			for(uint j(0); j < 4; ++j)
			{
				conKbit = getKmerRightward(args.conKmerSize, nodeVec, j, conStartPointKbit, parentPos);
				conKbitFreq = getFreq(args.conKmerTable, conKbit);
				if(conKbitFreq)
					if(args.conInvalidKmerTable != NULL)
						conKbitInvalid = getFreq(args.conInvalidKmerTable, conKbit);

				spaKbit = getKmerRightward(args.spaKmerSize+args.spaSize, nodeVec, j, spaStartPointKbit, parentPos);
//#ifdef _CORR_DEBUG
//				string spaSeq;
//				Bit2Seq(spaKbit, args.spaKmerSize+args.spaSize, spaSeq);
//#endif
				spaKbit = ((spaKbit & ((0x2llu << (2*(args.spaKmerSize/2)-1)) - 1)) |\
						( ( spaKbit & ( ( ( 0x2llu << ( 2* ( args.spaKmerSize- ( args.spaKmerSize/2 ) ) -1 ) ) - 1 ) << ( 2 * ( args.spaKmerSize/2+args.spaSize ) ) ) ) >> ( args.spaSize*2 ) ) );
				spaKbitFreq = getFreq(args.spaKmerTable, spaKbit);
				if(spaKbitFreq)
					if(args.spaInvalidKmerTable != NULL)
						spaKbitInvalid = getFreq(args.spaInvalidKmerTable, spaKbit);
//#ifdef _CORR_DEBUG
//				uint64_t tmpSpaKbit; string spaHeadSeq, spaTailSeq;
//				spaTailSeq = spaSeq.substr((args.spaKmerSize-(args.spaKmerSize/2))+args.spaSize, (args.spaKmerSize/2));
//				spaHeadSeq = spaSeq.substr(0, (args.spaKmerSize-(args.spaKmerSize/2)));
//				tmpSpaKbit = Seq2Bit(spaHeadSeq+spaTailSeq);
//				if(spaKbit != tmpSpaKbit)
//					mvnv(4, "spaKbit != tmpSpaKbit: %llu/%llu", spaKbit, tmpSpaKbit);
//#endif

				isEqualBase = ((bases[j] == read.read[cycleCur-1]) ? 1 : 0);
				changedCount = (isEqualBase ? (nodeVec[parentPos].change) : (nodeVec[parentPos].change + 1));
//#ifdef _CORR_DEBUG
//				mvnv(4, "nodeVec.push: parentPos/j/changedCount/isEqualBase/conFreq/spaFreq: %u/%c/%u/%u/%d/%d", parentPos, bases[j], changedCount, isEqualBase, conKbitFreq, spaKbitFreq);
//#endif
				if(conKbitFreq && spaKbitFreq && !conKbitInvalid && !spaKbitInvalid && (changedCount <= maxAllowChange))
				{
					if(read.altAllele.empty() || (read.altAllele[cycleCur-1] == 'N') || (bases[j] == read.altAllele[cycleCur-1]))
					{
						nodeVec.push_back(TreeNode_s(parentPos, j, changedCount, isEqualBase));
						nodeTmpIdx.push_back(++nodeVecPos);
					}
				}
			}
		}

		if(nodeTmpIdx.size() && nodeVecPos < args.maxNode)
		{
			nodeCurIdx.swap(nodeTmpIdx);
		}
		else
		{
			if(nodeVecPos >= args.maxNode)
				mvnv(2, "Node Vector Size %u Exceed Limit %u", nodeVecPos, args.maxNode);
			break;
		}
		++cycleCur;

		mvnvd(4, "BB Rightward: chkStart/chkEnd/cycleCur/nodeVecSize: %d/%d/%d/%lu", chkStart, chkEnd, cycleCur, nodeVec.size());
	}

	int minChangePos(nodeCurIdx[0]);
	int minChangeCount(nodeVec[minChangePos].change);
	int possiblePathCount(0);
	for(uint i(0); i < nodeCurIdx.size(); ++i)
	{
		uint & curPos(nodeCurIdx[i]);
		if(nodeVec[curPos].change < minChangeCount)
		{
			minChangeCount = nodeVec[curPos].change;
			minChangePos = curPos;
			possiblePathCount = 1;
		}
		else if(nodeVec[curPos].change == minChangeCount)
		{
			++possiblePathCount;
		}
	}

	lengthNeedTrim = (chkEnd - cycleCur + 1);
	if(args.minInterval)
	{
		if(possiblePathCount == 1)
		{
			uint pos(minChangePos), idx(cycleCur - 1);
			uint previousIdx(0), previousDeleted(0), previousPos;
			while(pos > 0)
			{
				if(!nodeVec[pos].same)
				{
					if(!previousIdx)
					{
						previousPos = pos;
						previousIdx = idx;
					}
					else if((previousIdx - idx) <= args.minInterval)
					{
						--minChangeCount;
						nodeVec[pos].base = (base2int(read.read[idx-1]) & 0x3);
						--nodeVec[pos].change;
						nodeVec[pos].same = 1;
						if(!previousDeleted)
						{
							--minChangeCount;
							nodeVec[previousPos].base = (base2int(read.read[previousIdx-1]) & 0x3);
							--nodeVec[previousPos].change;
							nodeVec[previousPos].same = 1;
							previousDeleted = 1;
						}
						if(!minChangeCount)
						{
							possiblePathCount = 0;
							break;
						}
						previousPos = pos;
						previousIdx = idx;
					}
				}
				pos = nodeVec[pos].pointer;
				--idx;
			}
		}
	}

	if((possiblePathCount == 1) && (!lengthNeedTrim || (lengthNeedTrim && canTrim)))
	{
		numCorrected = minChangeCount;
		uint pos(minChangePos), idx(cycleCur - 1);
		while(pos > 0)
		{
			if(!nodeVec[pos].same)
			{
				read.read[idx-1] = bases[nodeVec[pos].base];
				read.bbCorrectRecord.push_back(idx);
				read.bbCorrectCharRec.push_back(bases[nodeVec[pos].base]);
			}
			pos = nodeVec[pos].pointer;
			--idx;
		}
		read.trimRight = lengthNeedTrim;
		if(lengthNeedTrim)
			read.trimRight += (maxAllowChange - minChangeCount);
	}
	else
		numCorrected = 0;

	if(args.commentOn)
	{
		stringstream ss;
		ss<<"R"<<chkEnd<<"P"<<possiblePathCount<<"C"<<minChangeCount;
		read.comment += ss.str();
	}

	return numCorrected;
}


static int CorrectBasesLeftwardBB(const CorrGlobalVar_s & args, Read_s & read, const int chkStart, const int chkEnd, const int canTrim, int lengthNeedTrim)
{
	if(chkStart < 0 || chkEnd < 0)
	{
		mvnvd(4, "chkStart < 0 || chkEnd < 0");
		return 0;
	}
	int numCorrected(0);
	uint64_t conStartPointKbit(Seq2Bit(read.read.substr(chkStart, args.conKmerSize-1)));
	uint64_t spaStartPointKbit(Seq2Bit(read.read.substr(chkStart, args.spaKmerSize+args.spaSize-1)));

	vector<TreeNode_s> nodeVec;
	vector<uint> nodeCurIdx;
	nodeVec.push_back(TreeNode_s());
	nodeCurIdx.push_back(0);
	uint nodeVecPos(0);
	uint64_t conKbit, spaKbit;
	int conKbitFreq, spaKbitFreq;
	uint isEqualBase, changedCount;

	int cycleCur(chkStart);
	uint maxAllowChange(args.maxChangePerRead - read.fastCorrectCount - read.bbCorrectCount);
	if(maxAllowChange > 3)
		maxAllowChange = 3;
	while(cycleCur >= chkEnd)
	{
		vector<uint> nodeTmpIdx;

		for(uint i(0); i < nodeCurIdx.size(); ++i)
		{
			uint & parentPos(nodeCurIdx[i]);
			for(uint j(0); j < 4; ++j)
			{
				conKbit = getKmerLeftward(args.conKmerSize, nodeVec, j, conStartPointKbit, parentPos);
				conKbitFreq = getFreq(args.conKmerTable, conKbit);

				spaKbit = getKmerLeftward(args.spaKmerSize+args.spaSize, nodeVec, j, spaStartPointKbit, parentPos);
//#ifdef _CORR_DEBUG
//				string spaSeq;
//				Bit2Seq(spaKbit, args.spaKmerSize+args.spaSize, spaSeq);
//#endif
				spaKbit = ((spaKbit & ((0x2llu << (2*(args.spaKmerSize/2)-1)) - 1)) |\
						( ( spaKbit & ( ( ( 0x2llu << ( 2* ( args.spaKmerSize- ( args.spaKmerSize/2 ) ) -1 ) ) - 1 ) << ( 2 * ( args.spaKmerSize/2+args.spaSize ) ) ) ) >> ( args.spaSize*2 ) ) );
				spaKbitFreq = getFreq(args.spaKmerTable, spaKbit);
//#ifdef _CORR_DEBUG
//				uint64_t tmpSpaKbit; string spaHeadSeq, spaTailSeq;
//				spaTailSeq = spaSeq.substr((args.spaKmerSize-(args.spaKmerSize/2))+args.spaSize, (args.spaKmerSize/2));
//				spaHeadSeq = spaSeq.substr(0, (args.spaKmerSize-(args.spaKmerSize/2)));
//				tmpSpaKbit = Seq2Bit(spaHeadSeq+spaTailSeq);
//				if(spaKbit != tmpSpaKbit)
//					mvnv(4, "spaKbit != tmpSpaKbit: %llu/%llu", spaKbit, tmpSpaKbit);
//#endif

				isEqualBase = ((bases[j] == read.read[cycleCur-1]) ? 1 : 0);
				changedCount = (isEqualBase ? nodeVec[parentPos].change : (nodeVec[parentPos].change + 1));
//#ifdef _CORR_DEBUG
//				mvnv(4, "nodeVec.push: parentPos/j/changedCount/isEqualBase/conFreq/spaFreq: %u/%c/%u/%u/%d/%d", parentPos, bases[j], changedCount, isEqualBase, conKbitFreq, spaKbitFreq);
//#endif
				if(conKbitFreq && spaKbitFreq && (changedCount <= maxAllowChange))
				{
					if(read.altAllele.empty() || (read.altAllele[cycleCur-1] == 'N') || (bases[j] == read.altAllele[cycleCur-1]))
					{
						nodeVec.push_back(TreeNode_s(parentPos, j, changedCount, isEqualBase));
						nodeTmpIdx.push_back(++nodeVecPos);
					}
				}
			}
		}

		if(nodeTmpIdx.size() && nodeVecPos < args.maxNode)
		{
			nodeCurIdx.swap(nodeTmpIdx);
		}
		else
		{
			if(nodeVecPos >= args.maxNode)
				mvnv(2, "Node Vector Size %u Exceed Limit %u", nodeVecPos, args.maxNode);
			break;
		}
		--cycleCur;

		mvnvd(4, "BB Leftward: chkStart/chkEnd/cycleCur/nodeVecSize: %d/%d/%d/%lu", chkStart, chkEnd, cycleCur, nodeVec.size());

	}

	int minChangePos(nodeCurIdx[0]);
	int minChangeCount(nodeVec[minChangePos].change);
	int possiblePathCount(0);
	for(uint i(0); i < nodeCurIdx.size(); ++i)
	{
		uint & curPos(nodeCurIdx[i]);
		if(nodeVec[curPos].change < minChangeCount)
		{
			minChangeCount = nodeVec[curPos].change;
			minChangePos = curPos;
			possiblePathCount = 1;
		}
		else if(nodeVec[curPos].change == minChangeCount)
		{
			++possiblePathCount;
		}
	}

	lengthNeedTrim = (cycleCur - chkEnd + 1);
	if(args.minInterval)
	{
		if(possiblePathCount == 1)
		{
			uint pos(minChangePos), idx(cycleCur + 1);
			uint previousIdx(0), previousDeleted(0), previousPos;
			while(pos > 0)
			{
				if(!nodeVec[pos].same)
				{
					if(!previousIdx)
					{
						previousPos = pos;
						previousIdx = idx;
					}
					else if((idx - previousIdx) <= args.minInterval)
					{
						--minChangeCount;
						nodeVec[pos].base = (base2int(read.read[idx-1]) & 0x3);
						--nodeVec[pos].change;
						nodeVec[pos].same = 1;
						if(!previousDeleted)
						{
							--minChangeCount;
							nodeVec[previousPos].base = (base2int(read.read[previousIdx-1]) & 0x3);
							--nodeVec[previousPos].change;
							nodeVec[previousPos].same = 1;
							previousDeleted = 1;
						}
						if(!minChangeCount)
						{
							possiblePathCount = 0;
							break;
						}
						previousPos = pos;
						previousIdx = idx;
					}
				}
				pos = nodeVec[pos].pointer;
				++idx;
			}
		}
	}

	if((possiblePathCount == 1) && (!lengthNeedTrim || (lengthNeedTrim && canTrim)))
	{
		numCorrected = minChangeCount;
		uint pos(minChangePos), idx(cycleCur + 1);
		while(pos > 0)
		{
			if(!nodeVec[pos].same)
			{
				read.read[idx-1] = bases[nodeVec[pos].base];
				read.bbCorrectRecord.push_back(idx);
				read.bbCorrectCharRec.push_back(bases[nodeVec[pos].base]);
			}
			pos = nodeVec[pos].pointer;
			++idx;
		}
		read.trimLeft = lengthNeedTrim;
		if(lengthNeedTrim)
			read.trimLeft += (maxAllowChange - minChangeCount);
	}
	else
		numCorrected = 0;

	if(args.commentOn)
	{
		stringstream ss;
		ss<<"L"<<chkEnd<<"P"<<possiblePathCount<<"C"<<minChangeCount;
		read.comment += ss.str();
	}

	return numCorrected;
}


static uint64_t Correct1Read(const CorrGlobalVar_s & args, Read_s & read)
{
	mvnvd(4,"Marker 1.1");
	{
		vector<FreqContReg_s> fastFreqContReg;
		GetConKmerfreqBlock(args, read, fastFreqContReg);
		if(fastFreqContReg.size() > ((size_t)(read.read.size() * 0.1)))
		{
			read.shouldDelete = 1;
			return 0;
		}

		for(uint i(1); i < (fastFreqContReg.size() - 1); ++i)
		{
			if(fastFreqContReg[i].status != 0)
				continue;
			if(read.fastCorrectCount > args.maxChangePerRead)
				break;

			int isCorrected(0);
			if((fastFreqContReg[i].end - fastFreqContReg[i].start + 1) == args.conKmerSize)
				isCorrected = ConCorrect1BaseFast(args, read, fastFreqContReg[i]);
			if(isCorrected)
			{
				read.fastCorrectRecord.push_back(fastFreqContReg[i].end);
				read.fastCorrectCharRec.push_back(read.read[fastFreqContReg[i].end-1]);
				fastFreqContReg[i].status = 1;
				++read.fastCorrectCount;
			}
		}
	}
	mvnvd(4,"Marker 1.2");

	vector<FreqContReg_s> freqContReg;
	GetAllKmerfreqBlock(args, read, freqContReg);

	vector<FreqContReg_s> freqHighReg;
	GetHighFreqBlockFromCont(args, freqContReg, freqHighReg);

	mvnvd(4,"Marker 1.3");
	if(!freqHighReg.size())
	{
		read.shouldDelete = 1;
		return 0;
	}

	//Extend low freq edge
	{
		int edgeCutLen(args.minHighFreqLen / 2);
		int totalKmer(read.read.size() - args.spaKmerSize - args.spaSize +  1);
		for(uint i(0); i < freqHighReg.size(); ++i)
		{
			if(freqHighReg[i].start != 1)
			{
				freqHighReg[i].start += edgeCutLen;
				if(freqHighReg[i].start > totalKmer)
					freqHighReg[i].start -= edgeCutLen;
			}
			if(freqHighReg[i].end != totalKmer)
			{
				freqHighReg[i].end -= edgeCutLen;
				if(freqHighReg[i].end < 1)
					freqHighReg[i].end += edgeCutLen;
			}
		}
	}

	mvnvd(4,"Marker 1.4");
	if(freqHighReg.size() >= 2)
	{
		for(uint i(0); i < (freqHighReg.size()-1); ++i)
		{
			if((read.fastCorrectCount+read.fastCorrectCount) >= args.maxChangePerRead)
				break;

			int lengthNeedTrim(0);
			int numCorrect = CorrectBasesRightwardBB(args, read, (freqHighReg[i].end+args.spaKmerSize+args.spaSize), (freqHighReg[i+1].start+args.spaKmerSize+args.spaSize-2), 0, lengthNeedTrim);
			if(!lengthNeedTrim && numCorrect)
			{
				read.bbCorrectCount += numCorrect;
				freqContReg[freqHighReg[i].status-1].status = 1;
				continue;
			}
			else
			{
				lengthNeedTrim = 0;
				numCorrect = CorrectBasesLeftwardBB(args, read, (freqHighReg[i+1].start-1), (freqHighReg[i].end+1), 0, lengthNeedTrim);
				if(!lengthNeedTrim && numCorrect)
				{
					read.bbCorrectCount += numCorrect;
					freqContReg[freqHighReg[i].status-1].status = 1;
					continue;
				}
			}
		}
	}

	mvnvd(4,"Marker 1.5");
	FreqContReg_s freqMaxRegSingle;
	GetMaxHighfreqBlock(freqContReg, freqMaxRegSingle);
	mvnvd(4,"Marker 1.6");

	{
		uint highFreqStart(freqMaxRegSingle.start);
		int numCorrect(0);
		if(highFreqStart > 1)
		{
			numCorrect = CorrectBasesLeftwardBB(args, read, highFreqStart - 1, 1, args.trimAllow, read.trimLeft);
			if(numCorrect)
				read.bbCorrectCount += numCorrect;
		}
	}
	mvnvd(4,"Marker 1.7");
	{
		uint highFreqEnd(freqMaxRegSingle.end + args.spaKmerSize + args.spaSize - 1);
		int numCorrect(0);
		if(highFreqEnd < read.read.size())
		{
			numCorrect = CorrectBasesRightwardBB(args, read, highFreqEnd+1, read.read.size(), args.trimAllow, read.trimRight);
			if(numCorrect)
				read.bbCorrectCount += numCorrect;
		}
	}
	mvnvd(4,"Marker 1.8");

	if((read.trimLeft > args.trimAllow) || (read.trimRight > args.trimAllow))
		read.shouldDelete = 1;
	if((int)( read.read.size() - read.trimLeft - read.trimRight ) < args.minOutputReadLen)
		read.shouldDelete = 1;

	return (read.bbCorrectCount + read.fastCorrectCount);
}


static void FillAltAllele(Read_s & read1, Read_s & read2, const int connected)
{
	read1.altAllele.resize(read1.read.size(), 'N');
	read2.altAllele.resize(read2.read.size(), 'N');
	string::size_type pOrg(read1.read.size() - connected);
	string::size_type p1, p2;
	for(int i(0); i < connected; ++i)
	{
		p1 = pOrg + i;
		p2 = read2.read.size() - 1 - i;
		read1.altAllele[p1] = base2compbase(read2.read[p2]);
		read2.altAllele[p2] = base2compbase(read1.read[p1]);
	}
}


static uint64_t Correct2ReadsFile(const CorrGlobalVar_s & args, const int idx1, const int idx2)
{
	uint64_t inputReadsCount(0), inputBasesCount(0);
	uint64_t outputReadsCount(0), outputBasesCount(0);
	uint64_t trimReadsCount(0), trimBasesCount(0);
	uint64_t deleteReadsCount(0), deleteBasesCount(0);
	uint64_t fastCorrectAccum(0), bbCorrectAccum(0);
	int readLengthCache(0);

	ConnectArgs_s connArgs;
	if(args.joinRead || args.corrByOvl)
		ConnectInit(connArgs);
	igzstream I1, I2;
	OpenFileGZ(I1, args.readsFNs[idx1].c_str(), "rc");
	mvnv(1, "Processing file: %s", args.readsFNs[idx1].c_str());
	if(idx2 >= 0)
	{
		OpenFileGZ(I2, args.readsFNs[idx2].c_str(), "rc");
		mvnv(1, "Processing file: %s", args.readsFNs[idx2].c_str());
	}

	ogzstream fail1O, fail2O, new1O, new2O, connO;
	if(args.toDelete)
	{
		OpenFileGZ(fail1O, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx1] + ".fail.gz")).c_str(), "wc-");
		if(idx2 >= 0)
			OpenFileGZ(fail2O, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx2] + ".fail.gz")).c_str(), "wc-");
	}
	OpenFileGZ(new1O, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx1] + ".new.gz")).c_str(), "wc-");
	if(idx2 >= 0)
		OpenFileGZ(new2O, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx2] + ".new.gz")).c_str(), "wc-");
	if(args.joinRead)
		OpenFileGZ(connO, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx1] + ".conn.gz")).c_str(), "wc-");

	int readType(0);

//z for the switch of head reading
#define PRECESS1READ(x,y,z)\
{\
	if(readType == FQ)\
	{\
		if(z)\
			getline(y, x.head);\
		getline(y, x.read);\
		getline(y, x.qualHead);\
		getline(y, x.qual);\
	}\
	else if(readType == FA)\
	{\
		if(z)\
			getline(y, x.head);\
		getline(y, x.read);\
		x.qualHead = "+";\
		x.qual = string(x.read.size(), (char)(args.conn_qual_base + QUALMAX));\
	}\
	else\
		die("Should never reach here");\
	++inputReadsCount;\
	inputBasesCount += x.read.size();\
}
	//
	string tmpStr;
	while(true)
	{
		Read_s read1, read2;

		getline(I1, read1.head);
		if(read1.head.empty())
			break;
		if(!readType)
		{
			if(read1.head[0] == '@')
				readType = FQ;
			else if(read1.head[0] == '>')
				readType = FA;
			else
			{
				mvnv(0, "Irrecognized file type marker: '%c'", read1.head[0]);
				return 0;
			}
			mvnv(3, "Determined file type: %s", ((readType==FQ)?"FASTQ":"FASTA"));
		}

		PRECESS1READ(read1, I1, false);
		if(idx2 >= 0)
			PRECESS1READ(read2, I2, true);
#undef PRECESS1READ

		mvnvd(4,"Marker 1");

		int connected(0);
		if(idx2 >= 0 && (args.corrByOvl || args.joinRead))
		{
			Set1Read4Connect(connArgs, read1.head, read1.read, read1.qualHead, read1.qual, read2.head, read2.read, read2.qualHead, read2.qual);
			RevCompl(connArgs.read_b, connArgs.rread_b);
			connected = alignConnect(connArgs, connArgs.read_a, connArgs.rread_b, 0);
			mvnvd(4, "alignConnect: %d", connected);
		}
		mvnvd(4,"Marker 2");
		if(idx2 >= 0 && args.corrByOvl && connected >= connect::Min)
		{
			FillAltAllele(read1, read2, connected);
			mvnvd(4, "Filled AltAllele read1: %s", read1.altAllele.c_str());
			mvnvd(4, "Filled AltAllele read2: %s", read2.altAllele.c_str());
		}
		mvnvd(4, "Start correcting 1 read: %s", read1.head.c_str());
		mvnvd(4,"Marker 3");
		Correct1Read(args, read1);
		fastCorrectAccum += read1.fastCorrectCount;
		bbCorrectAccum += read1.bbCorrectCount;
		if(idx2 >= 0)
		{
			Correct1Read(args, read2);
			fastCorrectAccum += read2.fastCorrectCount;
			bbCorrectAccum += read2.bbCorrectCount;
		}
		mvnvd(4,"Marker 4");

		readLengthCache = read1.read.size();
		//Trim read
#define PROCESSTRIM(x) \
{\
	if(x.trimLeft || x.trimRight)\
	{\
		mvnvd(4, "Start trimming: %d/%d", x.trimLeft, x.trimRight);\
		++trimReadsCount;\
		trimBasesCount += (x.trimLeft + x.trimRight);\
		x.read = x.read.substr(x.trimLeft, x.read.size() - x.trimLeft - x.trimRight);\
		if(readType == FQ)\
			x.qual = x.qual.substr(x.trimLeft, x.qual.size() - x.trimLeft - x.trimRight);\
	}\
}
		//
		PROCESSTRIM(read1);
		if(idx2 >= 0)
			PROCESSTRIM(read2);
#undef PROCESSTRIM
		mvnvd(4,"Marker 5");

		//Process the header
#define ADDHEADER(x)\
{\
	stringstream ss;\
	ss<<x.head<<" ";\
	for(uint k(0); k < x.fastCorrectRecord.size(); ++k)\
		if((x.fastCorrectRecord[k] <= x.trimLeft) || (x.fastCorrectRecord[k] > (readLengthCache - x.trimRight)))\
			--x.fastCorrectCount;\
	for(uint k(0); k < x.bbCorrectRecord.size(); ++k)\
		if((x.bbCorrectRecord[k] <= x.trimLeft) || (x.bbCorrectRecord[k] > (readLengthCache - x.trimRight)))\
			--x.bbCorrectCount;\
	ss<<x.fastCorrectCount<<" "<<x.bbCorrectCount<<" "<<x.trimLeft<<" "<<x.trimRight<<" "<<x.shouldDelete<<" f:";\
	for(uint k(0); k < x.fastCorrectRecord.size(); ++k)\
	{\
		if(!((x.fastCorrectRecord[k] <= x.trimLeft) || (x.fastCorrectRecord[k] > (readLengthCache - x.trimRight))))\
		{\
			if(k)\
				ss<<"/";\
			ss<<x.fastCorrectRecord[k]<<x.fastCorrectCharRec[k];\
		}\
	}\
	ss<<" b:";\
	for(uint k(0); k < x.bbCorrectRecord.size(); ++k)\
	{\
		if(!((x.bbCorrectRecord[k] <= x.trimLeft) || (x.bbCorrectRecord[k] > (readLengthCache - x.trimRight))))\
		{\
			if(k)\
				ss<<"/";\
			ss<<x.bbCorrectRecord[k]<<x.bbCorrectCharRec[k];\
		}\
	}\
	if(args.commentOn)\
		ss<<" "<<x.comment;\
	x.head = ss.str();\
}
		//
		ADDHEADER(read1);
		mvnvd(4, "Processed header1: %s", read1.head.c_str());
		if(idx2 >= 0)
		{
			ADDHEADER(read2);
			mvnvd(4, "Processed header2: %s", read2.head.c_str());
		}
#undef ADDHEADER
		mvnvd(4,"Marker 6");

		//Output reads
#define PRTREAD(x,y) \
{\
	++outputReadsCount;\
	outputBasesCount += y.read.size();\
	if(readType == FQ)\
		x<<y.head<<"\n"<<y.read<<"\n"<<y.qualHead<<"\n"<<y.qual<<"\n";\
	else\
		x<<y.head<<"\n"<<y.read<<"\n";\
}
		//
		if(args.toDelete)
		{
			if(read1.shouldDelete || read2.shouldDelete)
			{
				++deleteReadsCount;
				deleteBasesCount += read1.read.size();
				PRTREAD(fail1O, read1);
				if(idx2 >= 0)
				{
					++deleteReadsCount;
					deleteBasesCount += read2.read.size();
					PRTREAD(fail2O, read2);
				}
			}
			else
			{
				PRTREAD(new1O, read1);
				if(idx2 >= 0)
					PRTREAD(new2O, read2);
			}
		}
		else
		{
			PRTREAD(new1O, read1);
			if(idx2 >= 0)
				PRTREAD(new2O, read2);
		}
#undef PRTREAD

		mvnvd(4,"Marker 7");
		//Join the read if read1 and read2 was not marked deleted and join read mode enabled
		if(idx2 >= 0 && args.joinRead && !read1.shouldDelete && !read2.shouldDelete)
		{
			++connArgs.total_pairs;
			if(read1.fastCorrectCount || read1.bbCorrectCount || read2.fastCorrectCount || read2.bbCorrectCount)
			{
				RevCompl(connArgs.read_b, connArgs.rread_b);
				connected = alignConnect(connArgs, connArgs.read_a, connArgs.rread_b, 1);
			}
			Analyze4Correction(connArgs, connected, connO);
		}
		mvnvd(4,"Marker 8");
	}

	//Output Statistics
	fstream statO;
	OpenFile(statO, RegisterFile(GetCWD() + "/" + GetFilenamePart(args.readsFNs[idx1]) + ".stat").c_str(), "wc-");
	{
		mvnvd(4, "Outputing stats");
		char str[4096];
		str[0] = '\0';
		sprintf(str, "%sInput_read_count\t%lu\n", str, inputReadsCount);
		sprintf(str, "%sInput_base_count\t%lu\n", str, inputBasesCount);
		sprintf(str, "%sOutput_read_count\t%lu\n", str, outputReadsCount);
		sprintf(str, "%sOutput_base_count\t%lu\n", str, outputBasesCount);
		sprintf(str, "%sTrim_read_count\t%lu\n", str, trimReadsCount);
		sprintf(str, "%sTrim_base_count\t%lu\n", str, trimBasesCount);
		sprintf(str, "%sDelete_read_count\t%lu\n", str, deleteReadsCount);
		sprintf(str, "%sDelete_base_count\t%lu\n", str, deleteBasesCount);
		sprintf(str, "%sCorrect_by_fast\t%lu\n", str, fastCorrectAccum);
		sprintf(str, "%sCorrect_by_bb\t%lu\n", str, bbCorrectAccum);
		sprintf(str, "%sOutput_base_percentage\t%f\n", str, ((double)outputBasesCount / inputBasesCount) * 100);
		sprintf(str, "%sTrim_base_percentage\t%f\n", str, ((double)trimBasesCount / inputBasesCount) * 100);
		sprintf(str, "%sDelete_base_percentage\t%f\n", str, ((double)deleteBasesCount / inputBasesCount) * 100);
		sprintf(str, "%sOutput_read_percentage\t%f\n", str, ((double)outputReadsCount / inputReadsCount) * 100);
		sprintf(str, "%sTrim_read_percentage\t%f\n", str, ((double)trimReadsCount / inputReadsCount) * 100);
		sprintf(str, "%sDelete_read_percentage\t%f\n", str, ((double)deleteReadsCount / inputReadsCount) * 100);
		sprintf(str, "%sCorrect_fast_percentage\t%f\n", str, ((double)fastCorrectAccum / inputBasesCount) * 100);
		sprintf(str, "%sCorrect_bb_percentage\t%f\n", str, ((double)bbCorrectAccum / inputBasesCount) * 100);
		sprintf(str, "%sCorrect_percentage\t%f\n", str, (((double)fastCorrectAccum + (double)bbCorrectAccum) / inputBasesCount) * 100);
		statO << str << "\n" << flush;
	}
	if(idx2 >= 0 && args.joinRead)
		OutputConnectStat(connArgs, statO);

	if(args.joinRead || args.corrByOvl)
		ConnectClean(connArgs);

	return outputReadsCount;
}


static void usage(char ** argv, const CorrGlobalVar_s & args)
{
	cerr<<boolalpha
	<<argv[0]<<" "<<argv[1]<<" [PARAMETERS] <Consecutive Kmer Frequency Table> <Space Kmer Frequency Table> <Reads Filenames List>"<<'\n'
	<<"    -k <int>     Consecutive Kmer Size, Default: "<< args.conKmerSize << '\n'
	<<"    -K <int>     Space Kmer Solid Part Size, Default: "<< args.spaKmerSize << '\n'
	<<"    -S <int>     Space Kmer Space Part Size, Default: "<< args.spaSize << '\n'
	<<"    -l <int>     Consecutive Kmer Frequency lower limit, Default: "<< args.conLowFreqCutoff << '\n'
	<<"    -e <int>     Consecutive Kmer Frequency upper limit, Default: "<< args.conHighFreqCutoff << '\n'
	<<"    -L <int>     Space Kmer Frequency lower limit, Default: "<< args.spaLowFreqCutoff << '\n'
	<<"    -E <int>     Space Kmer Frequency upper limit, Default: "<< args.spaHighFreqCutoff << '\n'
	<<"    -m <int>     Minimum consecutive high frequency kmers per high-freq block, Default: "<< args.minHighFreqLen << '\n'
	<<"    -c <int>     Maximum base alternation per read, Default: "<< args.maxChangePerRead << '\n'
	<<"    -t <int>     Maximum trimmed bases allowed on each side, Default: "<< args.trimAllow << '\n'
	<<"    -r <int>     Minimum read length after correction, Default: "<< args.minOutputReadLen << '\n'
	<<"    -I <int>     Minimum base interval between two correction, Default: "<< args.minInterval << '\n'
	<<"    -x           Delete correction failed reads, Default: "<< (bool)args.toDelete << '\n'
	<<"    -j           Connect reads after correction, use 2 times to output quality, Default: "<< (bool)args.joinRead << '\n'
	<<"    -y           Correct reads end by overlapping information, Default: "<< (bool)args.corrByOvl << '\n'
	<<"    -v           Increase Verbosity, 3 times max, Default: "<< (int)log2(verbosity) << '\n'
	<<"    -F           Append additional information to read ID (EXPERIMENTAL), Default: "<< (bool)args.commentOn << '\n'
	<<"    -h           This help"<< '\n'
	<<'\n'
	<<"Reads overlapping related parameters:"<< '\n'
    <<"    -q  <int>    Quality ASCII base (Default: 33)\n"
    <<"    -Q  <int>    Quality start range offset (Default: 2)\n"
	<<"    -a <int>     Minimal align length allowed, Default: "<< args.conn_Min << '\n'
	<<"    -A <int>     Maximal align length allowed, Default: "<< args.conn_Max << '\n'
	<<"    -i <int>     Mean align length induced alignment, Default: "<< args.conn_Mean << '\n'
	<<"    -u <float>   Identity cutoff of aligned region, Default: "<< args.conn_MATCH_CUTOFF << '\n'
	<<"    -d <float>   Maximum simularity between best and second best result, Default: "<< args.conn_DIF << '\n'
	<<'\n';

	exit(EXIT_FAILURE);
}


int Correct(int argc, char ** argv)
{
	int status;
	CorrGlobalVar_s args;
	if(argc < 3)
		usage(argv, args);

	//Parameter occupation: AEFIKLQS acdehijklmqrtuvxy
	while((status = getopt(argc-1, argv+1, "k:K:s:l:L:m:c:t:r:jyvxhq:FQ:a:A:i:u:d:E:e:I:")) != -1)
	{
		switch(status)
		{
			case 'k': args.conKmerSize = atoi(optarg); break;
			case 'K': args.spaKmerSize = atoi(optarg); break;
			case 'S': args.spaSize = atoi(optarg); break;
			case 'l': args.conLowFreqCutoff = atoi(optarg); break;
			case 'L': args.spaLowFreqCutoff = atoi(optarg); break;
			case 'e': args.conHighFreqCutoff = atoi(optarg); break;
			case 'E': args.spaHighFreqCutoff = atoi(optarg); break;
			case 'm': args.minHighFreqLen = atoi(optarg); break;
			case 'c': args.maxChangePerRead = atoi(optarg); break;
			case 't': args.trimAllow = atoi(optarg); break;
			case 'r': args.minOutputReadLen = atoi(optarg); break;
			case 'I': args.minInterval = atoi(optarg); break;
			case 'j': if(args.joinRead){++connect::QualityOutput;}; args.joinRead = 1; break;
			case 'y': args.corrByOvl = 1; break;
			case 'x': args.toDelete = 1; break;
			case 'F': args.commentOn = 1; break;
			case 'q': args.conn_qual_base = atoi (optarg); break;
			case 'Q': args.conn_qual_offset = atoi (optarg); break;
			case 'a': args.conn_Min = atoi(optarg); break;
			case 'A': args.conn_Max = atoi(optarg); break;
			case 'i': args.conn_Mean = atoi(optarg); break;
			case 'u': args.conn_MATCH_CUTOFF = atof(optarg); break;
			case 'd': args.conn_DIF = atof(optarg); break;
			case 'v': ModifyVerbosity(); break;
			case 'h': usage(argv, args); break;
			default: usage(argv, args); break;
		}
	}
	if(optind < (argc - 1))
	{
		args.conKmerTableFN = (argv+1)[optind++];
		args.spaKmerTableFN = (argv+1)[optind++];
		args.readsListFN = (argv+1)[optind++];
	}
	else
		usage(argv, args);

	if( (args.spaKmerSize + args.spaSize) % 2 == 0 )
		die("Space kmer length is not even number!");
	if( ((args.spaKmerSize / 2) + args.spaSize) > args.conKmerSize )
		die("Consecutive kmer cannot cover the space of space kmer");

	ImportToVector(args.readsListFN.c_str(), args.readsFNs, '\n', true);
	args.print("User defined parameters:");
	mvnv(1, "Start");

	uint64_t conKmerAmount, spaKmerAmount;
	args.conKmerMasker = ((0x2llu << (args.conKmerSize * 2 - 1)) - 1);
	args.spaHeadOffset = (args.spaKmerSize / 2);
	args.spaKmerTailMasker = ((0x2llu << ((args.spaKmerSize/2) * 2 - 1)) - 1 - 3);
	args.spaKmerHeadMasker = (((0x2llu << ((args.spaKmerSize - (args.spaKmerSize/2)) * 2 - 1)) - 1 - 3) << (args.spaHeadOffset * 2));
	{
		mvnv(1, "Importing consecutive kmer table...");
		args.conKmerTable = MapKmerfreq2Mem(args.conKmerTableFN, args.conKmerSize, conKmerAmount, args.conLowFreqCutoff, args.conHighFreqCutoff, args.conInvalidKmerTable);
//		RecordOneTime();
//		uint64_t theoryTotal(0), count(0), effectiveTotal(0); double lowFreqPercent(0.);
//		args.conKmerTable = ReadKmerfreq2Mem(args.conKmerTableFN, args.conKmerSize, theoryTotal, count, effectiveTotal, args.conLowFreqCutoff, args.conHighFreqCutoff, args.conInvalidKmerTable, lowFreqPercent);
//		conKmerAmount = theoryTotal;
//		mvnv(1, " Table theory maximum: %llu", theoryTotal);
//		mvnv(1, "     Total kmer count: %llu", count);
//		mvnv(1, "      Total kmer hits: %llu", effectiveTotal);
//		mvnv(1, "Low Freq Kmer Percent: %f", lowFreqPercent*100);
//		PrintDuration();
	}
	{
		mvnv(1, "Importing space kmer table...");
		args.spaKmerTable = MapKmerfreq2Mem(args.spaKmerTableFN, args.spaKmerSize, spaKmerAmount, args.spaLowFreqCutoff, args.spaHighFreqCutoff, args.spaInvalidKmerTable);
//		RecordOneTime();
//		uint64_t theoryTotal(0), count(0), effectiveTotal(0); double lowFreqPercent(0.);
//		args.spaKmerTable = ReadKmerfreq2Mem(args.spaKmerTableFN, args.spaKmerSize, theoryTotal, count, effectiveTotal, args.spaLowFreqCutoff, args.spaHighFreqCutoff, args.spaInvalidKmerTable, lowFreqPercent);
//		spaKmerAmount = theoryTotal;
//		mvnv(1, " Table theory maximum: %llu", theoryTotal);
//		mvnv(1, "     Total kmer count: %llu", count);
//		mvnv(1, "      Total kmer hits: %llu", effectiveTotal);
//		mvnv(1, "Low Freq Kmer Percent: %f", lowFreqPercent*100);
//		PrintDuration();
	}

	if(args.corrByOvl || args.joinRead)
	{
		if(args.readsFNs.size() % 2 != 0)
			die("File list %s include odd number of files.", args.readsListFN.c_str());
		for(uint i(0); i < args.readsFNs.size(); i+=2)
			Correct2ReadsFile(args, i, i+1);
	}
	else
	{
		for(uint i(0); i < args.readsFNs.size(); ++i)
			Correct2ReadsFile(args, i, -1);
	}

	mvnv(1, "End correcting reads");
	munmap(args.conKmerTable, sizeof(uint8_t) * (conKmerAmount/8+1));
	munmap(args.spaKmerTable, sizeof(uint8_t) * (spaKmerAmount/8+1));
	if(args.conInvalidKmerTable != NULL)
		munmap(args.conInvalidKmerTable, sizeof(uint8_t) * (conKmerAmount/8+1));
	if(args.spaInvalidKmerTable != NULL)
		munmap(args.spaInvalidKmerTable, sizeof(uint8_t) * (spaKmerAmount/8+1));
	mvnv(1, "Memory released");

	return status = 0;
}
