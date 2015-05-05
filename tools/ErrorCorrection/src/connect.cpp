#include "connect.h"

//#define DEBUG
//#define DEBUG_VERBOSE

//#define QD (cerr<< __FUNCTION__ << "\t" << __LINE__ << "\n")
#define QD

using namespace std;
using namespace connect;

//Global define
char cwd[1024];

double connect::MATCH_CUTOFF(0.9), connect::DIF(0.7), connect::B_CUTOFF(0.9);
int connect::Min(8), connect::Max(70), connect::Mean(0), connect::Nfilter(-1), connect::QualityOutput(0);
int connect::qual_base(33), connect::qual_offset(2);
int connect::tHead(0), connect::tTail(0);

long long Pairs(0);
string a_file, b_file, c_file, out_connect_file, fq1, fq2;
igzstream fin_a, fin_b, fin_c;
ofstream fout, fout_fq1, fout_fq2;
//Global define end

void Set1Read4Connect(ConnectArgs_s & args, const string & a1, const string & a2, const string & a3, const string & a4,\
											const string & b1, const string & b2, const string & b3, const string & b4)
{
	args.read_head_a = a1;
	args.a_id = a1.substr(1);
	args.read_a = a2;
	args.a_s = a3;
	args.a_q = a4;
	args.read_head_b = b1;
	args.b_id = b1.substr(1);
	args.read_b = b2;
	args.b_s = b3;
	args.b_q = b4;
}

void ConnectInit(ConnectArgs_s & args)
{
	int array_size = Max - Min + 1;
	args.Mismatch_array = new int[array_size];

	//Initialize mismatch penalty array
	for ( int i(0); i < array_size; ++i )
		args.Mismatch_array[i] = int(floor((double ( Min + i ) * double ( 1 - MATCH_CUTOFF * DIF )) + 0.5));

	//Initialize mismatch distribution
	args.Max_mismatch_num = int(floor((double ( Max ) * double( 1 - MATCH_CUTOFF )) + 0.5));
	args.mismatch_num_distr = new uint64_t[args.Max_mismatch_num + 1];

	bzero((void*)args.mismatch_num_distr, sizeof(uint64_t) * (args.Max_mismatch_num + 1));
}

void ConnectClean(ConnectArgs_s & args)
{
	delete[] args.Mismatch_array;
	delete[] args.mismatch_num_distr;
}

void usage (char ** argv)
{
	cerr << argv[0]<<" "<<argv[1]<<" [PARAMETERS]"
		 << "\n"
		 << "Options:\n"
		 << "    -a  <str>   Query a file, *.fq, *.fa\n"
		 << "    -b  <str>   Query b file, the same as a\n"
		 << "    -o  <str>   Output connected file in *.fa\n"
		 << "    -2  <str>   Output fail connected read1.fq\n"
		 << "    -3  <str>   Output fail connected read2.fq\n"
		 << "    -r  <int>   Trim x bp at head before use [" << tHead <<"]\n"
		 << "    -R  <int>   Trim x bp at tail before use [" << tTail << "]\n"
		 << "    -l  <int>   Minimal align length allowed [" << Min << "]\n"
		 << "    -m  <int>   Mean align length induced alignment (Default: No) [" << Mean << "]\n"
		 << "    -u  <int>   Maximal align length allowed [" << Max << "]\n"
		 << "    -c  <float> Identity cutoff of aligned region [" << MATCH_CUTOFF << "]\n"
		 << "    -d  <float> Maximum simularity between best and second best result. [" << DIF << "]\n"
		 << "    -B  <float> b threshold (Maximum %/100 of qual '0' base per alignment) [" << B_CUTOFF << "]\n"
		 << "    -t  <int>   Pairs to connect (Default: All) [" << Pairs << "]\n"
		 << "    -N  <int>   Maximum 'N' allowance for filtering (Default: Do not filter) [" << Nfilter << "]\n"
		 << "    -q          Output connecting quality [" << QualityOutput << "]\n"
		 << "    -h          Help\n"
		 << "\n"
		 << "Advance Options (Use with caution or consult your parents)\n"
		 << "    -x  <int>   Quality ASCII base (Default: 33)\n"
		 << "    -X  <int>   Quality start range offset (Default: 2)\n"
		 << "\n"
		 << "Example:\n"
		 << "    " << cwd << argv[0] << " -a in_1.fq -b in_2.fq -o out.fa -l 8 -u 70 -2 fail_1.fq -3 fail_2.fq\n"
		 << "\n"
		 << flush;
	exit ( 1 );
}

int getOptions ( int argc, char ** argv )
{
	int c;
	
	while((c = getopt(argc, argv, "a:b:e:o:2:3:l:m:u:c:d:B:t:nN:qx:X:r:R:h")) != -1)
	{
		switch ( c )
		{
			case 'a':
				a_file = optarg;
				break;
			case 'b':
				b_file = optarg;
				break;
			case 'e':
				c_file = optarg;
				break;
			case 'o':
				out_connect_file = optarg;
				break;
			case '2':
				fq1 = optarg;
				break;
			case '3':
				fq2 = optarg;
				break;
			case 'l':
				Min = atoi ( optarg );
				break;
			case 'm':
				Mean = atoi ( optarg );
				break;
			case 'u':
				Max = atoi ( optarg );
				break;
			case 'c':
				MATCH_CUTOFF = atof ( optarg );
				break;
			case 'd':
				DIF = atof ( optarg );
				break;
			case 'B':
				B_CUTOFF = atof ( optarg );
				break;
			case 't':
				Pairs = atoi ( optarg );
				break;
			case 'N':
				Nfilter = atoi ( optarg );
				break;
			case 'r':
				tHead = atoi ( optarg );
				break;
			case 'R':
				tTail = atoi ( optarg );
				break;
			case 'q':
				QualityOutput = 1;
				break;
			case 'x':
				qual_base = atoi ( optarg );
				break;
			case 'X':
				qual_offset = atoi ( optarg );
				break;
			case 'h':
				usage ( argv );
				break;
			default:
				return (optind - 1);
				break;
		}
	}

	if ( Max <= Min )
	{
		cerr << "Error Max or Min value!\n";
		exit ( 1 );
	}
	
	return argc;
};

#define base2int(base)			(char)(((base)&0x06)>>1)
#define int2base(seq)			"ACTG"[seq]			//0123
#define int2compbase(seq)		"TGAC"[seq]			//2301
#define int_comp(seq)			(char)(seq^0x02)	//(char)((0x4E>>((seq)<<1))&0x03)
//0123<-->2301
/*
ASCII:
A	00	1000001
C	01	1000011
T	10	1010100
G	11	1000111
			||
0x06	0000110
*/

inline void revCompl ( const string & seq, string & nseq )
{
	nseq.resize(seq.size());
	
	for ( int i(seq.size() - 1), j(0); i >= 0; --i, ++j )
		nseq[j] = int2compbase((int)base2int(seq[i]));
}

int alignConnect ( ConnectArgs_s & args, const string & read_a, const string & read_b, const int resubmisson)
{
	int match(0);

	double value[Max - Min + 1];
	bzero((void *) &value, sizeof(double) * (Max - Min + 1));

	int len = read_a.size();

	if ( len <= Min )
	{
		cerr << "Read too short: " << read_a.size() << " " << Min << endl << read_a << endl;
		return 0;
	}

	{
		string sa;
		string sb;
		int sum;
		int mismatch_sum;

		for ( int i(Min); i <= Max; ++i )
		{
			value[i - Min] = 0;

			if ( i > (int)read_a.length() || i > (int)read_b.length() ) {break;}

			sa = read_a.substr ( len - i );
			sb = read_b.substr ( 0, i );
			string::size_type a_len = sa.size();
			sum = 0;
			mismatch_sum = 0;

#ifdef	DEBUG
			if ( a_len != sb.size() ) {cerr << "Length error!" << endl; exit ( 1 );}
#endif 

			for ( uint j(0); j < a_len; ++j )
			{
				if ( sa[j] == 'N' || sa[j] != sb[j]) {++mismatch_sum;}
				else {++sum;}

				if ( mismatch_sum > args.Mismatch_array[i - Min] ) {sum = 0; break;}
			}

			if(sum)
				value[i - Min] = double ( double ( sum ) / double ( a_len ) * 100.0 );
			else
				value[i - Min] = 0.0;
		}
	}

	//Find max match
	double max_value(0.), max2_value(0.);
	int max_len(0), max2_len(0);
	

	for ( int i(0); i < Max - Min + 1; ++i )
	{
		if ( value[i] == 0.0 ) {continue;}

		if ( max_value <= value[i] )
		{
			max2_value = max_value;
			max2_len = max_len;
			max_value = value[i];
			max_len = Min + i;
			continue;
		}

		if ( max2_value < value[i] && max_value > value[i] )
		{
			max2_value = value[i];
			max2_len = Min + i;
		}
	}

	//Determine connect or not
	if ( max_value >= MATCH_CUTOFF * 100.0 )
	{
		if ( double ( max2_value ) < double ( max_value ) * DIF ) {return max_len;}
		else if(!resubmisson) {++args.confuse_pairs;}

		//Mean overlap size induced connection
		if ( (Mean > 0) && (max2_value > (MATCH_CUTOFF * 100.0)))
		{
			match = ( abs ( max_len - Mean ) > abs ( max2_len - Mean ) ) ? max2_len : max_len;
#ifdef DEBUG_VERBOSE
			cerr<<"Mean overlap size induced connection:\n "<<max_value<<"\t"<<max2_value<<"\t"<<match<<"\t"<<value[match-Min]<<"\n";
#endif
			return match;
		}
	}
	else if(!resubmisson) {++args.low_match_rate_pairs;}

	//Use the longer one as overlap length, negative to mark failure.
	match = ( max_len > max2_len ) ? -max_len : -max2_len;
	return match;
}

template <typename N>
void chkOpen(const string & fn, N & fh)
{
	if (!fh)
	{
		cerr << "Failed opening: " << cwd << fn << endl;
		exit ( 1 );
	}
}

int checkFile()
{
	fout.open ( out_connect_file.c_str() );
	chkOpen<ofstream>(out_connect_file, fout);
	
	fout_fq1.open ( fq1.c_str() );
	chkOpen<ofstream>(fq1, fout_fq1);
	
	fout_fq2.open ( fq2.c_str() );
	chkOpen<ofstream>(fq2, fout_fq2);
	
	if ( ( !a_file.empty() ) && ( !b_file.empty() ) && ( c_file.empty() ) )
	{
		cerr << "Pair-end connect:\n"
			 << "         Query 1: " << cwd << ((a_file.find("/")==string::npos)?(a_file):(a_file.substr(a_file.rfind("/")+1, string::npos))) << "\n"
			 << "         Query 2: " << cwd << ((b_file.find("/")==string::npos)?(b_file):(b_file.substr(b_file.rfind("/")+1, string::npos))) << "\n"
			 << "Output Connected: " << cwd << out_connect_file << "\n" 
			 << "   Output Fail 1: " << cwd << fq1 << "\n" 
			 << "   Output Fail 2: " << cwd << fq2 << "\n";
		
		fin_a.open ( a_file.c_str() );
		chkOpen<igzstream>(a_file, fin_a);

		fin_b.open ( b_file.c_str() );
		chkOpen<igzstream>(b_file, fin_b);
		
		return 1;
	}
	else
	{
		if ( !c_file.empty() && ( a_file.empty() ) && ( b_file.empty() ) )
		{
			cout << "Pair-end connect:\n";
			cout << "           Query: " << c_file << "\n"
				 << "Output Connected: " << out_connect_file << " \n" 
				 << "   Output Fail 1: " << fq1 << "\n" 
				 << "   Output Fail 2: " << fq2 << "\n" << flush;

			fin_c.open ( c_file.c_str() );
			chkOpen<igzstream>(c_file, fin_c);

			return 2;
		}
		else
		{
			cerr << "Input or output parameters error.\n";
			exit ( 1 );
		}
	}
}

//Filter by N.
inline int NFilter ( ConnectArgs_s & args )
{
	if ( count(args.read_a.begin(), args.read_a.end(), 'N') > Nfilter || count(args.read_b.begin(), args.read_b.end(), 'N') > Nfilter )
	{
		++args.low_quality_pairs;

		if ( !args.a_q.empty() )
		{
			fout_fq1 
				<< "@" << args.a_id << "\n"
				<< args.read_a << "\n"
				<< args.a_s << "\n"
				<< args.a_q << "\n";
			
			fout_fq2
				<< "@" << args.b_id << "\n"
				<< args.read_b << "\n"
				<< args.b_s << "\n"
				<< args.b_q << "\n";
		}
		else
		{
			fout_fq1
				<< ">" << args.a_id << "\n"
				<< args.read_a << "\n";
			
			fout_fq2
				<< ">" << args.b_id << "\n"
				<< args.read_b << "\n";
		}

		return 1;
	}

	return 0;
}

int Analyze4Correction (ConnectArgs_s & args, int connected, ogzstream & outFH )
{
	string seq;

	if ( connected == 0 )
		connected = Min - 1;

	if ( connected < Min )
	{
		connected = abs ( connected );
		++args.unconnect_pairs;
	}
	else
	{
		++args.connect_pairs;
		int mismatch_num(0);

		string overlap_part ( abs ( connected ), 'N' );
		string overlap_quality_part ( abs ( connected ), '-' );

		if ( !args.a_q.empty() )
		{
			{
				string::size_type pOrg(args.read_a.size() - connected);
				string::size_type p;
				for ( int i(0); i < connected; ++i )
				{
					p = pOrg + i;

					if ( args.read_a[p] != args.rread_b[i] )
					{
						++mismatch_num;
					}

					if ( args.a_q[p] >= args.b_q[args.read_a.size() - i - 1 ] )
					{
						overlap_part[i] = args.read_a[p];
						overlap_quality_part[i] = args.a_q[p];
					}
					else
					{
						overlap_part[i] = args.rread_b[i];
						overlap_quality_part[i] = args.b_q[args.read_a.size() - i - 1 ];
					}
				}
			}

			++args.mismatch_num_distr[mismatch_num];
			seq = args.read_a.substr ( 0, args.read_a.size() - connected ) + overlap_part + args.rread_b.substr ( connected );

			if ( QualityOutput )
			{
				string b_q_left = args.b_q.substr ( 0, args.read_a.size() - connected );
				string rev_b_q_left(b_q_left.size(), '\0');

				for ( int i(b_q_left.size() - 1), j(0); i >= 0; --i, ++j )
					rev_b_q_left[j] = b_q_left[i];

				string quality = args.a_q.substr ( 0, args.read_a.size() - connected ) + overlap_quality_part + rev_b_q_left;

				outFH << "@" << args.connect_pairs << "-'" << args.a_id << "'-'" << args.b_id << "'-" << connected << "-" << args.read_b.size() + args.read_a.size() - connected << "\n"
					 << seq << "\n"
					 <<"+\n"
					 << quality << "\n";
			}
			else
			{
				outFH << ">" << args.connect_pairs << "-'" << args.a_id << "'-'" << args.b_id << "'-" << connected << "-" << args.read_b.size() + args.read_a.size() - connected << "\n"
					 << seq << "\n";
			}
		}
		else
		{
			string::size_type pOrg = args.read_a.size() - connected;
			for ( int i(0); i < connected; ++i )
				if ( args.read_a[pOrg + i] != args.rread_b[i] )
					++mismatch_num;


			++args.mismatch_num_distr[mismatch_num];
			seq = args.read_a + args.rread_b.substr ( connected );
			outFH << ">" << args.connect_pairs << "\t" << args.a_id << "_" << args.b_id << "\t" << connected << "\t" << args.read_b.size() + args.read_a.size() - connected << "\n"
				 << seq << "\n";
		}

		++args.InsertSize_distr[args.read_b.size() + args.read_a.size() - connected];

	}

	return 0;
}

//Analysis connect result.
inline int analysisResult ( ConnectArgs_s & args, int connected )
{
	string seq;

	if ( connected == 0 )
		connected = Min - 1;

	if ( connected < Min )
	{
		connected = abs ( connected );
		++args.unconnect_pairs;
		
		if ( !args.a_q.empty() )
		{
			fout_fq1 
				<< "@" << args.a_id << "\n"
				<< args.read_a << "\n"
				<<"+\n"
				<< args.a_q << "\n";
			
			fout_fq2
				<< "@" << args.b_id << "\n"
				<< args.read_b << "\n"
				<<"+\n"
				<< args.b_q << "\n";
		}
		else
		{
			fout_fq1
				<< ">" << args.a_id << "\n"
				<< args.read_a << "\n";
			
			fout_fq2
				<< ">" << args.b_id << "\n"
				<< args.read_b << "\n";
		}
	}
	else
	{
		++args.connect_pairs;
		int mismatch_num(0);

		string overlap_part ( abs ( connected ), 'N' );
		string overlap_quality_part ( abs ( connected ), '-' );

		if ( !args.a_q.empty() )
		{
			{
				string::size_type pOrg(args.read_a.size() - connected);
				string::size_type p;
				for ( int i(0); i < connected; ++i )
				{
					p = pOrg + i;

					if ( args.read_a[p] != args.rread_b[i] )
					{
						++mismatch_num;
					}

					if ( args.a_q[p] >= args.b_q[args.read_a.size() - i - 1 ] )
					{
						overlap_part[i] = args.read_a[p];
						overlap_quality_part[i] = args.a_q[p];
					}
					else
					{
						overlap_part[i] = args.rread_b[i];
						overlap_quality_part[i] = args.b_q[args.read_a.size() - i - 1 ];
					}

#ifdef DEBUG_VERBOSE
					cout<<i<<"\t"<<p<<"\t"<<args.read_a[p]<<"\t"<<args.a_q[p]<<"\t"<<args.rread_b[i]<<"\t"<<args.read_b[args.read_a.size() -i - 1]<<"\t"<<args.b_q[args.read_a.size()- i - 1]<<"\t"<<overlap_part[i]<<endl;
#endif
				}
			}

			QD;

#ifdef DEBUG_VERBOSE
			cerr << mismatch_num << "\n";
#endif
			++args.mismatch_num_distr[mismatch_num];
			QD;
			seq = args.read_a.substr ( 0, args.read_a.size() - connected ) + overlap_part + args.rread_b.substr ( connected );

			QD;

			if ( QualityOutput )
			{
				string b_q_left = args.b_q.substr ( 0, args.read_a.size() - connected );
				string rev_b_q_left(b_q_left.size(), '\0');

				for ( int i(b_q_left.size() - 1), j(0); i >= 0; --i, ++j )
					rev_b_q_left[j] = b_q_left[i];

				string quality = args.a_q.substr ( 0, args.read_a.size() - connected ) + overlap_quality_part + rev_b_q_left;
#ifdef DEBUG_VERBOSE
				cout<<args.read_a.substr(0, args.read_a.size()- connected)<<"\n"<<overlap_part<<"\n"<<args.rread_b.substr(connected)<<endl;
#endif
				fout << "@" << args.connect_pairs << "-" << args.a_id << "-" << args.b_id << "-" << connected << "-" << args.read_b.size() + args.read_a.size() - connected << "\n"
					 << seq << "\n"
					 <<"+\n"
					 << quality << "\n";
			}
			else
			{
				fout << ">" << args.connect_pairs << "\t" << args.a_id << "-" << args.b_id << "\t" << connected << "\t" << args.read_b.size() + args.read_a.size() - connected << "\n"
					 << seq << "\n";
			}
		}
		else
		{
			string::size_type pOrg = args.read_a.size() - connected;
			for ( int i(0); i < connected; ++i )
			{
				if ( args.read_a[pOrg + i] != args.rread_b[i] )
				{
					++mismatch_num;
				}
			}

			++args.mismatch_num_distr[mismatch_num];
			seq = args.read_a + args.rread_b.substr ( connected );
			fout << ">" << args.connect_pairs << "\t" << args.a_id << "_" << args.b_id << "\t" << connected << "\t" << args.read_b.size() + args.read_a.size() - connected << "\n"
				 << seq << "\n";
		}

		++args.InsertSize_distr[args.read_b.size() + args.read_a.size() - connected];

	}

	if ( Pairs > 0 && args.connect_pairs == Pairs ) {cerr << "Connected pairs: " << args.connect_pairs << endl; return 2;}

	return 0;
}

const int FQ(1), FA(2);

void Routine1 (ConnectArgs_s & args)
{
	cerr << "Work routine 1 engaged...\n";

	int connected;
	int iType;
	while ( getline ( fin_a, args.read_head_a) )
	{
		connected = 0;

		//Read id of a read
		if ( args.read_head_a[0] == '@')
		{
			iType = FQ;
		}
		else if(args.read_head_a[0] == '>')
		{
			iType = FA;
		}
		else
		{
			cerr << "Weird record header a : "<< args.read_head_a << "\n";
			exit (1);
		}
		args.a_id = args.read_head_a.substr ( 1 );

		//Read strand and sequence quality
		if ( iType == FQ )              //for *.fq
		{
			getline ( fin_a, args.read_a );
			getline ( fin_a, args.a_s );
			getline ( fin_a, args.a_q );
#define TRIM(x) (x = x.substr(tHead, (x.size() - tHead - tTail)))
			TRIM(args.read_a);
			TRIM(args.a_q);

			getline ( fin_b, args.read_head_b );

			if ( args.read_head_b[0] == '@' )
			{
				args.b_id = args.read_head_b.substr ( 1 );
			}
			else
			{
				cerr << "Weird record header b: "<< args.read_head_b << "\n";
				exit (1);
			}

			getline ( fin_b, args.read_b );
			getline ( fin_b, args.b_s );
			getline ( fin_b, args.b_q );
			TRIM(args.read_b);
			TRIM(args.b_q);
		}
		else if (iType == FA)                            //for *.fa
		{
			getline ( fin_a, args.read_a );

			getline ( fin_b, args.read_head_b );

			if ( args.read_head_b[0] == '>' )
			{
				args.b_id = args.read_head_b.substr ( 1 );
			}
			else
			{
				cerr << "Weird record header b: "<< args.read_head_b << "\n";
				exit (1);
			}

			getline ( fin_b, args.read_b );
			TRIM(args.read_a);
			TRIM(args.read_b);
#undef TRIM
		}
		else
		{
			cerr<< "Will never reach here: "<< __LINE__ << "\n";
			exit (255);
		}

#ifdef DEBUG
		if(iType == FQ && args.a_id.substr(0, args.a_id.size() - 2).compare(args.b_id.substr(0, args.b_id.size() - 2)) != 0)
		{
		  cerr<<"ID of two reads mismatched:\n"<<args.a_id<<"\n"<<args.b_id<<"\n";
		  exit(1);
		}
#endif
		++args.total_pairs;

		//Filter by N.
		if ( Nfilter > -1 )
			if ( NFilter(args) )
				continue;

		//Filter by the quality value of read
		int ab_value(0), bb_value(0);

		if ( !args.a_q.empty() && B_CUTOFF > 0.0 )
		{
			for ( uint i(args.a_q.size() - 1); i >= (args.a_q.size() - Max); --i )
			{
				if ( args.a_q[i] == char(qual_base + qual_offset) ) { ++ab_value; }
			}
			for ( uint i(args.b_q.size() - 1); i >= (args.b_q.size() - Max); --i )
			{
				if ( args.b_q[i] == char(qual_base + qual_offset) ) { ++bb_value; }
			}

			if ( (B_CUTOFF < double ( ab_value ) / double ( Max )) || (B_CUTOFF < double ( bb_value ) / double ( Max )) )
			{
#ifdef DEBUG_VERBOSE
				cerr<<"B_CUTOFF test fail on: " << args.read_a << "\t" << args.read_b << "\n" << args.a_q << "\t" << args.b_q << "\n";
#endif
				++args.low_quality_pairs;

				fout_fq1
					<< "@" << args.a_id << "\n"
					<< args.read_a << "\n"
					<< args.a_s << "\n"
					//<<"+\n"
					<< args.a_q << "\n";

				fout_fq2
					<< "@" << args.b_id << "\n"
					<< args.read_b << "\n"
					<< args.b_s << "\n"
					//<<"+\n"
					<< args.b_q << "\n";

				continue;
			}
		}

		//Connect
		QD;
		revCompl ( args.read_b, args.rread_b );
		QD;
		connected = alignConnect ( args, args.read_a, args.rread_b );
		QD;

		//Analyze result
		int result = analysisResult ( args, connected );
		QD;

		if ( result == 0 || result == 1 ) {continue;}
		if ( result == 2 ) {break;}
	}

#define chg2per(x) (( ( double ) x ) / ( ( double ) args.total_pairs )  * 100.0)
	cerr << "Routine finished.\n"
		 << "       Total Pairs: " << args.total_pairs << "\n"
		 << "Undetermined Pairs: " << args.confuse_pairs << "  (" << chg2per(args.confuse_pairs) << "%)\n"
		 << "      Align Failed: " << args.low_match_rate_pairs << "  (" << chg2per(args.low_match_rate_pairs) << "%)\n"
		 << "         Connected: " << args.connect_pairs << "  (" << chg2per(args.connect_pairs) << "%)\n"
		 << "     Not Connected: " << args.unconnect_pairs << "  (" << chg2per(args.unconnect_pairs) << "%)\n"
		 << "       Low Quality: " << args.low_quality_pairs << "  (" << chg2per(args.low_quality_pairs) << "%)\n"
		 << flush;
#undef chg2per

}

void Routine2 (ConnectArgs_s & args)
{
	cerr << "Work routine 2 engaged...\n";

	int connected;
	while ( getline ( fin_c, args.read_head_a ) )
	{
		connected = 0;

		//Read id
		if ( args.read_head_a[0] == '>' )
		{
			args.a_id = args.read_head_a.substr ( 1 );
		}
		else
		{
			if ( args.read_head_a[0] == '@' )
			{
				cerr << "Parameter '-e' is only for interleaved *.fa file.\n";
				exit ( 1 );
			}
		}

		//Read strand
		getline ( fin_c, args.read_a );
		getline ( fin_c, args.read_head_b );

		if ( args.read_head_b[0] == '>' )
		{
			args.b_id = args.read_head_b.substr ( 1 );
		}
		else
		{
			cerr << "Weird record header b: "<< args.read_head_b << "\n";
			exit (1);
		}

		getline ( fin_c, args.read_b );

		++args.total_pairs;

		//Filter by N.
		if ( Nfilter > -1 )
			if ( NFilter(args) )
				continue;

		//Connect
		revCompl ( args.read_b, args.rread_b );
		connected = alignConnect ( args, args.read_a, args.rread_b );
		//Analyze result
		int result = analysisResult ( args, connected );

		if ( result == 0 || result == 1 ) {continue;}
		if ( result == 2 ) {break;}
	}

#define chg2per(x) (( ( double ) x ) / ( ( double ) args.total_pairs )  * 100.0)
	cerr << "Routine finished.\n"
	<< "       Total Pairs: " << args.total_pairs << "\n"
	<< "Undetermined Pairs: " << args.confuse_pairs << "  (" << chg2per(args.confuse_pairs) << "%)\n"
	<< "      Align Failed: " << args.low_match_rate_pairs << "  (" << chg2per(args.low_match_rate_pairs) << "%)\n"
	<< "         Connected: " << args.connect_pairs << "  (" << chg2per(args.connect_pairs) << "%)\n"
	<< "     Not Connected: " << args.unconnect_pairs << "  (" << chg2per(args.unconnect_pairs) << "%)\n"
	<< "       Low Quality: " << args.low_quality_pairs << "  (" << chg2per(args.low_quality_pairs) << "%)\n"
	<< flush;
#undef chg2per
}

void OutputConnectStat(ConnectArgs_s & args, fstream & O)
{
#define chg2per(x) (( ( double ) x ) / ( ( double ) args.total_pairs )  * 100.0)
	O << "Routine finished.\n"
	<< "       Total Pairs: " << args.total_pairs << "\n"
	<< "Undetermined Pairs: " << args.confuse_pairs << "  (" << chg2per(args.confuse_pairs) << "%)\n"
	<< "      Align Failed: " << args.low_match_rate_pairs << "  (" << chg2per(args.low_match_rate_pairs) << "%)\n"
	<< "         Connected: " << args.connect_pairs << "  (" << chg2per(args.connect_pairs) << "%)\n"
	<< "     Not Connected: " << args.unconnect_pairs << "  (" << chg2per(args.unconnect_pairs) << "%)\n"
	<< flush;
#undef chg2per

	//Output insert size distribution
	O << "insertSize" << "\t" << "#\n";
	for(map<int, uint64_t>::iterator it = args.InsertSize_distr.begin(); it != args.InsertSize_distr.end(); ++it)
	{
		O << it->first << "\t" << it->second << "\n";
	}
	O << "\n";

	//Output real data mismatch distribution
	O << "mismatches" << "\t" << "#" << endl;
	for ( int i(0); i <= args.Max_mismatch_num; ++i )
	{
		O << i << "\t" << args.mismatch_num_distr[i] << "\n";
	}
	O << "\n";
}

int Connect ( int argc, char ** argv )
{
	char* status;
	status = getcwd(cwd, sizeof(cwd)-1);
	strcat(cwd, "/");

	if ( argc <= 2 )
		usage(argv);

	ConnectArgs_s args;
	getOptions ( argc-1, argv+1 );
	int array_size = Max - Min + 1;
	ConnectInit(args);

	QD;
	//Start routine
	int check = checkFile();

	if ( check == 1 ) {Routine1(args);}
	if ( check == 2 ) {Routine2(args);}

	QD;
	cerr << endl;
	//Output Parameters
	cerr << "Command line:" << "\n";
	for(int i(0); i < argc; ++i)
		cerr << argv[i] << " ";
	cerr << "\n";

   
	//Output mismatch array
	int Max_allow_mismatch[array_size];

	for ( int i(0); i < array_size; ++i )
		Max_allow_mismatch[i] = int ( floor((double ( Min + i ) * double ( 1 - MATCH_CUTOFF )) + 0.5));
	cerr << "\n";

	cerr << "overlapRegionLength\tmaximumMismatches\n";

	for ( int i(0); i < array_size; ++i )
	{
		cerr << Min + i << "\t" << Max_allow_mismatch[i] << "\n";
	}
	cerr << "\n";
	
	//Output insert size distribution
	cerr << "insertSize" << "\t" << "#\n";
	for(map<int, uint64_t>::iterator it = args.InsertSize_distr.begin(); it != args.InsertSize_distr.end(); ++it)
	{
		cerr << it->first << "\t" << it->second << "\n";
	}
	cerr << "\n";

	//Output real data mismatch distribution
	cerr << "mismatches" << "\t" << "#" << endl;
	for ( int i(0); i <= args.Max_mismatch_num; ++i )
	{
		cerr << i << "\t" << args.mismatch_num_distr[i] << "\n";
	}
	cerr << "\n";
	
	ConnectClean(args);

	return 0;
}
