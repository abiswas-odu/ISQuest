#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<cstring>
#include<algorithm>
#include<iomanip>
#include<omp.h>
#include <stdlib.h>
#include <unistd.h>
 
#define MAX_TRANSHIT_LEN 100
#define MAX_SEED_LENGTH 3500
#define MAX_DIFF_HAMMING_DIST 10
#define MIN_CLASS_SIZE 2
#define MIN_IR_SCORE 13
#define MAX_IR_LEN 50
#define MIN_IR_LEN 15
#define IR_SEARCH_WINDOW 250
#define MIN_LEN_IR_SEARCH 500
#define MIN_GLOBAL_MATCHES_THRESH 0.80
#define MIN_GLOBAL_MATCHES_THRESH_STRICT 0.95
#define MIN_CLASS_LENGTH 50
#define MIN_CLASS_UNIQ_STRENGTH 1
#define MIN_READ_SUPPORT 10
#define CONSENSUS_THRESH 0.75
#define MIN_TRANSPOSASE_MATCH 50

#define IR_GAP 2.0
#define IR_MATCH 1.0
#define IR_MISMATCH 2.0

#define LONGEST_SEQ_GAP 10.0
#define LONGEST_SEQ_MATCH 1.0
#define LONGEST_SEQ_MISMATCH 2.0
#define LONGEST_SEQ_START 3

#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

using std::istringstream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;
using std::cin;
using std::cerr;
using std::map;
using std::multimap;
using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
const std::string currentDateTime();
void readFasta(ifstream& multiFile, map<string,string> &theMap);
string getRevComp(string s);

class IRepeat
{
	private:
	double score;
	int irlstart,irlend,irrstart,irrend;
	string seqName,seq;
	public:
	void setScore(double s);
	void setIrlstart(int s);
	void setIrlend(int s);
	void setIrrstart(int s);
	void setIrrend(int s);
	void setSeqName(string name);
	void setSeq(string s);
	double getScore() const;
	int getIrlstart() const;
	int getIrlend() const;
	int getIrrstart() const;
	int getIrrend() const;
	string getSeq() const;
	string getSeqName() const;
};
void generateFASTA(string fastaFileName, vector<string> &seqs,string baseName);
void generateSAM(string samFileName, vector<string> &beginSeqs, vector<string> &endSeqs,string refSeqID,string refSeq);


class ISeq
{
	private:
	string seq,baseSeq,seqName;
	IRepeat ir;
	vector<string> uniqBeginMatches,uniqEndMatches;
	int seed_start,seed_end;
	bool foundIR;
	int generatedIterCtr,beginCopyCtr, endCopyCtr;
	bool isInvertedRepeatFoundSW(string seq, string seqName, IRepeat &ir);
	
	public:
	int getSeedStart();
	int getSeedEnd();
	IRepeat getIRepeat();
	bool isInvertedRepeatFound();
	bool operator < (const ISeq& str) const;
	void addISeq(string hitSeqName,string hitSeq, string baseHitSeq, int iterCtr,vector<string> &uniqBeginMatches,vector<string> &uniqEndMatches,int thisSeedStart,int thisSeedEnd,int beginCopies,int endCopies);
	string getISSeq() const;
	string getSeqName() const;
	bool isNOtPrevFound(vector<ISeq> &allISSeqs,int posCtr);
	void writeFinalSeq(string outputPath,string baseName,vector<ISeq> &allISSeqs,int posCtr);
	void generateGKBFile(string fname, string seqName,string seq);
};

