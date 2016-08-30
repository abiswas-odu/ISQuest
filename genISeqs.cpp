#include "genISeqs.h"
#include "nwaln.h"
#include "swaln.h"

extern "C" {
	#include "gbfp.h"
}

void create_directory(string dirName)
{
	string cmdStr = "mkdir " + dirName;
    int ret = system(cmdStr.c_str());
}

void IRepeat::setScore(double s){
	score=s;}
void IRepeat::setIrlstart(int s){
	irlstart=s;}
void IRepeat::setIrlend(int s){
	irlend=s;}
void IRepeat::setIrrstart(int s){
	irrstart=s;}
void IRepeat::setIrrend(int s){
	irrend=s;}
void IRepeat::setSeqName(string name){
	seqName=name;}
void IRepeat::setSeq(string s){
	seq=s;}
double IRepeat::getScore() const {
	return score;}
int IRepeat::getIrlstart() const {
	return irlstart;}
int IRepeat::getIrlend() const {
	return irlend;}
int IRepeat::getIrrstart() const{
	return irrstart;}
int IRepeat::getIrrend() const {
	return irrend;}
string IRepeat::getSeq() const {
	return seq;}
string IRepeat::getSeqName() const {
	return seqName;}

string ISeq::getSeqName() const
{
	return seqName;
}
bool ISeq::isInvertedRepeatFound()
{
	return foundIR;
}
int ISeq::getSeedStart()
{
	return seed_start;
}
int ISeq::getSeedEnd()
{
	return seed_end;
}
IRepeat ISeq::getIRepeat()
{
	return ir;
}
bool ISeq::isInvertedRepeatFoundSW(string seqName, string seq, IRepeat &ir)
{
	if(seq.size()>MIN_LEN_IR_SEARCH)
	{
		string upStream = seq.substr(0,IR_SEARCH_WINDOW);
		string downStream = seq.substr(seq.size()-IR_SEARCH_WINDOW);
		int max_qstart = 0, max_qend = 0,max_sstart = 0, max_send = 0, irLen=0;
		downStream=getRevComp(downStream);
		double max_score = swalign(upStream,downStream,IR_GAP,IR_MATCH,IR_MISMATCH,max_qstart,max_qend,max_sstart,max_send,irLen);
		if(irLen >= MIN_IR_LEN && irLen <= MAX_IR_LEN && max_score>=MIN_IR_SCORE)
		{
			cout<<seqName<<"\n"<<seq<<"\n";
			cout<<"Max Score:"<<max_score<<",Downstream:("<<max_qstart<<","<<max_qend<<"),Upstream("<<seq.size()-max_sstart<<","<<seq.size()-max_send<<")\n";
			cout.flush();
			ir.setScore(max_score);
			ir.setSeqName(seqName);
			ir.setSeq(seq);
			ir.setIrlstart(max_qstart);
			ir.setIrlend(max_qend);
			ir.setIrrstart(seq.size()-max_send);
			ir.setIrrend(seq.size()-max_sstart);
			return true;
		}
	}
	return false;
}
void ISeq::addISeq(string hitSeqName,string hitSeq, string baseHitSeq,int iterCtr,vector<string> &uniqBegin,vector<string> &uniqEnd,int thisSeedStart,int thisSeedEnd,int beginCopies,int endCopies)
{
	seqName=hitSeqName;
	seq=hitSeq;
	baseSeq=baseHitSeq;
	foundIR = isInvertedRepeatFoundSW(hitSeqName,hitSeq, ir);
	generatedIterCtr=iterCtr;
	uniqBeginMatches=uniqBegin;
	uniqEndMatches=uniqEnd;
	seed_start=thisSeedStart;
	seed_end=thisSeedEnd;
	beginCopyCtr=beginCopies;
	endCopyCtr=endCopies;
}
string ISeq::getISSeq() const
{
	if(foundIR)
		return seq.substr(ir.getIrlstart(),(ir.getIrrend()-ir.getIrlstart()+1));
	else
		return seq;
}	
bool ISeq::operator < (const ISeq& str) const
{
    return (this->getISSeq().size() < str.getISSeq().size());
}
bool ISeq::isNOtPrevFound(vector<ISeq> &allISSeqs, int posCtr)
{
    for(int i=posCtr;i<allISSeqs.size();i++) {
		if(this->getSeqName() != allISSeqs[i].getSeqName()) {
			cout<<"GAlign:"<<this->getSeqName()<<" with "<<allISSeqs[i].getSeqName()<<"\n";
			int matches=0,mismatches=0,gaps=0,hitLen=0;
			double score = nwalign(allISSeqs[i].getISSeq(),this->getISSeq(),matches,mismatches,gaps,hitLen);
			
			int rev_matches=0,rev_mismatches=0,rev_gaps=0,rev_hitLen=0;
			double rev_score = nwalign(allISSeqs[i].getISSeq(),getRevComp(this->getISSeq()),rev_matches,rev_mismatches,rev_gaps,rev_hitLen);
			
			double minMatchesCtr1 = 0,minMatchesCtr2=0;
			minMatchesCtr1 = this->getISSeq().size()*MIN_GLOBAL_MATCHES_THRESH;
			minMatchesCtr2 = allISSeqs[i].getISSeq().size()*MIN_GLOBAL_MATCHES_THRESH;
			//Collapse repeating reads
			//if((matches >= minMatchesCtr1 && matches >= minMatchesCtr2) || (rev_matches >= minMatchesCtr1 && rev_matches >= minMatchesCtr2))
			if(matches >= minMatchesCtr1 || matches >= minMatchesCtr2 || rev_matches >= minMatchesCtr1 || rev_matches >= minMatchesCtr2)
			{
				cout<<"Global alignment found\n";
				cout.flush();
				return false;
			}
		}
   }
   return true;
}
void ISeq::writeFinalSeq(string outputPath,string baseName,vector<ISeq> &allISSeqs,int posCtr)
{
	bool isIsDupSeq = true;
		
	if(foundIR && this->isNOtPrevFound(allISSeqs, posCtr))
		isIsDupSeq = false;
	else if(this->isNOtPrevFound(allISSeqs, posCtr))
		isIsDupSeq = false;
	
	string fastaFileName = outputPath + "FinalResult/" + baseName +"_allFinalIS.fasta";
	ofstream fastaFinalFile;
	fastaFinalFile.open(fastaFileName.c_str(), std::ofstream::out | std::ofstream::app);
	
	string gbkFileName = outputPath + "FinalResult/" + baseName +"_allFinalIS.gb";
	ofstream gbkFinalFile;
	gbkFinalFile.open(gbkFileName.c_str(), std::ofstream::out | std::ofstream::app);
	
	string metaFileName = outputPath + "FinalResult/" + baseName +"_finalStat.csv";
	ofstream metaFinalFile;
	metaFinalFile.open(metaFileName.c_str(), std::ofstream::out | std::ofstream::app);
	
	if(!isIsDupSeq)
	{
		//Write fasta sequence for sam file
		cout<<"Writing Result SAM File\n";
		fastaFileName = outputPath + "FinalResult/" + baseName+"_"+SSTR(generatedIterCtr)+"_"+seqName+".fasta";
		ofstream fastaFile;
		fastaFile.open(fastaFileName.c_str());
		fastaFile<<">"<<seqName<<"\n";
		fastaFile<<baseSeq<<"\n";
		fastaFile.close();
		//Write stat file
		metaFinalFile<<seqName<<","<<beginCopyCtr<<","<<endCopyCtr<<"\n";
		//write SAM file
		generateSAM(outputPath+"FinalResult/" + baseName +"_"+SSTR(generatedIterCtr)+"_"+ seqName +".sam", uniqBeginMatches, uniqEndMatches,seqName,baseSeq);
		//Write genbank file if inverted repeat found...

		generateGKBFile(outputPath+"FinalResult/" + baseName+"_"+SSTR(generatedIterCtr)+"_"+ seqName +".gbk",seqName, seq);
		
		string gbk_str;
		string thisGBKFileName = outputPath+"FinalResult/" + baseName+"_"+SSTR(generatedIterCtr)+"_"+ seqName +".gbk";
		ifstream readGBKFile(thisGBKFileName.c_str());
		while(getline(readGBKFile,gbk_str))
		{
			gbkFinalFile<<gbk_str<<"\n";
		}
		gbkFinalFile<<"\n\n";

		fastaFinalFile<<">"<<this->getSeqName()<<"\n";
		fastaFinalFile<<this->getISSeq()<<"\n";
	}
	fastaFinalFile.close();
	gbkFinalFile.close();
	metaFinalFile.close();
}
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%d-%m-%Y", &tstruct);
    return buf;
}

void readFasta(ifstream& multiFile, map<string,string> &theMap)
{
	string str,id,seq;
	getline(multiFile,str);
	id=str.substr(1);
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			theMap.insert (pair<string,string>(id,seq) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	theMap.insert (pair<string,string>(id,seq) );
}
string getRevComp(string s)
{
	int ctgLen=s.length();
	char * cstr = new char [ctgLen+1];
	std::strcpy(cstr, s.c_str());
	char* retChars = new char[ctgLen+1];
	#pragma omp parallel for num_threads(2)
	for(int i=0; i<ctgLen; i++)
	{
		if(cstr[ctgLen-i-1]=='A')
			retChars[i] = 'T';
		else if(cstr[ctgLen-i-1]=='C')
			retChars[i] = 'G';
		else if(cstr[ctgLen-i-1]=='T')
			retChars[i] = 'A';
		else if(cstr[ctgLen-i-1]=='G')
			retChars[i] = 'C';
		else
			retChars[i] = cstr[ctgLen-i-1];
	}
	retChars[ctgLen]='\0';
	string retStr(retChars);
	delete[] retChars;
	delete[] cstr;
	return retStr;
}
bool isCutOverlap(int s,int e,int s1, int e1)
{
	if(s>=s1 && e<=e1)
		return true;
	else if(s1>=s && e1<=e)
		return true;
	return false;
}
int rangeOverlapLen(int s,int e,int s1, int e1,int &ds,int &de)
{
	int retLen=0;
	if(s >= s1 && s <= e1)
	{
		de=e<e1?e:e1;
		ds=s;
		retLen = de-ds;;
	}
	if(s1 >= s && s1 <= e)
	{
		de=e1<e?e1:e;
		ds=s1;
		retLen = de-ds;
	}
	return retLen;
}
map<string,string> getTransHitSeqs(vector<string>& isRelatedCtgBLASTResult, ifstream& CtgFile, string thisPath)
{
    multimap<string,string> ctg_cuts;
	map<string,string> ctg_list;
	map<string,string> transHit_list;
    map<string,string> retSeqs;
	map<string,string> transHitContigCoordMap;
    readFasta(CtgFile,ctg_list);
	for(int i=0;i<isRelatedCtgBLASTResult.size();i++)
	{
		std::vector<std::string> tok = split(isRelatedCtgBLASTResult[i], ',');
		string ctgID =  tok[0];
		int endIndx = tok.size();
		int qstart = atoi(tok[endIndx-2].c_str());
		int qend = atoi(tok[endIndx-1].c_str());

		std::multimap<string,string>::iterator iter=ctg_cuts.find(ctgID);

		if(iter != ctg_cuts.end())
		{
			pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
			ppp = ctg_cuts.equal_range(ctgID);
			int overlapFound=0;
			for (multimap<string,string>::iterator it2 = ppp.first; it2 != ppp.second; ++it2)
			{
				string cutRange = (*it2).second;
				std::vector<std::string> tok2 = split(cutRange, ',');
				int prev_qstart = atoi(tok2[0].c_str());
				int prev_qend = atoi(tok2[1].c_str());
				//Check if there is overlap with prev contig cut ranges
				if(isCutOverlap(qstart,qend,prev_qstart,prev_qend))
				{
					overlapFound=1;
					cout<<"Overlap Found..."<<ctgID<<","<<cutRange<<"\n";
					cout.flush();
					int cur_qsrart, cur_qend;
					if(prev_qstart < qstart)
						cur_qsrart = prev_qstart;
					else
						cur_qsrart = qstart;

					if(prev_qend < qend)
						cur_qend = qend;
					else
						cur_qend = prev_qend;

					cutRange =  SSTR(cur_qsrart) + "," + SSTR(cur_qend);
					it2->second = cutRange;
					break;
				}
			}
			//If the range is a another new one for this contig...
			if(overlapFound == 0)
			{
				string cutRange =  SSTR(qstart) + "," + SSTR(qend);
				ctg_cuts.insert(pair<string,string>(ctgID,cutRange));
			}
		}
		else
		{
			string cutRange =  SSTR(qstart) + "," + SSTR(qend);
			ctg_cuts.insert(pair<string,string>(ctgID,cutRange));
		}
	}
	cout<<"Overlap Calc Complete\n";
	cout.flush();
	int hitSeqCtr=0;
	for(std::map<string,string>::iterator it=ctg_cuts.begin(); it!=ctg_cuts.end(); ++it)
	{
		string key = (*it).first;
		string cutRange = (*it).second;
		std::vector<std::string> tok = split(cutRange, ',');
		//REDUCE TO 0 INDEX 
		int qstart = atoi(tok[0].c_str())-1;
		int qend = atoi(tok[1].c_str())-1;

		std::map<string,string>::iterator iter=ctg_list.find(key);
		if(iter != ctg_list.end())
		{
			//Check error point
			cout<<"Check:"<<key<<"\n";
			cout.flush();
			string seq = (*iter).second;
			string transHit = seq.substr(qstart,qend-qstart);

			if(transHit.size() > MAX_TRANSHIT_LEN)
			{
				//avoid duplicate tranposase hits even from different contigs...
				std::map<string,string>::iterator iter2=transHit_list.find(transHit);
				if(iter2 == transHit_list.end())
				{
					retSeqs.insert(pair<string,string>("transposaseHit_"+SSTR(hitSeqCtr),transHit));
					transHit_list.insert(pair<string,string>(transHit,key));
					transHitContigCoordMap.insert(pair<string,string>("transposaseHit_"+SSTR(hitSeqCtr),key+","+ SSTR(qstart) +","+SSTR(qend)));
					hitSeqCtr++;
				}
			}
		}
		else
		{
			cout<<"Not Found:"<<key<<"\n";
			cout.flush();
		}
    }
	//Write out the mapping of tranposase hit to contig ID
	ofstream coordFile;
	string coordFileName = thisPath+"IS_Contig_Map.txt";
	coordFile.open(coordFileName.c_str());
	for(std::map<string,string>::iterator it=transHitContigCoordMap.begin(); it!=transHitContigCoordMap.end(); ++it)
	{
		coordFile<<it->first<<","<<it->second<<"\n";
	}
	coordFile.close();
    return retSeqs;
}

void getRepeatSeqs(string hitSeqName,string hitSeq, map<string,string> &rawReads, vector<string> &blastHits,vector<string> &beginMatches,vector<string> &endMatches)
{
	for(int i=0;i<blastHits.size();i++)
	{
		string blast_str = blastHits[i];
		std::vector<std::string> tok = split(blast_str, ',');
		string readID =  tok[0];
		int qstart = atoi(tok[2].c_str());
		int qend = atoi(tok[3].c_str());

		int hstart = atoi(tok[4].c_str());
		int hend = atoi(tok[5].c_str());

		std::map<string,string>::iterator iter=rawReads.find(readID);
		string readSeq = (*iter).second;

		//Typical begin hit...
		if(hstart == 1){
			if(qstart>1 && qend==readSeq.size()){
				//isBorderHit = isFwd = isBeginBorder = 1;
				beginMatches.push_back(readSeq.substr(0,qstart));
			}
		}
		//Reverse hit at the rear
		else if(hstart==hitSeq.size()){
			if(qstart>1 && qend==readSeq.size()){
				//isBorderHit = 1;
				endMatches.push_back(getRevComp(readSeq.substr(0,qstart)));
			}
		}
		//
		else if(hend==1){
			if(qstart==1 && qend<readSeq.size()){
				//isBorderHit = isBeginBorder = 1;
				beginMatches.push_back(getRevComp(readSeq.substr(qend)));
			}
		}
		//Typical end hit...
		else if(hend==hitSeq.size()){
			if(qstart==1 && qend<readSeq.size()){
				//isBorderHit = isFwd = 1;
				endMatches.push_back(readSeq.substr(qend));
			}
		}
    }
}

void generateFASTA(string fastaFileName, vector<string> &seqs,string baseName)
{
	ofstream fastaFile;
	fastaFile.open(fastaFileName.c_str());
	for(int i=0;i<seqs.size();i++)
	{
		fastaFile<<">"<<baseName<<"_"<<i<<"\n";
		fastaFile<<seqs[i]<<"\n";
	}
	fastaFile.close();
}

void generateSAM(string samFileName, vector<string> &beginSeqs, vector<string> &endSeqs,string refSeqID,string refSeq)
{
	ofstream samFile;
	samFile.open(samFileName.c_str());
	samFile<<"@RG\tID:"<<refSeqID<<"\n";
	for(int i=0;i<beginSeqs.size();i++)
	{
		string revStr = beginSeqs[i];
		std::reverse(revStr.begin(),revStr.end());
		string fullSeq = revStr+refSeq;
		string cigarStr = SSTR(revStr.size()) + std::string("I") + SSTR(refSeq.size()) + std::string("M");
		string retStr = "Begin_" + SSTR(i) + std::string("\t") + "0" + std::string("\t") + refSeqID + std::string("\t") + "1" + std::string("\t255\t") + cigarStr + std::string("\t*\t0\t0\t") + fullSeq + std::string("\t*\tRG:Z:Unpaired reads assembled against ") + refSeqID + std::string(" \n");
		samFile<<retStr;
	}
	for(int i=0;i<endSeqs.size();i++)
	{
		string fullSeq = refSeq + endSeqs[i];
		string cigarStr = SSTR(refSeq.size()) + std::string("M") + SSTR(endSeqs[i].size()) + std::string("I");
		string retStr = "End_" + SSTR(i) + std::string("\t") + "0" + std::string("\t") + refSeqID + std::string("\t") + "1" + std::string("\t255\t") + cigarStr + std::string("\t*\t0\t0\t") + fullSeq + std::string("\t*\tRG:Z:Unpaired reads assembled against ") + refSeqID + std::string(" \n");
		samFile<<retStr;
	}
	samFile.close();
}

void ISeq::generateGKBFile(string fname, string seqName,string seq)
{
	ofstream gkbFile;
	std::transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
	gkbFile.open(fname.c_str());
	gkbFile<<"LOCUS       "<<seqName<<" "<<seq.size()<<" bp "<<"DNA linear UNA "<<currentDateTime()<<"\n";
	gkbFile<<"DEFINITION  .\n";
	gkbFile<<"ACCESSION   "<<seqName<<"\n";
	gkbFile<<"VERSION     "<<seqName<<"\n";
	gkbFile<<"KEYWORDS    .\n";
	gkbFile<<"SOURCE      \n";
	gkbFile<<"  ORGANISM  \n";
	gkbFile<<"            .\n";
	gkbFile<<"FEATURES             Location/Qualifiers\n";
	if(isInvertedRepeatFound())
	{
		gkbFile<<"     repeat_unit     "<<ir.getIrlstart()+1<<".."<<ir.getIrlend()+1<<"\n";
		gkbFile<<"                     /label=\"IRL\"\n";
		gkbFile<<"     repeat_unit     "<<ir.getIrrstart()+1<<".."<<ir.getIrrend()+1<<"\n";
		gkbFile<<"                     /label=\"IRR\"\n";
	}
	gkbFile<<"     seed_unit     "<<getSeedStart()+1<<".."<<getSeedEnd()+1<<"\n";
	gkbFile<<"                     /label=\"BLAST Transposase Seed\"\n";
	gkbFile<<"ORIGIN      \n";
	int isLineCtr=1, usedStr=0;
	for(int i=0;i<seq.size();i+=60)
	{
		gkbFile<<setw(9);
		gkbFile<<isLineCtr;

		isLineCtr+=60;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);

		usedStr+=10;
		gkbFile<<"\n";
	}
	gkbFile<<"//\n\n";
}

int hammingDist(string str1,string str2)
{
	string shortString,longString;
	if(str1.length()>str2.length())
	{
		shortString=str2;
		longString=str1;
	}
	else
	{
		shortString=str1;
		longString=str2;
	}
	int diffCtr=0;
	for(int i=0;i<shortString.length();i++)
	{
		if(shortString[i]!=longString[i])
			diffCtr++;
	}
	return diffCtr;
}
bool greatLength(const string& s1, const string& s2)
{
    return s1.length() > s2.length();
}
bool lessLength (const string& s1, const string& s2)
{
    return s1.length() < s2.length();
}

vector<string> myriadClasses(vector<string> seqs, map<string,int> readSupport)
{
	vector<string> classes,retClasses;
	vector<int> classCtr;
	sort(seqs.begin(),seqs.end(), greatLength);
	if(seqs.size()>0)
	{
		classes.push_back(seqs[0]);
		map<string,int>::iterator iter = readSupport.find(seqs[0]);
		classCtr.push_back(iter->second);
		cout<<"Unique sequence read 0"<<" support:"<<iter->second<<"\n";
		for(int i=1;i<seqs.size();i++)
		{
			iter = readSupport.find(seqs[i]);	
			cout<<"Unique sequence read "<<i<<" support:"<<iter->second<<"\n";
					
			if(seqs[i].size()>=MIN_CLASS_LENGTH)
			{
				int belongsToClass=-1;
				for(int j=0;j<classes.size();j++)
				{
					int matches=0,mismatches=0,gaps=0,hitLen=0;
					double score = nwalign(seqs[i],classes[j],matches,mismatches,gaps,hitLen);
					double minMatchesCtr1 = 0,minMatchesCtr2=0;
					minMatchesCtr1 = seqs[i].size()*MIN_GLOBAL_MATCHES_THRESH;
					if(matches >= minMatchesCtr1)
					{
						belongsToClass=j;
						break;
					}
				}
				if(belongsToClass==-1)
				{
					classes.push_back(seqs[i]);
					classCtr.push_back(iter->second);
				}
				else
					classCtr[belongsToClass]=classCtr[belongsToClass]+iter->second;
			}
		}
	
		for(int j=0;j<classCtr.size();j++)
		{
			cout<<"Class "<<j<<" support:"<<classCtr[j]<<"\n"<<classes[j]<<'\n';
			if(classCtr[j]>=MIN_READ_SUPPORT)
				retClasses.push_back(classes[j]);
		}		
	}
	return retClasses;
}
string longestCommonSeq(vector<string> seqs)
{
	string retStr="";
	if(seqs.size()>0)
	{
		std::sort(seqs.begin(),seqs.end(),lessLength);
		string base_seq = seqs[0];
		int max_base_start=0, min_base_end=seqs[0].size();
		for(int i=1;i<seqs.size();i++)
		{
			int base_start = 0, base_end = 0,curr_start = 0, curr_end = 0, matchLen=0;
			double max_score = swalign(base_seq,seqs[i],LONGEST_SEQ_GAP,LONGEST_SEQ_MATCH,LONGEST_SEQ_MISMATCH,base_start,base_end,curr_start,curr_end,matchLen);
			if(base_start>max_base_start)
				max_base_start=base_start;
			if(base_end<min_base_end)
				min_base_end=base_end;
		}
		if(max_base_start <= LONGEST_SEQ_START)
			retStr = base_seq.substr(0,min_base_end);
	}
	return retStr;
}
bool isInvertedRepeatFoundBLAST(string seq, string seqName, string blast_path,string thisPath, IRepeat &ir)
{
	if(seq.size()>1000)
	{
		string upStream = seq.substr(0,500);
		string downStream = seq.substr(seq.size()-500);

		ofstream irHitFile;
		string irHitFileName = thisPath+"upstream.fasta";
		irHitFile.open(irHitFileName.c_str());
		irHitFile<<">upstream"<<"\n"<<upStream<<"\n";
		irHitFile.close();

		string cmdStr = blast_path + "makeblastdb -in "+thisPath+"upstream.fasta -dbtype nucl -title upstream -out "+thisPath+"upstream";
		int ret = system(cmdStr.c_str());

		irHitFileName = thisPath+"downstream.fasta";
		irHitFile.open(irHitFileName.c_str());
		irHitFile<<">downstream"<<"\n"<<downStream<<"\n";
		irHitFile.close();

		cmdStr = blast_path + "blastn -task blastn -db "+thisPath+"upstream -query "+thisPath+"downstream.fasta -gapopen 4 -gapextend 4 -penalty -3 -reward 1 -out "+thisPath+"irHitResult.txt -outfmt \"10 qseqid sseqid score qstart qend sstart send\" ";
		ret = system(cmdStr.c_str());

		irHitFileName=thisPath+"irHitResult.txt";
		ifstream readResultFile(irHitFileName.c_str());
		if(readResultFile.fail())
		{
			cerr<<"\nThe "+thisPath+"irHitResult.txt file could not be opened.";
			return false;
		}
		string blast_str;
		/*Get best IR Hit*/
		int max_score = MIN_IR_SCORE, max_qstart = 0, max_qend = 0,max_sstart = 0, max_send = 0, isHitFound = 0;
		int max_transLen = 0;
		while(getline(readResultFile,blast_str))
		{
			std::vector<std::string> tok = split(blast_str, ',');
			int score = atoi(tok[2].c_str());
			int qstart = atoi(tok[3].c_str());
			int qend = atoi(tok[4].c_str());
			int sstart = atoi(tok[5].c_str());
			int send = atoi(tok[6].c_str());
			int transLen = seq.size()-500 + send - qend;
			if((qend-qstart) >= MIN_IR_LEN && (qend-qstart) <= MAX_IR_LEN && score >= max_score && transLen > max_transLen && (sstart>send))
			{
				isHitFound=1;
				max_score = score;
				max_qstart = qstart;
				max_qend = qend;
				max_sstart = sstart;
				max_send = send;
				max_transLen = seq.size()-500 + send - qend;
			}
		}
		if(isHitFound)
		{
			cout<<seqName<<"\n"<<seq<<"\n";
			cout<<"Max Score:"<<max_score<<",Downstream:("<<max_qstart<<","<<max_qend<<"),Upstream("<<seq.size()-500+max_sstart<<","<<seq.size()-500+max_send<<")\n";
			ir.setScore(max_score);
			ir.setSeqName(seqName);
			ir.setSeq(seq);
			ir.setIrlstart(max_qstart);
			ir.setIrlend(max_qend);
			ir.setIrrstart(seq.size()-500+max_send);
			ir.setIrrend(seq.size()-500+max_sstart);
			return true;
		}
	}
	return false;
}
void doBLAST(string readFilePath,string blast_path, string mpiPath, string thisPath, int splitSize)
{
	ifstream myfile(readFilePath.c_str());   
    unsigned seq_count =0;
	string str;
	while(getline(myfile,str))
	{
		if (str[0]=='>')
			seq_count++;
	}
	cout << "Sequences: " << seq_count << "\n";
	myfile.close();
	//Split the File into 128 pieces
	unsigned fileCtr=0,seqCtr=0;
	unsigned seqPerFile = seq_count/splitSize;
	myfile.open(readFilePath.c_str(),ifstream::in);
	string seqId,seq;
	ofstream blastQueryFile;
	string blastQueryFileName = thisPath+"blastQuery_"+SSTR(fileCtr)+".fasta";
	blastQueryFile.open(blastQueryFileName.c_str(),ofstream::out);
	getline(myfile,str);
	seqId=str;
	while(getline(myfile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			blastQueryFile<<seqId<<"\n"<<seq<<'\n';
			seqCtr++;
			if(seqCtr>=seqPerFile)
			{
				fileCtr++;
				blastQueryFile.close();
				blastQueryFileName = thisPath+"blastQuery_"+SSTR(fileCtr)+".fasta"; 
				blastQueryFile.open(blastQueryFileName.c_str(),ofstream::out);
				seqCtr=0;
			}
			seq="";
			seqId = str;
		}
		else
		{
			seq=seq+str;
		}
	}
	blastQueryFile<<seqId<<"\n"<<seq<<'\n';
	blastQueryFile.close();
	myfile.close();
			
	string cmdStr = "rm "+thisPath+"transHitsBLAST*";
	int ret = system(cmdStr.c_str());
	
	//Read file split into n parts...blastQuery_0...n.fasta
	cmdStr = "mpirun -n " + SSTR(splitSize) + " " + mpiPath + "mpiblast "+blast_path+ " " + thisPath;
	ret = system(cmdStr.c_str());
	
	cmdStr = "cat "+thisPath+"transHitsBLAST* > "+thisPath+"transHitsResult.txt";
	ret = system(cmdStr.c_str());

	cmdStr = "rm "+thisPath+"blastQuery_*";
	ret = system(cmdStr.c_str());
}
map<string,string> expandTransposaseSeqs(map<string,string> &hitSeqs,map<string,int> &hitStarts, map<string,int> &hitEnds, map<string,string> &rawReads, vector<ISeq> &allISSeqs,string readFilePath, string blast_path, string mpiPath,int mpiThreads, string thisPath, int iterCtr, string baseName,string writeIntermideateFiles)
{
	map<string,string> newHitSeqs;
	map<string,int> newHitStarts,newHitEnds;
	ofstream transHitFile;
	string transFileName = thisPath+"transHits.fasta";
    transHitFile.open(transFileName.c_str());

	for(std::map<string,string>::iterator iter = hitSeqs.begin(); iter != hitSeqs.end(); iter++) {
		transHitFile<<">"<<(*iter).first<<"\n"<<(*iter).second<<"\n";
    }
    transHitFile.close();

    string cmdStr = blast_path + "makeblastdb -in "+thisPath+"transHits.fasta -dbtype nucl -title transHits -out "+thisPath+"transHits";
    int ret = system(cmdStr.c_str());

	if(mpiPath=="" || thisPath=="")
	{
		cmdStr = blast_path + "blastn -task megablast -db "+thisPath+"transHits -evalue 1e-01 -ungapped -num_threads 32  -query " + readFilePath + " -out "+thisPath+"transHitResult.txt -outfmt \"10 qseqid sseqid qstart qend sstart send\" ";
		ret = system(cmdStr.c_str());
	}
	else
		doBLAST(readFilePath,blast_path,mpiPath,thisPath,mpiThreads);

    cout<<"Iteration "<<iterCtr<<" BLAST complete...\n";
	for(std::map<string,string>::iterator iter = hitSeqs.begin(); iter != hitSeqs.end(); iter++) {
		string thisHitSeq = (*iter).second;
		string thisHitSeqName = (*iter).first;
		vector<string> beginMatches, endMatches;

        string fname = thisPath + "transHitsResult.txt";
        ifstream transHitsBlastResult(fname.c_str());

        vector<string> blastHits;
        string blast_str;
        while(getline(transHitsBlastResult,blast_str))
        {
            std::vector<std::string> tok = split(blast_str, ',');
            string hitSeqName =  tok[1];
            if(hitSeqName == thisHitSeqName)
                blastHits.push_back(blast_str);
        }
        transHitsBlastResult.close();


		getRepeatSeqs(thisHitSeqName,thisHitSeq, rawReads, blastHits,beginMatches,endMatches);
		cout<<"\n\n\n\nTransposase Hit Sequence: "<<thisHitSeqName<<"\n";
		cout<<"Begin Matches:"<<beginMatches.size()<<"\n";
		cout<<"End Matches:"<<endMatches.size()<<"\n";

		for(int j=0;j<beginMatches.size();j++)
		{
			string revStr = beginMatches[j];
			std::reverse(revStr.begin(),revStr.end());
			beginMatches[j] = revStr;
		}
		vector<string> uniqBeginMatches;
		map<string,int> uniqBeginReadSupport;
		for(int j=0;j<beginMatches.size();j++)
		{
			int found=0;
			for(int k=0;k<beginMatches.size();k++)
			{
				if(j!=k){
					if (beginMatches[k].find(beginMatches[j]) != std::string::npos) {
						found++;
						break;
					}
				}
			}
			if(found==0)
			{
				uniqBeginMatches.push_back(beginMatches[j]);
				uniqBeginReadSupport.insert(pair<string,int>(beginMatches[j],0));
			}
		}
		std::sort(uniqBeginMatches.begin(),uniqBeginMatches.end(),lessLength);
		//cout<<"\nCount Unique Read Strength:"<<"\n";
		for(int j=0;j<beginMatches.size();j++)
		{
			int copyCtr=0;
			string uniqKey;
			for(map<string,int>::iterator iter = uniqBeginReadSupport.begin(); iter != uniqBeginReadSupport.end(); iter++) {
				if (iter->first.find(beginMatches[j]) != std::string::npos) {
					copyCtr++;
					uniqKey=iter->first;
				}
			}
			if(copyCtr==1)
			{
				map<string,int>::iterator iter = uniqBeginReadSupport.find(uniqKey);
				iter->second = iter->second + 1; 
			}
		}
		
		vector<string> beginClasses = myriadClasses(uniqBeginMatches,uniqBeginReadSupport);
		
		//Get Longest Common Left Match
		cout<<"Longest Common Left Sequence:"<<'\n';
		string commonBeginSeq = longestCommonSeq(beginClasses);
		cout<<commonBeginSeq<<'\n';
		cout<<"Class Count:"<<beginClasses.size()<<"\n";

		//The other side...
		vector<string> uniqEndMatches;
		map<string,int> uniqEndReadSupport;
		for(int j=0;j<endMatches.size();j++)
		{
			int found=0;
			for(int k=0;k<endMatches.size();k++)
			{
				if(j!=k){
					if (endMatches[k].find(endMatches[j]) != std::string::npos) {
						found=1;
						break;
					}
				}
			}
			if(found==0) {
				uniqEndMatches.push_back(endMatches[j]);
				uniqEndReadSupport.insert(pair<string,int>(endMatches[j],0));
			}
		}
		std::sort(uniqEndMatches.begin(),uniqEndMatches.end(),lessLength);
		//cout<<"\nCount Unique End Read Strength:"<<"\n";
		for(int j=0;j<endMatches.size();j++)
		{
			int copyCtr=0;
			string uniqKey;
			for(map<string,int>::iterator iter = uniqEndReadSupport.begin(); iter != uniqEndReadSupport.end(); iter++) {
				if (iter->first.find(endMatches[j]) != std::string::npos) {
					copyCtr++;
					uniqKey=iter->first;
				}
			}
			if(copyCtr==1)
			{
				map<string,int>::iterator iter = uniqEndReadSupport.find(uniqKey);
				iter->second = iter->second + 1; 
			}
		}
		//cout<<"\nRight Match Options:"<<"\n";

		vector<string> endClasses = myriadClasses(uniqEndMatches,uniqEndReadSupport);
		//Get Longest Common Right Match
		cout<<"\n\nLongest Common Right Sequence:"<<'\n';
		string commonEndSeq = longestCommonSeq(endClasses);
		cout<<commonEndSeq<<'\n';
		cout<<"Class Count:"<<endClasses.size()<<"\n";
		
		bool isFinalISForm = false;
		/*Conditions for IS Final Result*/
		if((beginClasses.size()>=MIN_CLASS_SIZE && endClasses.size()>=MIN_CLASS_SIZE) || 
			thisHitSeq.size() >= MAX_SEED_LENGTH || 
			(beginClasses.size()==0 && endClasses.size()==0))
		{
			isFinalISForm = true;
			std::reverse(commonBeginSeq.begin(),commonBeginSeq.end());
			std::map<string,int>::iterator iter=hitStarts.find(thisHitSeqName);
			int thisSeedStart = (*iter).second + commonBeginSeq.size();
			iter=hitEnds.find(thisHitSeqName);			
			int thisSeedEnd = (*iter).second + commonBeginSeq.size();
			
			ISeq newFinalSeq;
			newFinalSeq.addISeq(thisHitSeqName,commonBeginSeq + thisHitSeq + commonEndSeq,thisHitSeq,iterCtr,uniqBeginMatches,uniqEndMatches,thisSeedStart,thisSeedEnd,beginClasses.size(),endClasses.size());
			allISSeqs.push_back(newFinalSeq);
		}
		else if(uniqBeginMatches.size()>0 && uniqEndMatches.size()>0 && writeIntermideateFiles == "YES")
		{
			//Write fasta sequence for sam file
			string fastaFileName = thisPath+"level"+SSTR(iterCtr) + "/" + baseName+"_"+SSTR(iterCtr)+"_"+thisHitSeqName+".fasta";
			ofstream fastaFile;
			fastaFile.open(fastaFileName.c_str());
			fastaFile<<">"<<thisHitSeqName<<"\n";
			fastaFile<<thisHitSeq<<"\n";
			fastaFile.close();
			//write SAM file
			generateSAM(thisPath+"level"+SSTR(iterCtr) + "/" + baseName+"_"+SSTR(iterCtr)+"_"+ thisHitSeqName +".sam", uniqBeginMatches, uniqEndMatches,thisHitSeqName,thisHitSeq);
		}
		int childSeqID=0;
		if(!isFinalISForm)
		{
			for(int j=0;j<beginClasses.size();j++)
			{
				for(int k=0;k<endClasses.size();k++) 
				{
					string revStr = beginClasses[j];
					std::reverse(revStr.begin(),revStr.end());
					if(beginClasses.size()<MIN_CLASS_SIZE && endClasses.size()<MIN_CLASS_SIZE)
					{
						newHitSeqs.insert(pair<string,string>(revStr + thisHitSeq + endClasses[k],thisHitSeqName + "_" + SSTR(childSeqID)));
						std::map<string,int>::iterator iter=hitStarts.find(thisHitSeqName);
						newHitStarts.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second + revStr.size()));
						iter=hitEnds.find(thisHitSeqName);
						newHitEnds.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second + revStr.size()));
						childSeqID++;
					}
					else if(beginClasses.size()<MIN_CLASS_SIZE)
					{
						newHitSeqs.insert(pair<string,string>(revStr + thisHitSeq,thisHitSeqName + "_" + SSTR(childSeqID)));
						std::map<string,int>::iterator iter=hitStarts.find(thisHitSeqName);
						newHitStarts.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second + revStr.size()));
						iter=hitEnds.find(thisHitSeqName);
						newHitEnds.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second + revStr.size()));
						childSeqID++;
					}
					else if(endClasses.size()<3)
					{
						newHitSeqs.insert(pair<string,string>(thisHitSeq + endClasses[k],thisHitSeqName + "_" + SSTR(childSeqID)));
						std::map<string,int>::iterator iter=hitStarts.find(thisHitSeqName);
						newHitStarts.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second));
						iter=hitEnds.find(thisHitSeqName);
						newHitEnds.insert(pair<string,int>(thisHitSeqName + "_" + SSTR(childSeqID),(*iter).second));
						childSeqID++;
					}
				}
			}
		}
	}
	cout<<"Working ISeqs:"<<newHitSeqs.size()<<"\n";
	cout.flush();
	map<string,string> newReturnHitSeqs;
	for(std::map<string,string>::iterator iter = newHitSeqs.begin(); iter != newHitSeqs.end(); iter++) {
		string thisHitSeq = (*iter).first;
		string thisHitSeqName = (*iter).second;
		newReturnHitSeqs.insert(pair<string,string>(thisHitSeqName,thisHitSeq));
	}
	//update seed boundaries
	hitStarts = newHitStarts;
	hitEnds = newHitEnds;
	return newReturnHitSeqs;
}
static char *searchQualValue(const char *sWord, gb_feature *ptFeature) {
    gb_qualifier *i;
	char noteField[] = "note";
    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++)
        if (strcmp(i->sQualifier,noteField) != 0  && strstr(i->sValue,sWord) != NULL)
            return i->sValue;
    return NULL;
}
bool isAnnotationHit(vector<string> sWord, gb_feature *ptFeature)
{
	char *sValue = NULL;
	for(int i=0;i<sWord.size();i++)
	{
		if ((sValue = searchQualValue(sWord[i].c_str(), ptFeature)) != NULL) {
			return true;
		}
	}
	return false;
}
/* Some accession files cause errors. Black list them*/
vector<string> CreateBlackList(string blackListFile)
{
	vector<string> blackListAcc;
	string str;
	ifstream bListFile(blackListFile.c_str());

	if(bListFile.fail())
	{
		cerr<<"\nThe black list files could not be opened.";
	}
	else
	{
		while(getline(bListFile,str))
		{
			blackListAcc.push_back(str);
		}
	}
	return blackListAcc;
}
bool IsBlackList(string accCode,vector<string> &blackListAcc)
{
	if(std::find(blackListAcc.begin(), blackListAcc.end(), accCode)==blackListAcc.end() && accCode.find("BX") == std::string::npos)   //Check black listed acc codes
    {
		return 1;
	}
	return 0;
}
vector<string> findTransposaseHits(ifstream& ctgBlastFile, vector<string> &blackListAcc, string keywords,string genbank_path, string wget_path)
{
    
	std::vector<std::string> searchKeys = split(keywords, ',');

    multimap<string,string> blastStrings;
    vector<string> isRelatedCtgBLASTResult;
    string blast_str;
    while(getline(ctgBlastFile,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
        string accCode =  tok[1];
        blastStrings.insert(pair<string,string>(accCode,blast_str));
	}

	std::multimap<string,string>::iterator m_it, s_it;

	for (m_it = blastStrings.begin();  m_it != blastStrings.end();  m_it = s_it)
    {
        string accCode = (*m_it).first;
        if(IsBlackList(accCode,blackListAcc))   //Check black listed acc codes
        {
            // Load genbank File
            string accFilePath = genbank_path + accCode + ".gbk";
            ifstream gbkFile(accFilePath.c_str());
            gb_data **pptSeqData, *ptSeqData;
            gb_feature *ptFeature;
            if(gbkFile.fail())
            {
                cout<<"\nThe genebank files could not be found:"<<accFilePath<<"\n";
                cout<<"Downloading..."<<"\n";
				cout.flush();
                string urlStr = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=" + accCode;
                string cmdStr = wget_path + "wget \"" + urlStr +"\" -O " + accFilePath;
                int ret = system(cmdStr.c_str());
                if(ret==0)
                {
                    cout<<"\nThe genebank file downloaded:"<<accFilePath<<"\n";
                    char * cstr = new char [accFilePath.length()+1];
                    std::strcpy(cstr, accFilePath.c_str());
                    pptSeqData = parseGBFF(cstr);
                    delete[] cstr;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                cout<<"\nThe genebank file found:"<<accFilePath<<"\n";
				cout.flush();
                char * cstr = new char [accFilePath.length()+1];
                std::strcpy(cstr, accFilePath.c_str());
                pptSeqData = parseGBFF(cstr);
                delete[] cstr;
            }


            // Iterate over all map elements with given accCode
            pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
            ppp = blastStrings.equal_range(accCode);
			int ctr=0;
            for (s_it = ppp.first;  s_it != ppp.second;  ++s_it)
            {
				ctr++;
                blast_str = (*s_it).second;
                std::vector<std::string> tok = split(blast_str, ',');
                int endIndx = tok.size();
				int ori = atoi(tok[endIndx-2].c_str())<atoi(tok[endIndx-1].c_str())?0:1;
				int qs=atoi(tok[endIndx-4].c_str());
				int qe=atoi(tok[endIndx-3].c_str());
	            int sstart = atoi(tok[endIndx-2].c_str())<atoi(tok[endIndx-1].c_str())?atoi(tok[endIndx-2].c_str()):atoi(tok[endIndx-1].c_str());
                int send = atoi(tok[endIndx-2].c_str())<atoi(tok[endIndx-1].c_str())?atoi(tok[endIndx-1].c_str()):atoi(tok[endIndx-2].c_str());
                for (int i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
					for (int j = 0; j < ptSeqData->iFeatureNum; j++) {
						ptFeature = (ptSeqData->ptFeatures + j);
						int ds=0,de=0;
						if(rangeOverlapLen(sstart,send,ptFeature->lStart,ptFeature->lEnd,ds,de) > 50 &&
							isAnnotationHit(searchKeys,ptFeature))
						{
							if(ori==0)
							{
								int tempqs = qs + (ds-sstart);
								int tempqe = qe - (send-de);
								cout<<blast_str<<"\n"<<ptFeature->sFeature<<","<<ptFeature->lStart<<","<<ptFeature->cDirection<<","<<ptFeature->lEnd<<","<<"\n";
								cout<<blast_str + "," + SSTR(tempqs) + "," + SSTR(tempqe)<<"\n";
								isRelatedCtgBLASTResult.push_back(blast_str + "," + SSTR(tempqs) + "," + SSTR(tempqe));
							}
							else
							{
								int tempqs = qs + (send-de);
								int tempqe = qe - (ds-sstart);							
								cout<<blast_str<<"\n"<<ptFeature->sFeature<<","<<ptFeature->lStart<<","<<ptFeature->cDirection<<","<<ptFeature->lEnd<<","<<"\n";
								cout<<blast_str + "," + SSTR(tempqs) + "," + SSTR(tempqe)<<"\n";
								isRelatedCtgBLASTResult.push_back(blast_str + "," + SSTR(tempqs) + "," + SSTR(tempqe));
							}
						}
					}
				}
            }
			cout<<"Ctr:"<<ctr<<"\n";
            freeGBData(pptSeqData);

        }
        else       //just general search for keywords in blast string
        {
            pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
            ppp = blastStrings.equal_range(accCode);

            // Iterate over all map elements with given accCode

            for (s_it = ppp.first;  s_it != ppp.second;  ++s_it)
            {
                blast_str = (*s_it).second;
                for(int i=0;i<searchKeys.size();i++)
                {
                    if (blast_str.find(searchKeys[i]) != std::string::npos)
                        isRelatedCtgBLASTResult.push_back(blast_str);
                }
            }

        }
    }
    return isRelatedCtgBLASTResult;
}

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

/*function that show the help information*/
void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument1] [argument2] [argument3] [argument4]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-v  show version information"<<endl;
  cout<<"         "<<"[argument1]  Parameter file."<<endl;
  cout<<"         "<<"[argument2]  Path to the output directory."<<endl;
  cout<<"         "<<"[argument3]  Contigs in multi FASTA format."<<endl;
  cout<<"         "<<"[argument4]  Raw reads in multi FASTA format.."<<endl;
  cout<<"         "<<"[argument5]  Output file name prefix."<<endl;
  cout<<"example: "<<s<<" param.conf /home/xxx/324-958/ Contigs324-958.fasta 324-958_A_GGACTCCT-GTAAGGAG_L001_R_TRIM_PAIRED.fasta 324-958"<<endl;
}

int main(int argc, char **argv)
{
	//ISQuest
	//arg 1: Conf File
	//arg 2: O/P Directory
	//arg 3: Contigs
    //arg 4: Raw Reads
	//arg 5: Out Name
	/*if the program is ran without options ,it will show the usage and exit*/
	if(argc == 1)
	{
		showhelpinfo(argv[0]);
		exit(1);
	}
	/*use function getopt to get the arguments with option."hv" indicate 
	that option h,v are the options without arguments while u,p,s are the
	options with arguments*/
	char tmp;
	while((tmp=getopt(argc,argv,"hv"))!=-1)
	{
		switch(tmp)
		{
			/*option h show the help information*/
			case 'h':
				showhelpinfo(argv[0]);
				exit(1);
			break;
			/*option v show the version information*/
			case 'v':
				cout<<"The current version is 1.4.1"<<endl;
				exit(1);
			break;
			/*invaild input will get the help information*/
			default:
				showhelpinfo(argv[0]);
				exit(1);
			break;
		}
	}
	ifstream confFile(argv[1]);
    ifstream ctgRes(argv[3]);
	ifstream readFile(argv[4]);

	if(confFile.fail())
	{
		cerr<<"\nThe configuration file could not be opened:"<<string(argv[1]);
		return -1;
	}
	if(ctgRes.fail())
	{
		cerr<<"\nThe contig file could not be opened:"<<string(argv[3]);
		return -1;
	}
	if(readFile.fail())
	{
		cerr<<"\nThe read file could not be opened:"<<string(argv[4]);
		return -1;
	}
	string usePrevBLAST="",blast_str,conf_str,blast_path="",blast_db_path="",genbank_path="",wget_path="",mpiBLASTPath="",currentPath="",writeIntermideateFiles="NO",blackListFile="";
	currentPath=string(argv[2]);
	int maxLoopCtr=0,mpiThreads=4;
	while(getline(confFile,conf_str))
	{
	    if(conf_str.find("BLASTPATH") != std::string::npos)
	    {
            blast_path = conf_str.substr(conf_str.find("=")+1);
	    }
		 if(conf_str.find("MPIBLAST") != std::string::npos)
	    {
            mpiBLASTPath = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("BLASTDB") != std::string::npos)
	    {
            blast_db_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("LOOPCOUNTER") != std::string::npos)
	    {
            maxLoopCtr = atoi(conf_str.substr(conf_str.find("=")+1).c_str());
	    }
		if(conf_str.find("MPIHOSTCOUNT") != std::string::npos)
	    {
            mpiThreads = atoi(conf_str.substr(conf_str.find("=")+1).c_str());
	    }
		if(conf_str.find("GENBANKPATH") != std::string::npos)
	    {
            genbank_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("WGETPATH") != std::string::npos)
	    {
            wget_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("WRITE_INTERMIDEATE_FILES") != std::string::npos)
	    {
            writeIntermideateFiles = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("USEPREVBLAST") != std::string::npos)
	    {
            usePrevBLAST = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("BLACKLISTFILE") != std::string::npos)
	    {
            blackListFile = conf_str.substr(conf_str.find("=")+1);
	    }
	}
	if (access((blast_path + "blastn").c_str(), X_OK) == -1)
	{
		cerr<<"\nThe BLAST executable could not be found at:" + blast_path + "\n";
		return -1;
	}
	currentPath=currentPath+"/";
	map<string,string> rawReads;
	readFasta(readFile,rawReads);
	vector<string> blackAccCodes = CreateBlackList(blackListFile);
	string fname = currentPath + "ctgBLASTResult.txt";
	if(usePrevBLAST!="1" || !is_file_exist(fname.c_str()))
	{//string cmdStr = blast_path + "blastx -db "+ blast_db_path +"nr -evalue 10 -max_target_seqs 500 -threshold 99 -window_size 4 -gapopen 12 -gapextend 2 -num_threads 16 -query " + string(argv[2]) + " -out ctgBLASTResult.txt -outfmt \"10 qseqid sseqid stitle salltitles qstart qend sstart send\"";
		string cmdStr = blast_path + "blastn -task megablast -db "+ blast_db_path +"nt -evalue 1e-01 -max_target_seqs 100 -num_threads 32 -query " + string(argv[3]) + " -out "+currentPath+"ctgBLASTResult.txt -outfmt \"10 qseqid sacc stitle salltitles qstart qend sstart send\"";
		int ret = system(cmdStr.c_str());
	}
    
    ifstream ctgBlastFile(fname.c_str());
	vector<string> isRelatedCtgBLASTResult;
	string annotationSearchKeywords="transposase,intergrase,insertion sequence";
	isRelatedCtgBLASTResult = findTransposaseHits(ctgBlastFile,blackAccCodes,annotationSearchKeywords,genbank_path,wget_path);
	map<string,string> transHitSeqs = getTransHitSeqs(isRelatedCtgBLASTResult,ctgRes,currentPath);
	map<string,string> hitSeqs;
	map<string,int> hitStarts,hitEnds;
	for(std::map<string,string>::iterator iter = transHitSeqs.begin(); iter != transHitSeqs.end(); iter++) {
		string thisHitSeqName = (*iter).first;
		string thisHitSeq = (*iter).second;
		bool isSubString = false;
		for(std::map<string,string>::iterator iter2 = transHitSeqs.begin(); iter2 != transHitSeqs.end(); iter2++) {
			if((*iter2).second.find(thisHitSeq) != std::string::npos && thisHitSeqName != (*iter2).first) {
				isSubString=true;
				break;
			}
		}
		if(!isSubString)
		{
			hitSeqs.insert(pair<string,string>(thisHitSeqName,thisHitSeq));
			hitStarts.insert(pair<string,int>(thisHitSeqName,0));
			hitEnds.insert(pair<string,int>(thisHitSeqName,thisHitSeq.size()-1));
		}
			
	}
	int loopCtr=0;
	create_directory(currentPath + "FinalResult");
    vector<ISeq> allISSeqs;
    do{
		cout<<"Exe loop:"<<loopCtr<<"\n";
		cout.flush();
		if(writeIntermideateFiles == "YES")
			create_directory(currentPath + "level"+SSTR(loopCtr));
		hitSeqs = expandTransposaseSeqs(hitSeqs,hitStarts,hitEnds,rawReads, allISSeqs, string(argv[4]), blast_path, mpiBLASTPath,mpiThreads, currentPath, loopCtr, string(argv[5]),writeIntermideateFiles);
		loopCtr++;
	}while(hitSeqs.size()>0 && loopCtr < maxLoopCtr);
	sort(allISSeqs.begin(),allISSeqs.end());
	cout<<"Writing Final Result File with "<<allISSeqs.size()<<" sequences. \n";
	for(int i=0;i<allISSeqs.size();i++){
		allISSeqs[i].writeFinalSeq(currentPath,string(argv[5]),allISSeqs,i);
	}
	readFile.close();
    ctgBlastFile.close();
	confFile.close();
	ctgRes.close();
	return 0;
}
