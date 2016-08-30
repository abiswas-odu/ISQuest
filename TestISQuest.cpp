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

#include "nwaln.h"
#define MIN_GLOBAL_MATCHES_THRESH 0.95

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

int main(int argc, char **argv)
{
	//ISQuest Testing Module
	//arg 1: ISQuest Final Result 
	//arg 2: Known IS 
	ifstream isQuestFile(argv[1]);
	ifstream knownISFile(argv[2]);

	if(isQuestFile.fail() || knownISFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	map<string,string> isQuest;
	readFasta(isQuestFile,isQuest);

	map<string,string> knownIS;
	readFasta(knownISFile,knownIS);
	int alignFound=0, alignNotFound=0; 
	for(std::map<string,string>::iterator iter = knownIS.begin(); iter != knownIS.end(); iter++) {
		string thisISName = (*iter).first;
		string thisISSeq = (*iter).second;
		bool isAligned = false;
		for(std::map<string,string>::iterator iter2 = isQuest.begin(); iter2 != isQuest.end(); iter2++) {
			string thisFoundName = (*iter2).first;
			string thisFoundSeq = (*iter2).second;
			int matches=0,mismatches=0,gaps=0,hitLen=0;
			cout<<"GAlign:"<<thisISName<<" with "<<thisFoundName<<"\n";
			double score = nwalign(thisISSeq,thisFoundSeq,matches,mismatches,gaps,hitLen);
			double minMatchesCtr1 = 0,minMatchesCtr2=0;
			minMatchesCtr1 = thisISSeq.size()*MIN_GLOBAL_MATCHES_THRESH;
			minMatchesCtr2 = thisFoundSeq.size()*MIN_GLOBAL_MATCHES_THRESH;
			if(matches >= minMatchesCtr1 && matches >= minMatchesCtr2)
			{
				cout<<"Global alignment found\n";
				alignFound++;
				isAligned=true;
				break;
			}
		}
		if(!isAligned)
		{
			alignNotFound++;
		}
	}
	cout<<"IS Found:"<<alignFound<<"\nIS Not Found:"<<alignNotFound;
	isQuestFile.close();
	knownISFile.close();
	return 0;
}