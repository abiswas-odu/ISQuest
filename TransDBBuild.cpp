#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<cstring>
#include<dirent.h>
#include<omp.h>
#include<ctime>
#include <unistd.h>

extern "C" {
	#include "gbfp.h"
}

using std::istringstream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;
using std::cin;
using std::cerr;
using namespace std;

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

void initBlkList(vector<string> &blackListAcc)
{
	//blackListAcc.push_back("CP000896");
	blackListAcc.push_back("AP011136");
	blackListAcc.push_back("CP005080");
	blackListAcc.push_back("CP001331");
	blackListAcc.push_back("CP001398");
	
	blackListAcc.push_back("NR_041897");
	blackListAcc.push_back("NR_041901");
	blackListAcc.push_back("NR_041905");
	
	blackListAcc.push_back("NR_043588");
	blackListAcc.push_back("NR_043904");
	blackListAcc.push_back("NR_043905");
	blackListAcc.push_back("M65227");
}
static char *getQualValue(char *sQualifier, gb_feature *ptFeature) {
    static char *null = NULL;
    gb_qualifier *i;

    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++) 
        if (strcmp(sQualifier, i->sQualifier) == 0)
            return i->sValue;
    return null;
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
/*function that show the help information*/
void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument1]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-v  show version information"<<endl;
  cout<<"         "<<"[argument1]  Path to Genbank file"<<endl;
  cout<<"example: "<<s<<" /research/gbkFiles"<<endl;
}
int main(int argc, char **argv)
{
	//ISQuest
	//arg 1: *.gbk file path
	
	if(argc == 1)
	{
		showhelpinfo(argv[0]);
		exit(1);
	}
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
				cout<<"The current version is 1.0"<<endl;
				exit(1);
			break;
			/*invaild input will get the help information*/
			default:
				showhelpinfo(argv[0]);
				exit(1);
			break;
		}
	}
	
	string locPath = string(argv[1]);
	vector<string> blackListAcc;
	initBlkList(blackListAcc);
	DIR *pDIR;
	char locusField[] = "locus_tag";
	char productField[] = "product";
	struct dirent *entry;
	string annotationSearchKeywords="transposase,intergrase,insertion,mobile";
	std::vector<std::string> searchKeys = split(annotationSearchKeywords, ',');
	vector<string> fileNames;
	
	if(pDIR=opendir(locPath.c_str())){
		while(entry = readdir(pDIR)){
			if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){
				string fname = string(entry->d_name);
				fileNames.push_back(fname);
			}
		}
	}
	closedir(pDIR);
	const clock_t begin_time = clock();
	#pragma omp parallel for num_threads(10)
	for(int i=0;i<100;i++)
	{
		string baseCode = fileNames[i].substr(0,fileNames[i].find("."));
		string filePath = locPath + "/" + fileNames[i];
		string cmdStr= "gzcat " + filePath + " >" + baseCode;  
		int ret = system(cmdStr.c_str());
	
		ofstream ofile;
		string oFileName =  "output/"+baseCode+".fasta";
		ofile.open(oFileName.c_str(), ios::out | ios::trunc);
		
		gb_data **pptSeqData, *ptSeqData;
		gb_feature *ptFeature;
		
		char * cstr = new char [baseCode.length()+1];
		std::strcpy(cstr, baseCode.c_str());
		pptSeqData = parseGBFF(cstr);
		for (int i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
			for (int j = 0; j < ptSeqData->iFeatureNum; j++) {
				ptFeature = (ptSeqData->ptFeatures + j);
				if(isAnnotationHit(searchKeys,ptFeature) && ptSeqData->sSequence)
				{
					ofile<<">"<<ptSeqData->sAccession<<"|"<<ptSeqData->sOrganism<<"|"<<ptFeature->sFeature<<"|"<<ptFeature->cDirection<<"|"<<ptFeature->lStart<<"|"<<ptFeature->lEnd;
					if(getQualValue(locusField,ptFeature)!=NULL)
						ofile<<"|"<<getQualValue(locusField,ptFeature);
					if(getQualValue(productField,ptFeature)!=NULL)
						ofile<<"|"<<getQualValue(productField,ptFeature);
					ofile<<"\n";
					ofile<<getSequence(ptSeqData->sSequence,ptFeature)<<"\n";
				}
			}
		}
		freeGBData(pptSeqData);
		delete[] cstr;
		cmdStr= "rm " + baseCode;  
		ret = system(cmdStr.c_str()); 
		ofile.close();
	}
	cout <<"Total time taken:"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<"s\n\n";
	return 0;
}