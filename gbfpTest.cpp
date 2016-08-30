#include<cstring>
#include<string>
#include<fstream>
extern "C" {
	#include "gbfp.h"
}
#include "gtest/gtest.h"

namespace {

using namespace std;

// The fixture for testing. 
class GbfpTesting : public ::testing::Test {  
 public:
	gb_data **pptSeqData, *ptSeqData;
    gb_feature *ptFeature;
	virtual void SetUp() {
	}

	virtual void TearDown() {
	}
};

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

TEST_F (GbfpTesting, parse) {
	string accFilePath="testData/gbcon145";
	ifstream gbkFile(accFilePath.c_str());

	if(gbkFile)
    {
        char * cstr = new char [accFilePath.length()+1];
        std::strcpy(cstr, accFilePath.c_str());
        pptSeqData = parseGBFF(cstr);
        delete[] cstr;
    }
	char locusField[] = "locus_tag";
	char productField[] = "product";
	string annotationSearchKeywords="transposase,intergrase,insertion,mobile";
	std::vector<std::string> searchKeys = split(annotationSearchKeywords, ',');
	for (int i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
		for (int j = 0; j < ptSeqData->iFeatureNum; j++) {
			ptFeature = (ptSeqData->ptFeatures + j);
			if(isAnnotationHit(searchKeys,ptFeature) && ptSeqData->sSequence)
			{
				cout<<">"<<ptSeqData->sAccession<<"|"<<ptSeqData->sOrganism<<"|"<<ptFeature->sFeature<<"|"<<ptFeature->cDirection<<"|"<<ptFeature->lStart<<"|"<<ptFeature->lEnd;
				if(getQualValue(locusField,ptFeature)!=NULL)
					cout<<"|"<<getQualValue(locusField,ptFeature);
				if(getQualValue(productField,ptFeature)!=NULL)
					cout<<"|"<<getQualValue(productField,ptFeature);
				cout<<"\n";
				getSequence(ptSeqData->sSequence,ptFeature);
			}
		}
	}
	freeGBData(pptSeqData);
}

}  // namespace
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}