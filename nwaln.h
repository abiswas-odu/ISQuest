#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>

#define GLOBAL_GAP 2.0
#define GLOBAL_MATCH 5.0
#define GLOBAL_MISMATCH 9.026

using namespace std;

double global_similarity_score(char a,char b);
double global_find_array_max(double array[],int length, int &ind);
double nwalign(string seq_a, string seq_b, int &matchesCtr, int &mismatchesCtr,int &gapCtr, int &hitLen);
