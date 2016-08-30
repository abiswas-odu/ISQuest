#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
 
using namespace std;

double similarity_score(char a,char b,double match, double mismatch);
double find_array_max(double array[],int length, int &ind);
double swalign(string seq_a, string seq_b, double gap, double match, double mismatch, int &seq_a_start, int &seq_a_end,int &seq_b_start,int &seq_b_end, int &hitLen);