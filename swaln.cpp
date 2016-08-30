#include "swaln.h"

double swalign(string seq_a, string seq_b, double gap, double match, double mismatch, int &seq_a_start, int &seq_a_end,int &seq_b_start,int &seq_b_end, int &hitLen){

  // string s_a=seq_a,s_b=seq_b;
  int N_a = seq_a.length();                     // get the actual lengths of the sequences
  int N_b = seq_b.length();
  int rows = N_a+1;
  int cols = N_b+1;
  ////////////////////////////////////////////////

  // initialize H
  double **H = (double **) malloc(sizeof(double *)*rows);
  for(int i=0; i<rows; i++){
    H[i] = (double *) malloc(sizeof(double)*cols);
  }
  // Populate H
  for(int i=0;i<=N_a;i++){
    for(int j=0;j<=N_b;j++){
      H[i][j]=0.;
    }
  }
  double temp[4];
  // Index matrices to remember the 'path' for backtracking
  int **I_i = (int **) malloc(sizeof(int *)*rows);
  for(int i=0; i<rows; i++){
    I_i[i] = (int *) malloc(sizeof(int)*cols);
  }

  int **I_j = (int **) malloc(sizeof(int *)*rows);
  for(int i=0; i<rows; i++){
    I_j[i] = (int *) malloc(sizeof(int)*cols);
  }

  // here comes the actual algorithm
  int ind = 0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1],match,mismatch);
      temp[1] = H[i-1][j]-gap;
      temp[2] = H[i][j-1]-gap;
      temp[3] = 0.;
      H[i][j] = find_array_max(temp,4,ind);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
   	I_i[i][j] = i-1;
	I_j[i][j] = j-1;
	break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
     	I_i[i][j] = i-1;
	I_j[i][j] = j;
	break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
      	I_i[i][j] = i;
	I_j[i][j] = j-1;
	break;
      case 3:                                  // (i,j) is the beginning of a subsequence
      	I_i[i][j] = i;
	I_j[i][j] = j;
	break;
      }
    }
  }

  // search H for the maximal score
  double H_max = 0.;
  int i_max=0,j_max=0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      if(H[i][j]>H_max){
	H_max = H[i][j];
	i_max = i;
	j_max = j;
      }
    }
  }

  // Backtracking from H_max
  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;

  char *consensus_a = (char *) malloc(sizeof(char)*(N_a+N_b+2));
  char *consensus_b = (char *) malloc(sizeof(char)*(N_a+N_b+2));

  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){

    if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
    else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

    if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
    else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick++;
    }
	//cout<<"\nA,"<<i_max<<","<<current_i;
	//cout<<"\nB,"<<j_max<<","<<current_j;
	//cout<<"\nTick:"<<tick<<",Score:"<<H_max;
  // Output of the consensus motif to the console
  //cout<<endl<<"***********************************************"<<endl;
  //cout<<"The alignment of the sequences"<<endl<<endl;
  //for(int i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<endl;
  //for(int i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<endl<<endl;
  //cout<<"is for the parameters  mu = "<<mismatch<<" and delta = "<<gap<<" given by"<<endl<<endl;
  //for(int i=tick-1;i>=0;i--) cout<<consensus_a[i];
  //cout<<endl;
  //for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
  //cout<<endl;

  seq_a_start = current_i;
  seq_a_end = i_max;
  seq_b_start = current_j;
  seq_b_end = j_max;
  hitLen = tick;

  free(consensus_a);
  free(consensus_b);

  for(int i=0;i<rows;i++)
    free(H[i]);
  free(H);

  for(int i=0;i<rows;i++)
    free(I_i[i]);
  free(I_i);

  for(int i=0;i<rows;i++)
    free(I_j[i]);
  free(I_j);

  return H_max;
} // END of main

/////////////////////////////////////////////////////////////////////////////

double similarity_score(char a,char b,double match, double mismatch){
  double result;
  if(a==b){
      result=match;
    }
  else{
      result=-mismatch;
    }
  return result;
}

/////////////////////////////////////////////////////////////////////////////

double find_array_max(double array[],int length, int &ind){
  double max = array[0];            // start with max = first element
  ind = 0;
  for(int i = 1; i<length; i++){
      if(array[i] > max){
	max = array[i];
	ind = i;
      }
  }
  return max;                    // return highest value in array
}
