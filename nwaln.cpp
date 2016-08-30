#include "nwaln.h"

double nwalign(string seq_a, string seq_b, int &matchesCtr, int &mismatchesCtr,int &gapCtr, int &hitLen){

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
  double temp[3];
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
      temp[0] = H[i-1][j-1]+global_similarity_score(seq_a[i-1],seq_b[j-1]);
      temp[1] = H[i-1][j]-GLOBAL_GAP;
      temp[2] = H[i][j-1]-GLOBAL_GAP;
      H[i][j] = global_find_array_max(temp,3,ind);
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
      }
    }
  }

  // search last row and column of H for the maximal score
  double H_max = H[N_a][N_b];
  int i_max=N_a,j_max=N_b;
  for(int i=1;i<=N_a;i++){
        if(H[i][N_b]>H_max){
            H_max = H[i][N_b];
            i_max = i;
            j_max = N_b;
        }
  }
  for(int j=1;j<=N_b;j++){
        if(H[N_a][j]>H_max){
            H_max = H[N_a][j];
            i_max = N_a;
            j_max = j;
        }
  }


  // Backtracking from H_max
  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;
  matchesCtr=0;
  mismatchesCtr=0;
  gapCtr=0;

  char *consensus_a = (char *) malloc(sizeof(char)*(N_a+N_b+2));
  char *consensus_b = (char *) malloc(sizeof(char)*(N_a+N_b+2));

  while((next_j!=0) && (next_i!=0)){

    if(next_i==current_i)
    {
        consensus_a[tick] = '-';                  // deletion in A
        gapCtr++;
    }
    else
    {
        consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
        if(seq_a[current_i-1] == seq_b[current_j-1]) matchesCtr++;
        else   mismatchesCtr++;
    }

    if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
    else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick++;
    }
	//cout<<"\nTick:"<<tick<<",Score:"<<H_max;
  // Output of the consensus motif to the console
  //cout<<endl<<"***********************************************"<<endl;
  //cout<<"The global alignment of the sequences"<<endl<<endl;
  //for(int i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<endl;
  //for(int i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<endl<<endl;
  //cout<<"is for the parameters  mu = "<<GLOBAL_MISMATCH<<" and delta = "<<GLOBAL_GAP<<" given by"<<endl<<endl;
  //for(int i=tick-1;i>=0;i--) cout<<consensus_a[i];
  //cout<<endl;
  //for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
  //cout<<endl<<"Match:"<<matchesCtr<<",Mismatch:"<<mismatchesCtr<<",Gaps:"<<gapCtr<<endl;
  //cout<<endl;

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

double global_similarity_score(char a,char b){
  double result;
  if(a==b){
      result=GLOBAL_MATCH;
    }
  else{
      result=-GLOBAL_MISMATCH;
    }
  return result;
}

/////////////////////////////////////////////////////////////////////////////

double global_find_array_max(double array[],int length, int &ind){
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
