#include <cmath>
#include <iostream>
#include <stdlib.h>     /* atoi */
#include <vector>
#include <cstdlib>
#include <time.h>       /* time */
#include <chrono>

using namespace std;
using namespace std::chrono;

void print(vector< vector<double> > mat);
void fillRandom(vector < vector<double> >& mat);
void multiply(vector < vector<double> > matA, vector < vector<double> > matB,  vector < vector<double> >& result);

int main(int argc, char** argv)
{

  //Arguments handeling
  if(argc<=1){
    cout << "\nNot enough arguments\nUsage\n---------------------------\n1-number of samples per experiment\n";
    return 0;
  }
  int iterations = atoi(argv[1]);

  cout << "n\t\tsamples\t\tavgTime(ms)\t\tstd(ms)\n";
  for(int n=200; n<=2000; n+=200)
  {

    cout << n << "\t\t" << iterations << "\t\t";

    double sum = 0, square_sum=0;

    for(int iterI=0; iterI<iterations; iterI++){
      // initialize matrixes
      srand(time(NULL));
      vector< vector<double> > matA(n, vector<double>(n,0));
      vector< vector<double> > matB(n, vector<double>(n,0));
      vector< vector<double> > matC(n, vector<double>(n,0));

      // Fill Arrays with ramdom values
      fillRandom(matA);
      fillRandom(matB);

      // Call Multiply
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      multiply(matA, matB, matC);
      high_resolution_clock::time_point t2 = high_resolution_clock::now();

      // Calculate time difference
      auto duration = duration_cast<microseconds>( t2 - t1 ).count();

      //Add to sum and square sum
      sum+=duration;
      square_sum+= duration*duration;
    }

    double mean = sum/iterations ;
    double mean_square = square_sum/iterations ;
    printf ("%4.2f\t\t\t%4.2f\n", mean/1000, sqrt(mean_square - mean*mean)/1000 );
    // print(matC);

  }
	return 0;
}

void print(vector< vector<double> > mat){
  int n = mat.size();
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++)
     cout << mat[i][j] << "\t";
    cout << endl;
  }
}

void fillRandom(vector < vector<double> >& mat){
  int n = mat.size();
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        mat[i][j] = rand()%10;
    }
  }
}

// Multiplication Function
void multiply(vector < vector<double> > matA,
  vector < vector<double> > matB,  vector < vector<double> >& result){
  int n = matA.size();
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        double sum = 0;
        for(int k=0; k<n; k++){
          sum+= matA[i][k]*matB[k][j];
        }
        result[i][j] = sum;
    }
  }
}
