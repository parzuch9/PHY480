/*
PHY 480 - Computational Physics - Project 1

Solve the one-dimensional Poisson equation with Dirichlet boundary conditions 
by rewritting it as a set of linear equations
*/

#include<iostream>
using std::cout; using std::endl;
#include<cstdio>
using std::atoi;

int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);

  double ** AllocateMatrix(int, int);
  void DeallocateMatrix(double **, int, int); 

  cout<<"n: "<<n<<endl;
  return 0;
}

// Allocate memory for a matrix and initialize the elements to zero

double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}

// Free memory

void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}

