//PROJECT 2

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <cstdlib>
#include <math.h>
#include "time.h"
//#include "jacobi.h"
// use namespace for output and input
using namespace std;
using namespace arma;

const double kPi = 3.1415926535897;

ofstream ofile;

/*
The function int comp()
is a utility function for the library function qsort()
to sort double numbers after increasing values.
*/
int comp(const double *val_1, const double *val_2)
{
if((*val_1) <= (*val_2)) return -1;
else if((*val_1) > (*val_2)) return +1;
else
return 0;
} // End: function comp()

double maxoffdiag ( mat A, int * k, int * l, int n ) {
  double max = 0.0;
  for ( int i = 0; i < n; i++ ) {
    for ( int j = i + 1; j < n; j++ ) {
      if ( fabs(A(i,j)) > max ) {
	max = fabs(A(i,j));
	*l = i;
	*k = j;
      }
    }
  }
  return max;
}

// Function to find the values of cos and sin
mat rotate ( mat A, mat &R, int k, int l, int n ) {
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau =( A(l,l) - A(k,k))/(2*A(k,l));
    if ( tau >= 0 ) {
      t = 1/( tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1/(-tau + sqrt(1.0 + tau*tau));
    }
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  // changing the matrix elements with indices k and l
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0; // hard-coding of the zeros
  A(l,k) = 0.0;
  // and then we change the remaining elements
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
    // Finally, we compute the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
  return A;
}

vec* jacobi_method ( mat A, mat &R, int n ) {
  // Setting up the eigenvector matrix
  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < n; j++ ) {
      if ( i == j ) {
	R(i,j) = 1.0;
      } else {
	R(i,j) = 0.0;
      }
    }
  }
  int k, l;
  double epsilon = 1.0e-10;
  double max_number_iterations = (double) n * (double) n * (double) n;
  int iterations = 0;
  double max_offdiag = maxoffdiag ( A, &k, &l, n );
  //cout<<"fabs (max_offdiag): "<<fabs(max_offdiag)<<" > epsilon"<<epsilon<<endl;
  //cout<<"Iterations: "<<iterations<<" < max_number_iterations: "<<max_number_iterations<<endl;
  while ( fabs (max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max_offdiag = maxoffdiag( A, &k, &l, n );
    A = rotate ( A, R, k, l, n );
    iterations++;
  }
  cout << "Number of iterations: " << iterations << "\n";
  //cout<<"A: "<<A<<endl;

  vec eigval[n];
  for (int i=0;i<n;i++){
    for (int j=0;i,n;i++){
      if (i==j) {
	eigval[i] = A(i,j);
      }
    }
  }

  qsort (eigval, n, sizeof(double), (int(*)(const void *,const void *))comp);

  return eigval;
}


void rotation(mat A, int n) {

    vec identity1(n);
    for (int i=0; i<n; i++) {
      double one = cdot(A.col(i),A.col(i));
      identity1[i] = one;
    }
    
    cout<<"v^{T}*v: "<<identity1<<endl;

    mat rot = zeros<mat>(n,n);

    double theta = kPi/2;
    
    rot(0,0) = cos(theta);
    rot(0,1) = -sin(theta);
    rot(1,0) = sin(theta);
    rot(1,1) = cos(theta);
    for (int i=2;i<n;i++)  rot(i,i) = 1;
    
    mat B = zeros<mat>(n,n);
    
    for (int i=0; i<n; i++) {
      vec C = A.col(i);
      B.col(i) = solve(rot,A.col(i));
    }
   
    vec identity2(n);
    for (int i=0; i<n; i++) {
      identity2[i] = cdot(B.col(i),B.col(i));
    }
 
    cout<<"w^{T}*w: "<<identity2<<endl;
    
    return;
}


int main(int argc, char *argv[]){
  int n; //number of mesh points
    string filename;

    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and number of mesh points" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        n = atoi(argv[2]);
    }
    
    clock_t start, finish;
    
    vec p(n); 
    
    p[0] = 0; p[n] = 10;
    mat A = zeros<mat>(n,n); 
    double h = (p[n]-p[0])/n;
    double hh = h*h;
    cout<<"h: "<<h<<endl;


    for (int i=1;i<n;i++) p[i] = h*i;
 
    A(0,0) = 2/hh;
    A(0,1) = -1/hh;
    A(n-1,n-2) = -1/hh;
    A(n-1,n-1) = 2/hh+p[n-1]*p[n-1];

    for (int i=1;i<n-1;i++) {
      A(i,i-1) = -1/hh;
      A(i,i) = 2/hh+p[i]*p[i];
      A(i,i+1) = -1/hh;
    }
    
    vec eigval;
    mat eigvec;
    
    //call armadillo function to find eigenvalues of symmetric matricies
    //start = clock();
    //eig_sym(eigval,eigvec,A);
    //finish = clock();

    //cout<< "Time for armadillo eig_sym, non-interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;

    //cout<<"A: "<<A<<endl;
    //cout<<"Eigenvalues: "<<eigval<<endl;
    // cout<<"Eigenvectors: "<<eigvec<<endl;
    
    //  rotation(eigvec, n);

    mat R = zeros<mat>(n,n);

    //start = clock();
    //mat B = jacobi_method(A,R,n);
    //finish = clock();
    
    //cout<< "Time for jacobi_method, non-interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;

    mat R_test = zeros<mat>(3,3);

    mat test = zeros<mat>(3,3);

    test(0,0) = 1;
    test(1,0) = test(0,1) = 2;
    test(1,1) = 1;
    test(2,0) = test(0,2) = 3;
    test(2,1) = test(1,2) = 4;
    test(2,2) = 5;
    
    cout<<"Test, "<<test<<" has eigenvalues of: "<<jacobi_method(test,R_test,3);
    
    int count = 0;
    int maxdim [4] = {50,10,5,5};
    int j=0;
    

    //WHAT IS WRONG
    //
    //

    double w = 1;
    
    for (int i=1;i<n;i++) p[i] = h*i;
    
      cout<<"frequency w = "<<w<<endl;
      for (int i=1;i<n-1;i++) {
	A(i,i-1) = -1/hh;
	A(i,i) = 2/hh+w*w*p[i]*p[i]+1/p[i];
	A(i,i+1) = -1/hh;
      }

      start = clock();
      vec* D = jacobi_method(A,R,n);
      finish = clock();

      // cout<<"Eigenvectors: "<<endl;
      //cout<<R<<endl;
      cout<<"Eigenvalues: "<<endl;
      cout<<D<<endl;
      //cout<<"Eigvalue 1: "<<D(0,0)<<endl;
      //cout<<"2: "<<D(1,1)<<endl;
      //cout<<"3: "<<D(2,2)<<endl;

      //pg 230, qsort ch 7.8










    /*
    vec p(n); 

    p[0] = 0; p[n] = 50;
    mat A = zeros<mat>(n,n); 
    double h = (p[n]-p[0])/n;
    double hh = h*h;
    cout<<"h: "<<h<<endl;
    
    for (int i=1;i<n;i++) p[i] = h*i;
    */
    /*
    for (double w : {0.01,0.25,1.0,5.0}) {

      /*
    p.zeros();
    p[0] = 0; p[n] = maxdim[j];
    //mat A = zeros<mat>(n,n); 
    double h = (p[n]-p[0])/n;
    double hh = h*h;
    cout<<"h: "<<h<<endl;
    j++;

      
    for (int i=1;i<n;i++) p[i] = h*i;
      */
    /*
      cout<<"frequency w = "<<w<<endl;
      for (int i=1;i<n-1;i++) {
	A(i,i-1) = -1/hh;
	A(i,i) = 2/hh+w*w*p[i]*p[i];
	A(i,i+1) = -1/hh;
      }
      /*
      start = clock();
      eig_sym(eigval,eigvec,A);
      finish = clock();
      
      //cout<<"Eigenvalues, non-interacting: "<<eigval<<endl;
      //cout<<"Eigenvectors, interacting: "<<eigvec<<endl;

      cout<< "Time for armadillo eig_sym, non interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
     */
    /*
      start = clock();
      mat C = jacobi_method(A,R,n);
      finish = clock();
      
      cout<< "Time for jacobi_method, non interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
      
      cout<<"Eigvalue 1: "<<C(0,0)<<endl;
      cout<<"2: "<<C(1,1)<<endl;
      cout<<"3: "<<C(2,2)<<endl;
      
      ofile.open(filename+to_string(count));
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i = 0; i < n;i++) {
        ofile << setw(15) << setprecision(8) << p[i];
        ofile << setw(15) << setprecision(8) << R(i,0)<<endl;//eigvec(i,0) << endl;
	
      }
      ofile.close();
      count++;

      for (int i=1;i<n-1;i++) {
	A(i,i) += 1/p[i];
      }
      /*
      start = clock();
      eig_sym(eigval,eigvec,A);
      finish = clock();

      cout<< "Time for armadillo eig_sym, interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
    *//*
      start = clock();
      mat D = jacobi_method(A,R,n);
      finish = clock();
    
      cout<< "Time for jacobi_method, interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
      
      cout<<"Eigvalue 1: "<<D(0,0)<<endl;
      cout<<"2: "<<D(1,1)<<endl;
      cout<<"3: "<<D(2,2)<<endl;

      ofile.open(filename+to_string(count));
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i = 0; i < n;i++) {
        ofile << setw(15) << setprecision(8) << p[i];
        ofile << setw(15) << setprecision(8) << R(i,0)<<endl;//eigvec(i,0) << endl;
      }
      ofile.close();
      count++;
      
    }
*/
  return 0;
}
