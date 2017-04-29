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


//functions for the potential term
double potential(double p, double w) {
  return w*w*p*p;
}

double potentialIn(double p, double w) {
  return w*w*p*p+1/p;
}

//find max off diag element
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

vec jacobi_method ( mat A, mat &R, int n ) {
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

  //rotate while an off diagonal element is greater than the tolerance
  while ( fabs (max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max_offdiag = maxoffdiag( A, &k, &l, n );
    A = rotate ( A, R, k, l, n );
    iterations++;
  }
  cout << "Number of iterations: " << iterations << "\n";
  
  vec eigval = zeros<vec>(n);
 
  for (int i=0;i<n;i++){
    eigval(i) = A(i,i);
  }
  //sort eigenvalues
  vec eigval2(n);
  eigval2 = sort(eigval);
  vec index(n);
  //find order of eigenvectors
  for (int i=0;i<n;i++) {
    uvec test(n);
    test =  find (eigval2 == eigval(i));
    index(i) = test(0);
  }
  //fix order of eigenvectors
  mat eigvec = zeros<mat>(n,n); 
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      eigvec(j,i) = R(j,index(i));
    }
  }
  R = eigvec;
  return eigval;
}

//function that tests unitary transformations - orthoganility and dot product check
void rotation(mat A, int n) {

    vec identity1(n);
    for (int i=0; i<n; i++) {
      double one = cdot(A.col(i),A.col(i));
      identity1[i] = one;
    }
    
    cout<<"v^{T}*v: "<<identity1<<endl;
    mat rot = zeros<mat>(n,n);

    double theta = kPi/2;
    
    //create a unitary transformation matrix
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
    
    vec p = zeros<vec>(n);
    vec V = zeros<vec>(n);
    
    //Set max r value and calculate step size
    double Rmin = 0; double Rmax = 10;
    mat A = zeros<mat>(n,n); 
    double h = (Rmax-Rmin)/n;
    double hh = h*h;
    cout<<"h: "<<h<<endl;


    //frequency
    double w = 1;
  

    vec eigval = zeros<vec>(n);   
    mat eigvec = zeros<vec>(n);
    
    mat R = zeros<mat>(n,n);

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
    for (double w : {0.01,0.25,1.0,5.0}) {
      cout<<"w = "<<w<<endl;
   
      //Initialize the potential
    for (int i=0;i<n;i++){
      p(i) = Rmin + h*(i+1);
      V(i) = potential(p(i),w);
    }

    //Initialize the hamiltonian matrix
    A(0,0) = 2/hh+V(0);
    A(0,1) = -1/hh;
    A(n-1,n-2) = -1/hh;
    A(n-1,n-1) = 2/hh+V(n-1);

    for (int i=1;i<n-1;i++) {
      A(i,i-1) = -1/hh;
      A(i,i) = 2/hh+V(i);
      A(i,i+1) = -1/hh;
    }
    

    //run for non-interacting with armadillo eig_sym
    /*
      start = clock();
      eig_sym(eigval,eigvec,A);
      finish = clock();
      
      cout<< "Time for armadillo eig_sym, non interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
      cout<<"Eigenvalues: "<<endl;
      cout<<eigval(0)<<endl;
      cout<<eigval(1)<<endl;
      cout<<eigval(2)<<endl;
    */

      //run for non-interacting with jacobi rotation algorithm
      start = clock();
      vec C = jacobi_method(A,R,n);
      finish = clock();
      
      cout<< "Time for jacobi_method, non interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
         
      cout<<"Eigvalue 1: "<<C(0)<<endl;
      cout<<"2: "<<C(1)<<endl;
      cout<<"3: "<<C(2)<<endl;
      
    
      //uncomment for output to files
      /*
      ofile.open(filename+to_string(count));
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i = 0; i < n;i++) {
        ofile << setw(15) << setprecision(8) << p(i);
	//R(i,0) for jacobi, eigvec(i,0) for arma
        ofile << setw(15) << setprecision(8) << R(i,0)<<endl;//eigvec(i,0) << endl;	
      }
      ofile.close();      
      count++;
      */

      // change potential to interacting case
      for (int i=0;i<n;i++){
        p(i) = Rmin + h*(i+1);
        V(i) = potentialIn(p(i),w);
      }
   
      A(0,0) = 2/hh+V(0);
      A(n-1,n-1) = 2/hh+V(n-1);
      for (int i=1;i<n-1;i++) {
        A(i,i) = 2/hh+V(i);
      }
    
         
     //run for interacting with armadillo eig_sym
      /*
      start = clock();
      eig_sym(eigval,eigvec,A);
      finish = clock();
           
      cout<< "Time for armadillo eig_sym, interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
      
      cout<<"Eigenvalues: "<<endl;
      cout<<eigval(0)<<endl;
      cout<<eigval(1)<<endl;
      cout<<eigval(2)<<endl;
      */
      
     //run for interacting with jacobi rotation algorithm
      start = clock();
      vec D = jacobi_method(A,R,n);
      finish = clock();
    
      cout<< "Time for jacobi_method, interacting: "<< (double)(finish - start)/((double)CLOCKS_PER_SEC)<<endl;
      
      
      cout<<"Eigvalue 1: "<<D(0)<<endl;
      cout<<"2: "<<D(1)<<endl;
      cout<<"3: "<<D(2)<<endl;
      
      //uncomment for output to files
      /*
      ofile.open(filename+to_string(count));
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i = 0; i < n;i++) {
        ofile << setw(15) << setprecision(8) << p(i);
	//R(i,0) for jacobi, eigvec(i,0) for arma
        ofile << setw(15) << setprecision(8) << R(i,0)<<endl;//eigvec(i,0) << endl;
      }   
      ofile.close();
      count++;
      */
      
    }

  return 0;
}
