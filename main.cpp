#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <algorithm>   // std:: max_element
#include <armadillo>

using namespace std;
using namespace arma;

// function declarations

double off_diagonal(mat &, int *, int *, int);
void rotation (mat &, int, int, int n);

int main(int argc, char* argv[])
{
    int n = atoi(argv[1]);
    int i, k, l;

    double rho_min = 0.0;
    double rho_max = 5.0;
    double h = (rho_max - rho_min)/(n+1);
    double epsilon = 1.0e-8;
    double aij;

    double* rho = new double[n];
    double* V = new double [n];


    mat A = zeros(n,n);

    for(k=0; k<n; k++) {

        rho[k] = rho_min + (k+1)*h;
        V[k] = rho[k]*rho[k];

        for (l=0; l<n; l++){

            if (k==l) {
                A(k,l) = (2./(h*h)) + V[k];
            }
            else if (k==l-1) {
                A(k,l) = -1/(h*h);
            }
            else if (k ==l+1) {
                A(k,l) = -1/ (h*h);
            }
            else {
                A(k,l) = 0.0;
            }
        }
    }


    aij = off_diagonal(A, &k, &l, n);

    while (epsilon < aij) {
        rotation(A, k, l, n);
        aij = off_diagonal(A, &k, &l, n);
    }

    vec lambda = zeros(n);
    for (i=0; i<n; i++) {
        lambda[i] =(A(i,i));
    }
    lambda = sort(lambda);
    cout << "Three lowest eigenvalues: " << lambda(0) << "  "  << lambda(1) << "   "  <<  lambda(2) << endl;
    return 0;
}

// max value upper off-diagonal

double off_diagonal (mat & A, int * k, int * l, int n ) {

    double max_value = 0.0;
    int i, j;

    for (i=0; i<n-1; i++) {
        for (j=i+1; j<n-1; j++) {
            if (fabs(A(i,j)) > max_value) {
                max_value = fabs(A(i,j));
                *k = i;
                *l = j;
            }
        }
    }

    return max_value;
}

void rotation (mat & A, int k, int l, int n){
    double t, tau, s, c;
    int i;

    // Finding tan, cos, sin

    tau = (A(l,l) - A(k,k)) / (2.0*A(k,l));

    if (tau < 0) {
        t = -1.0 / (-tau + sqrt(1+tau*tau));
    }

    if (tau > 0) {
        t = 1.0 / (tau + sqrt(1+tau*tau));
    }

    c = 1.0 / sqrt(1+t*t);
    s = t*c;


    // Performing the rotation

    double a_ll, a_kk, a_ik, a_il;
    a_kk = A(k,k);
    a_ll = A(l,l);

    A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2.0*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    for (i=0; i<n; i++) {
        if ( i != k &&  i!= l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);
        }

    }
}

