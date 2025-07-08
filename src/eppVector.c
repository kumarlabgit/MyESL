#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <armadillo>
#include <math.h>
//#include "matrix.h"
#include "epph.h" /* This is the head file that contains the implementation of the used functions*/


/*
 Lp Norm Regularized Euclidean Projection

        min  1/2 ||x- v||_2^2 + rho * ||x||_p

 Usage (in Matlab):
 [x, c, iter_step]=epp(v, n, rho, p, c0);

 Usage in C:
 epp(x, c, iter_step, v, n, rho, p, c0);

 The function epp implements the following three functions
 epp1(x, v, n, rho) for p=1
 epp2(x, v, n, rho) for p=2
 eppInf(x, c, iter_step, v,  n, rho, c0) for p=inf
 eppO(x, c, iter_step, v,   n, rho, p) for other p

------------------------------------------------------------

  Here, the input and output are of Vector form.


 Written by Jun Liu, May 18th, 2009
 For any problem, please contact: j.liu@asu.edu

 */

//arma::vec eppVector(const arma::vec& v, const arma::vec& ind, int k, int n, const arma::vec& lambda_over_L_times_gWeight, double q);
//void eppVector(double * x, double * v, double * ind, int k, int n, double * rho, double p)
arma::vec eppVector(arma::vec& v_in, const arma::uvec& ind, int k, int n, arma::vec rho, double p){
	double *x;
	x = (double*) malloc(n*sizeof(double));
	double* v = v_in.memptr();
    int i, *iter_step;
    double c0, c;
	double *px, *pv;

    iter_step=(int *)malloc(sizeof(int)*2);

    c0=0;
    for(i=0; i<k; i++){

        px=x+(int)ind[i];
		pv=v+(int)ind[i];

        epp(px, &c, iter_step, pv, (int)(ind[i+1]-ind[i]), rho[i], p, c0);

    }
    arma::vec x_vec(&x[0], n);
    free(x);
    free(iter_step);
    return x_vec;
}


