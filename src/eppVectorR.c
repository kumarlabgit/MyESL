#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <armadillo>
#include <math.h>
//#include "matrix.h"

/*
   min  1/2 ( ||x- u||_2^2 + ||t-v||_2^2 )
   s.t.  ||x_j||_2 <= t_j

 */

//std::pair<arma::vec, arma::vec> eppVectorR(const arma::vec& u, const arma::vec& v, const arma::vec& ind, int n, int k);
//void eppVectorR(double *x, double * t, double * u, double * v, double * ind, int n, int k)
std::pair<arma::vec, arma::vec> eppVectorR(arma::vec& u_in, const arma::vec& v_in, const arma::uvec& ind_in, int n, int k){
	double *x;
	x = (double*) malloc(n*sizeof(double));
	double *t;
	t = (double*) malloc(n*sizeof(double));

	const double* u = u_in.memptr();
	const double* v = v_in.memptr();
	const unsigned long long int* ind = ind_in.memptr();

    int i, j;
    double temp;

	/* compute the 2 norm of each group
	*/

	for(j=0;j<k;j++){
		temp=0;
		for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
			temp+= u[i]* u[i];
        temp=sqrt(temp);
        /*temp contains the 2-norm of of each row of u*/

        if(temp > fabs(v[j])){
           t[j]=(temp + v[j])/2;

           for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
               x[i]= t[j] / temp * u[i];
        }
        else
           if(temp <= v[j]){
               t[j]=v[j];

               for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
                   x[i]= u[i];
            }
            else{
                t[j]=0;

               for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
                   x[i]=0;
            }

	}
//    std::pair<arma::vec, arma::vec> vec_pair = {arma::vec x_vec(&x[0], n), arma::vec t_vec(&t[0], n)};
    std::pair<arma::vec, arma::vec> vec_pair = {arma::vec(&x[0], n), arma::vec(&t[0], n)};
    free(x);
    free(t);
    return vec_pair;
}