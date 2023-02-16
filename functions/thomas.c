#include "mex.h"
#include <stdlib.h>
#include <complex.h>

void thomas(const mxArray * A,const mxArray * B,const mxArray * C,const mxArray * D, mxArray * X, size_t N){

    mwSize i; 

    # Typecasting to perform complex operations nativelly
	double complex *a = (double complex *)mxGetComplexDoubles(A);
    double complex *b = (double complex *)mxGetComplexDoubles(B);
    double complex *c = (double complex *)mxGetComplexDoubles(C);
    double complex *d = (double complex *)mxGetComplexDoubles(D);
    double complex *x = (double complex *)mxGetComplexDoubles(X);

	double complex *cprime = mxCalloc(N, sizeof(double complex));

	cprime[0] = c[0]/b[0];
	x[0] = d[0]/b[0];

	double complex m;

	for (i = 1; i < N; i++) {
		x[i] = d[i];
		m = 1.0 / (b[i] - a[i] * cprime[i - 1]);
		cprime[i] = c[i] * m;
		x[i] = (x[i] - a[i] * x[i - 1]) * m;
	}

	for (i = N - 1; i-- > 0; ){
		x[i] = (x[i] - cprime[i] * x[i + 1]);
	}

	mxFree(cprime);

};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){


    size_t ncols;  


    ncols = mxGetM(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix( 1, (mwSize)ncols, mxCOMPLEX);

    // mexPrintf("\n%d\n",ncols);
    
    thomas(prhs[0], prhs[1], prhs[2], prhs[3], plhs[0], (size_t)ncols);
}

