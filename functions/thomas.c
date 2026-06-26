// Compile it with mex -R2018a CFLAGS='$CFLAGS -fopenmp -ffast-math' LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS='-O3 -DNDEBUG' thomas.c
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "mex.h"

void thomas(
	const mxArray * A,
	const mxArray * B, 
	const mxArray * C, 
	const mxArray * Tminus, 
	const mxArray * GammaMinus, 
	const mxArray * Tplus, 
	const mxArray * V, 
	mxArray * X, size_t N, size_t M){

    mwSize i, i0, i1, k, n;

    // Bi1                                 + GammaMinus.*v    +  Typlus(M+1:end).*v(M+1:end)
    // Tyminus(1:end-M).*v(1:end-M)        + GammaMinus.*v    +  BiN
 
    // b = complex(vertcat(Bi1,Tyminus(1:end-M).*v(1:end-M)) + ...
    // GammaMinus.*v + ...
    // vertcat(Typlus(M+1:end).*v(M+1:end),BiN));


    // d[i] = Tyminus[i]*v[(i-M)%N] + GammaMinus[i]*v[i] + Typlus[i]*v[(i+M)%N]

    // Typecasting to perform complex operations nativelly
	float complex *a = (float complex *)mxGetComplexSingles(A);
    float complex *b = (float complex *)mxGetComplexSingles(B);
    float complex *c = (float complex *)mxGetComplexSingles(C);

	float complex *tminus = (float complex *)mxGetComplexSingles(Tminus);
    float complex *gammaMinus = (float complex *)mxGetComplexSingles(GammaMinus);
    float complex *tplus = (float complex *)mxGetComplexSingles(Tplus);
    float complex *v = (float complex *)mxGetComplexSingles(V);

    float complex *x = (float complex *)mxGetComplexSingles(X);

	float complex m;
    float complex d;

	float complex *cprime;
	cprime = mxCalloc(N, sizeof(float complex));



    // Set number of threads
    omp_set_num_threads(omp_get_max_threads());


    // Thomas algorithm. Exploiting natural parallelization of this large linear system
	#pragma omp parallel for private(m, d)
	for (n = 0; n < M; n++){

		i0 = n*M;
    	
		d = tminus[i0]*v[(i0-M)%N] + gammaMinus[i0]*v[i0] + tplus[i0]*v[(i0+M)%N];

		// First solution
		cprime[i0] = c[i0]/b[i0];
		x[i0] = d/b[i0];
 
		// First point is n*M
		// Start this loop on the next point (n*M+1)
		// Go until the last point of this line at index (n+1)*M-1
		for (i = i0 + 1; i < (n+1)*M; i++) {
			x[i] = tminus[i]*v[(i-M)%N] + gammaMinus[i]*v[i] + tplus[i]*v[(i+M)%N];;
			m = 1.0 / (b[i] - a[i] * cprime[i - 1]);
			cprime[i] = c[i] * m;
			x[i] = (x[i] - a[i] * x[i - 1]) * m;
		}
 
		for (i = (n+1)*M - 1; i-- > i0; ){
			x[i] = (x[i] - cprime[i] * x[i + 1]);
		}
	}
	mxFree(cprime);

};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

    size_t ncols, mcols;  

    ncols = mxGetM(prhs[0]);
    mcols = (size_t) mxGetScalar(prhs[7]);
    // (mwSize m, mwSize n, mxClassID classid, mxComplexity ComplexFlag);
    plhs[0] = mxCreateNumericMatrix( 1, (mwSize)ncols, mxSINGLE_CLASS, mxCOMPLEX);
    
    thomas(
    	prhs[0], // A
    	prhs[1], // B
    	prhs[2], // C
    	prhs[3], // Tminus
    	prhs[4], // GammaMinus
    	prhs[5], // Tplus
    	prhs[6], // V
    	plhs[0], // X
    	(size_t)ncols,  // M*M
    	(size_t)mcols); // M
}



