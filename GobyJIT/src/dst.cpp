//***************************************************************************************
//	Using fftw3 package of complex fftw to create dst (discrete sine transformation) 
//***************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "cblas.h"

double error(double *x, double *y, int N);
void initial(double *in, int N);
void dst(double *in, double *out, int N);
void idst(double *in, double *out, int N);

int main()
{
	printf("\n");
	int N, p, q, r;
	printf(" Input size N = 2^p * 3^q * 5^r, (p, q, r) = ");
	scanf("%d %d %d", &p, &q, &r);
	N = pow(2, p) * pow(3, q) * pow(5, r) - 1;
	printf(" N = %d \n\n", N);

	// Input vector *in, and propose to get the output with *out = dst(*in)
	double *in, *out, *out_idst;
	in = (double *) malloc(N*sizeof(double));
	out = (double *) malloc(N*sizeof(double));
	out_idst = (double *) malloc(N*sizeof(double));

	initial(in, N);

	dst(in, out, N);
	idst(out, out_idst, N);

	printf(" The max error of input *in between idst(dst(*in)) \n error = %e \n", error(in, out_idst, N));
	return 0;
}

double error(double *x, double *y, int N)
{
	int i;
	double e, temp;
	e = 0.0;
	for (i=0; i<N; i++)
	{
		temp = fabs(x[i] - y[i]);
		if (temp > e)	e = temp;
	}
	return e;
}

void initial(double *in, int N)
{
	int i;
	for (i=0; i<N; i++)	in[i] = 1.0*i;
}

void dst(double *in, double *out, int N)
{
	int i, L;

	// expand data size to 2*N + 2
	L = 2*N + 2;

	fftw_complex *ex_in, *ex_out;
	ex_in = (fftw_complex *) malloc(L*sizeof(fftw_complex));
	ex_out = (fftw_complex *) malloc(L*sizeof(fftw_complex));

	ex_in[0][0] = ex_in[0][1] = ex_in[N+1][0] = ex_in[N+1][1] = 0.0;
	for (i=0; i<N; i++)
	{
		ex_in[i+1][0] = in[i];
		ex_in[i+1][1] = 0.0;
		ex_in[N+i+2][0] = -1.0*in[N-1-i];
		ex_in[N+i+2][1] = 0.0;
	}

	fftw_plan plan;
	plan = fftw_plan_dft_1d(L, ex_in, ex_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	// After fft(y[k]), Y[k] = fft(y[k]), Sx[k] = i*Y[k+1]/2
	for (i=0; i<N; i++)	out[i] = -1.0*ex_out[i+1][1]/2.0;

	free(ex_in);
	free(ex_out);
}

void idst(double *in, double *out, int N)
{
	dst(in, out, N);

	double s = 2.0/(N+1);
	cblas_dscal(N, s, out, 1);
}

