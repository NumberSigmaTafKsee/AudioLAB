/***************************************************************
*	Author: Mona Abid    									   *
*	KIT IPE  												   *
*	Date: 19/03/2018 ; 16h32		     					   *
****************************************************************/
#include "XcorrFFT.h"
#include "mex.h"


/*------------------------------------------ MEX OF THE CROSS-CORRELATION --------------------------------------------*/
int main()
{
	return 0;
}
/*The getway function that replaces the main*/
//-----------------Working Version-------------------------
//void mexFunction(int nlhs, mxArray *plhs[],
//	int nrhs, const mxArray*prhs[]) {
//
//	int visualize = 0;
//
//	 //Input arguments
//	float * tab1_r, *tab2_r;/*input arrays 1xn1 and 1xn1*/
//	int n1, m; /*arrays size*/
//	int nthreads, RESAMPLE_FACTOR;
//
//	//Output arguments
//	float* Result_r, *Result_i; /*Output array*/
//
//					 /* create a pointer to the real data in the input array  */
//	tab1_r = (float *)mxGetPr(prhs[0]);
//	tab2_r = (float *)mxGetPr(prhs[1]);
//
//	n1 = mxGetM(prhs[0]);
//	m = mxGetN(prhs[0]);
//
//	Result_r = (float *)mxGetPr(prhs[2]);
//	RESAMPLE_FACTOR= (int)mxGetScalar(prhs[3]);
//
//
//	nthreads = 0;
//	if (nrhs > 4) {
//		nthreads = (int)mxGetScalar(prhs[4]);
//	}
//
//
//	if (nthreads == -1) {
//		//void fftw_cleanup_threads(void);
//		fftw_cleanup();
//		mexPrintf("Clean up FFTW ressources and threads - no leaking of ressources (on recompile etc) ");
//		return;
//	}
//	if (visualize == 1) {
//		mexPrintf(" nthreads:  %i\n", nthreads); //1 
//		mexPrintf(" n:  %i\n", n1);//10000
//		mexPrintf(" m:  %i\n", m);//10000
//	}
//
//	Myinterp_fft_Debug(tab1_r, tab2_r, n1, m, Result_r, nthreads, RESAMPLE_FACTOR);
//
//
//}
#include<time.h>
/***************************************************************
*	Author: Mona Abid    									   *
*	KIT IPE  												   *
*	Last Version: xcorr + upsampling + max_index               *
*	Last modified: 11/06/2018 ; 13h02		     			   *
****************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[]) {

	int visualize = 0;

	//Input arguments
	float * tab1_r, *tab2_r;/*input arrays 1xn1 and 1xn1*/
	int n1, m,nResampled; /*arrays size*/
	int nthreads, RESAMPLE_FACTOR;
	mxArray *timeExec[2];

	int TimeComputation = 1;
	//Output arguments
	float* Result_r, *Result_i; /*Output array*/
	int *lag_max_i;
	int * start, *end;
	double time_spent;
/* create a pointer to the real data in the input array  */
	tab1_r = (float *)mxGetPr(prhs[0]);
	tab2_r = (float *)mxGetPr(prhs[1]);
	Result_r = (float *)mxGetPr(prhs[2]);
	n1 = mxGetM(prhs[0]);
	m = mxGetN(prhs[0]);
	nResampled = mxGetM(prhs[2]);

	RESAMPLE_FACTOR = (int)mxGetScalar(prhs[3]);
	start = (int *)mxGetPr(prhs[4]);
	end = (int *)mxGetPr(prhs[5]);
	lag_max_i = (int *)mxGetPr(prhs[6]);

	//time_spent = (double)mxGetScalar(prhs[7]);
	nthreads = 0;
	if (nrhs > 7) {
		nthreads = (int)mxGetScalar(prhs[7]);
	}


	//Myinterp_fft_DebugTHREADS(tab1_r, tab2_r, n1, m, Result_r, nthreads, RESAMPLE_FACTOR);
	
	if (visualize == 1) {
		mexPrintf(" nthreads:  %i\n", nthreads); //1 
		mexPrintf(" n:  %i\n", n1);//10000
		mexPrintf(" m:  %i\n", m);//10000
	}
	//for (int i; i < 100000; i++) {};
	//if (TimeComputation == 1) mexCallMATLAB(0, NULL, 0, NULL, "tic");
	//clock_t begintime = clock();
	Myinterp_fft_Debug(tab1_r, tab2_r, n1, m, Result_r, nthreads, RESAMPLE_FACTOR);
	//clock_t endtime = clock();
	//time_spent = (double)(endtime - begintime) / CLOCKS_PER_SEC;
	//if (TimeComputation == 1) {
		//mexCallMATLAB(0, &timeExec[0], 0, NULL, "toc");
		//for (int i; i < 10000; i++) {};
		//mexCallMATLAB(0, NULL, 0, NULL, "tic");
	//}
	max_index_2Darray_according_1D(Result_r, nResampled, m, start, end, lag_max_i);
	//if (TimeComputation == 1) mexCallMATLAB(0, &timeExec[1], 0, NULL, "toc");


	//float max_v;
	//	for (int j = 0; j < m; j++) {
	//		max_v = Result_r[ j * n1 ];
	//			for (int i = 0; i < n1; i++) {// Change the limits
	//				if (Result_r[j * n1 + i] > max_v) {
	//					max_v = Result_r[j * n1 + i];
	//					lag_max_i[j] = i + 1;	// +1 For matlab ... To take it away after test
	//				}

	//			}
	//
	//	}


}

