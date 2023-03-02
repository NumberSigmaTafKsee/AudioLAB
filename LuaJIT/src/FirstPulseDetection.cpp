/***************************************************************
*	Author: Mona Abid    									   *
*	KIT IPE  												   *
*	Date: 04/05/2018 ; 15h01     	     					   *
****************************************************************/


#include "XcorrFFT.h"
#include"mex.h"

void max_index_2Darray_according_1D(float *a, int n, int m, int * start, int *end, int * max_i)
{
	//int TimeComputation = 1;
	mxArray *timeExec[2];
	int visualize = 0;
	float max_v;

	
	for (int j = 0; j < m; j++) {
		

		max_i[j] = start[j] ;// +1 For matlab ... To take it away after test
		max_v = a[start[j] - 1 + j * n];
		
		/*if (visualize == 1) {
			mexPrintf(" j iteration :  %i\n", j);
			mexPrintf(" start[j] :  %i\n", start[j]);
			mexPrintf(" end[j] :  %i\n", end[j]);
			mexPrintf(" max_i intialised :  %i\n", start[j] - 1 + j * n);
			mexPrintf(" max_v intialised :  %f\n", max_v);
		}*/
		//for (int k; k < 100000; k++) {};
		//if (TimeComputation == 1) mexCallMATLAB(0, NULL, 0, NULL, "tic");

		for (int i = start[j]-1 ; i < end[j] ; i++) {// Change the limits


			if (a[j * n + i] > max_v) {
				max_v = a[j * n + i];
				max_i[j] = i+1;	// +1 For matlab ... To take it away after test
			}
	
		}
		//mexCallMATLAB(0, &timeExec[1], 0, NULL, "toc");

		/*if (visualize == 1) {
			mexPrintf(" max_i :  %i\n", max_i[j]);
			mexPrintf(" Ptr max_i :  %p\n", max_i + j);
		}*/
	}

}
//
//
//void max_index_2Darray_according_1D(float *a, int n, int m, float * start, float *end, int * max_i)
//{
//	float max_v;
//	for (int j = 0; j < m; j++) {
//		max_i[j] = (int)start[j];
//		max_v = a[(int)(start[j] * (j + 1))];
//		mexPrintf(" start[j] :  %f\n", start[j]);
//		mexPrintf(" j :  %i\n", j);
//		mexPrintf(" index :  %i\n", (int)(start[j] * (j + 1)));
//		mexPrintf(" end[j] :  %f\n", end[j]);
//		mexPrintf(" max_v :  %f\n", max_v);
//		for (int i = (int)start[j]; i < (int)end[j] + 1; i++) {// Change the limits
//			if (a[j * n + i] > max_v) {
//				max_v = a[j * n + i];
//				max_i[j] = i;
//				mexPrintf(" max_i :  %i\n", max_i[j]);
//				mexPrintf(" Ptr max_i :  %p\n", max_i + j);
//			}
//		}
//	}
//}