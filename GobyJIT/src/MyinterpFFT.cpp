/***************************************************************
*	Author: Mona Abid    									   *
*	KIT IPE  												   *
*	Date: 17/05/2018 ; 13h47		     					   *
*	Last Modified : 29/06/2018	; 17h00						   *
****************************************************************/
#include "XcorrFFT.h"
#include <windows.h> 
#include <fstream>
#define REAL 0
#define IMAG 1



/*------------------------ Interpolation using Fourier domain + xcorr ------------------------------------------------------
* This function calculates the cross correlation between 2 input data arrays tab1 and tab2 by considering their    *
* length respectively n1 and n2. The output data is stored in result.										        *
* NOTE: The code is limited to the case the arrays tab1 and tab2 have the same size and containg at least 2 values *
* The expected format of the input is a LINE ARRAY VECTOR
--------------------------------------------------------------------------------------------------------------------*/
#include"mex.h"


void Myinterp_fft_Debug(float *tab1_r, float *tab2_r, int n1, int m, float *Result_r, int nthreads, int RESAMPLE_FACTOR)

{
	const fftwf_r2r_kind kind1[2] = { (fftwf_r2r_kind)0, (fftwf_r2r_kind)1 }; // R2HC -- cf fftw3.h :  FFTW_R2HC=0, FFTW_HC2R=1


	FILE *pFile, *fid, *fid_;
	int wisdom_type = FFTW_MEASURE; // FFTW_MEASURE; //FFTW_EXHAUSTIVE; // , FFTW_PATIENT;
	char filename[500] = "";
	int visualize = 0;
	int benchmark = 1;
	static int Index;
	static int old_n1 = -1;
	static int nthreads_old = -1;
	static int old_m = -1;
	static int fftw_threads_initalized = -1;
	static char name[35] = " ";
	char namenew_[500] = "";
	float powerFactor = 1;
	
	unsigned int wisdomtable[4] = { 0, 8, 32,64 }; //correspondes to MEASURE, Exhaustive, patient, estimate, 


	int n = n1;
	int nresample_r = n1*RESAMPLE_FACTOR;///for real from hermetian symmetry for upsample 2 -> stays the same (!)

	std::complex<float> ComplauxResult;
	if ((nthreads != nthreads_old) & (fftw_threads_initalized != -1)) {
		fftwf_cleanup_threads(); if (visualize == 1) mexPrintf("#threads cleaned\n");
	}// ------ thread number changed -> clean up -------
	if ((nthreads > 0) & (fftw_threads_initalized == -1)) {
		fftwf_init_threads(); fftw_threads_initalized = nthreads; if (visualize == 1) mexPrintf("#threads initalized\n");
		if (fftwf_init_threads() == 0) {
			mexPrintf("NO ENOUGH MEMORY"); return;
			//Could think of a recursive function that allows to reduce the ressources
		}
	} // -------- first thread init -------
	  // ------- Check the inputData pointer ---------
	if (tab1_r == NULL) { mexPrintf("NULL pointer -1!"); return; }
	if (tab2_r == NULL) { mexPrintf("NULL pointer 0!"); return; }
	// ------- Check the Alignement ---------
	if (((int)tab1_r & 0xF) != 0) { mexPrintf(" nonaligned 1\n"); return; }
	if (((int)tab2_r & 0xF) != 0) { mexPrintf(" nonaligned 2\n"); return; }
	if (((int)Result_r & 0xF) != 0) { mexPrintf(" nonaligned 3\n"); return; }

	static float *  tab3_r = NULL;
	static float * output1_r = NULL;
	static float * output2_r = NULL;

	static fftwf_plan plan1 = NULL;
	static fftwf_plan plan2 = NULL;
	static fftwf_plan plan3 = NULL;

	

	if ((n1 != old_n1) | (nthreads_old != nthreads) | (m != old_m))
		//reset cache
	{
		if (visualize == 1) mexPrintf("oldm:%d mi:%d oldn1:%d n1:%d nhread_old:%d nthread:%d", old_m, m, old_n1, n1, nthreads_old, nthreads);

		old_n1 = n1;
		nthreads_old = nthreads;
		old_m = m;
		if (visualize == 1) mexPrintf("mi:%d n1:%d nthread:%d", m, n1, nthreads);

		if (tab1_r == NULL) { mexPrintf("NULL pointer 6!"); return; }
		if (tab2_r == NULL) { mexPrintf("NULL pointer 7!"); return; }

		if (output1_r != NULL) free(output1_r);
		if (output2_r != NULL) free(output2_r);
		if (tab3_r != NULL) free(tab3_r);
	

		if (plan1 != NULL) fftwf_destroy_plan(plan1);
		if (plan2 != NULL) fftwf_destroy_plan(plan2);
		if (plan3 != NULL) fftwf_destroy_plan(plan3);

		tab3_r = (float *)calloc(m *  nresample_r, sizeof(float));
		output1_r = (float *)calloc(m *  nresample_r, sizeof(float));// Data is stored in adifferent way //Hermitien Symmetry => Mirror frequencies not allocated -- +1 for the even case : considering the sampling half frequency 
		output2_r = (float *)calloc(m *  nresample_r, sizeof(float)); // Hermitien Symmetry=> Mirror frequencies not allocated -- +1 for the even case : considering the sampling half frequency 


		if (output1_r == NULL) { mexPrintf("NULL pointer 8!"); return; }
		if (output2_r == NULL) { mexPrintf("NULL pointer 9!"); return; }


		if (nthreads > 0) {
			fftwf_plan_with_nthreads(nthreads); if (visualize == 1) mexPrintf("#threads = %i\n", nthreads);
		} // Called only once 

		if (name[0] == ' ')
		{
			system("hostname > hostname.txt");
			fid_ = fopen("hostname.txt", "r");//W:\\abid\\USCT\\Source_code\\TestSaving
			fscanf(fid_, "%s", name);
			fclose(fid_);

		}



		//--------------CREATING AND SAVING INTO LOCAL SUBFOLDER_WISDOM_FILES--------------!!
		int l=n;
		//if(RESAMPLE_FACTOR == 1) l= 0;
		int k;
		do {
			for (int kk = 3; kk >= 0; kk--) {// Evalute the best wisdom available in the folder
				k = wisdomtable[kk];
				if ((benchmark != 1) | (k == wisdom_type)) {  // Enforce the usage of the selected widom type (no better wisdom even if existing for benchmarking)
					sprintf(namenew_, "W:\\abid\\USCT\\Source_code\\TestWisdom\\%s_fftw_wisdom_%d_%d_%d_%d.fftw", name, l, m, nthreads, k);

					fid = fopen(namenew_, "r"); //check if wisdom for this problem already exist
					if (fid != NULL) {
						if (visualize == 1) mexPrintf("reuse wisdom from file");
						fclose(fid);
						if (fftwf_import_wisdom_from_filename(namenew_) == 0) {
							mexPrintf("wisdom import error"); std::remove(namenew_);
						}
						else {
							wisdom_type = k;  // import succeded, updated wisdom considered to create the plans
							break; //end loop End search because best is found
						}
					}
				}
			}

			l = l + (nresample_r - n); // first run l=n, second run l=nresample_r, third would be more

		} while (l == nresample_r);


		//if (visualize == 1) mexPrintf("wisdom_type= %d\n", wisdom_type);
		if (visualize == 1) mexPrintf("n:%d m:%d tab1_r_Ptr:%p output1_Ptr:%p nresample_r:%d wisdom_type:%d", n, m, tab1_r, output1_r, nresample_r, wisdom_type);

mexPrintf("First Element Adress: %p\n Second Element Adress: %p\n First Element Data: %f\n Second Element Data: %f\n", &tab1_r[0],&tab1_r[1],tab1_r[0],tab1_r[1]);
		plan1 = fftwf_plan_many_r2r(1, &n, m, tab1_r, NULL, 1, n, output1_r, NULL,  1, nresample_r, kind1, wisdom_type);//|| FFTW_UNALIGNED|| FFTW_PRESERVE_INPUT
		plan2 = fftwf_plan_many_r2r(1, &n, m, tab2_r, NULL, 1, n, output2_r, NULL,  1, nresample_r, kind1, wisdom_type); // FFTW_WISDOM_ONLY |
		plan3 = fftwf_plan_many_r2r(1, &nresample_r, m, tab3_r, NULL, 1, nresample_r, Result_r, NULL, 1, nresample_r, kind1+1, wisdom_type); //nresample_r or n for the 7th arguement, same result
		

		if (visualize == 1) mexPrintf("nresample_r : %i Ptr tab3:%o  Ptr Result_r :  %p", nresample_r, tab3_r, Result_r);
		l = n;
		//if (RESAMPLE_FACTOR == 1) l = 0;
		do {
			sprintf(namenew_, "W:\\abid\\USCT\\Source_code\\TestWisdom\\%s_fftw_wisdom_%d_%d_%d_%d.fftw", name, l, m, nthreads, wisdom_type);
			fid = fopen(namenew_, "r"); //check if wisdom for this problem already exist
										//Add one second for saving the plans on the hardware !!!!!!!!!!!!!!!!!!!!!!!!
			if (fid == NULL) {
				fftwf_export_wisdom_to_filename(namenew_);
				Sleep(1000);
				if (fftwf_export_wisdom_to_filename(namenew_) == 0) { mexPrintf("wisdom save error"); }
				fftwf_forget_wisdom();
				//mexPrintf("kill wisdom after saving!");
			}
			else
				fclose(fid);
			l = l + (nresample_r - n); // first run l=n, second run l=nresample_r, third would be more
			if (visualize == 1) mexPrintf("L=%d",l);
		} while (l == nresample_r);


		if (plan1 == NULL) { mexPrintf("NULL pointer 10!"); return; }
		if (plan2 == NULL) { mexPrintf("NULL pointer 11!"); return; }
		if (plan3 == NULL) { mexPrintf("NULL pointer 12!"); return; }
	}


	//Create a DFT plan
	fftwf_execute(plan1);
	
	fftwf_execute(plan2);

	/*Calculus of  fft(tab1)* conj(fft(tab2)) */
	float Normalization;
	int indexREAL = 0;
	int indexIMAG = 0;
	int indexIMAGShifted = 0;
	//DataLayout: http://www.fftw.org/fftw3_doc/The-Halfcomplex_002dformat-DFT.html
	for (int i = 0; i < m; i++) {
		
			for (int j = 1; j < n / 2; j++) //running through the frequency bins (beside 0 and SF/2)
			{
				
				indexREAL = j + i* (nresample_r); 
				indexIMAG = indexREAL + n - 1 - j - j + 1;
				indexIMAGShifted = nresample_r - n + 1 + indexIMAG - 1;
					
				Normalization = (sqrt(output1_r[indexREAL] * output1_r[indexREAL] + output1_r[indexIMAG] * output1_r[indexIMAG]));
				tab3_r[indexREAL] = (output1_r[indexREAL] * output2_r[indexREAL] + output1_r[indexIMAG] * output2_r[indexIMAG]);// / Normalization;
				tab3_r[indexIMAGShifted] = (output1_r[indexIMAG] * output2_r[indexREAL] - output1_r[indexREAL] * output2_r[indexIMAG]);// / Normalization;
				}
				
	
	//j= 0; 0 --Frequency
	indexREAL = i* (nresample_r);
	Normalization = (sqrt(output1_r[indexREAL] * output1_r[indexREAL]));
	tab3_r[indexREAL] = (output1_r[indexREAL] * output2_r[indexREAL]);// / Normalization;
	
	//j = n/2; FS/2--Frequency
	indexREAL = n / 2 + i* (nresample_r);
	Normalization = (sqrt(output1_r[indexREAL] * output1_r[indexREAL]));
	tab3_r[indexREAL] = ((output1_r[indexREAL] * output2_r[indexREAL])) / 2;// *Normalization);
	if (RESAMPLE_FACTOR == 1)
		tab3_r[indexREAL] = tab3_r[indexREAL] * 2; 

	}



	float scale = float(1.0 / nresample_r);

	//Create a IDFT plan
	fftwf_execute(plan3);


	for (int i = 0; i < nresample_r * m; i++)
		Result_r[i] *= scale;


}


