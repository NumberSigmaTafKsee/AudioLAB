// *===========================================================================* //
//
//  __::((xodVAFilter.cpp))::__
//
//  ___::((XODMK Programming Industries))::___
//  ___::((XODMK:CGBW:BarutanBreaks:djoto:2020))::___
//
//
//	Purpose: C++ implementation of Virtual Analog Filters
//			 IIR Toplogy-Preserving Transform (TPT) Filters
//			 "The Art Of VA Filter Design" - Vadim Zavalishin
//			 VA filters: Moog Ladder 4-pole
//
//	Revision History: Feb 08, 2017 - initial
//	Revision History: Mar 10, 2020 - current
//
// *===========================================================================* //


#include <iostream>
#include <math.h>

#include "xodVAFilter_base.h"
#include "xodVAFilter.h"



// *--------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass with Feedback ---* //

void onePoleTPTFB_LP::setAlpha_LP(float alpha) {
	fAlpha = alpha;
}

// TPT 1-pole filter - Zavalishin p46 (the Art of VA Design)
void onePoleTPTFB_LP::doFilterStage_LP(float xn, float& z1fb, float& ynLP) {
	float v = (xn - z1)*fAlpha;
	ynLP = v + z1;
	z1 = ynLP + v;
	z1fb = z1;

	//std::cout<<"xn_LPR = "<<xn<<",      v_LPR = "<<v<<",      LP_LPR = "<<ynLP<<",      z1_LPR = "<<z1<<std::endl;

}



// *--------------------------------------------------------* //
// *--- Stereo TPT Moog Half Ladder Low-Pass filter ---* //


void xodMoogLadder4P::initialize(float newSampleRate) {

	fAlpha0 = 0;
	fBeta1 = 0;
	fBeta2 = 0;
	fBeta3 = 0;
	fBeta4 = 0;
	G = 0;
	K = 0;

	LPF1.initialize_LP(newSampleRate);
	LPF2.initialize_LP(newSampleRate);
	LPF3.initialize_LP(newSampleRate);
	LPF4.initialize_LP(newSampleRate);

}

void xodMoogLadder4P::setFcAndRes(float cutoff, float resonance, float sampleRate) {

	// prewarp for BZT
	double wd = 2*pi*cutoff;
	double T  = 1/(double)sampleRate;
	double wa = (2/T)*tan(wd*T/2);
	double g  = wa*T/2;

	// G - the feedforward coeff in the VA One Pole
	float G = g/(1.0 + g);

	LPF1.setAlpha_LP(G);
	LPF2.setAlpha_LP(G);
	LPF3.setAlpha_LP(G);
	LPF4.setAlpha_LP(G);


	fBeta1 = (G*G*G) / (1.0 + g);
	fBeta2 = (G*G) / (1.0 + g);
	fBeta3 = G / (1.0 + g);
	fBeta4 = 1 / (1.0 + g);


	// calculate alpha0
	// for 2nd order, K = 2 is max so limit it there
	// ** fixed-point implementation, K=2 requires many integer bits to prevent overflow
	// (future enhancement -> use internal data type to handle bit-growth)
	// currently limit K resonance to less than 2.0 to prevent overflow:
	K = resonance;
	if(K > 2.0)
		K = 2.0;

	fAlpha0 = 1.0 / (1.0 + K*G*G*G*G);


	std::cout<<std::endl<<"TB_G: "<<G<<std::endl;
	std::cout<<"TB_K: "<<K<<std::endl;
	//std::cout<<"fAlpha0 = "<<fAlpha0<<std::endl;
	//std::cout<<"setFcAndRes_Ref: g = "<<g<<",        G = "<<G<<",        K = "<<K<<",        beta1 = "<<fBeta1<<",        beta2 = "<<fBeta2<<",        beta3 = "<<fBeta3<<std::endl;

}


void xodMoogLadder4P::advance(float xn, float& yn) {

	static float SM = 0;
	//static float SM_R = 0;

	static float yn_LP1;
	static float yn_LP2;
	static float yn_LP3;
//	static float ynR_LP2;

	static float z1fb_1;
	static float z1fb_2;
	static float z1fb_3;
	static float z1fb_4;

//	static float z1fb_R1;
//	static float z1fb_R2;
//	static float z1fb_R3;

	SM = fBeta1*z1fb_1 + fBeta2*z1fb_2 + fBeta3*z1fb_3 + fBeta4*z1fb_4;		// sum 4 internal Z1 states

	float un = fAlpha0*(xn - K*SM);

	//std::cout<<"xnL = "<<xnL<<",	fb1 = "<<z1fb_L1<<",	fb2 = "<<z1fb_L2<<",	fb3 = "<<z1fb_L3<<",	fAlpha0 = "<<fAlpha0<<",	SM_L = "<<SM_L<<",	unL_REF = "<<unL<<std::endl;

	LPF1.doFilterStage_LP(un, z1fb_1, yn_LP1);
	LPF2.doFilterStage_LP(yn_LP1, z1fb_2, yn_LP2);
	LPF3.doFilterStage_LP(yn_LP2, z1fb_3, yn_LP3);
	LPF4.doFilterStage_LP(yn_LP3, z1fb_4, yn);


	//std::cout<<"SM_REF = "<<SM_L<<",        K_REF = "<<K<<",        unL_REF = "<<unL<<std::endl;

}

// *--------------------------------------------------------* //