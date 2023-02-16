// *===========================================================================* //
//
//  __::((xodVAFilter.h))::__
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


#ifndef __XODVAFILTER_H__
#define __XODVAFILTER_H__


#include <iostream>
#include <math.h>

#include "xodVAFilter_base.h"


// *--------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass with Feedback ---* //

class onePoleTPTFB_LP {
public:

protected:
	// controls
	float fAlpha;			// cutoff
	float sampleRate;	// fs
	float z1;			// z-1 register

public:
	inline void initialize_LP(float newSampleRate) {
		sampleRate = newSampleRate;
		z1 = 0;
	}

	void setAlpha_LP(float alpha);
	void doFilterStage_LP(float xn, float& z1fb, float& ynLP);
};


// *--------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass Reference Model ---* //

// class onePoleTPT_LPRef2 {
// public:

// protected:
// 	// controls
// 	float fAlpha;			// cutoff
// 	float sampleRate;	// fs
// 	float z1;			// z-1 register

// public:
// 	inline void initialize_LPRef(float newSampleRate) {
// 		sampleRate = newSampleRate;
// 		z1 = 0;
// 	}

// 	//float getSampleRate_LPRef(){return sampleRate;}
// 	//float getZ1regValue_LPRef(){return z1;}
// 	void setAlpha_LPRef(float alpha);
// 	void doFilterStage_LPRef(float xn, float& z1fb, float& ynLP);
// };


// *--------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass Reference Model ---* //

// class onePoleTPT_LPRef3 {
// public:

// protected:
// 	// controls
// 	float fAlpha;			// cutoff
// 	float sampleRate;	// fs
// 	float z1;			// z-1 register

// public:
// 	inline void initialize_LPRef(float newSampleRate) {
// 		sampleRate = newSampleRate;
// 		z1 = 0;
// 	}

// 	//float getSampleRate_LPRef(){return sampleRate;}
// 	//float getZ1regValue_LPRef(){return z1;}
// 	void setAlpha_LPRef(float alpha);
// 	void doFilterStage_LPRef(float xn, float& z1fb, float& ynLP);
// };


// *--------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass Reference Model ---* //

// class onePoleTPT_LPRef4 {
// public:

// protected:
// 	// controls
// 	float fAlpha;			// cutoff
// 	float sampleRate;	// fs
// 	float z1;			// z-1 register

// public:
// 	inline void initialize_LPRef(float newSampleRate) {
// 		sampleRate = newSampleRate;
// 		z1 = 0;
// 	}

// 	//float getSampleRate_LPRef(){return sampleRate;}
// 	//float getZ1regValue_LPRef(){return z1;}
// 	void setAlpha_LPRef(float alpha);
// 	void doFilterStage_LPRef(float xn, float& z1fb, float& ynLP);
// };




// *--------------------------------------------------------* //



// *--------------------------------------------------------* //
// *--- Stereo Moog Ladder 4-pole Filter ---* //

class xodMoogLadder4P {
public:

	onePoleTPTFB_LP LPF1;
	onePoleTPTFB_LP LPF2;
	onePoleTPTFB_LP LPF3;
	onePoleTPTFB_LP LPF4;

protected:
	// controls
	float sampleRate;	// fs

	// for loop scalar
	float fAlpha0;

	float fBeta1;
	float fBeta2;
	float fBeta3;
	float fBeta4;

	float G;

	float K;

public:
	void initialize(float newSampleRate);
	void setFcAndRes(float cutoff, float resonance, float sampleRate);
	void advance(float xn, float& yn);
};

// *--------------------------------------------------------* //



#endif // __XODVAFILTER_H__
