// *===========================================================================* //
//
//  __::((xodVAFilter_base.h))::__
//
//  ___::((XODMK Programming Industries))::___
//  ___::((XODMK:CGBW:BarutanBreaks:djoto:2020))::___
//
//
//	Purpose: C++ header for Virtual Analog Filters
//			 IIR Toplogy-Preserving Transform (TPT) Filters
//			 "The Art Of VA Filter Design" - Vadim Zavalishin
//
//	Revision History: Feb 08, 2017 - initial
//	Revision History: Mar 10, 2020 - current
//
// *===========================================================================* //

#ifndef __XODVAFILTER_BASE_H__
#define __XODVAFILTER_BASE_H__


#include <math.h>



const double pi = 3.141595926536;
const double Fs = 48000;						// Fixed clock frequency

// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass Model ---* //

class onePoleTPT_LP {
public:

protected:
	// controls
	float G;			// cutoff
	float sampleRate;	// fs
	float z1;			// z-1 register

public:
	inline void initialize_LP(float newSampleRate) {
		sampleRate = newSampleRate;
		z1 = 0;
	}

	float getSampleRate_LP(){return sampleRate;}
	float getZ1regValue_LP(){return z1;}
	void setFc_LP(float fc);
	void doFilterStage_LP(float xn, float& ynLP);
};



// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT High-Pass Model ---* //

class onePoleTPT_HP {
public:

protected:
	// controls
	float G;			// cutoff
	float sampleRate;	// fs
	float z1;			// z-1 register

public:
	inline void initialize_HP(float newSampleRate) {
		sampleRate = newSampleRate;
		z1 = 0;
	}

	float getSampleRate_HP(){return sampleRate;}
	float getZ1regValue_HP(){return z1;}
	void setFc_HP(float fc);
	void doFilterStage_HP(float xn, float& ynHP);
};


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass + High-Pass Model ---* //

class onePoleTPT_LPHP {
public:

protected:
	// controls
	float G;			// cutoff
	float sampleRate;	// fs
	float z1;			// z-1 register

public:
	inline void initialize_LPHP(float newSampleRate) {
		sampleRate = newSampleRate;
		z1 = 0;
	}

	float getSampleRate_LPHP(){return sampleRate;}
	float getZ1regValue_LPHP(){return z1;}
	void setFc_LPHP(float fc);
	void doFilterStage_LPHP(float xn, float& ynLP, float& ynHP);
};


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT All-Pass Model ---* //

class onePoleTPT_AP {
public:

protected:
	// controls
	float G;			// cutoff
	float sampleRate;	// fs
	float z1;			// z-1 register

public:
	inline void initialize_AP(float newSampleRate) {
		sampleRate = newSampleRate;
		z1 = 0;
	}

	float getSampleRate_AP(){return sampleRate;}
	float getZ1regValue_AP(){return z1;}
	void setFc_AP(float fc);
	void doFilterStage_AP(float xn, float& ynAP);
};

// *---------------------------------------------------------------------------* //




#endif // __XODVAFILTER_BASE_H__
