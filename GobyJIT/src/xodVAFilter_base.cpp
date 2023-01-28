// *===========================================================================* //
//
//  __::((xodVAFilter_base.cpp))::__
//
//  ___::((XODMK Programming Industries))::___
//  ___::((XODMK:CGBW:BarutanBreaks:djoto:2020))::___
//
//
//	Purpose: C++ implementation of Virtual Analog Filters
//			 IIR Toplogy-Preserving Transform (TPT) Filters
//			 "The Art Of VA Filter Design" - Vadim Zavalishin
//			 VA 1-pole base filters: LP, HP, LP+HP, AP	
//
//	Revision History: Feb 08, 2017 - initial
//	Revision History: Mar 10, 2020 - current
//
// *===========================================================================* //

#include <iostream>
#include <math.h>

#include "xodVAFilter_base.h"


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass Model ---* //

void onePoleTPT_LP::setFc_LP(float fc) {
	// pre-warp the cutoff - bilinear-transform filters
	float wd = 2*pi*fc;
	float T = 1/sampleRate;
	float wa = (2/T)*tan(wd*T/2);
	float g = wa*T/2;

	// calculate big G value - Zavalishin p46 (the Art of VA Design)
	// 0.0 < G < 1.0 (from simulation)
	G = g/(1.0 + g);

	std::cout<<"TPT G = "<<G<<std::endl;
}

void onePoleTPT_LP::doFilterStage_LP(float xn, float& ynLP) {
	float v = (xn - z1)*G;
	ynLP = v + z1;
	z1 = ynLP + v;

}


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT High-Pass Model ---* //

void onePoleTPT_HP::setFc_HP(float fc) {
	// prewarp the cutoff - bilinear-transform filters
	float wd = 2*pi*fc;
	float T = 1/sampleRate;
	float wa = (2/T)*tan(wd*T/2);
	float g = wa*T/2;

	// calculate big G value - Zavalishin p46 (the Art of VA Design)
	// 0.0 < G < 1.0 (from simulation)
	G = g/(1.0 + g);

	std::cout<<"TPT G = "<<G<<std::endl;
}

void onePoleTPT_HP::doFilterStage_HP(float xn, float& ynHP) {
	float v = (xn - z1)*G;
	float ynLP = v + z1;
	ynHP = xn - ynLP;
	z1 = ynLP + v;

}


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT Low-Pass + High-Pass Model ---* //

void onePoleTPT_LPHP::setFc_LPHP(float fc) {
	// prewarp the cutoff - bilinear-transform filters
	float wd = 2*pi*fc;
	float T = 1/sampleRate;
	float wa = (2/T)*tan(wd*T/2);
	float g = wa*T/2;

	// calculate big G value - Zavalishin p46 (the Art of VA Design)
	// 0.0 < G < 1.0 (from simulation)
	G = g/(1.0 + g);

	std::cout<<"TPT G = "<<G<<std::endl;
}

void onePoleTPT_LPHP::doFilterStage_LPHP(float xn, float& ynLP, float& ynHP) {
	float v = (xn - z1)*G;
	ynLP = v + z1;
	ynHP = xn - ynLP;
	z1 = ynLP + v;

}


// *---------------------------------------------------------------------------* //
// *--- 1-pole TPT All-Pass Model ---* //

void onePoleTPT_AP::setFc_AP(float fc) {
	// prewarp the cutoff - bilinear-transform filters
	float wd = 2*pi*fc;
	float T = 1/sampleRate;
	float wa = (2/T)*tan(wd*T/2);
	float g = wa*T/2;

	// calculate big G value - Zavalishin p46 (the Art of VA Design)
	// 0.0 < G < 1.0 (from simulation)
	G = g/(1.0 + g);

	std::cout<<"TPT G = "<<G<<std::endl;
}

void onePoleTPT_AP::doFilterStage_AP(float xn, float& ynAP) {
	float v = (xn - z1)*G;
	float LP = v + z1;
	float HP = xn - LP;
	ynAP = LP - HP;
	z1 = LP + v;

}

// *---------------------------------------------------------------------------* //
