// https://github.com/antoniograzioli/Autodafe/blob/master/src/Biquad.cpp

//
//  Biquad.h
//
//  Created by Nigel Redmon on 11/24/12
//  EarLevel Engineering: earlevel.com
//  Copyright 2012 Nigel Redmon
//
//  For a complete explanation of the Biquad code:
//  http://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
//
//  License:
//
//  This source code is provided as is, without warranty.
//  You may copy and distribute verbatim copies of this document.
//  You may modify and use this source code to create binary code
//  for your own purposes, free or commercial.
//

#pragma once

namespace Filters::AutoDafe
{
	enum {
		bq_type_lowpass = 0,
		bq_type_highpass,
		bq_type_bandpass,
		bq_type_notch,
		bq_type_peak,
		bq_type_lowshelf,
		bq_type_highshelf
	};

	class Biquad : FilterProcessor
	{
	public:
		Biquad();
		Biquad(int type, double Fc, double Q, double peakGainDB);
		~Biquad();
		void setType(int type);
		void setQ(double Q);
		void setFc(double Fc);
		void setPeakGain(double peakGainDB);
		void setBiquad(int type, double Fc, double Q, double peakGain);
		float process(float in);

	protected:
		void calcBiquad(void);

		int type;
		double a0, a1, a2, b1, b2;
		double Fc, Q, peakGain;
		double z1, z2;
	};

	inline float Biquad::process(float in) {
		double out = in * a0 + z1;
		z1 = in * a1 + z2 - b1 * out;
		z2 = in * a2 - b2 * out;
		return out;
	}

	//
	//  Biquad.cpp
	//
	//  Created by Nigel Redmon on 11/24/12
	//  EarLevel Engineering: earlevel.com
	//  Copyright 2012 Nigel Redmon
	//
	//  For a complete explanation of the Biquad code:
	//  http://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
	//
	//  License:
	//
	//  This source code is provided as is, without warranty.
	//  You may copy and distribute verbatim copies of this document.
	//  You may modify and use this source code to create binary code
	//  for your own purposes, free or commercial.
	//

	#include <math.h>
	#include "Biquad.h"

	Biquad::Biquad() : FilterProcessor()
	{
		type = bq_type_lowpass;
		a0 = 1.0;
		a1 = a2 = b1 = b2 = 0.0;
		Fc = 0.50;
		Q = 0.707;
		peakGain = 0.0;
		z1 = z2 = 0.0;
	}

	Biquad::Biquad(int type, double Fc, double Q, double peakGainDB) 
	: FilterProcessor()
	{
		setBiquad(type, Fc, Q, peakGainDB);
		z1 = z2 = 0.0;
	}

	Biquad::~Biquad() {
	}

	void Biquad::setType(int type) {
		this->type = type;
		calcBiquad();
	}

	void Biquad::setQ(double Q) {
		this->Q = Q;
		calcBiquad();
	}

	void Biquad::setFc(double Fc) {
		this->Fc = Fc;
		calcBiquad();
	}

	void Biquad::setPeakGain(double peakGainDB) {
		this->peakGain = peakGainDB;
		calcBiquad();
	}

	void Biquad::setBiquad(int type, double Fc, double Q, double peakGainDB) {
		this->type = type;
		this->Q = Q;
		this->Fc = Fc;
		setPeakGain(peakGainDB);
	}

	void Biquad::calcBiquad(void) {
		double norm;
		double V = pow(10, fabs(peakGain) / 20.0);
		double K = tan(M_PI * Fc);
		switch (this->type) {
		case bq_type_lowpass:
			norm = 1 / (1 + K / Q + K * K);
			a0 = K * K * norm;
			a1 = 2 * a0;
			a2 = a0;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - K / Q + K * K) * norm;
			break;

		case bq_type_highpass:
			norm = 1 / (1 + K / Q + K * K);
			a0 = 1 * norm;
			a1 = -2 * a0;
			a2 = a0;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - K / Q + K * K) * norm;
			break;

		case bq_type_bandpass:
			norm = 1 / (1 + K / Q + K * K);
			a0 = K / Q * norm;
			a1 = 0;
			a2 = -a0;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - K / Q + K * K) * norm;
			break;

		case bq_type_notch:
			norm = 1 / (1 + K / Q + K * K);
			a0 = (1 + K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = a0;
			b1 = a1;
			b2 = (1 - K / Q + K * K) * norm;
			break;

		case bq_type_peak:
			if (peakGain >= 0) {    // boost
				norm = 1 / (1 + 1 / Q * K + K * K);
				a0 = (1 + V / Q * K + K * K) * norm;
				a1 = 2 * (K * K - 1) * norm;
				a2 = (1 - V / Q * K + K * K) * norm;
				b1 = a1;
				b2 = (1 - 1 / Q * K + K * K) * norm;
			}
			else {    // cut
				norm = 1 / (1 + V / Q * K + K * K);
				a0 = (1 + 1 / Q * K + K * K) * norm;
				a1 = 2 * (K * K - 1) * norm;
				a2 = (1 - 1 / Q * K + K * K) * norm;
				b1 = a1;
				b2 = (1 - V / Q * K + K * K) * norm;
			}
			break;
		case bq_type_lowshelf:
			if (peakGain >= 0) {    // boost
				norm = 1 / (1 + sqrt(2) * K + K * K);
				a0 = (1 + sqrt(2 * V) * K + V * K * K) * norm;
				a1 = 2 * (V * K * K - 1) * norm;
				a2 = (1 - sqrt(2 * V) * K + V * K * K) * norm;
				b1 = 2 * (K * K - 1) * norm;
				b2 = (1 - sqrt(2) * K + K * K) * norm;
			}
			else {    // cut
				norm = 1 / (1 + sqrt(2 * V) * K + V * K * K);
				a0 = (1 + sqrt(2) * K + K * K) * norm;
				a1 = 2 * (K * K - 1) * norm;
				a2 = (1 - sqrt(2) * K + K * K) * norm;
				b1 = 2 * (V * K * K - 1) * norm;
				b2 = (1 - sqrt(2 * V) * K + V * K * K) * norm;
			}
			break;
		case bq_type_highshelf:
			if (peakGain >= 0) {    // boost
				norm = 1 / (1 + sqrt(2) * K + K * K);
				a0 = (V + sqrt(2 * V) * K + K * K) * norm;
				a1 = 2 * (K * K - V) * norm;
				a2 = (V - sqrt(2 * V) * K + K * K) * norm;
				b1 = 2 * (K * K - 1) * norm;
				b2 = (1 - sqrt(2) * K + K * K) * norm;
			}
			else {    // cut
				norm = 1 / (V + sqrt(2 * V) * K + K * K);
				a0 = (1 + sqrt(2) * K + K * K) * norm;
				a1 = 2 * (K * K - 1) * norm;
				a2 = (1 - sqrt(2) * K + K * K) * norm;
				b1 = 2 * (K * K - V) * norm;
				b2 = (V - sqrt(2 * V) * K + K * K) * norm;
			}
			break;
		}

		return;
	}


	struct FixedFilter : public FilterProcessor
	{
		
		FixedFilter() {
			eq[0] = 1.0;
			eq[1] = 1.0;
			eq[2] = 1.0;
			eq[3] = 1.0;
			eq[4] = 1.0;
			eq[5] = 1.0;
			eq[6] = 1.0;
			eq[7] = 1.0;

			bq1 = new Biquad();
			bq2 = new Biquad();
			bq3 = new Biquad();
			bq4 = new Biquad();
			bq5 = new Biquad();
			bq6 = new Biquad();
			bq7 = new Biquad();
			bq8 = new Biquad();
		}
		~FixedFilter() 
		{
			delete bq1;
			delete bq2;
			delete bq3;
			delete bq4;
			delete bq5;
			delete bq6;
			delete bq7;
			delete bq8;
		}

		Biquad *bq1;
		Biquad *bq2;
		Biquad *bq3;
		Biquad *bq4;
		Biquad *bq5;
		Biquad *bq6;
		Biquad *bq7;
		Biquad *bq8;
		
		
		void step();
	};

	FixedFilter::SeqEQ(double EQ[8]) {
		memcpy(eq,EQ,8*sizeof(double));		
	};

	void FixedFilter::step(double input) 
	{	
		bq1->setBiquad(bq_type_peak, 75.0 / engineGetSampleRate(), 5,eq[0]);
		bq2->setBiquad(bq_type_peak, 125.0 / engineGetSampleRate(), 5, eq[1]);
		bq3->setBiquad(bq_type_peak, 250.0 / engineGetSampleRate(), 5, eq[2]);
		bq4->setBiquad(bq_type_peak, 500.0 / engineGetSampleRate(), 5, eq[3]);
		bq5->setBiquad(bq_type_peak, 1000.0 / engineGetSampleRate(), 5, eq[4]);
		bq6->setBiquad(bq_type_peak, 2000.0 / engineGetSampleRate(), 5, eq[5]);
		bq7->setBiquad(bq_type_peak, 4000.0 / engineGetSampleRate(), 5, eq[6]);
		bq8->setBiquad(bq_type_peak, 8000.0 / engineGetSampleRate(), 5, eq[7]);
			
		double out = bq1->process(input/5.0);		
		out = bq2->process(out);
		out = bq3->process(out);
		out = bq4->process(out);
		out = bq5->process(out);
		out = bq6->process(out);
		out = bq7->process(out);
		out = bq8->process(out);		
		return out*5.0;
	}
}