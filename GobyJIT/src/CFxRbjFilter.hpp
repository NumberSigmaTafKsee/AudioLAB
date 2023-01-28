#pragma once

#include <cmath>
#include "SoundObject.hpp"

namespace Filters::RBJ
{
	class RBJFilter : public FilterProcessor
	{
	public:
		
		RBJFilter() : FilterProcessor()
		{
			// reset filter coeffs
			b0a0=b1a0=b2a0=a1a0=a2a0=0.0;
			
			// reset in/out history
			ou1=ou2=in1=in2=0.0f;	
		};

		double filter(double in0)
		{
			// filter
			double const yn = b0a0*in0 + b1a0*in1 + b2a0*in2 - a1a0*ou1 - a2a0*ou2;

			// push in/out buffers
			in2=in1;
			in1=in0;
			ou2=ou1;
			ou1=yn;

			// return output
			return yn;
		};

		double Tick(double I, double A = 1, double F = 0, double Q = 0) {
			double cut = freq;
			double r   = this->q;
			calc_filter_coeffs(type,freq + F*freq, sr, this->q+Q, dbGain,bandwidth);
			double x   = std::tanh(A*filter(A*I));
			calc_filter_coeffs(type,cut, sr, r, dbGain,bandwidth);
			return x;
		}
		void calc_filter_coeffs(int const type,double const frequency,double const sample_rate,double const q,double const db_gain,bool q_is_bandwidth)
		{
			// temp pi
			double const temp_pi=3.1415926535897932384626433832795;
			
			if(frequency < 0) frequency = 0;
			if(frequency > (sample_rate/2-1)) frequency = sample_rate/2-1;
			if(q < 0) q = 0;
			if(db_gain < 0) db_gain = 0;

			// temp coef vars
			double alpha,a0,a1,a2,b0,b1,b2;
			freq = frequency;
			sr   = sample_rate;
			Q    = q;
			dbGain = db_gain;
			bandwidth = q_is_bandwidth
			// peaking, lowshelf and hishelf
			if(type>=6)
			{
				double const A		=	pow(10.0,(db_gain/40.0));
				double const omega	=	2.0*temp_pi*frequency/sample_rate;
				double const tsin	=	sin(omega);
				double const tcos	=	cos(omega);
				
				if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*q*omega/tsin);
				else
				alpha=tsin/(2.0*q);

				double const beta	=	sqrt(A)/q;
				
				// peaking
				if(type==6)
				{
					b0=double(1.0+alpha*A);
					b1=double(-2.0*tcos);
					b2=double(1.0-alpha*A);
					a0=double(1.0+alpha/A);
					a1=double(-2.0*tcos);
					a2=double(1.0-alpha/A);
				}
				
				// lowshelf
				if(type==7)
				{
					b0=double(A*((A+1.0)-(A-1.0)*tcos+beta*tsin));
					b1=double(2.0*A*((A-1.0)-(A+1.0)*tcos));
					b2=double(A*((A+1.0)-(A-1.0)*tcos-beta*tsin));
					a0=double((A+1.0)+(A-1.0)*tcos+beta*tsin);
					a1=double(-2.0*((A-1.0)+(A+1.0)*tcos));
					a2=double((A+1.0)+(A-1.0)*tcos-beta*tsin);
				}

				// hishelf
				if(type==8)
				{
					b0=double(A*((A+1.0)+(A-1.0)*tcos+beta*tsin));
					b1=double(-2.0*A*((A-1.0)+(A+1.0)*tcos));
					b2=double(A*((A+1.0)+(A-1.0)*tcos-beta*tsin));
					a0=double((A+1.0)-(A-1.0)*tcos+beta*tsin);
					a1=double(2.0*((A-1.0)-(A+1.0)*tcos));
					a2=double((A+1.0)-(A-1.0)*tcos-beta*tsin);
				}
			}
			else
			{
				// other filters
				double const omega	=	2.0*temp_pi*frequency/sample_rate;
				double const tsin	=	sin(omega);
				double const tcos	=	cos(omega);

				if(q_is_bandwidth)
				alpha=tsin*sinh(log(2.0)/2.0*q*omega/tsin);
				else
				alpha=tsin/(2.0*q);

				
				// lowpass
				if(type==0)
				{
					b0=(1.0-tcos)/2.0;
					b1=1.0-tcos;
					b2=(1.0-tcos)/2.0;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// hipass
				if(type==1)
				{
					b0=(1.0+tcos)/2.0;
					b1=-(1.0+tcos);
					b2=(1.0+tcos)/2.0;
					a0=1.0+ alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// bandpass csg
				if(type==2)
				{
					b0=tsin/2.0;
					b1=0.0;
					b2=-tsin/2;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// bandpass czpg
				if(type==3)
				{
					b0=alpha;
					b1=0.0;
					b2=-alpha;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// notch
				if(type==4)
				{
					b0=1.0;
					b1=-2.0*tcos;
					b2=1.0;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// allpass
				if(type==5)
				{
					b0=1.0-alpha;
					b1=-2.0*tcos;
					b2=1.0+alpha;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}
			}

			// set filter coeffs
			b0a0=double(b0/a0);
			b1a0=double(b1/a0);
			b2a0=double(b2/a0);
			a1a0=double(a1/a0);
			a2a0=double(a2/a0);
		};

		void setCutoff(double f) {
			freq = f;
		}
		void setQ(double Q) {
			q = Q;
		}
		void setGain(double db) {
			dbGain = db;
		}
	private:

		double freq,sr,q,dbGain;
		bool bandwidth;
		// filter coeffs
		double b0a0,b1a0,b2a0,a1a0,a2a0;

		// in/out history
		double ou1,ou2,in1,in2;
	};

	struct StereoRBJFilter : public StereoFXProcessor
	{
		RBJFilter left,right;

		StereoRBJFilter() :
		StereoFXProcessr()
		{

		}
		void setCutoff(double c) {
			left.setCutoff(c);
			right.setCutoff(c);
		}
		void setQ(double q) {
			left.setQ(q);
			right.setQ(q);
		}
		void setGain(double gain) {
			left.setGain(gain);
			right.setGain(gain);
		}
		double Tick(double iL, double iR, double &L, double &R, double A=1, double X=0, double Y=0)
		{
			L = left.Tick(iL,A,X,Y);
			R = right.Tick(iR,A,X,Y);
			return 0.5*(L+R);
		}
	};
}