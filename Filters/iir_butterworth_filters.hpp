#pragma once

namespace iir_filters
{
	std::complex<DspFloatType> butterworthpole(DspFloatType k, DspFloatType n)
	{
		DspFloatType p = M_PI * ((2 * k + n - 1) / (2 * n));
		return std::complex<DspFloatType>(-std::cos(p), -std::sin(p));
	}
	std::complex<DspFloatType> ButterworthPoles(DspFloatType K, DspFloatType N)
	{
		DspFloatType theta = ((2*K+N-1)/(2*N))*M_PI;
		DspFloatType sk = cos(theta);
		DspFloatType omegak = -sin(theta);
		std::complex<DspFloatType> p(sk,omegak);
		return p;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// Butterworth
	/////////////////////////////////////////////////////////////////////////////////////////////
		BiquadSOS butterlp(int order, double Q=1.0)
		{
			BiquadSOS sos;    
			size_t n = 1;
			if(order %2 != 0) {
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
				DspFloatType x1 = abs(p1);
				DspFloatType x2 = 0;
						
				// (s-p1)        
				c.z[0] = 1.0/x1;
				c.z[1] = 0.0;
				c.z[2] = 0.0;
				c.p[0] = 1.0;      
				c.p[1] = 1/x1;
				c.p[2] = 0.0;
				sos.push_back(c);        
				n++;
			}
				
			for(size_t i = n; i < order; i += 2)
			{
				std::complex<DspFloatType> p1  = ButterworthPoles(i,order);
				std::complex<DspFloatType> p2  = ButterworthPoles(i+1,order);
				
				DspFloatType x1 = abs(p1*p2);
				DspFloatType x2 = abs(-p1-p2);
				
				// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
				BiquadSection c;
				c.z[0] = 1.0/x1;
				c.z[1] = 0.0;
				c.z[2] = 0.0;
				c.p[0] = 1;        
				c.p[1] = ((1.0/Q)*x2)/x1;
				c.p[2] = 1/x1;    
				sos.push_back(c);
			}
			return sos;
		}
		BiquadSOS butterhp(int order, double Q=1.0)
		{
			BiquadSOS sos;    
			size_t n = 1;
			
			if(order %2 != 0) {
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
				DspFloatType x1 = abs(p1);
				

				c.z[0] = 0.0;
				c.z[1] = 1.0/x1;
				c.z[2] = 0.0;
				c.p[0] = 1.0;
				c.p[1] = 1/x1;
				c.p[2] = 0.0;

				sos.push_back(c);        
				n++;
			}
				
			for(size_t i = n; i < order; i += 2)
			{
				std::complex<DspFloatType> p1  = ButterworthPoles(i,order);
				std::complex<DspFloatType> p2  = ButterworthPoles(i+1,order);
				
				DspFloatType x1 = abs(p1*p2);
				DspFloatType x2 = abs(-p1-p2);
				
				// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
				BiquadSection c;
				c.z[0] = 0.0;
				c.z[1] = 0.0;
				c.z[2] = 1.0/x1;
				c.p[0] = 1;        
				c.p[1] = (1.0/Q)*x2/x1;
				c.p[2] = 1.0/x1;    
				sos.push_back(c);
			}
			
			return sos;
		}

		BiquadSOS butterlp2hp(int order, double Q=1.0)
		{
			BiquadSOS sos;    
			size_t n = 1;
			
			if(order %2 != 0) {
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
				DspFloatType x1 = abs(p1);
				

				c.z[0] = 0.0;
				c.z[1] = 1.0/x1;
				c.z[2] = 0.0;
				c.p[0] = 1.0;
				c.p[1] = 1/x1;
				c.p[2] = 0.0;

				sos.push_back(c);        
				n++;
			}
				
			for(size_t i = n; i < order; i += 2)
			{
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
				std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
				
				std::vector<std::complex<DspFloatType>> zeros,poles;
				poles.push_back(p1);
				poles.push_back(p2);
				DspFloatType gain;
				lp2hp(zeros,poles,1.0,gain);
				// (1-z)(1-z) = 1 -z1-z2 +z1z2
				
				DspFloatType x1 = abs(poles[0]*poles[1]);
				c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
				c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
				c.z[2] = 1.0/x1;
				c.p[0] = 1.0;
				c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
				c.p[2] = 1.0/x1;

				sos.push_back(c);
			}
			
			return sos;
		}
		BiquadSOS butterlp2bp(int order, double Q=1.0)
		{
			BiquadSOS sos;    
			size_t n = 1;
						
			for(size_t i = n; i < order; i += 2)
			{
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
				std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
				std::vector<std::complex<DspFloatType>> zeros,poles;
				poles.push_back(p1);
				poles.push_back(p2);
				DspFloatType gain;
				lp2bp(zeros,poles,1.0,0.5,gain);
				// (1-z)(1-z) = 1 -z1-z2 +z1z2
				
				DspFloatType x1 = abs(poles[0]*poles[1]);
				c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
				c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
				c.z[2] = 1.0/x1;
				c.p[0] = 1.0;
				c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
				c.p[2] = 1.0/x1;

				sos.push_back(c);
			}
			
			return sos;
		}
		BiquadSOS butterlp2bs(int order, double Q=1.0)
		{
			BiquadSOS sos;    
			size_t n = 1;
						
			for(size_t i = n; i < order; i += 2)
			{
				BiquadSection c;
				std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
				std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
				std::vector<std::complex<DspFloatType>> zeros,poles;
				poles.push_back(p1);
				poles.push_back(p2);
				DspFloatType gain;
				lp2bs(zeros,poles,1.0,0.5,gain);
				// (1-z)(1-z) = 1 -z1-z2 +z1z2
				
				DspFloatType x1 = abs(poles[0]*poles[1]);
				c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
				c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
				c.z[2] = 1.0/x1;
				c.p[0] = 1.0;
				c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
				c.p[2] = 1.0/x1;

				sos.push_back(c);
			}
			
			return sos;
		}
		
		BiquadSection butter2lp(double Q=1.0)
		{    
			std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
			std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
		
			DspFloatType x1 = abs(p1*p2);
			DspFloatType x2 = abs(-p1-p2);
			//std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
			// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
			BiquadSection c;
			c.z[0] = 1.0/x1;
			c.z[1] = 0.0;
			c.z[2] = 0.0;
			c.p[0] = 1;    
			c.p[1] = (1.0/Q)*x2/x1;
			c.p[2] = 1.0/x1;        

			return c;
		}
		BiquadSection butter2hp(double Q=1.0)
		{    
			std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
			std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
		
			DspFloatType x1 = abs(p1*p2);
			DspFloatType x2 = abs(-p1-p2);
			//std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
			// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
			BiquadSection c;
			c.z[0] = 0.0;
			c.z[1] = 0.0;
			c.z[2] = 1.0/x1;
			c.p[0] = 1;    
			c.p[1] = (1.0/Q)*x2/x1;
			c.p[2] = 1.0/x1;        

			return c;
		}
		BiquadSection butterlp2hp2(double Q=1.0)
		{    
			std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
			std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
		
			std::vector<std::complex<DspFloatType>> zeros,poles;
			poles.push_back(p1);
			poles.push_back(p2);
			DspFloatType gain;
			lp2hp(zeros,poles,1.0,gain);
			// (1-z)(1-z) = 1 -z1-z2 +z1z2
			BiquadSection c;        
			DspFloatType x1 = abs(poles[0]*poles[1]);
			c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
			c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
			c.z[2] = 1.0/x1;
			c.p[0] = 1.0;
			c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
			c.p[2] = 1.0/x1;

			return c;
		}
		BiquadSection butterlp2bp2(double Q=1.0)
		{    
			std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
			std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
		
			std::vector<std::complex<DspFloatType>> zeros,poles;
			poles.push_back(p1);
			poles.push_back(p2);
			DspFloatType gain;
			// i dont not really know what this should be normalized 1.0,0 or 1.0,0.5?
			lp2bp(zeros,poles,1.0,0.5,gain);
			// (1-z)(1-z) = 1 -z1-z2 +z1z2
			BiquadSection c;        
			DspFloatType x1 = abs(poles[0]*poles[1]);
			c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
			c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
			c.z[2] = 1.0/x1;
			c.p[0] = 1.0;
			c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
			c.p[2] = 1.0/x1;
			return c;
		}
		BiquadSection butterlp2bs2(double Q=1.0)
		{    
			std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
			std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
		
			std::vector<std::complex<DspFloatType>> zeros,poles;
			poles.push_back(p1);
			poles.push_back(p2);
			DspFloatType gain;
			lp2bs(zeros,poles,1.0,0.5,gain);
			// (1-z)(1-z) = 1 -z1-z2 +z1z2
			BiquadSection c;        
			DspFloatType x1 = abs(poles[0]*poles[1]);
			c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
			c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
			c.z[2] = 1.0/x1;
			c.p[0] = 1.0;
			c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
			c.p[2] = 1.0/x1;
			return c;
		}
};				

