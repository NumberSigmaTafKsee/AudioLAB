#pragma once
#include <Eigen/Core>
#include <unsupported/Eigen/Polynomials>

namespace iir_filters
{
	// Hn(s) = H0/d * PROD( (s^2 + ai) / (s^2 + bi*s + ci))
	// d = s+p0 n=odd
	// d = 1    n=even

	DspFloatType factorial(DspFloatType x)
	{
		if (x == 0)
			return 1;
		return x * factorial(x - 1);
	}
	DspFloatType binomial(DspFloatType n, DspFloatType k)
	{
		if(n == 0 || k == 0) return 1;
		return factorial(n) / (factorial(k) * factorial(fabs(n - k)));
	}

	DspFloatType bk(DspFloatType n, DspFloatType k) {
		DspFloatType num = factorial(2*n-k);
		DspFloatType den = pow(2.0,n-k)*factorial(k)*factorial(n-k);
		return num/den;
	}

	// calculate bessel filter coefficients qn(s) order n 
	// it will return the polynomial coefficients starting at b0
	// b0,b1*s,b2*s^2 + ... bn*s^n
	// these use the polynomial root solver to find the poles (complex)
	std::vector<DspFloatType> qn(size_t n) 
	{
		std::vector<DspFloatType> r(n+1);
		r[0] = bk(n,0);
		for(size_t k=1; k <= n; k++) {
			DspFloatType b = bk(n,k);
			r[k] = b;
		}
		return r;
	}

	std::vector<std::complex<DspFloatType>> bessel_roots(uint32_t order)
	{
		std::vector<DspFloatType> bessel = qn(order);	
		Eigen::VectorXd coeff(8);
		for(size_t i = 0; i < bessel.size(); i++) coeff[i] = bessel[i];
		Eigen::PolynomialSolver<DspFloatType, Eigen::Dynamic> solver;
		solver.compute(coeff);
		const Eigen::PolynomialSolver<DspFloatType, Eigen::Dynamic>::RootsType &r = solver.roots();
		std::vector<std::complex<DspFloatType>> out;
		
		// real must go first
		for(size_t i = 0; i < r.size(); i++) {
			if(r[i].imag() == 0) {
				out.push_back(r[i]);
				break;
			}
		}
		for(size_t i = 0; i < r.size(); i++) {
			if(r[i].imag() != 0) {
				out.push_back(r[i]);
			}
		}
		return out;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// Bessel Filter
	/////////////////////////////////////////////////////////////////////////////////////////////    
	BiquadSOS bessellp(int order, double Q=1.0)
	{
		
		BiquadSOS sos;    
		std::vector<std::complex<DspFloatType>> bessel_poles = bessel_roots(order);		
		size_t n = 1;
		

		if(order %2 != 0) {
			BiquadSection c;
			std::complex<DspFloatType> p1 = bessel_poles[0];
			DspFloatType x1 = abs(p1);
			DspFloatType x2 = 0;
					
			// (s-p1)        
			c.z[0] = 1.0/x1;
			c.z[1] = 0.0;
			c.z[2] = 0.0;
			c.p[0] = 1.0;;        
			c.p[1] = 1/x1;
			c.p[2] = 0.0;
			sos.push_back(c);        
			n++;
		}
			
		for(size_t i = n; i < order; i += 2)
		{
			std::complex<DspFloatType> p1 = bessel_poles[i];
			std::complex<DspFloatType> p2 = bessel_poles[i+1];
			
			DspFloatType x1 = abs(p1*p2);
			DspFloatType x2 = abs(-p1-p2);
			
			// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
			BiquadSection c;
			c.z[0] = 1.0/x1;
			c.z[1] = 0.0;
			c.z[2] = 0.0;
			c.p[0] = 1;        
			c.p[1] = (1.0/Q)*x2/x1;
			c.p[2] = 1/x1;    
			sos.push_back(c);
		}
		return sos;
	}
	
	BiquadSOS besselhp(int order, double Q=1.0)
	{		
		BiquadSOS sos;    
		std::vector<std::complex<DspFloatType>> bessel_poles = bessel_roots(order);
		size_t n = 1;
		if(order % 2 != 0) {
			BiquadSection c;
			std::complex<DspFloatType> p1  = bessel_poles[0];
			DspFloatType x1 = abs(p1);
			DspFloatType x2 = 0;
					
			// (s-p1)        
			c.z[0] = 0.0;
			c.z[1] = 1.0;
			c.z[2] = 0.0;
			c.p[0] = 1.0;;        
			c.p[1] = 1/x1;
			c.p[2] = 0.0;
			sos.push_back(c);        
			n++;
		}
			
		for(size_t i = n; i < order; i += 2)
		{
			std::complex<DspFloatType> p1  = bessel_poles[i];
			std::complex<DspFloatType> p2  = bessel_poles[i+1];
			
			DspFloatType x1 = abs(p1*p2);
			DspFloatType x2 = abs(-p1-p2);
			
			// (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
			BiquadSection c;
			
			c.z[0] = 0.0;
			c.z[1] = 0.0;
			c.z[2] = 1.0/x1;
			c.p[0] = 1;        
			c.p[1] = (1.0/Q)*x2/x1;
			c.p[2] = 1/x1;  

			sos.push_back(c);
		}
		return sos;
	}
}
