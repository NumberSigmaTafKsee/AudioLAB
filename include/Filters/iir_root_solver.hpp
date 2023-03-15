#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/Polynomials>


namespace iir_filters
{
	struct RootSolver
	{
		// will solve any polynomial transfer function into a filter
		// H(s) = Y(s) / X(s)
		// zeros = Y(s)
		// poles = X(s)		
		
		std::vector<DspFloatType> z,p,xd,yd;
		
		RootSolver(size_t n) {
			z.resize(n);
			p.resize(n);
			xd.resize(n);
			yd.resize(n);
			clear();
		}
		void setZero(size_t i, DspFloatType v) {
			z[i] = v;
		}
		void setPole(size_t i, DspFloatType v) {
			p[i] = v;
		}
		void setZeros(const std::vector<DspFloatType> & v) {
			z = v;
		}
		void setPoles(const std::vector<DspFloatType> & v) {
			p = v;
		}
		void clear() {
			std::fill(xd.begin(),xd.end(),0.0);
			std::fill(yd.begin(),yd.end(),0.0);			
		}
		DspFloatType Tick(DspFloatType I)
		{
			T sum = 0;
			sum = xd[0] * z[0];			
			for(size_t i = 1; i < z.size(); i++)
				sum += xd[i] * z[i];
			for(size_t i = 1; i < p.size(); i++)
				sum -= yd[i] * p[i];
			for(size_t i = 1; i < xd.size(); i++)
				xd[i] = xd[i-1];
			for(size_t i = 1; i < yd.size(); i++)
				yd[i] = yd[i-1];
			xd[0] = I;
			xd[1] = sum;
		}
		std::vector<std::complex<DspFloatType>> solveRoots(const std::vector<double> & coeffs)
		{
			for(size_t i = 0; i < bessel.size(); i++) coeff[i] = coeffs[i];
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
			
	};
}
