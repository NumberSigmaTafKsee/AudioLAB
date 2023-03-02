#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "PolynomialRoots.hpp"
#include "LaguerrePolynomialRoots.hpp"
#include <Eigen/Core>
#include <unsupported/Eigen/Polynomials>

typedef double DspFloatType;

DspFloatType factorial(DspFloatType x)
{
    if (x == 0) return 1;
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

// doesn't work
int TestLaguerre()
{
  std::vector<DspFloatType> bessel = qn(7);

  deque<complex<DspFloatType>> P,R;
  for(size_t i = 0; i < bessel.size(); i++) P.push_back(std::complex<DspFloatType>(bessel[i],0));
  LaguerreMethod L(P);
  R = L.solve_roots();

  // Display the equation to solve
  for (int i = P.size()-1; i >=1 ; i--) cout << P[i] <<"*x^" << i << " + ";
  cout << P[0] << "= 0"<< endl;

  //Display roots
  cout << "--------- ROOTS ---------" << endl;
  for (int i = 0; i < R.size(); i++) cout << R[i] << endl;
  
  //In order to gather the quality of the roots found, the lines below evaluate the polynom at every
  //root, R_i, and takes the maximum and the mean deviation to zero 
  cout << endl << endl;
  complex<DspFloatType> P_x(0.,0.);
  DspFloatType mean_err = 0, max_err = -1e12, abs_err;
  cout << "--------- Root analysis ---------" << endl;
  for (int j = 0; j < R.size(); j++)
  { 
    P_x = complex<DspFloatType>(0,0);
    for (int i = 0; i < P.size(); i++) P_x += P[i]*pow(R[j], i);
    abs_err = abs(P_x);
    if (abs(P_x) > max_err) max_err = abs_err;
    mean_err += abs_err;
    cout << "Error{Root[" << j << "]}= " << abs_err << endl;
  }
  cout << endl;
  cout << "Mean error = sum(|P(R_i)|, i={1, N}/N = " << mean_err/R.size() << endl;
  cout << "Maximum error = max(|P(R_i)|), i={1, N}/N = " << max_err << endl;
  return 0;
}

// works
void TestEigen()
{
	std::vector<DspFloatType> bessel = qn(7);	
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
	for(size_t i = 0; i < out.size(); i++) {
		std::cout << out[i] << std::endl;
	}
}

int main()
{
	TestEigen();
}
