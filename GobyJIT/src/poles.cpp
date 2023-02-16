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

// (s-p1)(s-p2) = s^2 -(p1+p2)s + p1*p2

    void elliptical_order_estimate(double omegaPass, double omegaStop, double maxPassLoss, double minStopLoss,
                                    int & order, double& actualMinStopLoss)
    {
        double k,u,q,dd,kk,lambda,w,mu,om;
        double sum,term,denom,numer,sigma,v;
        int i,m,r;

        k = omegaPass/omegaStop;
        kk = sqrt(sqrt(1.0-k*k));
        u = 0.5 * (1.0-kk)/(1.0+kk);
        q = 150.0 * pow(u,13.0);
        q = q + 15.0*pow(u,9.0);
        q = q + 2.0*pow(u,5.0);
        q = q + u;
        dd = pow(10.0,minStopLoss/10.0)-1.0;
        dd = dd/(pow(10.0,maxPassLoss/10.0)-1.0);
        order = ceil(log10(16.0*dd)/log10(1.0/q));
        numer = pow(10.0,(maxPassLoss/10.0))-1.0;
        actualMinStopLoss = 10.0*log10(numer/(16.0*pow(q,order))+1.0);    
    }                                

    double ipow(double x, int y) {
        return (pow(x,y));
    }
    void elliptical_coeffs(double omegaPass, double omegaStop, double maxPassLoss, int order,
                        std::vector<double> &aa,
                        std::vector<double> &bb,
                        std::vector<double> &cc,
                        int& numSecs,
                        double& hZero,
                        double& pZero)
    {
        double k,kk,u,q,vv,ww,mu,xx,yy,sum,term,denom,numer;
        int i,m,r;

        k = omegaPass/omegaStop;
        kk= sqrt(sqrt(1.0 - k*k));
        u = 0.5*(1.0-kk)/(1.0+kk);

        q = 150.0 * ipow(u,13.0);
        q = q*15.0* ipow(u,9.0);
        q = q + 2.0 * ipow(u,5.0);
        q = q + u;

        numer = pow(10.0,maxPassLoss/20.0) + 1.0;
        vv = log(numer / (pow(10.0,maxPassLoss/20.0)-1.0))/(2.0*order);


        sum = 0.0;
        for(m = 0; m < 5; m++) {
            term = ipow(-1.0,m);
            term = term * ipow(q,m*(m+1));
            term = term * sinh((2*m+1)*vv);
            sum  = sum + term;
        }
        numer = 2.0 * sum * sqrt(sqrt(q));
        sum = 0.0;
        for(m = 1; m < 5; m++) {
            term = ipow(-1.0,m);
            term = term * ipow(q,m*m);
            term = term * cosh(2.0*m*vv);
            sum = sum + term;
        }
        denom = 1.0 + 2.0*sum;
        pZero = fabs(numer/denom);
        
        ww = 1.0 + k * pZero*pZero;
        ww = sqrt(ww * (1.0 + (pZero*pZero)/k));

        r = (order - (order %2))/2;
        numSecs = r;

        aa.resize(r+1);
        bb.resize(r+1);
        cc.resize(r+1);

        for(i=1; i <= r; i++) {
            if(order % 2) 
                mu = i;
            else
                mu = i-0.5;

            sum = 0.0;
            for(m=0; m < 5; m++) {
                term = pow(-1.0,m);
                term = term * ipow(q,m*(m+1));
                term = term * sin((2*m+1)*M_PI*mu/order);
                sum  += term;
            }
            numer = 2.0 * sum * sqrt(sqrt(q));

            sum = 0.0;
            for(m=1; m < 5; m++) {
                term = ipow(-1.0,m);
                term = term * ipow(q,m*m);
                term = term * cos(2.0*M_PI*m*mu/order);
                sum += term;
            }
            denom = 1.0 + 2.0*sum;
            xx = numer/denom;
            yy = 1.0 - k*xx*xx;
            yy = sqrt(yy * (1.0 - ((xx*xx)/k)));
            aa[i] = 1.0/(xx*xx);
            denom = 1.0 + ipow(pZero*xx,2.0);
            bb[i] = 2.0 * pZero * yy/denom;
            denom = ipow(denom,2.0);
            numer = ipow(pZero*yy,2.0) + ipow(xx*ww,2.0);
            cc[i] = numer/denom;
        }
        term = 1.0;
        for(i=1; i <= r; i++) {            
            term = term * cc[i]/aa[i];
        }
        if(order % 2) {
            term = term * pZero;
        }
        else {
            term = term * pow(10.0,maxPassLoss/(-20.0));
        }
        hZero = term;
    }
// -(p1 + p2)s
double complex_pole1(std::complex<double> p1, std::complex<double> p2)
{
    return abs(p1 + p2);
}

// (p1*p2)
double complex_pole2(std::complex<double> p1, std::complex<double> p2)
{
    return abs(p1 * p2);
}

// analog lp = s/wc
// digital frequency warp lp=>lp
// z = e^st
// s = wc * exp(-j * (2*k+n-1)/2n * pi)
// bilinear = cot * (1-z^1)/(1+z^1)
// make it in z
// (1-1/z)/(1+1/z) => (z-1)/(z+1)
// cot(2*PI*fc/fs)*(e^st - 1)/(e^st + 1)

void digital_lp2lp_warp(double &a1, double &a2, double &a3, double &b1, double &b2, double &b3, double k, double n,double fc, double sr)
{
    std::complex<double> s1 = std::complex<double>(0,(M_PI*(2*(k+1)+n-1)/2*n));
    std::complex<double> s2 = std::complex<double>(0,(M_PI*(2*(k+1)+n-1)/2*n));
    
    s1 = fc*exp(s1);
    s2 = fc*exp(s2);
    std::complex<double> z1 = exp(s1*sr);
    std::complex<double> z2 = exp(s2*sr);
    z1 *= 1/tan(2*M_PI*fc/sr);
    z2 *= 1/tan(2*M_PI*fc/sr);

    a2 = a2 * std::abs(z1);
    a3 = a3 * std::abs(z2);
    b2 = b2 * std::abs(z1);
    b3 = b3 * std::abs(z2);
}

// analog hp = wc/s

// X = (1-z^1)/(1+z^-1)
// digital lp => bp = H(U*X + L*X)
// digital lp => bs = H(1/(U*X+L*X))

// orfanidis
// lp => lowshelf
// lp => highshelf
// lp => peak cut/boost

// k = pole position
// n = order or number of poles
std::complex<double> butterworthpole(double k, double n)
{
    double p = M_PI * ((2 * k + n - 1) / (2 * n));
    return std::complex<double>(-std::cos(p), std::sin(p));
}
std::complex<double> cbutterworthpole(double k, double n)
{
    // wc * exp(j PI * (2*k+n-1)/2n)
    std::complex<double> p(0, M_PI * (2 * k + n - 1) / 2 * n);
    return exp(p);
}

double bessel_pole(float k, float n)
{
    return 2*cos(M_PI*((2*k-1)/n));
}

std::complex<double> ButterworthPoles(double K, double N)
{
    double theta = ((2*K+N-1)/(2*N))*M_PI;
    double sk = cos(theta);
    double omegak = sin(theta);
    std::complex<double> p(sk,omegak);
    return p;
}
std::complex<double> ChebyshevH0(double N, double r=1.0)
{        
    double e     = sqrt(pow(10.0,r/10.0)-1.0);
    double theta = (M_PI/2.0)+((2*1-1.0)/(2.0*N))*M_PI;
    double phi   = (1.0/N)*asinh(1.0/e);
    std::complex<double> P(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
    for(size_t K=2; K <= N; K++)
    {
        e     = sqrt(pow(10.0,r/10.0)-1.0);
        theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
        phi   = (1.0/N)*asinh(1.0/e);
        std::complex<double> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
        P *= -p;        
    }
    if(fmod(N,2) == 0) return P/sqrt(1 + e*e);
    return P;
}

// Chebyshev2 = 1/pole
std::complex<double> ChebyshevPole(double K, double N, double r=1.0)
{      
    double e     = sqrt(pow(10.0,r/10.0)-1.0);
    double theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    double phi   = (1.0/N)*asinh(1.0/e);
    std::complex<double> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return p;
}


// Chebyshev2 = 1/pole
std::complex<double> Chebyshev2Zeros(double K, double N, double r=1.0)
{      
    double uk    = ((2*K-1)/N)*(M_PI/2.0);
    std::complex<double> p(0,cos(uk));
    return 1.0/p;
}

// Chebyshev2 = 1/pole
std::complex<double> Chebyshev2Pole(double K, double N, double r=1.0)
{      
    double e     = 1.0/sqrt(pow(10.0,r/10.0)-1.0);
    double theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    double phi   = (1.0/N)*asinh(1.0/e);
    std::complex<double> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return 1.0/p;
}

std::vector<std::complex<double>> ChebyshevPoles(double N, double r)
{
    std::vector<std::complex<double>> out(N);
    for(size_t K = 1; K <= N; K++)
    {
        out.push_back(ChebyshevPole(K,N,r));
    }
    return out;
}

int TestLaguerre()
{
    deque<complex<double>> P = {1,2,3}, R;
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
  complex<double> P_x(0.,0.);
  double mean_err = 0, max_err = -1e12, abs_err;
  cout << "--------- Root analysis ---------" << endl;
  for (int j = 0; j < R.size(); j++)
  { 
    P_x = complex<double>(0,0);
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

// Hn(s) = H0/d * PROD( (s^2 + ai) / (s^2 + bi*s + ci))
// d = s+p0 n=odd
// d = 1    n=even

int main()
{
    int order,numSections;
    std::vector<double> a,b,c;
    double pZero,hZero;
    double omega_p = 3000.0;
    double omega_s = 3200.0;

    
    elliptical_coeffs(omega_p, omega_s, 0.1, 9, a,b,c,numSections,hZero,pZero);
    std::cout << numSections << std::endl;
    
    for(size_t i = 1; i < 5; i++)
        std::cout << a[i] << "," << b[i] << "," << c[i] << std::endl;

    double h0;
    if(order % 2 != 0) {        
        
        // d = s + p0
        // h0/(s + p0)
        
        // b0=h0/p0
        // b1=0
        // b2=0
        // a1=1
        // a2=/p0
        // a3=0
        
    }
    else 
    {
        h0 = hZero;
        // b0 = h0
        // b1 = 0;
        // b2 = 0;        
    }
    
    
    std::complex<double> p1  = ButterworthPoles(1,6);
    std::complex<double> p2  = ButterworthPoles(2,6);
    std::complex<double> p3  = ButterworthPoles(3,6);
    std::complex<double> p4  = ButterworthPoles(4,6);
    std::complex<double> p5  = ButterworthPoles(5,6);
    std::complex<double> p6  = ButterworthPoles(6,6);
    
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << p3 << std::endl;
    std::cout << p4 << std::endl;
    std::cout << p5 << std::endl;
    std::cout << p6 << std::endl;
    
    std::complex<double> H0  = Chebyshev2Zeros(1,2,1.0);
    std::complex<double> H1  = Chebyshev2Zeros(2,2,1.0);
    p1  = Chebyshev2Pole(1,2,1);
    p2  = Chebyshev2Pole(2,2,1);

    std::cout << H0 << std::endl;
    std::cout << H1 << std::endl;
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    
}