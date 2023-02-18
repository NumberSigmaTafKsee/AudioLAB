#pragma once

// Chebyshev2 = 1/pole
std::complex<DspFloatType> Chebyshev2Zeros(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType uk    = ((2*K-1)/N)*(M_PI/2.0);
    std::complex<DspFloatType> p(0,cos(uk));
    return 1.0/p;
}

// Chebyshev2 = 1/pole
std::complex<DspFloatType> Chebyshev2Pole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType e     = 1.0/sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return 1.0/p;
}
    

/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 2
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby2lp(int order, double Q=1.0,double rips=1.0)
    {        
        BiquadSOS sos;            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,order,1.0);
            std::complex<DspFloatType> p1  = Chebyshev2Pole(1,order);

            DspFloatType x1 = abs(p1);
            DspFloatType x2 = abs(H0);
                    
            // (s-p1)        
            c.z[0] = x2/x1;
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
            std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,2,1.0);
            std::complex<DspFloatType> H1  = Chebyshev2Zeros(2,2,1.0);
            std::complex<DspFloatType> p1  = Chebyshev2Pole(1,2,1);
            std::complex<DspFloatType> p2  = Chebyshev2Pole(2,2,1);
    
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            DspFloatType z1 = abs(H0*H1);
            DspFloatType z2 = abs(-H0-H1);

            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;

            c.z[0] = z1/x1;
            c.z[1] = z2/x1;
            c.z[2] = 0;
            c.p[0] = 1;
            // radius is the same thing but goes from 0..1
            // 0 = most resonant
            // 1 = least resonant
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            
            sos.push_back(c);
        }
        return sos;
    }    
    

