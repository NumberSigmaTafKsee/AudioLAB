std::complex<DspFloatType> ChebyshevH0(DspFloatType N, DspFloatType r=1.0)
{        
    DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*1-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> P(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
    for(size_t K=2; K <= N; K++)
    {
        e     = sqrt(pow(10.0,r/10.0)-1.0);
        theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
        phi   = (1.0/N)*asinh(1.0/e);
        std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
        P *= -p;        
    }
    if(fmod(N,2) == 0) return P/sqrt(1 + e*e);
    return P;
}

// Chebyshev2 = 1/pole
std::complex<DspFloatType> ChebyshevPole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return p;
}

std::vector<std::complex<DspFloatType>> ChebyshevPoles(DspFloatType N, DspFloatType r)
{
    std::vector<std::complex<DspFloatType>> out(N);
    for(size_t K = 1; K <= N; K++)
    {
        out.push_back(ChebyshevPole(K,N,r));
    }
    return out;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 1
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby1lp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;            
        std::complex<DspFloatType> H0  = ChebyshevH0(order);            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ChebyshevPole(1,order);
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = abs(H0)/x1;
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
            std::complex<DspFloatType> p1  = ChebyshevPole(i,order);
            std::complex<DspFloatType> p2  = ChebyshevPole(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = abs(H0)/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }

    BiquadSOS cheby1hp(int order, double Q=1.0)
    {
        std::complex<DspFloatType> H0  = ChebyshevH0(order);            
        BiquadSOS sos;            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ChebyshevPole(1,order);
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = 0.0;
            c.z[1] = abs(H0)/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = ChebyshevPole(i,order);
            std::complex<DspFloatType> p2  = ChebyshevPole(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 0.0;
            c.z[1] = 0.0;
            c.z[2] = abs(H0)/x1;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }



