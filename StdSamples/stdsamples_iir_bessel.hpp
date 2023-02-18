#pragma once

#include "audio_iir_biquads.hpp"

/*
function genfilter(order)
  [Z,P,K] = besself(order,1,1,'s')
  for i = 1:size(P,1)
   printf("std::complex<DspFloatType>(%f,%f),\n",real(P(i)),imag(P(i)));
  endfor
endfunction
*/
namespace IIRFilters::Bessel
{
///////////////////////////////////
// Bessel filter pole/zero tables
///////////////////////////////////

    std::complex<DspFloatType> bessel_poles_2[] = { std::complex<DspFloatType>(-0.5001,0.8660), 
                                            std::complex<DspFloatType>(-0.5001,-0.8660)};
    std::complex<DspFloatType> bessel_poles_3[] = { std::complex<DspFloatType>(-0.9416,0),
                                            std::complex<DspFloatType>(-0.7456,0.7114),
                                            std::complex<DspFloatType>(-0.7456,-0.7114)};                                          
    std::complex<DspFloatType> bessel_poles_4[] = { std::complex<DspFloatType>(-0.6572,0.0302), 
                                            std::complex<DspFloatType>(-0.6572,0.0302), 
                                            std::complex<DspFloatType>(-0.9048,0.2709),                                          
                                            std::complex<DspFloatType>(-0.9048,-0.2709)};
    std::complex<DspFloatType> bessel_poles_5[] = { std::complex<DspFloatType>(-0.9264,0),
                                            std::complex<DspFloatType>(-0.5906,0.9072),
                                            std::complex<DspFloatType>(-0.5906,-0.9072),
                                            std::complex<DspFloatType>(-0.8516,0.4427),
                                            std::complex<DspFloatType>(-0.8515,-0.4427)};
    std::complex<DspFloatType> bessel_poles_6[] = { 
                                            std::complex<DspFloatType>(-0.5386,0.9617),
                                            std::complex<DspFloatType>(-0.5386,-0.9617),
                                            std::complex<DspFloatType>(-0.7997,0.5622),
                                            std::complex<DspFloatType>(-0.7997,-0.5622),
                                            std::complex<DspFloatType>(-0.9094,0.1857),
                                            std::complex<DspFloatType>(-0.9094,-0.1857)};
    std::complex<DspFloatType> bessel_poles_7[] = { std::complex<DspFloatType>(-0.9195,0),
                                            std::complex<DspFloatType>(-0.4967,1.0025), 
                                            std::complex<DspFloatType>(-0.4967,-1.0025),
                                            std::complex<DspFloatType>(-0.7527,0.6505),
                                            std::complex<DspFloatType>(-0.7527,-0.6505),
                                            std::complex<DspFloatType>(-0.8800,0.3217),
                                            std::complex<DspFloatType>(-0.8800,-0.3217)};
    std::complex<DspFloatType> bessel_poles_8[] = { 
                                            std::complex<DspFloatType>(-0.4622,1.0344),
                                            std::complex<DspFloatType>(-0.4622,-1.0344),
                                            std::complex<DspFloatType>(-0.7111,0.7187),
                                            std::complex<DspFloatType>(-0.7111,-0.7187),
                                            std::complex<DspFloatType>(-0.8473,0.4259),
                                            std::complex<DspFloatType>(-0.8473,-0.4259),
                                            std::complex<DspFloatType>(-0.9097,0.1412),
                                            std::complex<DspFloatType>(-0.9097,-0.1412) };
    


    std::complex<DspFloatType> *bessel_poles[] = {
        bessel_poles_2,
        bessel_poles_3,
        bessel_poles_4,
        bessel_poles_5,
        bessel_poles_6,
        bessel_poles_7,
        bessel_poles_8,
    };

/////////////////////////////////////////////////////////////////////////////////////////////
// Bessel Filter
/////////////////////////////////////////////////////////////////////////////////////////////    
    BiquadSOS bessellp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;    
        
        if(order <  2) order = 2;
        if(order >  8) order = 8;
        
        size_t n = 0;
        std::complex<DspFloatType> *poles = bessel_poles[order-2];

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1=poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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

    // todo polysolver
    BiquadSOS besselhp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;    
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        size_t n = 1;
        std::complex<DspFloatType> *poles = bessel_poles[order-2];

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = bessel_poles_2[0];
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
            std::complex<DspFloatType> p1  = bessel_poles_2[i];
            std::complex<DspFloatType> p2  = bessel_poles_2[i+1];
            
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
    BiquadSOS butterlp2hp(int order, double Q=1.0)
    {
        BiquadSOS sos;
        if(order <  2) order = 2;
        if(order >  8) order = 8;        
        size_t n = 1;
        std::complex<DspFloatType> *poles = bessel_poles[order-2];

                
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[0];
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
            std::complex<DspFloatType> p1  = poles[i];
            std::complex<DspFloatType> p2  = poles[i+1];
            
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
        if(order <  2) order = 2;
        if(order >  8) order = 8;
        
        size_t n = 1;
        std::complex<DspFloatType> *bpoles = bessel_poles[order-2];

                
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = bpoles[0];
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
            std::complex<DspFloatType> p1  = bpoles[i];
            std::complex<DspFloatType> p2  = bpoles[i+1];
            
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
        if(order <  2) order = 2;
        if(order >  8) order = 8;
        
        size_t n = 1;
        std::complex<DspFloatType> *bpoles = bessel_poles[order-2];
        
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = bpoles[0];
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
            std::complex<DspFloatType> p1  = bpoles[i];
            std::complex<DspFloatType> p2  = bpoles[i+1];
            
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2bs(zeros,poles,1.0,0.5,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain * abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain * abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q) * abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;

            sos.push_back(c);
        }        
        return sos;
    }
}