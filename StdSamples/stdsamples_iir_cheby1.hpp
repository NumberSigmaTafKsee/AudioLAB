#pragma once
#include "audio_iir_biquads.hpp"

/*
function genfilter(order)
  [Z,P,K] = cheby1(order,1,1,'s')
  for i = 1:size(P,1)
   printf("std::complex<DspFloatType>(%f,%f),\n",real(P(i)),imag(P(i)));
  endfor
endfunction
*/
namespace IIRFilters::Chebyshev1
{
    std::complex<DSPFloatType> cheby1_2_gain(0.9826,0);
    std::complex<DspFloatType> cheby1_2_poles[] = {
        std::complex<DspFloatType>(-0.548867,-0.895129),
        std::complex<DspFloatType>(-0.548867,0.895129),
    };
    std::complex<DSPFloatType> cheby1_3_gain(0.4913,0);
    std::complex<DspFloatType> cheby1_3_poles[] = {
        std::complex<DspFloatType>(-0.494171,0.000000),
        std::complex<DspFloatType>(-0.247085,-0.965999),        
        std::complex<DspFloatType>(-0.247085,0.965999),
    };
    std::complex<DSPFloatType> cheby1_4_gain(2.4565e-01,6.1843e-18);
    std::complex<DspFloatType> cheby1_4_poles[] = {
        std::complex<DspFloatType>(-0.139536,-0.983379),
        std::complex<DspFloatType>(-0.336870,-0.407329),
        std::complex<DspFloatType>(-0.336870,0.407329),
        std::complex<DspFloatType>(-0.139536,0.983379),
    };
    std::complex<DSPFloatType> cheby1_5_gain(1.2283e-01,8.6736e-18);
    std::complex<DspFloatType> cheby1_5_poles[] = {
        std::complex<DspFloatType>(-0.089458,-0.990107),
        std::complex<DspFloatType>(-0.234205,-0.611920),
        std::complex<DspFloatType>(-0.289493,0.000000),
        std::complex<DspFloatType>(-0.234205,0.611920),
        std::complex<DspFloatType>(-0.089458,0.990107),
    };
    std::complex<DSPFloatType> cheby1_6_gain(6.1413e-02,-4.6382e-18);
    std::complex<DspFloatType> cheby1_6_poles[] = {
        std::complex<DspFloatType>(-0.062181,-0.993411),
        std::complex<DspFloatType>(-0.169882,-0.727227),
        std::complex<DspFloatType>(-0.232063,-0.266184),
        std::complex<DspFloatType>(-0.232063,0.266184),
        std::complex<DspFloatType>(-0.169882,0.727227),
        std::complex<DspFloatType>(-0.062181,0.993411),
    };
    std::complex<DSPFloatType> cheby1_7_gain(3.0707e-02, 1.3010e-18);
    std::complex<DspFloatType> cheby1_7_poles[] = {
        std::complex<DspFloatType>(-0.205414,0.000000),
        std::complex<DspFloatType>(-0.045709,-0.995284),
        std::complex<DspFloatType>(-0.128074,-0.798156),
        std::complex<DspFloatType>(-0.185072,-0.442943),        
        std::complex<DspFloatType>(-0.185072,0.442943),
        std::complex<DspFloatType>(-0.128074,0.798156),
        std::complex<DspFloatType>(-0.045709,0.995284),
    };
    std::complex<DSPFloatType> cheby1_8_gain(1.5353e-02,8.6967e-19);
    std::complex<DspFloatType> cheby1_8_poles[] = {
        std::complex<DspFloatType>(-0.035008,-0.996451),
        std::complex<DspFloatType>(-0.099695,-0.844751),
        std::complex<DspFloatType>(-0.149204,-0.564444),
        std::complex<DspFloatType>(-0.175998,-0.198206),
        std::complex<DspFloatType>(-0.175998,0.198206),
        std::complex<DspFloatType>(-0.149204,0.564444),
    };
    std::complex<DSPFloatType> cheby1_gains[] = {
        cheby1_2_gain,
        cheby1_3_gain,
        cheby1_4_gain,
        cheby1_5_gain,
        cheby1_6_gain,
        cheby1_7_gain,
        cheby1_8_gain,
    };
    std::complex<DSPFloatType> *cheby1_poles[] = {
        cheby1_2_poles,
        cheby1_3_poles,
        cheby1_4_poles,
        cheby1_5_poles,
        cheby1_6_poles,
        cheby1_7_poles,
        cheby1_8_poles,
    };

/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 1
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby1lp(int order, double Q=1.0)
    {        
        BiquadSOS sos;
        if(order < 2) order = 2;            
        if(order > 7) order = 8;
        std::complex<DspFloatType> H0  = cheby1_gains[order-2];            
        std::complex<DspFloatType> * poles = cheby1_poles[order-2];
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
        BiquadSOS sos;
        if(order < 2) order = 2;            
        if(order > 7) order = 8;
        std::complex<DspFloatType> H0  = cheby1_gains[order-2];            
        std::complex<DspFloatType> * poles = cheby1_poles[order-2];
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];            
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
