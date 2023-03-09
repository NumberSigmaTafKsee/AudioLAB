#pragma once

#include "stdsamples_iir_biquads.hpp"

/*
function genfilter(order)
  [Z,P,K] = cheby2(order,1,1,'s')
  for i = 1:size(Z,1)
   printf("std::complex<DspFloatType>(%f,%f),\n",real(Z(i)),imag(Z(i)));
  endfor
  for i = 1:size(P,1)
   printf("std::complex<DspFloatType>(%f,%f),\n",real(P(i)),imag(P(i)));
  endfor
endfunction
*/
namespace IIRFilters::Chebyshev2
{
    std::complex<DspFloatType> cheby2_2_gain = std::complex<DspFloatType>(0.8913,0);
    std::complex<DspFloatType> cheby2_2_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.414214),
        std::complex<DspFloatType>(-0.000000,-1.414214),
    };
    std::complex<DspFloatType> cheby2_2_poles[] = {
        std::complex<DspFloatType>(-0.311324,-1.298299),
        std::complex<DspFloatType>(-0.311324,1.298299),
    };
    std::complex<DspFloatType> cheby2_3_gain = std::complex<DspFloatType>(5.8957,0);
    std::complex<DspFloatType> cheby2_3_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.154701),
        std::complex<DspFloatType>(-0.000000,-1.154701),
    };
    std::complex<DspFloatType> cheby2_3_poles[] = {
        std::complex<DspFloatType>(-6.106489,-0.000000),
        std::complex<DspFloatType>(-0.105405,-1.129687),        
        std::complex<DspFloatType>(-0.105405,1.129687),
    };
    std::complex<DspFloatType> cheby2_4_gain = std::complex<DspFloatType>(0.8913,0);
    std::complex<DspFloatType> cheby2_4_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.082392),
        std::complex<DspFloatType>(0.000000,2.613126),
        std::complex<DspFloatType>(-0.000000,-2.613126),
        std::complex<DspFloatType>(-0.000000,-1.082392),
    };
    std::complex<DspFloatType> cheby2_4_poles[] = {
        std::complex<DspFloatType>(-0.054008,-1.071629),
        std::complex<DspFloatType>(-0.701365,-2.387691),
        std::complex<DspFloatType>(-0.701365,2.387691),
        std::complex<DspFloatType>(-0.054008,1.071629),
    };
    std::complex<DspFloatType> cheby2_5_gain = std::complex<DspFloatType>(9.8261,0);
    std::complex<DspFloatType> cheby2_5_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.051462),
        std::complex<DspFloatType>(0.000000,1.701302),
        std::complex<DspFloatType>(-0.000000,-1.701302),
        std::complex<DspFloatType>(-0.000000,-1.051462),
    };
    std::complex<DspFloatType> cheby2_5_poles[] = {
        std::complex<DspFloatType>(-10.206345,-0.000000),
        std::complex<DspFloatType>(-0.033122,-1.045402),
        std::complex<DspFloatType>(-0.223227,-1.663234),        
        std::complex<DspFloatType>(-0.223227,1.663234),
        std::complex<DspFloatType>(-0.033122,1.045402),
    };
    std::complex<DspFloatType> cheby2_6_gain = std::complex<DspFloatType>(0.8913,0);
    std::complex<DspFloatType> cheby2_6_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.035276),
        std::complex<DspFloatType>(0.000000,1.414214),
        std::complex<DspFloatType>(0.000000,3.863703),
        std::complex<DspFloatType>(-0.000000,-3.863703),
        std::complex<DspFloatType>(-0.000000,-1.414214),
        std::complex<DspFloatType>(-0.000000,-1.035276),
    };
    std::complex<DspFloatType> cheby2_6_poles[] = {
        std::complex<DspFloatType>(-0.022478,-1.031356),
        std::complex<DspFloatType>(-0.113895,-1.400264),
        std::complex<DspFloatType>(-1.070345,-3.525988),
        std::complex<DspFloatType>(-1.070345,3.525988),
        std::complex<DspFloatType>(-0.113895,1.400264),
        std::complex<DspFloatType>(-0.022478,1.031356),
    };
    std::complex<DspFloatType> cheby2_7_gain = std::complex<DspFloatType>(13.757,0);
    std::complex<DspFloatType> cheby2_7_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.025717),
        std::complex<DspFloatType>(0.000000,1.279048),
        std::complex<DspFloatType>(0.000000,2.304765),
        std::complex<DspFloatType>(-0.000000,-2.304765),
        std::complex<DspFloatType>(-0.000000,-1.279048),
        std::complex<DspFloatType>(-0.000000,-1.025717),
    };
    std::complex<DspFloatType> cheby2_7_poles[] = {
        std::complex<DspFloatType>(-14.300043,-0.000000),
        std::complex<DspFloatType>(-0.016288,-1.022959),
        std::complex<DspFloatType>(-0.070763,-1.271995),
        std::complex<DspFloatType>(-0.326203,-2.251897),        
        std::complex<DspFloatType>(-0.326203,2.251897),
        std::complex<DspFloatType>(-0.070763,1.271995),
        std::complex<DspFloatType>(-0.016288,1.022959),
    };
    std::complex<DspFloatType> cheby2_8_gain = std::complex<DspFloatType>(0.8913,0);
    std::complex<DspFloatType> cheby2_8_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.019591),
        std::complex<DspFloatType>(0.000000,1.202690),
        std::complex<DspFloatType>(0.000000,1.799952),
        std::complex<DspFloatType>(0.000000,5.125831),
        std::complex<DspFloatType>(-0.000000,-5.125831),
        std::complex<DspFloatType>(-0.000000,-1.799952),
        std::complex<DspFloatType>(-0.000000,-1.202690),
        std::complex<DspFloatType>(-0.000000,-1.019591),
    };
    std::complex<DspFloatType> cheby2_8_poles[] = {
        std::complex<DspFloatType>(-0.012359,-1.017538),
        std::complex<DspFloatType>(-0.048898,-1.198450),
        std::complex<DspFloatType>(-0.162825,-1.781713),
        std::complex<DspFloatType>(-1.435344,-4.675639),
        std::complex<DspFloatType>(-1.435344,4.675639),
        std::complex<DspFloatType>(-0.162825,1.781713),
        std::complex<DspFloatType>(-0.048898,1.198450),
        std::complex<DspFloatType>(-0.012359,1.017538),
    };
    std::complex<DspFloatType> cheby2_gains[] = {
        cheby2_2_gain,
        cheby2_3_gain,
        cheby2_4_gain,
        cheby2_5_gain,
        cheby2_6_gain,
        cheby2_7_gain,
        cheby2_8_gain,
    };
    std::complex<DspFloatType> *cheby2_zeros[] = {
        cheby2_2_zeros,
        cheby2_3_zeros,
        cheby2_4_zeros,
        cheby2_5_zeros,
        cheby2_6_zeros,
        cheby2_7_zeros,
        cheby2_8_zeros,
    };
    std::complex<DspFloatType> *cheby2_poles[] = {
        cheby2_2_poles,
        cheby2_3_poles,
        cheby2_4_poles,
        cheby2_5_poles,
        cheby2_6_poles,
        cheby2_7_poles,
        cheby2_8_poles,
    };


/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 2
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby2lp(int order, double Q=1.0,double rips=1.0)
    {        
        BiquadSOS sos;  
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        std::complex<DspFloatType> * czeros = cheby2_zeros[order-1];          
        std::complex<DspFloatType> * cpoles = cheby2_poles[order-1];          
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = cheby2_gains[order-1]*czeros[n];
            std::complex<DspFloatType> p1  = cpoles[n];

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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> H0  = cheby2_gains[order-1]*czeros[n]; 
            std::complex<DspFloatType> H1  = cheby2_gains[order-1]*czeros[n+1]; 
            std::complex<DspFloatType> p1  = cpoles[n++];
            std::complex<DspFloatType> p2  = cpoles[n++];
    
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
            
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            
            sos.push_back(c);
        }
        return sos;
    }    
}
