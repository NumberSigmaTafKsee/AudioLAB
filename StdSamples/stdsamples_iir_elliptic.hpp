#pragma once

#include "stdsamples_iir_biquads.hpp"

namespace IIRFilters::Ellliptic
{
    std::complex<DspFloatType> ellip_2_gain = std::complex<DspFloatType>(3.1620e-03,0);
    std::complex<DspFloatType> ellip_2_zeros[] = {
        std::complex<DspFloatType>(0.000000,7.294427),
        std::complex<DspFloatType>(-0.000000,-7.294427),
    };
    std::complex<DspFloatType> ellip_2_poles[] = {
        std::complex<DspFloatType>(-0.115703,0.720179),
        std::complex<DspFloatType>(-0.115703,-0.720179),
    };
    std::complex<DspFloatType> ellip_3_gain = std::complex<DspFloatType>(0.017775,0);
    std::complex<DspFloatType> ellip_3_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.280738),
        std::complex<DspFloatType>(-0.000000,-2.280738),
    };
    std::complex<DspFloatType> ellip_3_poles[] = {  
        std::complex<DspFloatType>(-0.117388,0.000000),
        std::complex<DspFloatType>(-0.049807,0.886094),
        std::complex<DspFloatType>(-0.049807,-0.886094),        
    };
    std::complex<DspFloatType> ellip_4_gain = std::complex<DspFloatType>(3.1620e-03,0);
    std::complex<DspFloatType> ellip_4_zeros[] = {
        std::complex<DspFloatType>(0.000000,3.035120),
        std::complex<DspFloatType>(0.000000,1.433499),
        std::complex<DspFloatType>(-0.000000,-3.035120),
        std::complex<DspFloatType>(-0.000000,-1.433499),
    };
    std::complex<DspFloatType> ellip_4_poles[] = {  
        std::complex<DspFloatType>(-0.083369,0.450273),
        std::complex<DspFloatType>(-0.022592,0.949851),
        std::complex<DspFloatType>(-0.083369,-0.450273),
        std::complex<DspFloatType>(-0.022592,-0.949851),
    };
    std::complex<DspFloatType> ellip_5_gain = std::complex<DspFloatType>(0.013033,0);
    std::complex<DspFloatType> ellip_5_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.594131),
        std::complex<DspFloatType>(0.000000,1.172108),
        std::complex<DspFloatType>(-0.000000,-1.594131),
        std::complex<DspFloatType>(-0.000000,-1.172108),
    };
    std::complex<DspFloatType> ellip_5_poles[] = {  
        std::complex<DspFloatType>(-0.091161,0.000000),
        std::complex<DspFloatType>(-0.049256,0.720872),
        std::complex<DspFloatType>(-0.010192,0.977717),
        std::complex<DspFloatType>(-0.049256,-0.720872),
        std::complex<DspFloatType>(-0.010192,-0.977717),
    };
    std::complex<DspFloatType> ellip_6_gain = std::complex<DspFloatType>(3.1628e-03,0);
    std::complex<DspFloatType> ellip_6_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.663135),
        std::complex<DspFloatType>(0.000000,1.226580),
        std::complex<DspFloatType>(0.000000,1.072668),
        std::complex<DspFloatType>(-0.000000,-2.663135),
        std::complex<DspFloatType>(-0.000000,-1.226580),
        std::complex<DspFloatType>(-0.000000,-1.072668),
    };
    std::complex<DspFloatType> ellip_6_poles[] = {  
        std::complex<DspFloatType>(-0.074612,0.401042),
        std::complex<DspFloatType>(-0.025369,0.867234),
        std::complex<DspFloatType>(-0.004554,0.990115),
        std::complex<DspFloatType>(-0.074612,-0.401042),
        std::complex<DspFloatType>(-0.025369,-0.867234),
        std::complex<DspFloatType>(-0.004554,-0.990115),
    };
    std::complex<DspFloatType> ellip_7_gain = std::complex<DspFloatType>(0.012329,0);
    std::complex<DspFloatType> ellip_7_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.505495),
        std::complex<DspFloatType>(0.000000,1.094314),
        std::complex<DspFloatType>(0.000000,1.031498),
        std::complex<DspFloatType>(-0.000000,-1.505495),
        std::complex<DspFloatType>(-0.000000,-1.094314),
        std::complex<DspFloatType>(-0.000000,-1.031498),
    };
    std::complex<DspFloatType> ellip_7_poles[] = {  
        std::complex<DspFloatType>(-0.086432,0.000000),
        std::complex<DspFloatType>(-0.047099,0.684688),
        std::complex<DspFloatType>(-0.012071,0.939200),
        std::complex<DspFloatType>(-0.002024,0.995621),
        std::complex<DspFloatType>(-0.047099,-0.684688),
        std::complex<DspFloatType>(-0.012071,-0.939200),
        std::complex<DspFloatType>(-0.002024,-0.995621),        
    };
    std::complex<DspFloatType> ellip_8_gain = std::complex<DspFloatType>(3.1628e-03,0);
    std::complex<DspFloatType> ellip_8_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.599225),
        std::complex<DspFloatType>(0.000000,1.196001),
        std::complex<DspFloatType>(0.000000,1.040615),
        std::complex<DspFloatType>(0.000000,1.013796),
        std::complex<DspFloatType>(-0.000000,-2.599225),
        std::complex<DspFloatType>(-0.000000,-1.196001),
        std::complex<DspFloatType>(-0.000000,-1.040615),
        std::complex<DspFloatType>(-0.000000,-1.013796),
    };
    std::complex<DspFloatType> ellip_8_poles[] = {  
        std::complex<DspFloatType>(-0.072877,0.391644),
        std::complex<DspFloatType>(-0.024974,0.847710),
        std::complex<DspFloatType>(-0.005518,0.972693),
        std::complex<DspFloatType>(-0.000896,0.998063),
        std::complex<DspFloatType>(-0.072877,-0.391644),
        std::complex<DspFloatType>(-0.024974,-0.847710),
        std::complex<DspFloatType>(-0.005518,-0.972693),
        std::complex<DspFloatType>(-0.000896,-0.998063),
    };
    std::complex<DspFloatType> ellip_gains[] = {
        ellip_2_gain,
        ellip_3_gain,
        ellip_4_gain,
        ellip_5_gain,
        ellip_6_gain,
        ellip_7_gain,
        ellip_8_gain,
    };
    std::complex<DspFloatType> *ellip_zeros[] = {
        ellip_2_zeros,
        ellip_3_zeros,
        ellip_4_zeros,
        ellip_5_zeros,
        ellip_6_zeros,
        ellip_7_zeros,
        ellip_8_zeros,
    };
    std::complex<DspFloatType> *ellip_poles[] = {
        ellip_2_poles,
        ellip_3_poles,
        ellip_4_poles,
        ellip_5_poles,
        ellip_6_poles,
        ellip_7_poles,
        ellip_8_poles,
    };


/////////////////////////////////////////////////////////////////////////////////////////////
// Elliptic Filter
/////////////////////////////////////////////////////////////////////////////////////////////    
    BiquadSOS elliplp(int order, double Q=1.0)
    {        
        BiquadSOS sos;  
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        std::complex<DspFloatType> * czeros = ellip_zeros[order-1];          
        std::complex<DspFloatType> * cpoles = ellip_poles[order-1];          
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = ellip_gains[order-1]*czeros[n];
            std::complex<DspFloatType> p1  = cpoles[n++];

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
            std::complex<DspFloatType> H0  = ellip_gains[order-1]*czeros[n]; 
            std::complex<DspFloatType> H1  = ellip_gains[order-1]*czeros[n+1]; 
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
            // radius is the same thing but goes from 0..1
            // 0 = most resonant
            // 1 = least resonant
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            
            sos.push_back(c);
        }
        return sos;
    }    
};
