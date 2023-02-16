#pragma once

/*
function genfilter(order)
  [Z,P,K] = butter(order,1,'s')
  for i = 1:size(P,1)
   printf("std::complex<DspFloatType>(%f,%f),\n",real(P(i)),imag(P(i)));
  endfor

endfunction
*/

#include "stdsamples_iir_biquads.hpp"
namespace IIRFilters::Butterworth
{
    std::complex<DspFloatType> butter_2_poles[] = {
        std::complex<DspFloatType>(-0.707107,0.707107),
        std::complex<DspFloatType>(-0.707107,-0.707107)
    };
    std::complex<DspFloatType> butter_3_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),
        std::complex<DspFloatType>(-0.500000,0.866025),   
        std::complex<DspFloatType>(-0.500000,-0.866025),
    };
    std::complex<DspFloatType> butter_4_poles[] = {
        std::complex<DspFloatType>(-0.382683,0.923880),
        std::complex<DspFloatType>(-0.923880,0.382683),
        std::complex<DspFloatType>(-0.923880,-0.382683),
        std::complex<DspFloatType>(-0.382683,-0.923880),
    };
    std::complex<DspFloatType> butter_5_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),
        std::complex<DspFloatType>(-0.309017,0.951057),
        std::complex<DspFloatType>(-0.809017,0.587785),
        std::complex<DspFloatType>(-0.809017,-0.587785),
        std::complex<DspFloatType>(-0.309017,-0.951057),
    };
    std::complex<DspFloatType> butter_6_poles[] = {
        std::complex<DspFloatType>(-0.258819,0.965926),
        std::complex<DspFloatType>(-0.707107,0.707107),
        std::complex<DspFloatType>(-0.965926,0.258819),
        std::complex<DspFloatType>(-0.965926,-0.258819),
        std::complex<DspFloatType>(-0.707107,-0.707107),
        std::complex<DspFloatType>(-0.258819,-0.965926),
    };
    std::complex<DspFloatType> butter_7_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),    
        std::complex<DspFloatType>(-0.222521,0.974928),
        std::complex<DspFloatType>(-0.623490,0.781831),
        std::complex<DspFloatType>(-0.900969,0.433884),    
        std::complex<DspFloatType>(-0.900969,-0.433884),
        std::complex<DspFloatType>(-0.623490,-0.781831),
        std::complex<DspFloatType>(-0.222521,-0.974928),
    };
    std::complex<DspFloatType> butter_8_poles[] = {
        std::complex<DspFloatType>(-0.195090,0.980785),
        std::complex<DspFloatType>(-0.555570,0.831470),
        std::complex<DspFloatType>(-0.831470,0.555570),
        std::complex<DspFloatType>(-0.980785,0.195090),
        std::complex<DspFloatType>(-0.980785,-0.195090),
        std::complex<DspFloatType>(-0.831470,-0.555570),
        std::complex<DspFloatType>(-0.555570,-0.831470),
        std::complex<DspFloatType>(-0.195090,-0.980785),
    };
    std::complex<DspFloatType> *butter_poles[] = {
        butter_2_poles,
        butter_3_poles,
        butter_4_poles,
        butter_5_poles,
        butter_6_poles,
        butter_7_poles,
        butter_8_poles,
    };

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Butterworth
    /////////////////////////////////////////////////////////////////////////////////////////////
        BiquadSOS butterlp(int order, double Q=1.0)
        {
            BiquadSOS sos; 
            if(order <= 2) order = 2;
            if(order >  8) order = 8;
            std::complex<DspFloatType> * poles = butter_poles[order-2];
            size_t n = 0;
            if(order %2 != 0) {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];
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
        BiquadSOS butterhp(int order, double Q=1.0)
        {
            BiquadSOS sos;    
            if(order <= 2) order = 2;
            if(order >  8) order = 8;
            std::complex<DspFloatType> * poles = butter_poles[order-2];
            size_t n = 0;

            if(order %2 != 0) {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];
                DspFloatType x1 = abs(p1);
                

                c.z[0] = 0.0;
                c.z[1] = 1.0/x1;
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
                c.z[2] = 1.0/x1;
                c.p[0] = 1;        
                c.p[1] = (1.0/Q)*x2/x1;
                c.p[2] = 1.0/x1;    
                sos.push_back(c);
            }
            
            return sos;
        }

        BiquadSOS butterlp2hp(int order, double Q=1.0)
        {
            BiquadSOS sos;    
            if(order <= 2) order = 2;
            if(order >  8) order = 8;
            std::complex<DspFloatType> * poles = butter_poles[order-2];
            size_t n = 0;

            if(order %2 != 0) {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];
                DspFloatType x1 = abs(p1);
                

                c.z[0] = 0.0;
                c.z[1] = 1.0/x1;
                c.z[2] = 0.0;
                c.p[0] = 1.0;
                c.p[1] = 1/x1;
                c.p[2] = 0.0;

                sos.push_back(c);                        
            }
                
            for(size_t i = n; i < order; i += 2)
            {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];
                std::complex<DspFloatType> p2  = poles[n++];
                
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
            if(order <= 2) order = 2;
            if(order >  8) order = 8;
            std::complex<DspFloatType> * poles = butter_poles[order-2];
            size_t n = 0;
            for(size_t i = 1; i < order; i += 2)
            {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];        
                std::complex<DspFloatType> p2  = poles[n++];
                std::vector<std::complex<DspFloatType>> zeros,poles;
                poles.push_back(p1);
                poles.push_back(p2);
                DspFloatType gain;
                lp2bp(zeros,poles,1.0,0.5,gain);
                // (1-z)(1-z) = 1 -z1-z2 +z1z2
                
                DspFloatType x1 = abs(poles[0]*poles[1]);
                c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
                c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
                c.z[2] = gain/x1;
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
            if(order <= 2) order = 2;
            if(order >  8) order = 8;
            std::complex<DspFloatType> * poles = butter_poles[order-2];
            size_t n = 0;
            for(size_t i = 1; i < order; i += 2)
            {
                BiquadSection c;
                std::complex<DspFloatType> p1  = poles[n++];
                std::complex<DspFloatType> p2  = poles[n++];
                std::vector<std::complex<DspFloatType>> zeros,poles;
                poles.push_back(p1);
                poles.push_back(p2);
                DspFloatType gain;
                lp2bs(zeros,poles,1.0,0.5,gain);
                // (1-z)(1-z) = 1 -z1-z2 +z1z2
                
                DspFloatType x1 = abs(poles[0]*poles[1]);
                c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
                c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
                c.z[2] = gain/x1;
                c.p[0] = 1.0;
                c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
                c.p[2] = 1.0/x1;

                sos.push_back(c);
            }
            
            return sos;
        }
        
        BiquadSection butter2lp(double Q=1.0)
        {   
            std::complex<DspFloatType> * poles = butter_poles[0];

            std::complex<DspFloatType> p1  = poles[0];
            std::complex<DspFloatType> p2  = poles[1];
        
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            //std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 1.0/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1;    
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1.0/x1;        

            return c;
        }
        BiquadSection butter2hp(double Q=1.0)
        {    
            std::complex<DspFloatType> * poles = butter_poles[0];

            std::complex<DspFloatType> p1  = poles[0];
            std::complex<DspFloatType> p2  = poles[1];
        
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            //std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 0.0;
            c.z[1] = 0.0;
            c.z[2] = 1.0/x1;
            c.p[0] = 1;    
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1.0/x1;        

            return c;
        }
        BiquadSection butterlp2hp2(double Q=1.0)
        {    
            std::complex<DspFloatType> * bpoles = butter_poles[0];

            std::complex<DspFloatType> p1  = bpoles[0];
            std::complex<DspFloatType> p2  = bpoles[1];
        
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2hp(zeros,poles,1.0,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            BiquadSection c;        
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;

            return c;
        }
        BiquadSection butterlp2bp2(double Q=1.0)
        {    
            std::complex<DspFloatType> * bpoles = butter_poles[0];

            std::complex<DspFloatType> p1  = bpoles[0];
            std::complex<DspFloatType> p2  = bpoles[1];
        
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            // i dont not really know what this should be normalized 1.0,0 or 1.0,0.5?
            lp2bp(zeros,poles,1.0,0.5,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            BiquadSection c;        
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;
            return c;
        }
        BiquadSection butterlp2bs2(double Q=1.0)
        {    
            std::complex<DspFloatType> * bpoles = butter_poles[0];

            std::complex<DspFloatType> p1  = bpoles[0];
            std::complex<DspFloatType> p2  = bpoles[1];
        
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2bs(zeros,poles,1.0,0.5,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            BiquadSection c;        
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;
            return c;
        }
}
