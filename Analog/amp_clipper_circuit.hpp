// Title:   clipperCircuit
// Author:  Ben Holmes
// Date:    2017/12/19
// License: GPL v3.0

// A class that implements a diode clipper circuits
#pragma once
#include <cmath>
#include "GenericSoundObject.hpp"

namespace FX::Distortion::ClipperCircuit
{
    template<class T>
    class ClipperCircuit : public GSSoundProcessor<T>
    {
    public:
        
        // Default public constructor
        ClipperCircuit() : ClipperCircuit(48e3) {};
        
        // Delegated constructor for use with Sample Rate
        ClipperCircuit(T fs)
        : ClipperCircuit(fs, 1., 1.) {};
        
        // Constructor for
        ClipperCircuit(T fs, T ideality, T asymmetry)
        : sFreq(fs), Is(1e-7), C1(1e-6), asym(asymmetry),
        state(0.), yn1(0.), GSSoundProcessor<T>()
        {
            setIdeality(ideality);
        }
        
        // Destructor
        ~ClipperCircuit(){};
        
        // Single sample calculation
        T run(const T uVal)
        {
            T residual = 10;
            unsigned int numIters = 0;
            T step = 0.;
            T funcValue = 0.;
            T derivativeValue = 0.;
            
            // Set last solution as current iterate.
            T yCur = yn1;
            
            // Combine constant terms.
            const T pConst = uVal/R1 + state;
            
            T ePos = 0.;
            T eNeg = 0.;
            // Find solution iteratively
            while ((residual>tolerance) && (numIters < maximumIters))
            {
                // Precalculate exponentials as they are used twice.
                ePos = std::exp(yCur/(NVt));
                eNeg = std::exp(-yCur/(asym*NVt));
                
                // Calculate step
                funcValue = pConst - (yCur/R1)
                - ((2.*C1*sFreq)*yCur) - (Is*(ePos - eNeg));
                
                derivativeValue = (-1./R1) - (2.*C1*sFreq)
                - (Is*((ePos/(NVt)) + (eNeg/(asym*NVt))));
                
                step = -funcValue/derivativeValue;
                
                // Residual is only step size.
                residual = std::fabs(step);
                
                yCur += step;
                
                numIters++;
            }
            
            // Fail-safe if iterating explodes.
            if(yCur!=yCur)
            {
                yCur = 0;
            }
            
            // Update state and past-sample
            state = (4.*C1*sFreq)*yCur - state;
            yn1 = yCur;
            
            return yCur * getForwardVoltage();
        }
        void ProcessSIMD(size_t n, T * in, T * out)
		{
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = run(in[i]);
		}
		void ProcessBlock(size_t n, T * in, T * out)
		{
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = run(in[i]);
		}
		void ProcessInplace(size_t n, T * buffer)
		{
			#pragma omp simd aligned(buffer)
			for(size_t i = 0; i < n; i++) buffer[i] = run(buffer[i]);
		}
	
        // Set ideality factor, combined with thermal voltage
        void setIdeality(const T idealityFactor)
        {
            NVt = idealityFactor*Vt;
        }
        
        // Set capacitor value.
        void setCapacitance(const T capacitance)
        {
            C1 = capacitance;
        }
        
        // Forward voltage generally indicates voltage amplitude of
        // the diode clipper, and can be used to scale the gain of
        // and effect.
        T getForwardVoltage() const
        {
            return NVt*std::log(1.e-3/(2*Is) + 1);
        }
        
        // Set asymmetry in the exponent of the diodes.
        void setAsymmetry(const T asymmetry)
        {
            asym = asymmetry;
        }
        
    protected:
        
        // Approximate exp as lim( 1 + x/n)^n as n tends to infinity
        // Using n = 1024;
        T fastExp(const T index)
        {
            
            double exp = 1. + index/1024.;
            exp *= exp; exp *= exp; exp *= exp;
            exp *= exp; exp *= exp; exp *= exp;
            exp *= exp; exp *= exp; exp *= exp;
            exp *= exp;
            
            return exp;
        }
        
        // Sampling frequency
        T sFreq;
        
        // Circuit parameters
        T Is;               // Saturation current
        T NVt;              // Ideality factor combined with thermal voltage
        T C1;               // Capacitor
        T asym = 1;         // Multiplier of one ideality
        
        // Circuit constants
        const T Vt = 25.83e-3; // Thermal voltage
        const T R1 = 2200;     // Resistor
        
        // Model state
        T state;   // capacitor state
        T yn1;     // previous iterative solution
        
        // Iterative solver parameters
        const T tolerance = 1e-6;
        const unsigned int maximumIters = 100;
        const unsigned int maxSubIters = 5;
    };
}
