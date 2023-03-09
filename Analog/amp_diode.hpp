#pragma once
#include "GenericSoundObject.hpp"

namespace FX::Distortion::Diode
{
    ///////////////////////////////////////////////////////////////////////////
    // Diode Saturation
    ///////////////////////////////////////////////////////////////////////////
    inline DspFloatType Diode(DspFloatType x, DspFloatType Vt = 0.0253,DspFloatType eta = 1.68,DspFloatType Is = .105)
    {
        return Is * (exp(0.1*x/(eta*Vt))-1);
    }

    template<typename T>
    struct DiodeClipperNR : public GSSoundProcessor<T>
    {
        T controlledR;
        T Id, C, Ve, Vp, R, err;
        T Fs, T;
        T vNom, vDenom;
        T vin;
        T vout, voutTemp, voutOld;
        T beta, betaM1;    
        int oversample;
        std::vector<T> blockInput, blockOutput, blockOutputDownsampled;
        T oldBlockOutput;
        //juce::dsp::ProcessorDuplicator< juce::dsp::IIR::Filter <T>, juce::dsp::IIR::Coefficients<T>> lowPassFilter;

        DiodeClipperNR(T sampleRate=44100.0) : GSSoundProcessor<T>()
        {
            // Use this method as the place to do any pre-playback
            // initialisation that you need..
            //oversampling->reset();
            //oversampling->initProcessing(static_cast<size_t> (samplesPerBlock));
            oversample = 16;
            // Set the constants
            Fs = sampleRate;
            T = 1 / Fs;
            C = 100e-9;
            R = 1e3;
            Vp = 0.17;
            Ve = 0.023;
            Id = 2.52e-9;
            err = 0.1e-3; // err for stopping iterations
            // set controlled values to starting values (redundant maybe delit later)
            controlledR = 1.0;
            voutOld = 0;
            beta = 0.125;
            betaM1 = 1 - beta;
        }
        T gdExp(T vc)
        {
            return Id * (std::exp(vc / (2.0 * Ve)) - 1.0);
        }
        T gdExpDiff(T vc)
        {
            return (Id * std::exp(vc / (2.0 * Ve))) / (2.0 * Ve);
        }
        T limiter(T val)
        {
            if (val < -1.0)
            {
                val = -1.0;
                return val;
            }
            else if (val > 1.0)
            {
                val = 1.0;
                return val;
            }
            return val;
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            voutTemp = 1;
            vout = 0;
            vin = controlledR * I;
            
            int itter = 0;
            while (std::abs(voutTemp - vout) > err && itter < 20) {
                voutTemp = vout;
                vNom = T * voutTemp * R * gdExpDiff(-voutTemp) + T * R * gdExp(-voutTemp) + voutOld * R * C + T * (gdExpDiff(voutTemp) * R * voutTemp - R * gdExp(voutTemp) + vin);
                vDenom = T * R * gdExpDiff(voutTemp) + T * R * gdExpDiff(-voutTemp) + R * C + T;
                vout = vNom / vDenom;
                itter++;
            }        
            voutOld = vout;
            
            return limiter(vout);
        }
	void ProcessSIMD(size_t n, T * in, T * out) {
		#pragma omp simd
		for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
	}
	void ProcessBlock(size_t n, T * in, T * out) {
		ProcessSIMD(n,in,out);
	}
	void ProcessInplace(size_t n, T * buffer) {
		ProcessSIMD(n,buffer,buffer);
	}
    };

    
	template<typename T>
    struct DiodeClipperFP : : public GSSoundProcessor<T>
    {
        T controlledR;
        T Id, C,  Ve, Vp, R, err;
        T Fs, T;
        
        T vin;
        T vout, voutTemp, voutOld;

        int upsamplingScale;

        //juce::dsp::ProcessorDuplicator< juce::dsp::IIR::Filter <T>, juce::dsp::IIR::Coefficients<T>> lowPassFilter;
        
        T gdExpDiff(T vc)
        {
            return (Id * std::exp(vc / (2.0 * Ve))) / (2.0 * Ve);
        }
        T gdExp(T vc) 
        {
            return Id * (std::exp(vc / (2 * Ve)) - 1);
        }
        T gdPoly(T vc)
        {
            return Vp * vc*vc*vc*vc * Heaviside(vc);
        }
        T Heaviside(T vc)
        {
            if (vc >= 0)
                return 1;
            else return 0;
        }
        T limiter(T val)
        {
            if (val < -1.0)
            {
                val = -1.0;
                return val;
            }
            else if (val > 1.0)
            {
                val = 1.0;
                return val;
            }
            return val;
        }
        DiodeClipperFP(T sampleRate=44100.0) : GSSoundProcessor<T>()
        {
            //oversampling->reset();
            //oversampling->initProcessing(static_cast<size_t> (samplesPerBlock));

            //juce::dsp::ProcessSpec spec;
            //spec.sampleRate = sampleRate * 16;
            //spec.maximumBlockSize = samplesPerBlock * 15;
            //spec.numChannels = getTotalNumOutputChannels();

            //lowPassFilter.prepare(spec);
            //lowPassFilter.reset();
            // Set the constants
            Fs = sampleRate;
            T = 1 / Fs;
            C = 100e-9;
            R = 1e3;
            Vp = 0.17;
            Ve = 0.023;
            Id = 2.52e-9;
            err = 0.1e-3; // err for stopping iterations
            // set controlled values to starting values (redundant maybe delit later)
            controlledR = 1.0;
            voutOld = 0;
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            int itter=0;
            voutTemp = 1;
            vout = 0;
            vin = controlledR * I;
            while (std::abs(voutTemp - vout) > err && itter++ < 20) {
                voutTemp = vout;
                vout = R * T / (R * C + T) * (gdExp(-voutTemp) - gdExp(voutTemp)) + T / (R * C + T) * vin + R * C / (R * C + T) * voutOld;
            }
            
            voutOld = vout;
            //vout = vin;
            vout = limiter(vout);
            return vout;            
        }
	void ProcessSIMD(size_t n, T * in, T * out) {
		#pragma omp simd
		for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
	}
	void ProcessBlock(size_t n, T * in, T * out) {
		ProcessSIMD(n,in,out);
	}
	void ProcessInplace(size_t n, T * buffer) {
		ProcessSIMD(n,buffer,buffer);
	}
    };
}

