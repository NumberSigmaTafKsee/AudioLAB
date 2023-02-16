#pragma once

#include <cmath>
#include "SoundObject.hpp"

#include "FX/Amplifiers.hpp"
#include "FX/Amplifier.hpp"
#include "FX/CSmoothFilters.hpp"
#include "FX/WDFDiodeClipper.hpp"
#include "FX/DiodeSimulator.hpp"
#include "FX/Diode.hpp"

namespace Analog
{
    // Amplifier = rtspice, wavy/rt-wdf, smartguitar, wavenet
    // Waveshapers
    // Clippers
    enum {
        RTSPICE,        
        RTWDF,          
        FXWDF,          // Will Pirkins FXObjects
        SMARTGUITAR,    // Neutron Network https://github.com/GuitarML/SmartGuitarAmp
        SMARTPEDAL,     // https://github.com/GuitarML/SmartGuitarPedal
        PEDALNET,       // https://github.com/GuitarML/pedalnet
        PEDALNETRT,     // https://github.com/GuitarML/PedalNetRT
        GUITARAMP,      // AudioTK        
        TOOBAMP,        // Neutron Network https://github.com/rerdavies/ToobAmp
        WAVENET,        // https://github.com/damskaggep/WaveNetVA
        WAVESHAPER,     // y = function<foo> 
    };

    struct VCA1 : AmplifierProcessor
    {
    private:
        WDFDiodeClipper clip;
        Random noise;
    public:
        DspFloatType gain;
        FX::Filters::Smoothers::CSmoothFilter gainSmooth,vtSmooth,etaSmooth,isSmooth;
        DspFloatType Vt = 0.0253;
        DspFloatType Eta = 1.68;
        DspFloatType Is = .105;    
        DspFloatType clipMin=-1.0,clipMax=1.0;
        DspFloatType dcBias = 0.0;
        uint32_t     oversample=4;
        

        VCA1(DspFloatType g = 1.0, DspFloatType sampleRate=44100.0) 
        :   AmplifierProcessor(),
            gainSmooth(sampleRate*oversample,1/0.1),
            vtSmooth(sampleRate*oversample,1/0.1),
            etaSmooth(sampleRate*oversample,1/0.1),
            isSmooth(sampleRate*oversample,1/0.1) 
        {
            gain = pow(10,g/20.0);
            clip.prepare(sampleRate*oversample);
            clip.setCircuitParams(1000.0f);
        }

        enum {
            PORT_GAIN,
            PORT_VT,
            PORT_ETA,
            PORT_IS,
            PORT_DCBIAS,
            PORT_RESET,
            PORT_CIRCUIT,
            PORT_RANDOMIZE,
        };
        void setPort(int port, DspFloatType v)
        {
            switch(port)
            {
                case PORT_GAIN: gain = pow(10.0,v/20.0); break;
                case PORT_VT: Vt = v; break;
                case PORT_ETA: Eta = v; break;
                case PORT_IS: Is = v; break;
                case PORT_DCBIAS: dcBias = v; break;
                case PORT_RESET: Reset(); break;
                case PORT_CIRCUIT: clip.setCircuitParams(v); break;
                case PORT_RANDOMIZE: randomize(); break;
            }
        }
        void Reset()
        {
            Vt = 0.0253;
            Eta = 1.68;
            Is = 0.105;
        }
        void Randomize()
        {
            int n = noise.randint(-1,1);
            Vt += n <= 0? -1.0:1.0f*noise.frand()*0.001;
            if(Vt < 0.01) Vt = 0.01;
            if(Vt > 0.05) Vt = 0.05;
            
            n = noise.randint(-1,1);
            Eta += n <= 0? -1.0:1.0f*noise.frand()*0.01;
            if(Eta < 1.5) Eta = 1.5;
            if(Eta > 1.7) Eta = 1.7;
            
            n = noise.randint(-1,1);        
            Is += n <= 0? -1.0:1.0f*noise.frand()*0.005;
            if(Is < 0.05) Is = 0.05;
            if(Is > 0.2) Is = 0.2;
            
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = -1, DspFloatType Y = 1)
        {       
            DspFloatType r = gain*(I+dcBias);                 
            if(r < X) r = X;
            if(r > Y) r = Y;            
            r = clip.processSample(A*r) - dcBias;            
            return clamp(r,-1.0,1.0);
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {            
            FX::Distortion::amp_vector(n,gain,in,out);
            if(dcBias != 0) FX::Distortion::bias_vector(n,dcBias,out);
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType I = out[i];
                clip.processSample(I);                
                out[i] = I;            
            }
            if(dcBias != 0) FX::Distortion::bias_vector(n,-dcBias,out);
            FX::Distortion::clamp_vector(n,out,clipMin,clipMax);            
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DspFloatType * in) {
            ProcessSIMD(n,in,in);
        }
    };

    struct VCA
    {
        FX::Distortion::Amplifier::Amplifier A,B;
        double morph = 0.5;
        double gain  = 1;
        double pregain = 1;
        VCA(double G = 0.0, double P = 0.0) {
            setGain(G);
            setPreGain(P);
        }
        void setGain(double G) {
            gain = std::pow(10.0,G/20.0);
        }
        void setPreGain(double P) {
            pregain = std::pow(10.0,P/20.0);
        }
        DspFloatType Tick(DspFloatType I,DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType x = this->A.Tick(I,A,X,Y);
            return x + morph*(B.Tick(I) - x);
        }
        void Randomize() {
            A.RandomClip();
            B.RandomClip();
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType x = A.Tick(in[i]);
                out[i] = x + morph*(B.Tick(in[i]) - x);
                out[i] *= gain;
            }
        }
        void ProcessInplace(size_t n, DspFloatType * in) {
            ProcessBlock(n,in,in);
        }
    };
}

