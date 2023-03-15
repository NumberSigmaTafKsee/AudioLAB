#pragma once
#include "FX/LiquidNeuron.hpp"

namespace Liquid
{
    struct LiquidMoog : public FilterProcessor
    {                        
        DspFloatType f,fs,k,p,scale,r,y1,y2,y3,y4,oldx,oldy1,oldy2,oldy3;
        DspFloatType cutoff,Q;
        DspFloatType F,R,P,SCALE,K;
        DspFloatType x;
        DspFloatType envScale = 0.5f;
        int oversample = 4;
        DspFloatType low,high,band,notch;
        DspFloatType machoScale = 1.0;
        DspFloatType preGain=1.0;
        DspFloatType postGain=1.0;
        DspFloatType dcBias=0.0;
        DspFloatType cMin=-1.0;
        DspFloatType cMax=1.0;

		std::array<float,1024> A_a,X_a,Y_a;
		
        bool sat1=false;
        bool sat2=false;
        bool sat3=false;
        bool sat4=false;
        bool sigd1=false;
        bool sigd2=false;
        bool invertP = false;
                
                        
        enum {
            ONEPOLE,
            TWOPOLE,
            THREEPOLE,
            FOURPOLE,
        };

        int poles = FOURPOLE;

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_SAT1,
            PORT_SAT2,
            PORT_SAT3,
            PORT_SAT4,
            PORT_SIGD1,
            PORT_SIGD2,
            PORT_INVERTP,
            PORT_ONEPOLE,
            PORT_TWOPOLE,
            PORT_THREEPOLE,
            PORT_FOURPOLE,
            PORT_DCBIAS,
            PORT_OVERSAMPLE,
            PORT_ENVSCALE,
            PORT_PREGAIN,
            PORT_POSTGAIN,
        };
       
		enum {
			N_P,
			N_SCALE,
			N_K,
			N_OLDX,
			N_OLDY1,
			N_OLDY2,
			N_OLDY3,
			N_CUTOFF,
			N_Q,
			N_LAST
		};
			
        std::array<LiquidNeuron*,N_LAST> neurons;
        

        LiquidMoog(DspFloatType sampleRate, DspFloatType cutoff, DspFloatType resonance)
        : FilterProcessor()
        {   
			for(size_t i = 0; i < N_LAST; i++) neurons[i] = new LiquidNeuron(sampleRate, 1.0f);             
            coefficients(sampleRate,cutoff,resonance);
            x=y1=y2=y3=y4=oldx=oldy1=oldy2=oldy3=0;
            low=high=band=notch=0;
            
        }
		~LiquidMoog() {
			for(size_t i = 0; i < N_LAST; i++) {
				if(neurons[i]) delete neurons[i];
			}
		}
        void coefficients(DspFloatType sampleRate,DspFloatType frequency, DspFloatType resonance) 
        {
            F = frequency;            
            R = resonance;
            fs = sampleRate;

			/*
			neurons[N_CUTOFF]->setTarget(F);
			neurons[N_Q]->setTarget(R);
			
            cutoff  = neurons[N_CUTOFF]->Tick(); 
            Q       = neurons[N_Q]->Tick(); 
            */
            cutoff  = F;
            Q       = R;
            
            f =  cutoff / (oversample*sampleRate); //[0 - 1]
            K = 3.6*f - 1.6*f*f -1; //(Empirical tuning)
                        
            P = (k+1)*0.5;
            if(invertP) P = -P;
                        
            SCALE = std::exp((1-P)*1.386249);            
            SCALE *= machoScale;
            
        }
        void setCutoff(DspFloatType c) {     
            c = clamp(c,0,fs/2.0);
            coefficients(fs,c,Q);
        }
        void setResonance(DspFloatType res) {
            res = clamp(res,0,1);
            coefficients(fs,cutoff,res);
        }
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_CUTOFF: setCutoff(v); break;
            case PORT_RESONANCE: setResonance(v); break;
            case PORT_SAT1: sat1 = !sat1; break;
            case PORT_SAT2: sat2 = !sat2; break;
            case PORT_SAT3: sat3 = !sat3; break;
            case PORT_SAT4: sat4 = !sat4; break;
            case PORT_SIGD1: sigd1 = !sigd1; break;
            case PORT_SIGD2: sigd2 = !sigd2; break;
            case PORT_INVERTP: invertP = !invertP; break;
            case PORT_ONEPOLE: poles = ONEPOLE; break;
            case PORT_TWOPOLE: poles = TWOPOLE; break;
            case PORT_THREEPOLE: poles = THREEPOLE; break;
            case PORT_FOURPOLE: poles = FOURPOLE; break;            
            case PORT_OVERSAMPLE: oversample = (int)v; break;
            case PORT_ENVSCALE: envScale = v; break;
            case PORT_PREGAIN: preGain = v; break;
            case PORT_POSTGAIN: postGain = v; break;
            }
        }
        DspFloatType saturate(DspFloatType x)
        {
            // changing it will change the resonance but too much and it only resonates
            return 1.03 * FX::Distortion::serpent_curve(x);
        }
        DspFloatType MoogI_Tick(DspFloatType input, DspFloatType A=1, DspFloatType X=2, DspFloatType Y=1)
        {
            Undenormal denormal;
            DspFloatType c = F;
            DspFloatType res = R;
                                
            X *= envScale;
            X = clamp(X,0,1);
            
            coefficients(fs,c*X,res*Y);
            
            neurons[N_P]->setTarget(P);
            neurons[N_SCALE]->setTarget(SCALE);
            neurons[N_K]->setTarget(K);
            
            p 		= neurons[N_P]->Tick(); 
            scale 	= neurons[N_SCALE]->Tick(); 
            k 		= neurons[N_K]->Tick(); 
            r 		= Q*scale;
                      
            input *= std::pow(10.0,preGain/20.0f);
            if(sigd1) input = FX::Distortion::sigmoid_function(input);

            for(size_t i = 0; i < oversample; i++)
            {
                //--Inverted feed back for corner peaking
                x = (input+dcBias) - r*y4;
                x *= A;
                //Four cascaded onepole filters (bilinear transform)
                y1=x*p + oldx*p - k*y1;      
                
                // saturate very badly messed the reso and cutoff
                if(sat1) y1 = saturate(y1);
                y2=y1*p + oldy1*p - k*y2;                        
                if(sat2) y2 = saturate(y2);
                y3=y2*p + oldy2*p - k*y3;                        
                if(sat3) y3 = saturate(y3);
                y4=y3*p + oldy3*p - k*y4;                        
                if(sat4) y4 = saturate(y4);

                // super but resonance needs to be tamed
                if(sigd2) y4 = FX::Distortion::serpent_curve(A*y4);
                                       
				neurons[N_OLDX]->setTarget(x);
				neurons[N_OLDY1]->setTarget(y1);
				neurons[N_OLDY2]->setTarget(y2);
				neurons[N_OLDY3]->setTarget(y3);
                
                oldx  = neurons[N_OLDX]->Tick(); 
                oldy1 = neurons[N_OLDY1]->Tick(); 
                oldy2 = neurons[N_OLDY2]->Tick(); 
                oldy3 = neurons[N_OLDY3]->Tick(); 
            }

            coefficients(fs,c,res);
            
            y4 -= dcBias;            
            DspFloatType out = y4;
            if(out < cMin) out = cMin;
            if(out > cMax) out = cMax;

            if(poles == ONEPOLE) out = y1;
            else if(poles == TWOPOLE) out = y2;
            else if(poles == THREEPOLE) out = y3;
            else if(poles == FOURPOLE) out = y4;
            
            out *= std::pow(10.0,postGain/20.0);
            if(sigd2) out = FX::Distortion::serpent_curve(out);
            return out;
        }

        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return MoogI_Tick(I,A,X,Y);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {            
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DspFloatType * in) {
            ProcessSIMD(n,in,in);
        }
    };
}
