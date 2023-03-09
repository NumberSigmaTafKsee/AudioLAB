#pragma once
#include <cmath>
#include "ClipFunctions.hpp"

namespace Filters
{
    struct StateVariableFilter
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        DspFloatType cutoff,scale,fs,low,high,band,notch;
        enum  {
			LP, HP, BP, BS
		};
		int type = LP;
        StateVariableFilter(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(DspFloatType F) { cutoff = F; }
        void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
        {
            Undenormal denormal;
            DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = Distortion::tanhify(I);
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            return low;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out){
			Undenormal denormal;			
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++){
				DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
				//--beginloop
				//I = Distortion::tanhify(I);
				const DspFloatType I = in[i];
				low = low + f * band;
				high = scale * I - low - scale*band;
				band = f * high + band;
				notch = high + low;
				switch(type)
				{
					default:
					case LP: out[i] = low; break;
					case HP: out[i] = high; break;
					case BP: out[i] = band; break;
					case BS: out[i] = notch; break;
				}
			}
		}
    };
    struct StateVariableFilterLP
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        DspFloatType cutoff,scale,fs,low,high,band,notch;
            
        StateVariableFilterLP(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(DspFloatType F) { cutoff = F; }
        void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
        {
            Undenormal denormal;
            DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = Distortion::tanhify(I);
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            return low;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out){
			Undenormal denormal;			
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++){
				DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
				//--beginloop
				//I = Distortion::tanhify(I);
				const DspFloatType I = in[i];
				low = low + f * band;
				high = scale * I - low - scale*band;
				band = f * high + band;
				notch = high + low;
				out[i] = low;
			}
		}
	};
    struct StateVariableFilterHP
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        DspFloatType cutoff,scale,fs,low,high,band,notch;
            
        StateVariableFilterHP(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(DspFloatType F) { cutoff = F; }
        void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
        {
            Undenormal denormal;
            DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = Distortion::tanhify(I);
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            return low;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out){
			Undenormal denormal;			
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++){
				DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
				//--beginloop
				//I = Distortion::tanhify(I);
				const DspFloatType I = in[i];
				low = low + f * band;
				high = scale * I - low - scale*band;
				band = f * high + band;
				notch = high + low;
				out[i] = high;
			}
		}
	};
	struct StateVariableFilterBP
	{
		/*
		cutoff = cutoff freq in Hz
		fs = sampling frequency //(e.g. 44100Hz)
		f = 2 sin (pi * cutoff / fs) //[approximately]
		q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
		low = lowpass output
		high = highpass output
		band = bandpass output
		notch = notch output
		*/
		DspFloatType cutoff,scale,fs,low,high,band,notch;
			
		StateVariableFilterBP(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
			scale = Q;
			cutoff= Fc;
			fs    = Fs;
			low=high=band=notch=0;
		}
		void setCutoff(DspFloatType F) { cutoff = F; }
		void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }
		DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
		{
			Undenormal denormal;
			DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
			//--beginloop
			//I = Distortion::tanhify(I);
			low = low + f * band;
			high = scale * I - low - scale*band;
			band = f * high + band;
			notch = high + low;
			return low;
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out){
			Undenormal denormal;			
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++){
				DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
				//--beginloop
				//I = Distortion::tanhify(I);
				const DspFloatType I = in[i];
				low = low + f * band;
				high = scale * I - low - scale*band;
				band = f * high + band;
				notch = high + low;
				out[i] = band;
			}
		}
	};
	struct StateVariableFilterBS
	{
		/*
		cutoff = cutoff freq in Hz
		fs = sampling frequency //(e.g. 44100Hz)
		f = 2 sin (pi * cutoff / fs) //[approximately]
		q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
		low = lowpass output
		high = highpass output
		band = bandpass output
		notch = notch output
		*/
		DspFloatType cutoff,scale,fs,low,high,band,notch;
			
		StateVariableFilterBS(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
			scale = Q;
			cutoff= Fc;
			fs    = Fs;
			low=high=band=notch=0;
		}
		void setCutoff(DspFloatType F) { cutoff = F; }
		void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }
		DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
		{
			Undenormal denormal;
			DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
			//--beginloop
			//I = Distortion::tanhify(I);
			low = low + f * band;
			high = scale * I - low - scale*band;
			band = f * high + band;
			notch = high + low;
			return low;
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out){
			Undenormal denormal;			
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++){
				DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
				//--beginloop
				//I = Distortion::tanhify(I);
				const DspFloatType I = in[i];
				low = low + f * band;
				high = scale * I - low - scale*band;
				band = f * high + band;
				notch = high + low;
				out[i] = notch;
			}
		}
	};
}
