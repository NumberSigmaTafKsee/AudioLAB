#pragma once
#include "StdSamples/stdsamples.hpp"

namespace Filters::FIR
{
    struct FIR
    {
        AudioDSP::sample_vector<DspFloatType> taps;
        AudioDSP::sample_vector<DspFloatType> sr;
        DspFloatType fs,f,fu;
        int    filter_type;
        DspFloatType lambda,phi;

        enum {
            Lowpass,
            Highpass,
            Bandpass,
            Bandstop,
            Peak,
            Notch,
            LowShelf,
            HighShelf,
        };
        FIR(size_t n, int type, DspFloatType sampleRate, DspFloatType fc, DspFloatType _fu = 0, DspFloatType g=1)
        {
            fs = sampleRate;
            f  = fc;
            fu = _fu;
            filter_type = type;
            lambda = M_PI * fc / (fs/2.0);        
            phi = M_PI * fu * (fs/2.0);
            taps.resize(n);
            sr.resize(n);
            memset(sr.data(),0,n*sizeof(DspFloatType));
            switch(type)
            {
                case Lowpass: DesignLowPass(); break;
                case Highpass: DesignHighPass(); break;
                case Bandpass: DesignBandPass(); break;
                case Bandstop: DesignBandStop(); break;
                case Peak: DesignPeak(g); break;
                case Notch: DesignNotch(g); break;
                case LowShelf: DesignLowShelf(g); break;
                case HighShelf: DesignHighShelf(g); break;
            }
        }
        void DesignLowPass()
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = lambda / M_PI;
                else t[n] = std::sin( mm * lambda ) / (mm * M_PI);
            }
        }
        void DesignHighPass()
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = 1.0 - lambda / M_PI;
                else t[n] = -std::sin( mm * lambda ) / (mm * M_PI);
            }    
        }
        void DesignBandPass()
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = (    std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);
            }
        }
        void DesignBandStop()
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = -(   std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);
            }
        }
        void DesignPeak(DspFloatType g)
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                DspFloatType bs = -(   std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);
                DspFloatType bp = (   std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);

                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = g*bp + bs;
            }
        }
        void DesignNotch(DspFloatType g)
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                DspFloatType bs = -(   std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);
                DspFloatType bp = (   std::sin( mm * phi ) -
                                    std::sin( mm * lambda )   ) / (mm * M_PI);

                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = bp + g*bs;
            }
        }
        void DesignLowShelf(DspFloatType g)
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = g*(std::sin( mm * lambda ) / (mm * M_PI)) - (std::sin( mm * lambda ) / (mm * M_PI)); 
            }
        }
        void DesignHighShelf(DspFloatType g)
        {
            int n;
            DspFloatType mm;
            DspFloatType * t = this->taps.data();
			#pragma omp simd aligned(t)
            for(n = 0; n < taps.size(); n++){
                mm = n - (taps.size() - 1.0) / 2.0;
                if( mm == 0.0 ) t[n] = (phi - lambda) / M_PI;
                else t[n] = (std::sin( mm * lambda ) / (mm * M_PI)) - g*(std::sin( mm * lambda ) / (mm * M_PI)); 
            }
        }
        DspFloatType Tick(DspFloatType data_sample)
        {
            int i;
            DspFloatType result;
            Undenormal denormals;
			DspFloatType * s = this->sr.data();
			DspFloatType * t = this->taps.data();
            #pragma omp simd aligned(s)
            for(i = taps.size() - 1; i >= 1; i--){
                s[i] = s[i-1];
            }	
            s[0] = data_sample;

            result = 0;            
            #pragma omp simd aligned(s,t)
            for(i = 0; i < taps.size(); i++) result += sr[i] * taps[i];

            return result;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			Undenormal denormals;
			DspFloatType * _s = sr.data();
			DspFloatType * _t = taps.data();
			#pragma omp simd aligned(in,out)
			for(size_t s = 0; s < n; s++)
			{
				int i;
				DspFloatType result;				

				
				const DspFloatType data_sample = in[s];
				for(i = taps.size() - 1; i >= 1; i--){
					_s[i] = _s[i-1];
				}	
				_s[0] = data_sample;

				result = 0;            
				
				for(i = 0; i < taps.size(); i++) result += _s[i] * _t[i];

				out[s] = result;
			}
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType  *out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessSIMD(n,out,out);
		}
        DspFloatType Tick2x(DspFloatType x) {
            Tick(x);
            return Tick(0);
        }
    };
}
