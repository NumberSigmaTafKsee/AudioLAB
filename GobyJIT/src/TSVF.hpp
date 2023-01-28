#pragma once

#include "SoundAlchemy.hpp"
    
namespace SoundAlchemy
{    
    template<typename T>
    struct TChamberlinSVF 
    {
    /*
    //Input/Output
        I - input sample
        L - lowpass output sample
        B - bandpass output sample
        H - highpass output sample
        N - notch output sample
        F1 - Frequency control parameter
        Q1 - Q control parameter
        D1 - delay associated with bandpass output
        D2 - delay associated with low-pass output
    */
        T L;
        T H;
        T B;
        T N;
        T Q;
        T Fs;
        T F;
        T D1,D2;

        TChamberlinSVF(T sr=44100.0f)
        {
            L = H = B = N = 0.0f;
            D1 = D2 = 0.0f;
            F = 1;
            Fs = sr;
        }
        void SetCutoff(T freq) { F = freq; }
        void SetQ(T q) { Q = 2*q; }

        T Tick(T I)
        {
            // parameters:
            T Q1 = 1/Q;
            // where Q1 goes from 2 to 0, ie Q goes from .5 to infinity

            // simple frequency tuning with error towards nyquist
            // F is the filter's center frequency, and Fs is the sampling rate
            T F1 = 2*M_PI*F/Fs;

            // ideal tuning:
            F1 = 2 * std::sin( M_PI * F / Fs );

            // algorithm
            // loop
            L = D2 + F1 * D1;
            H = I - L - Q1*D1;
            B = F1 * H + D1;
            N = H + L;

            // store delays
            D1 = B;
            D2 = L;       
            return L;
        }

        T Variabilize(T lp, T I)
        {
            // parameters:
            T Q1 = 1/Q;
            // where Q1 goes from 2 to 0, ie Q goes from .5 to infinity

            // simple frequency tuning with error towards nyquist
            // F is the filter's center frequency, and Fs is the sampling rate
            T F1 = 2*M_PI*F/Fs;

            // ideal tuning:
            F1 = 2 * std::sin( M_PI * F / Fs );

            // algorithm
            // loop
            L = lp*D2 + F1 * D1;
            H = I - L - Q1*D1;
            B = F1 * H + D1;
            N = H + L;

            // store delays
            D1 = B;
            D2 = L;       
            return L;
        }
        
        void Process(vector<T> & buffer, vector<T>  amp = {}, vector<T> cutoff = {}, vector<T> q = {})
        {
            T freq = F;
            T res  = Q;
            T A = 1.0f;
            for(size_t i = 0; i < buffer.size(); i++)
            {
                if(cutoff.size() > 0) {
                    T tc = freq*cutoff[i];            
                    SetCutoff(freq+tc);
                }
                if(q.size() > 0) {
                    T tq = res*q[i];
                    SetQ(res + tq);
                }
                if(amp.size() > 0) {
                    A = amp[i];
                }
                buffer[i] = A*Tick(buffer[i]);
            }
            SetCutoff(freq);
            SetQ(Q);        
        }

        T GetLowPassOut()  { return L; }
        T GetHighPassOut() { return H; }
        T GetBandPassOut() { return B; }
        T GetNotchOUt()    { return N; }
    };

    template<typename T>
    struct TStateVariable
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
        scale = q
        low=high=band=0;
        */
        
        T fs;
        T cutoff;
        T q;
        T low;
        T high;
        T band;
        T notch;
        

        TStateVariable(T sr, T cutoff, T q) {
            fs = sr;
            this->cutoff = clamp(cutoff,0,sr/2);
            this->q = clamp(q,0,1);        
            low=high=band=notch=0;
        }

        T Tick(T I, T A = 1, T F = 0, T Q = 0) {
            //--beginloop        
            T ft = 2 * std::sin(M_PI * (cutoff + F*cutoff)/fs);
            T scale = q + (Q*q);
            A = clamp(A,-1,1);
            F = clamp(F,-1,1);
            Q = clamp(Q,-1,1);
            low = low + ft * band;
            high = scale * I - low - q*band;
            band = ft * high + band;
            notch = high + low;
            low *= A;
            high *= A;
            band *= A;
            notch *= A;
            //--endloop
            return low;
        }        
    };
}