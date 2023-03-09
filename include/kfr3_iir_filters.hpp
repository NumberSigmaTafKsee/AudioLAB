#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <random> 

#include <Undenormal.hpp>

#include <kfr/kfr.h>
#include <kfr/dft.hpp>
#include <kfr/io.hpp>
#include <kfr/math.hpp>
#include <kfr/dsp.hpp>


namespace kfr3
{
    typedef float DspFloatType;

    struct FilterCoefficients
    {
        DspFloatType a[2];
        DspFloatType b[3];
    };

    struct BiquadSection
    {
        DspFloatType z[3];
        DspFloatType p[3];

        BiquadSection()
        {
            memset(z, 0, sizeof(z));
            memset(p, 0, sizeof(p));
        }
        BiquadSection(const FilterCoefficients &c)
        {
            z[0] = c.b[0];
            z[1] = c.b[1];
            z[2] = c.b[2];
            p[0] = c.a[0];
            p[1] = c.a[1];
        }
        BiquadSection(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
        {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[1] = p2;
        }
        BiquadSection(const BiquadSection &b)
        {
            memcpy(z, b.z, sizeof(z));
            memcpy(p, b.p, sizeof(p));
        }
        void setCoefficients(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
        {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[1] = p2;
        }
        void setCoefficients(DspFloatType n[3], DspFloatType d[2])
        {
            memcpy(z, n, sizeof(z));
            memcpy(p, d, sizeof(p));
        }
        void setCoefficients(const FilterCoefficients &c)
        {
            z[0] = c.b[0];
            z[1] = c.b[1];
            z[2] = c.b[2];
            p[0] = c.a[0];
            p[1] = c.a[1];
        }
        BiquadSection &operator=(const BiquadSection &b)
        {
            memcpy(z, b.z, sizeof(z));
            memcpy(p, b.p, sizeof(p));
            return *this;
        }

        void print()
        {
            std::cout << z[0] << " + " << z[1] << " z^-1 + " << z[2] << " z^-1\n";
            std::cout << "-------------------------------------------------------------\n";
            std::cout << " 1 + " << p[0] << +" z^-1 + " << p[1] << " z^-2\n";
        }
    };


    using BiquadSOS = std::vector<BiquadSection>;
    
    struct BiquadTransposedTypeII 
    {
        BiquadSection biquad;
        DspFloatType x, y, d1, d2;

        BiquadTransposedTypeII() 
        {
            x = y = 0;
            d1 = d2 = 0;
        }
        BiquadTransposedTypeII(const BiquadSection &b) : biquad(b)
        {
            x = y = 0;
            d1 = d2 = 0;
        }
        BiquadTransposedTypeII &operator=(const BiquadTransposedTypeII &b)
        {
            biquad = b.biquad;
            x = b.x;
            y = b.y;
            d1 = b.d1;
            d2 = b.d2;
            return *this;
        }
        void setCoefficients(const BiquadSection &b)
        {
            biquad = b;
        }        
        void setBiquad(const BiquadSection &b)
        {
            biquad = b;
        }

        // transposed is just flip - to +
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
        {
            Undenormal denormal;
            x  = I;
            y  = biquad.z[0] * x + d1;
            d1 = biquad.z[1] * x - biquad.p[0] * y + d2;
            d2 = biquad.z[2] * x - biquad.p[1] * y;
            return A * y;
        }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
	{
		Undenormal denormal;
		#pragma omp simd aligned(in,out)
		for(size_t i = 0; i < n; i++) {
			x  = in[i];
			y  = biquad.z[0] * x + d1;
			d1 = biquad.z[1] * x - biquad.p[0] * y + d2;
			d2 = biquad.z[2] * x - biquad.p[1] * y;
			out[i] = y;
		}
	}
	void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
		ProcessSIMD(n,in,out);
	}
	void ProcessInplace(size_t n, DspFloatType * buffer) {
		ProcessSIMD(n,buffer,buffer);
	}
    };

    struct BiquadFilterBase
    {        
        std::vector<BiquadTransposedTypeII> biquads;
        size_t Order;
        DspFloatType Fc,Fs,Q,rQ,G,ripple,bandstop,Fu,Fl;

        BiquadFilterBase() 
        {
        }
        BiquadFilterBase(const BiquadSOS &s)
        {
            setCoefficients(s);
        }
        void setCoefficients(const BiquadSOS &s)
        {            
            biquads.resize(s.size());
            for (size_t i = 0; i < s.size(); i++)
            {
                biquads[i].setCoefficients(s[i]);
            }
        }
        virtual void setCutoff(DspFloatType f) = 0;
        virtual void setQ(DspFloatType q) = 0;

        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 1, DspFloatType Y = 1)
        {
            DspFloatType o = biquads[0].Tick(I, A, X, Y);            
            DspFloatType f = Fc;
            DspFloatType q = rQ;

            setCutoff(f * fabs(X));
            setQ(q * fabs(Y));            
            for (size_t i = 1; i < biquads.size(); i++)
                o = biquads[i].Tick(o, A, X, Y);
            setCutoff(f);
            setQ(q);                
            return A * o;
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++)
                out[i] = Tick(in[i]);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A, DspFloatType * X, DspFloatType * Y) {
            for(size_t i = 0; i < n; i++)
                out[i] = Tick(in[i],A[i],X[i],Y[i]);
        }
    };

    struct BiquadLowpassFilter : BiquadFilterBase
    {        
        BiquadLowpassFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_lowpass(Fc/Fs,Q);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }                
    };


    struct BiquadAllpassFilter : BiquadFilterBase
    {        
        BiquadAllpassFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_allpass(Fc/Fs,Q);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }                
    };

    struct BiquadHighpassFilter : BiquadFilterBase
    {        
        BiquadHighpassFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_highpass(Fc/Fs,Q);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }                
    };


    struct BiquadBandpassFilter : BiquadFilterBase
    {        
        BiquadBandpassFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_bandpass(Fc/Fs,Q);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }                
    };

    struct BiquadNotchFilter : BiquadFilterBase
    {        
        BiquadNotchFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_notch(Fc/Fs,Q);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }                
    };

    struct BiquadPeakFilter : BiquadFilterBase
    {        
        BiquadPeakFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType q, DspFloatType gain)
        {             
            Order = order;
            Fs = fs;
            rQ = q;
            Q = pow(q,order);
            G = gain;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_peak(Fc/Fs,Q,G);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }              
        void setGain(DspFloatType gain) {
            G = gain;
            setCutoff(Fc);
        }  
    };

    struct BiquadLowshelfFilter : BiquadFilterBase
    {        
        BiquadLowshelfFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType gain)
        {             
            Order = order;
            Fs = fs;                        
            G = gain;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_lowshelf(Fc/Fs,G);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }              
        void setGain(DspFloatType gain) {
            G = gain;
            setCutoff(Fc);
        }  
    };

    struct BiquadHighshelfFilter : BiquadFilterBase
    {        
        BiquadHighshelfFilter(size_t order, DspFloatType fc, DspFloatType fs, DspFloatType gain)
        {             
            Order = order;
            Fs = fs;            
            G = gain;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            for(size_t i = 0; i < Order/2; i++)
            {
                kfr::biquad_params<DspFloatType> bq = kfr::biquad_highshelf(Fc/Fs,G);
                BiquadSection sec;
                sec.z[0] = bq.b0;
                sec.z[1] = bq.b1;
                sec.z[2] = bq.b2;
                sec.p[0] = bq.a1;
                sec.p[1] = bq.a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }
        void setQ(DspFloatType q) {
            rQ = q;
            Q = pow(q,Order);
            setCutoff(Fc);
        }              
        void setGain(DspFloatType gain) {
            G = gain;
            setCutoff(Fc);
        }  
    };


    struct BesselLowPassFilter : BiquadFilterBase
    {        
        BesselLowPassFilter(size_t order, DspFloatType fc, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_lowpass(kfr::bessel<DspFloatType>(Order),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct BesselHighPassFilter : BiquadFilterBase
    {        
        BesselHighPassFilter(size_t order, DspFloatType fc, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_highpass(kfr::bessel<DspFloatType>(Order),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct BesselBandPassFilter : BiquadFilterBase
    {        
        BesselBandPassFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandpass(kfr::bessel<DspFloatType>(Order),fl,fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct BesselBandStopFilter : BiquadFilterBase
    {        
        BesselBandStopFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandstop(kfr::butterworth<DspFloatType>(Order),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };


    struct ButterworthLowPassFilter : BiquadFilterBase
    {        
        ButterworthLowPassFilter(size_t order, DspFloatType fc, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_lowpass(kfr::butterworth<DspFloatType>(Order),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ButterworthHighPassFilter : BiquadFilterBase
    {        
        ButterworthHighPassFilter(size_t order, DspFloatType fc, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_highpass(kfr::butterworth<DspFloatType>(Order),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ButterworthBandPassFilter : BiquadFilterBase
    {        
        ButterworthBandPassFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandpass(kfr::butterworth<DspFloatType>(Order),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ButterworthBandStopFilter : BiquadFilterBase
    {        
        ButterworthBandStopFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandstop(kfr::butterworth<DspFloatType>(Order),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevILowPassFilter : BiquadFilterBase
    {                
        ChebyshevILowPassFilter(size_t order, DspFloatType fc, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_lowpass(kfr::chebyshev1<DspFloatType>(Order,ripple),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIHighPassFilter : BiquadFilterBase
    {        
        ChebyshevIHighPassFilter(size_t order, DspFloatType fc, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_highpass(kfr::chebyshev1<DspFloatType>(Order,ripple),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIBandPassFilter : BiquadFilterBase
    {        
        ChebyshevIBandPassFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandpass(kfr::chebyshev1<DspFloatType>(Order,ripple),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIBandStopFilter : BiquadFilterBase
    {        
        ChebyshevIBandStopFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl,DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandstop(kfr::chebyshev1<DspFloatType>(Order,ripple),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };


    struct ChebyshevIILowPassFilter : BiquadFilterBase
    {                
        ChebyshevIILowPassFilter(size_t order, DspFloatType fc, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_lowpass(kfr::chebyshev1<DspFloatType>(Order,ripple),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIIHighPassFilter : BiquadFilterBase
    {        
        ChebyshevIIHighPassFilter(size_t order, DspFloatType fc, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fc);            
        }
        void setCutoff(DspFloatType fc) {
            if(fc <= 0 || fc >= Fs/2.0) return;
            Fc = fc;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_highpass(kfr::chebyshev2<DspFloatType>(Order,ripple),Fc,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIIBandPassFilter : BiquadFilterBase
    {        
        ChebyshevIIBandPassFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl, DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandpass(kfr::chebyshev2<DspFloatType>(Order,ripple),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };

    struct ChebyshevIIBandStopFilter : BiquadFilterBase
    {        
        ChebyshevIIBandStopFilter(size_t order, DspFloatType fl, DspFloatType fu, DspFloatType rippleDb, DspFloatType fs)
        {             
            Order = order;
            Fs = fs;                        
            ripple = rippleDb;
            setCutoff(fl,fu);            
        }
        void setCutoff(DspFloatType fl,DspFloatType fu) {
            Fu = fu;
            Fl = fl;
            BiquadSOS sos;
            kfr::zpk<DspFloatType> zpk = kfr::iir_bandstop(kfr::chebyshev2<DspFloatType>(Order,ripple),Fl,Fu,Fs);
            std::vector<kfr::biquad_params<DspFloatType>> bqs = kfr::to_sos(zpk);
            for(auto i = bqs.begin(); i != bqs.end(); i++) {
                BiquadSection sec;
                sec.z[0] = i->b0;
                sec.z[1] = i->b1;
                sec.z[2] = i->b2;
                sec.p[0] = i->a1;
                sec.p[1] = i->a2;
                sec.p[2] = 0;
                sos.push_back(sec);
            }
            setCoefficients(sos);
        }        
    };
    
    template<typename T>
    void filter(BiquadFilterBase * filter, size_t n, DspFloatType * in, DspFloatType * out) {
		filter->ProcessSIMD(n,in,out);
	}
}
