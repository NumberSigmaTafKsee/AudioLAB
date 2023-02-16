#pragma once


#include "SoundObject.hpp"
#include "include/Undenormal.hpp"
#include "stdsamples.hpp"
#include "stdsamples_vectorize.hpp"
#include "stdsamples_fftw.hpp"
#include "stdsamples_iir_butterworth.hpp"
#include "FX/Amplifiers.hpp"


////////////////////
// Biquad
////////////////////
struct FilterBase
{
	virtual floatType Tick(floatType I, floatType A=1, floatType X=1, floatType Y=1) = 0;
};

struct Biquad : public FilterBase
{
    floatType z[3];
    floatType p[3];
    floatType x[2];
    floatType y[2];

    Biquad() {
        x[0] = x[1] = 0;
        y[0] = y[1] = 0;    
    }
    void setCoeffs(floatType Z[3], floatType P[3]) {
        memcpy(z,Z,sizeof(z));
        memcpy(p,P,sizeof(p));
    }
    void clear() {
        memset(x,0,sizeof(x));
        memset(y,0,sizeof(y));
    }
    floatType Tick(floatType I, floatType A=1, floatType X=1, floatType Y=1)
    {
        floatType r = I*z[0] + x[0]*z[1] + x[1] * z[2] - y[0]*p[0] - y[1]*p[1];
        x[1] = x[0];
        x[0] = I;
        y[1] = y[0];
        y[0] = r;
        return r;
    }
};

struct ButterworthLowpassFilter : public FilterBase
{
    int order;
    floatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthLowpassFilter(int order, floatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(floatType f, floatType Q) {
        auto x = IIRFilters::Butterworth::butterlp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(floatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(floatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    floatType Tick(floatType I, floatType A=1, floatType X=1, floatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, floatType * in, floatType * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<floatType> impulse_response(size_t n)
    {
        std::vector<floatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
    std::vector<floatType> frequency_response(int ir_size) {    
        AudioDSP::FFTPlanRealFloat fftPlan(ir_size);
        std::vector<std::complex<floatType>> fr(ir_size);
        std::vector<floatType> r(ir_size);            
        r[0] = Tick(1.0);
        for(size_t i = 1; i < ir_size; i++)
            r[i] = Tick(0);        
        AudioDSP::fft(fftPlan,r.data(),fr.data());                
        for(size_t i = 0; i < ir_size; i++)
		{
			fr[i] /= std::complex<floatType>(ir_size,ir_size);
		}  
        for(size_t i = 0; i < ir_size; i++) r[i] = std::abs(fr[i]);        
        return r;
    }
};

struct ButterworthHighpassFilter : public FilterBase
{
    int order;
    floatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthHighpassFilter(int order, floatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(floatType f, floatType Q) {
        auto x = IIRFilters::Butterworth::butterhp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(floatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(floatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    floatType Tick(floatType I, floatType A=1, floatType X=1, floatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, floatType * in, floatType * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<floatType> impulse_response(size_t n)
    {
        std::vector<floatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
    std::vector<floatType> frequency_response(int ir_size) {    
        AudioDSP::FFTPlanRealFloat fftPlan(ir_size);
        std::vector<std::complex<floatType>> fr(ir_size);
        std::vector<floatType> r(ir_size);            
        r[0] = Tick(1.0);
        for(size_t i = 1; i < ir_size; i++)
            r[i] = Tick(0);        
        AudioDSP::fft(fftPlan,r.data(),fr.data());                
        for(size_t i = 0; i < ir_size; i++)
		{
			fr[i] /= std::complex<floatType>(ir_size,ir_size);
		}  
        for(size_t i = 0; i < ir_size; i++) r[i] = std::abs(fr[i]);        
        return r;
    }
};

void RemoveDC(std::vector<floatType> & out)
{
	// remove dc offset
	const int sz = out.size();
	AudioDSP::FFTPlanRealFloat fftPlan(sz);
	std::vector<std::complex<floatType>> fr(sz);	
	AudioDSP::fft(fftPlan,out.data(),fr.data());                
	fr[0] = fr[fr.size()-1] = std::complex<floatType>(0,0);
	for(size_t i = 0; i < fr.size(); i++) fr[i] /= std::complex<floatType>(sz,sz);
	AudioDSP::ifft(fftPlan,fr.data(),out.data());	
}

void Norm(std::vector<floatType> & out)
{
	// normalize
	auto max = *std::max_element(out.begin(),out.end());
	auto min = *std::min_element(out.begin(),out.end());
	div(out.size(),out.data(),max,out.data());
}

void RemoveDCNorm(std::vector<floatType> & out)
{
	RemoveDC(out);
	Norm(out);
}
void RemoveDCClamp(std::vector<floatType> & out)
{	
	RemoveDC(out);
	FX::Distortion::clamp_vector(out.size(),out.data(),out.data(),-1,1);
}

