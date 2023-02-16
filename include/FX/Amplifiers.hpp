// these are waveshapers
#pragma once
#include <functional>
#include <random>
#include "ClipFunctions.hpp"
#include "SoundObject.hpp"
// this was a bad idea
DspFloatType pre_gain  = 1.0;
DspFloatType post_gain = 1.0;


namespace FX::Distortion
{
    struct AmplifierFunction1 : public AmplifierProcessor
    {
        std::function<DspFloatType (DspFloatType I)> func;
        AmplifierFunction1() : AmplifierProcessor() {

        }
        AmplifierFunction1(std::function<DspFloatType (DspFloatType I)> f) : AmplifierProcessor(),func(f)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*func(I);
        }
    };
    struct AmplifierFunction2 : public AmplifierProcessor
    {
        std::function<DspFloatType (DspFloatType I, DspFloatType X)> func;
        AmplifierFunction2() : AmplifierProcessor() {
            
        }
        AmplifierFunction2(std::function<DspFloatType (DspFloatType I, DspFloatType X)> f) : AmplifierProcessor(),func(f)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*func(I,X);
        }
    };
    struct AmplifierFunction3 : public AmplifierProcessor
    {
        std::function<DspFloatType (DspFloatType I, DspFloatType X, DspFloatType Y)> func;
        AmplifierFunction3() : AmplifierProcessor() {
            
        }
        AmplifierFunction3(std::function<DspFloatType (DspFloatType I, DspFloatType X, DspFloatType Y)> f) : AmplifierProcessor(),func(f)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*func(I,X,Y);
        }
    };
    struct AmplifierFunction4 : public AmplifierProcessor
    {
        std::function<DspFloatType (DspFloatType I, DspFloatType A, DspFloatType X, DspFloatType Y)> func;
        AmplifierFunction4() : AmplifierProcessor() {
            
        }
        AmplifierFunction4(std::function<DspFloatType (DspFloatType I, DspFloatType A, DspFloatType X, DspFloatType Y)> f) : AmplifierProcessor(),func(f)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return func(I,A,X,Y);
        }
    };

	// std::clamp = annoying template
	inline DspFloatType amp_clamp(DspFloatType x, DspFloatType a, DspFloatType b) {
        return x < a? a: x > b? b : x;
    }
    inline DspFloatType udo1(DspFloatType x, DspFloatType g = 1.0)
    {
        DspFloatType ax = fabs(x*pre_gain);
        if(ax == 0) return 0;
        return amp_clamp(post_gain*((x/ax)*(1-exp(g*(x*x)/ax))),-1.0,1.0);
    }
    inline void udo1_simd(size_t n, DspFloatType * in, DspFloatType g = 1.0)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x  = in[i];
            DspFloatType ax = fabs(x*pre_gain);
            if(ax == 0) in[i] = 0;
            else in[i] = amp_clamp(post_gain*((x/ax)*(1-exp(g*(x*x)/ax))),-1.0,1.0);
        }
    }
    
    inline DspFloatType Fold(DspFloatType x)
    {
        if( x > 1) return Fold(1-(x-1));
        if( x < -1) return Fold(-1-(x+1));
        return x;
    }
    inline void fold_vector(size_t n, DspFloatType * buffer) {
        // i dont think this will do any simd but why not
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            buffer[i] = Fold(buffer[i]);
    }
    inline DspFloatType Wrap(DspFloatType x) {
        if( x > 1) return fmod(x,1)-1;
        if( x < 1) return 1-fmod(-x,1);
        return x;
    }
    inline void wrap_vector(size_t n, DspFloatType * buffer) {
        // i dont think this will do any simd but why not
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            buffer[i] = Wrap(buffer[i]);
    }
    inline DspFloatType SinFold(DspFloatType x) {
        return sin(2*M_PI*x);
    }
    inline void sinfold_vector(size_t n, DspFloatType * buffer) {
        // i dont think this will do any simd but why not
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            buffer[i] = SinFold(buffer[i]);
    }
    inline DspFloatType cheby(int n, DspFloatType x)
    {
        if(n == 0) return 1;
        if(n == 1) return x;
        return 2*x*cheby(n-1,x) - cheby(n-2,x);
    }
    inline void cheby_vector(size_t n, int o, DspFloatType * buffer) {
        // i dont think this will do any simd but why not
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            buffer[i] = cheby(o,buffer[i]);
    }
    inline DspFloatType cheby_polynomial(int n, DspFloatType x)
    {
        if(n == 0) return 1.0;
        if(n == 1) return x;
        return 2*x*cheby_polynomial(n-1,x) - cheby_polynomial(n-2,x);
    }
    inline void cheby_polynomial_vector(size_t n, int o, DspFloatType * buffer) {
        // i dont think this will do any simd but why not
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            buffer[i] = cheby_polynomial(o,buffer[i]);
    }
    
    inline void clamp_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType a, DspFloatType b) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            DspFloatType x = in[i];
            out[i] = x < a? a: x > b? b : x;
        }
    }
    inline void bias_vector(size_t n, DspFloatType bias, DspFloatType *out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {            
            out[i] = out[i] + bias;
        }
    }
    inline void clamp_vector(size_t n, DspFloatType * out, DspFloatType a, DspFloatType b) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            DspFloatType x = out[i];
            out[i] = x < a? a: x > b? b : x;
        }
    }
    inline DspFloatType preamp(DspFloatType x) {
        return pre_gain * x;
    }
    inline void preamp_vector(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            DspFloatType x = in[i] * pre_gain;
            out[i] = x;
        }
    }
    inline void amp_vector(size_t n, DspFloatType gain, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {            
            out[i] *= gain;
        }
    }
    inline void amp_vector(size_t n, DspFloatType gain, DspFloatType *in, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {            
            out[i] = in[i]*gain;
        }
    }
    inline DspFloatType postamp(DspFloatType x) {
        return post_gain * x;
    }
    inline void postamp_vector(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            DspFloatType x = in[i] * post_gain;
            out[i] = x;
        }
    }
    inline DspFloatType tanh_normal(DspFloatType x, DspFloatType K=10,DspFloatType r = 1.0) {
        return amp_clamp(post_gain*std::tanh(pre_gain*K*x) / std::tanh(r),-1.0,1.0);        
    }
    inline void tanh_normal_vector(size_t n, DspFloatType *in, DspFloatType * out, DspFloatType K=10,DspFloatType r = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            DspFloatType x = in[i]*pre_gain;
            out[i] = amp_clamp(post_gain* std::tanh(pre_gain*K*x) / std::tanh(r),-1.0,1.0);
        }
    }
    inline DspFloatType positive_signal(DspFloatType x) {
        return (x+1)/2;
    }
    inline DspFloatType negative_signal(DspFloatType x) {
        return (x-1)/2;
    }
    

    inline DspFloatType sigmoid_function(DspFloatType x, DspFloatType K=10) {
        return amp_clamp(2*post_gain*(1.0 / (1.0 + std::exp(-K*pre_gain*x))-0.5),-1.0,1.0);
    }
    inline void sigmoid_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType K=10) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];
            out[i] = amp_clamp(2*post_gain*(1.0 / (1.0 + std::exp(-K*pre_gain*x))-0.5),-1.0,1.0);
        }        
    }
    inline DspFloatType sigmoid_minus(DspFloatType x, DspFloatType K=10) {
        return -sigmoid_function(x,K);
    }
    inline void sigmoid_minus_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType K=10) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];
            out[i] = -amp_clamp(2*post_gain*(1.0 / (1.0 + std::exp(-K*pre_gain*x))-0.5),-1.0,1.0);
        }        
    }
    inline DspFloatType bpsigmoid(DspFloatType x, DspFloatType g = 10) {
        return sigmoid_function(-x,g);
    }
    inline void bpsigmoid_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType K=10) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = -in[i];
            out[i] = amp_clamp(2*post_gain*(1.0 / (1.0 + std::exp(-K*pre_gain*x))-0.5),-1.0,1.0);
        }        
    }
    inline DspFloatType full_rectify(DspFloatType x) {
        return amp_clamp(x*std::abs(x),0,1);
    }
    inline void full_rectify_vector(size_t n, DspFloatType *in) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) in[i] = in[i]*std::abs(in[i]);
        clamp_vector(n,in,0,1);        
    }
    inline DspFloatType half_rectify(DspFloatType x) {
        return amp_clamp(x,0,1);
    }
    inline void half_rectify_vector(size_t n, DspFloatType *x) {
        clamp_vector(n,x,0,1);
    }
    inline DspFloatType modulated_signals(DspFloatType a, DspFloatType b) {
        return amp_clamp(post_gain*a*b,-1.0,1.0);
    }
    

    inline DspFloatType circular_modulated_signals(DspFloatType a, DspFloatType b) {
        return amp_clamp(post_gain*std::fmod(pre_gain*a,pre_gain*b),-1.0,1.0);
    }
    inline DspFloatType positive_modulated_signals(DspFloatType a, DspFloatType b) {
        return amp_clamp(positive_signal(a)*positive_signal(b),-1.0,1.0);
    }
    inline DspFloatType negative_modulated_signals(DspFloatType a, DspFloatType b) {
        return amp_clamp(negative_signal(a)*negative_signal(b),-1.0,1.0);
    }

    inline DspFloatType sigmoidDistortionFunction(DspFloatType x, DspFloatType gain,
                                            DspFloatType max, DspFloatType dc) {                                                        
        return amp_clamp(max * gain * x / sqrt(1 + (gain * std::pow(gain * x, 2))) + dc,-1.0,1.0);
    }
    inline void sigmoid_distortion_functionvector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType gain,DspFloatType max, DspFloatType dc) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];            
            out[i] = sigmoidDistortionFunction(x,gain,max,dc);
        }        
    }
    inline DspFloatType asymmetricSigmoidDistortionFunction(DspFloatType x) 
    {
        // Cutoff for chopping top
        static DspFloatType cutoff = 0.05;
        static DspFloatType slope = 0.1;
        static DspFloatType gain = 20;
        static DspFloatType max = 0.3;
        static DspFloatType dc = 0;
        // Calculate constant to add to linear region to make it join up with the
        // sigmoid function
        static DspFloatType b = sigmoidDistortionFunction(x, gain, max, dc) - slope * cutoff;
        if (x > cutoff) {
            return slope * x + b;
        } else {
            return sigmoidDistortionFunction(x, gain, max, dc);
        }
    }

    inline DspFloatType assymetric_sigmoid(DspFloatType I, DspFloatType A = 1, DspFloatType X = -1, DspFloatType Y = 1) {        
        DspFloatType x = asymmetricSigmoidDistortionFunction(A*pre_gain*I);
        if( std::isnan(x)) {
            return (I < 0)? -1.0:1.0;
        }
        x = x < X? X : x > Y? Y : x;
        return amp_clamp(post_gain*x,-1.0,1.0);
    }
    inline void assymetric_sigmoid_vector(size_t n, DspFloatType * in, DspFloatType *out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];
            out[i] = assymetric_sigmoid(in[i]);
        }        
    }
    
    inline DspFloatType asymmetricSigmoidDistortionFunction2(DspFloatType x) {
        // Cutoff for chopping top
        static DspFloatType cutoff = 0.05;
        static DspFloatType gain = 20;
        static DspFloatType max = 0.3;
        static DspFloatType dc = 0;
        
        if (x > cutoff) {
            return sigmoidDistortionFunction(sigmoidDistortionFunction(x, gain, max, dc), gain * 2, max, dc);
        } else {
            return sigmoidDistortionFunction(x, gain, max, dc);
        }
    }
    

    inline DspFloatType assymetric_sigmoid2(DspFloatType I, DspFloatType A = 1, DspFloatType X = -1, DspFloatType Y = 1) {
        DspFloatType x = asymmetricSigmoidDistortionFunction2(A*pre_gain*I);
        return amp_clamp(post_gain*x,X,Y);
    }
    inline void assymetric_sigmoid2_vector(size_t n, DspFloatType * in, DspFloatType *out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];
            out[i] = assymetric_sigmoid2(in[i]);
        }        
    }
    inline DspFloatType distortionFunction(DspFloatType x) {
        if (x < -0.08905) {
            // Assume x >= -1
            // Therefore, this first interval is actually -1 <= x < -0.08905
            return -(3 / 4) * (1 - std::pow((1 - (-x - 0.032847)), 12) +
                            (1 / 3) * (-x - 0.032847)) +
                0.01;
        } else if (x < 0.320018) {
            return -6.153 * std::pow(x, 2) + 3.9375 * x;
        } else {
            // Assume x <= 1
            // Therefore, this last interval is actually 0.320018 <= x <= 1
            return 0.630035;
        }
    }
    
    inline DspFloatType distortion_function(DspFloatType I, DspFloatType A=1,DspFloatType X = -1, DspFloatType Y=1) 
    {
        DspFloatType x = distortionFunction(A*pre_gain*I);
        return amp_clamp(post_gain*x,X,Y);
    }
    inline void distortion_function_vector(size_t n, DspFloatType * in, DspFloatType *out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = in[i];
            out[i] = distortion_function(in[i]);
        }        
    }

    inline DspFloatType cubic_distortion(DspFloatType in, DspFloatType A = 1, DspFloatType X = -1, DspFloatType Y = 1)
    {        
        DspFloatType r = A*(in - (1.0/3.0)*std::pow(pre_gain*in,3.0));    
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void cubic_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = G*(in[i] - (1.0/3.0)*std::pow(pre_gain*in[i],3.0));                
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }        
    }
    inline DspFloatType asin_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType r = (2.f / M_PI) * std::asin(pre_gain*in * A);
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void asin_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::asin(in[i] * pre_gain* G);
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType acos_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType r = (2.f / M_PI) * std::acos(in * pre_gain * A);
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,X,Y);
    }
    inline void acos_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::acos(in[i] * pre_gain * G);
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1,1);
        }
    }
    
    inline DspFloatType atan_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType r = (2.f / M_PI) * std::atan(in * pre_gain * A);
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void atan_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::atan(in[i] * pre_gain * G);
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType asinh_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType r = (2.f / M_PI) * std::asinh(in * pre_gain * A);    
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(A*r,-1.0,1.0);    
    }
    inline void asinh_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::asinh(in[i] * pre_gain * G);    
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);    
        }
    }
    inline DspFloatType acosh_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType r = (2.f / M_PI) * std::acosh(in * pre_gain * A);
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void acosh_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::acosh(in[i] * pre_gain * G);
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType atanh_distortion(DspFloatType in, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType r = (2.f / M_PI) * std::atanh(in * pre_gain * A);
        if(std::isnan(r)) {
            return in < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void atanh_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType r = (2.f / M_PI) * std::atanh(in[i] * pre_gain * G);
            if(std::isnan(r)) {
                out[i] = in < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType exp_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType sign = x < 0? -1.0:1.0;
        DspFloatType r = sign * (1.f - std::exp(-std::fabs(x * pre_gain * A)));
        if(std::isnan(r)) {
            return x < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        return amp_clamp(post_gain*r,-1.0,1.0);
    }
    inline void exp_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            DspFloatType sign = in[i] < 0? -1.0:1.0;
            DspFloatType r = sign * (1.f - std::exp(-std::fabs(in[i] * pre_gain * G)));
            if(std::isnan(r)) {
                out[i] = r < 0? -1.0 : 1.0;
            }            
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType dc_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType dc = (DspFloatType)rand() / (DspFloatType)(RAND_MAX);
        DspFloatType r = 0;
        if(x < 0) r = cubic_distortion(x - X*dc,A,X,Y);
        else r = cubic_distortion(x+Y*dc,A,X,Y);
        return amp_clamp(post_gain*r,X,Y);
    }
    inline void dc_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {        
        DspFloatType dc = (DspFloatType)rand() / (DspFloatType)(RAND_MAX);
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {            
            DspFloatType r = 0;
            DspFloatType x = in[i];
            if(x < 0) r = cubic_distortion(x - 0.01*dc,G);
            else r = cubic_distortion(x+0.01*dc,G);
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }

    inline DspFloatType bipolar_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType r = 0;
        if(x > 0) r = atan_distortion(x,A,X,Y);
        else r = cubic_distortion(x,A,X,Y);
        return amp_clamp(post_gain*r,X,Y);
    }
    inline void bipolar_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {            
            DspFloatType r = 0;
            DspFloatType x = in[i];
            if(x > 0) r = atan_distortion(x,G);
            else r = cubic_distortion(x,G);
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }
    inline DspFloatType quadratic_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {         
        DspFloatType r = x;
        if(x >= 0 && x < M_PI/2) r = atan_distortion(x,A,X,Y);
        else if(x >= (M_PI/2) && x < M_PI) r = atan_distortion(x,A,X,Y);
        else if(x >= M_PI && x < (3*M_PI/4)) r = exp_distortion(x,A,X,Y);
        else r= exp_distortion(x,A,X,Y);
        return amp_clamp(post_gain*r,X,Y);        
    }
    inline void quadratic_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {           
            DspFloatType x = in[i],r;
            if(x >= 0 && x < M_PI/2) r = atan_distortion(x,G);
            else if(x >= (M_PI/2) && x < M_PI) r = atan_distortion(x,G);
            else if(x >= M_PI && x < (3*M_PI/4)) r = exp_distortion(x,G);
            else r= exp_distortion(x,G);
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);         
        }
    }
    inline DspFloatType quadratic2_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {            
        DspFloatType r = x;
        if(x >= 0 && x < M_PI/2) r= atan_distortion(x,A,X,Y);
        else if(x >= (M_PI/2) && x < M_PI) r= cubic_distortion(x,A,X,Y);
        else if(x >= M_PI && x < (3*M_PI/4)) r= atan_distortion(x,A,X,Y);
        else r= cubic_distortion(x,A,X,Y);
        return amp_clamp(post_gain*r,X,Y);
    }
    inline void quadratic2_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType x = in[i],r;
            if(x >= 0 && x < M_PI/2) r= atan_distortion(x,G);
            else if(x >= (M_PI/2) && x < M_PI) r= cubic_distortion(x,G);
            else if(x >= M_PI && x < (3*M_PI/4)) r= atan_distortion(x,G);
            else r= cubic_distortion(x,G);
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);        
        }
    }
    inline DspFloatType quadratic3_distortion(DspFloatType x, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {            
        DspFloatType r = x;
        if(x >= 0 && x < M_PI/2) r=cubic_distortion(x,A,X,Y);
        else if(x >= (M_PI/2) && x < M_PI) r=cubic_distortion(x,A,X,Y);
        else if(x >= M_PI && x < (3*M_PI/4)) r=atan_distortion(x,A,X,Y);
        else r= atan_distortion(x,A,X,Y);
        return amp_clamp(post_gain*r,X,Y);
    }
    inline void quadratic3_distortion_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType r,x = in[i];
            if(x >= 0 && x < M_PI/2) r=cubic_distortion(x,G);
            else if(x >= (M_PI/2) && x < M_PI) r=cubic_distortion(x,G);
            else if(x >= M_PI && x < (3*M_PI/4)) r=atan_distortion(x,G);
            else r= atan_distortion(x,G);
            out[i] = amp_clamp(post_gain*r,-1.0,1.0);        
        }
    }

    inline DspFloatType parametric_clip(DspFloatType input, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType softClip = (2.f / M_PI) * std::atan(input * A);
        DspFloatType blend = input * (X*0.5 + softClip * Y*0.5);
        return amp_clamp(post_gain*blend,X,Y);
    }

    inline void parametric_clip_vector(size_t n, DspFloatType * in, DspFloatType *out, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType softClip = (2.f / M_PI) * std::atan(in[i] * G);            
            out[i] = amp_clamp(post_gain*softClip,-1.0,1.0);
        }
    }
    inline DspFloatType arcTanDistortion (DspFloatType input, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y = 1)
    {     
        DspFloatType gain = pre_gain*(A + 1.0);    
        DspFloatType out = 2.0 / M_PI * std::atan(gain *  input);    
        out = out / std::log(gain);
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,X,Y);
    }
    inline void arctandistortion_vector(size_t n, DspFloatType * in, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType gain = pre_gain*(G + 1.0);    
            DspFloatType out = 2.0 / M_PI * std::atan(gain * in[i]);    
            out = out / std::log(gain);
            if(std::isnan(out)) {
                output[i] =  in[i] < 0? -1.0 : 1.0;
            }            
            else output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }
    inline DspFloatType softClipper (DspFloatType input, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {     
        DspFloatType newInput = pre_gain*(input * A);
        DspFloatType out = 0.0;
        
        if (newInput >= 1.0)
            out = 1.0;
        else if ((newInput > -1) && (newInput < 1))
            out = (3.0 / 2.0) * (newInput - (std::pow(newInput, 3.0) / 3.0));
        else if (newInput <= -1)
            out = -1.0;
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,X,Y);    
    }
    inline void softclipper_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType newInput = pre_gain*(input[i] * G);
            DspFloatType out = 0.0;            
            if (newInput >= 1.0)
                out = 1.0;
            else if ((newInput > -1) && (newInput < 1))
                out = (3.0 / 2.0) * (newInput - (std::pow(newInput, 3.0) / 3.0));
            else if (newInput <= -1)
                out = -1.0;            
            output[i] = amp_clamp(post_gain*out,-1.0,1.0);    
        }        
    }

    inline DspFloatType errorf(DspFloatType x, DspFloatType K = 10, DspFloatType X =-1, DspFloatType Y=1) {
        DspFloatType r = std::erf(pre_gain*K*x);
        if(std::isnan(r)) {
            return x < 0? -1.0 : 1.0;
        }
        r = r < X? X : r > Y? Y : r;
        r *= post_gain;
        return amp_clamp(r,X,Y);
    }

    inline void errorf_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType r = std::erf(pre_gain*G*input[i]);
            if(std::isnan(r)) {
                output[i] = r < 0? -1.0 : 1.0;
            }                
            else output[i] = amp_clamp(post_gain*r,-1.0,1.0);
        }
    }

    inline DspFloatType sigmoided (DspFloatType input, DspFloatType A=1, DspFloatType X = -1, DspFloatType Y=1)
    {        
        DspFloatType gain = pre_gain * (A + 1.0);    
        DspFloatType out = (2.0 * (1.0 / (1.0 + std::exp(-gain * input)))) - 1;    
        out = (out) / (std::log(gain));
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }


    inline DspFloatType hardclip(DspFloatType input, DspFloatType A=1, DspFloatType X=-1, DspFloatType Y=1) {
        input *= A*pre_gain;
        input = input < X? X : input > Y? Y : input;
        return amp_clamp(post_gain*input,-1.0,1.0);
    }

    inline DspFloatType hyperbolicTangent (DspFloatType input, DspFloatType gain=1, DspFloatType X = -1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        DspFloatType out = (std::tanh(gain * pre_gain * input)) / (std::tanh(gain));
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }

    inline DspFloatType diodeClipping (DspFloatType input, DspFloatType gain=1, DspFloatType X=-1, DspFloatType Y=1)
    {
        DspFloatType diodeClippingAlgorithm = std::exp((0.1 * pre_gain * input) / (0.0253 * 1.68)) - 1.0;
        DspFloatType out = 2 / M_PI * std::atan(diodeClippingAlgorithm * (gain * 16));
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    inline void diode_clipping_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {   
            DspFloatType diodeClippingAlgorithm = std::exp((0.1 * pre_gain * input[i]) / (0.0253 * 1.68)) - 1.0;
            DspFloatType out = 2 / M_PI * std::atan(diodeClippingAlgorithm * (G * 16));
            if(std::isnan(out)) {
                out = input < 0? -1.0 : 1.0;
            }            
            output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    inline DspFloatType fuzzExponential (DspFloatType input, DspFloatType gain=1, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain *= pre_gain;
        DspFloatType newInput = input * gain;
        DspFloatType out;
        
        //Soft clipping
        if (newInput < 0.0)
            out = -1.0 *  (1.0 - std::exp(-abs(newInput)));
        else
            out = 1.0 * (1.0 - std::exp(-abs(newInput)));
    
        //Half Wave Rectifier
        out = 0.5 * (out + abs(out));
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    inline void fuzz_exponential_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType newInput = G * pre_gain * input[i];
            DspFloatType out;
        
            //Soft clipping
            if (newInput < 0.0)
                out = -1.0 *  (1.0 - std::exp(-abs(newInput)));
            else
                out = 1.0 * (1.0 - std::exp(-abs(newInput)));
        
            //Half Wave Rectifier
            out = 0.5 * (out + abs(out));
            if(std::isnan(out)) 
               output[i] = input < 0? -1.0 : 1.0;
            else
                output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    inline DspFloatType pieceWiseOverdrive (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        DspFloatType newInput = pre_gain*input * (gain) ;
        DspFloatType out = 0.0;
        
        if (abs(newInput) <= 1.0 / 3.0)
            out = 2.0 * newInput;
        else if (abs(newInput) > 2.0 / 3.0)
        {
            if (newInput > 0.0)
                out = newInput;
            if (newInput < 0.0)
                out = -newInput;
        } else
        {
            if (newInput > 0.0)
                out = (3.0 - std::pow((2.0 - newInput * 3.0), 2.0)) / 3.0;
            if (newInput < 0.0)
                out = -(3.0 - std::pow((2.0 - newInput * 3.0), 2.0)) / 3.0;
        }
        
        out = (out / std::log(gain + 1.0));
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);    
    }
    inline void piecewise_overdrive_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType newInput = pre_gain*input[i]*G;
            DspFloatType out = 0.0;
            
            if (abs(newInput) <= 1.0 / 3.0)
                out = 2.0 * newInput;
            else if (abs(newInput) > 2.0 / 3.0)
            {
                if (newInput > 0.0)
                    out = newInput;
                if (newInput < 0.0)
                    out = -newInput;
            } else
            {
                if (newInput > 0.0)
                    out = (3.0 - std::pow((2.0 - newInput * 3.0), 2.0)) / 3.0;
                if (newInput < 0.0)
                    out = -(3.0 - std::pow((2.0 - newInput * 3.0), 2.0)) / 3.0;
            }            
            out = (out / std::log(G+ 1.0));            
            if(std::isnan(out))
                output[i] = input < 0? -1.0 : 1.0;
            else
                output[i] =  amp_clamp(post_gain*out,-1.0,1.0);    
        }
    }

    inline DspFloatType tube (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        DspFloatType Q = -1.5; //more negative = more linear
        DspFloatType distortion = 5; //higher number = higher distortion
        DspFloatType out;
        
        DspFloatType newInput = pre_gain * input * (gain / 10);
        
        if (Q == 0) {
            out = newInput / (1 - std::exp(-distortion * newInput));
            if (newInput == Q)
            {
                out = 1 / distortion;
            }
        } else {            
            if (newInput == Q)            
                out = (1 / distortion) + (Q / (1 - std::exp(distortion * Q)));
            else
                out = ((newInput - Q) / (1 - std::exp(-distortion * (newInput - Q)))) + (Q / (1 - std::exp(distortion * Q)));
        }
        
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    inline void tube_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType gain = (G + 1.0);
            DspFloatType Q = -1.5; //more negative = more linear
            DspFloatType distortion = 5; //higher number = higher distortion
            DspFloatType out;
            
            DspFloatType newInput = pre_gain * input[i] * (gain / 10);
            
            if (Q == 0) {                
                if (newInput == Q)                
                    out = 1 / distortion;                
                else
                    out = newInput / (1 - std::exp(-distortion * newInput));
            } else {            
                if (newInput == Q)            
                    out = (1 / distortion) + (Q / (1 - std::exp(distortion * Q)));
                else
                    out = ((newInput - Q) / (1 - std::exp(-distortion * (newInput - Q)))) + (Q / (1 - std::exp(distortion * Q)));
            }
            
            if(std::isnan(out)) {
                output[i] =  input[i] < 0? -1.0 : 1.0;
            }
            else output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    inline DspFloatType arraya (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        auto newInput = pre_gain * input;
        
        //Arraya

        auto out = ((3.0 * newInput) / 2.0) * (1.0 - (std::pow(newInput, 2.0) / 3.0));
        
    //    Fuzz Exponential
        if (out < 0.0)
            out = 1.0 * ((1.0 - std::exp(abs(out)) / (std::exp(1.0) - 1.0)));
        else
            out = -1.0 * ((1.0 - std::exp(abs(out)) / (std::exp(1.0) - 1.0)));
        
        //Exponential 2
    //    out = (std::exp(1.0) - std::exp(1.0 - out)) / (std::exp(1.0) - 1.0);
        
    //    out = 0.5 * (out + abs(out));
    //    out = abs(out);
        
        if (gain >= 10.0)
            out = out * (gain / 100.0);
        else
            out = out * (0.1);
        
        //Arraya
        out = ((3.0 * out) / 2.0) * (1.0 - (std::pow(out, 2.0) / 3.0));
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    inline void arraya_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType gain = (G + 1.0);           
            auto newInput = pre_gain * input[i] * gain;
        
        //Arraya
            auto out = ((3.0 * newInput) / 2.0) * (1.0 - (std::pow(newInput, 2.0) / 3.0));
            
        //    Fuzz Exponential
            if (out < 0.0)
                out = 1.0 * ((1.0 - std::exp(abs(out)) / (std::exp(1.0) - 1.0)));
            else
                out = -1.0 * ((1.0 - std::exp(abs(out)) / (std::exp(1.0) - 1.0)));
            
            //Exponential 2
        //    out = (std::exp(1.0) - std::exp(1.0 - out)) / (std::exp(1.0) - 1.0);
            
        //    out = 0.5 * (out + abs(out));
        //    out = abs(out);
            
            if (gain >= 10.0)
                out = out * (gain / 100.0);
            else
                out = out * (0.1);
            
            //Arraya
            out = ((3.0 * out) / 2.0) * (1.0 - (std::pow(out, 2.0) / 3.0));
            if(std::isnan(out)) {
                output[i] = input < 0? -1.0 : 1.0;
            }
            else 
                output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    DspFloatType gallo (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        DspFloatType a = -0.01;
        DspFloatType b = 0.7;
        DspFloatType k1 = std::pow(a, 2.0);
        DspFloatType k2 = 1 + (2 * a);
        DspFloatType k3 = std::pow(b, 2.0);
        DspFloatType k4 = 1 - (2 * b);
        DspFloatType out_1 = 0.0;
        
        auto newInput = pre_gain * input * gain;
        
        if (newInput < a)
            out_1 = (k1 + newInput) / (k2 - newInput);
        if (newInput >= a && newInput <= b)
            out_1 = newInput;
        if (newInput > b)
            out_1 = (newInput - k3) / (newInput + k4);
        
        return amp_clamp(post_gain*out_1,X,Y);
    }
    inline void gallo_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType gain = (G + 1.0);
            DspFloatType a = -0.01;
            DspFloatType b = 0.7;
            DspFloatType k1 = std::pow(a, 2.0);
            DspFloatType k2 = 1 + (2 * a);
            DspFloatType k3 = std::pow(b, 2.0);
            DspFloatType k4 = 1 - (2 * b);
            DspFloatType out_1 = 0.0;
            
            auto newInput = pre_gain * input[i] * gain;
            
            if (newInput < a)
                out_1 = (k1 + newInput) / (k2 - newInput);
            if (newInput >= a && newInput <= b)
                out_1 = newInput;
            if (newInput > b)
                out_1 = (newInput - k3) / (newInput + k4);

            if(std::isnan(out_1)) 
                output[i] = input < 0? -1.0 : 1.0;
            else        
                output[i] = amp_clamp(post_gain*out_1,-1,1);
        }
    }

    DspFloatType doubleSoftClipper (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = (gain + 1.0);
        auto slope = 2.0;
        auto upperLim = 0.8;
        auto lowerLim = -1.0;
        auto upperSkew = 1.0;
        auto lowerSkew = 1.0;
        auto xOffFactor = 0.0;
        auto out = 0.0;
        
        auto xOff = (1.0 / slope) * std::pow(slope, xOffFactor);
        
        input *= (gain / 10.0);
        
        if (input > 0.0)
        {
            input = (input - xOff) * upperSkew;
            
            if (input >= 1.0 / slope)
                out = upperLim;
            else if (input <= -1.0 / slope)
                out = 0.0;
            else
                out = (3.0 / 2.0) * upperLim * (slope * input - std::pow(slope * input, 3.0) / 3.0) / 2.0 + (upperLim / 2.0);
        } else
        {
            input = (input + xOff) * lowerSkew;
            
            if (input >= 1.0 / slope)
                out = 0.0;
            else if (input <= -1.0 / slope)
                out = lowerLim;
            else
                out = (3.0 / 2.0) * -lowerLim * (slope * input - std::pow(slope * input, 3.0) / 3.0) / 2.0 + (lowerLim / 2.0);
        }
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    inline void double_soft_clipper_vector(size_t n, DspFloatType * in, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType gain = (G + 1.0);
            auto slope = 2.0;
            auto upperLim = 0.8;
            auto lowerLim = -1.0;
            auto upperSkew = 1.0;
            auto lowerSkew = 1.0;
            auto xOffFactor = 0.0;
            auto out = 0.0;
            
            auto xOff = (1.0 / slope) * std::pow(slope, xOffFactor);
            
            DspFloatType input = pre_gain * in[i] * (gain / 10.0);
            
            if (input > 0.0)
            {
                input = (input - xOff) * upperSkew;
                
                if (input >= 1.0 / slope)
                    out = upperLim;
                else if (input <= -1.0 / slope)
                    out = 0.0;
                else
                    out = (3.0 / 2.0) * upperLim * (slope * input - std::pow(slope * input, 3.0) / 3.0) / 2.0 + (upperLim / 2.0);
            } else
            {
                input = (input + xOff) * lowerSkew;
                
                if (input >= 1.0 / slope)
                    out = 0.0;
                else if (input <= -1.0 / slope)
                    out = lowerLim;
                else
                    out = (3.0 / 2.0) * -lowerLim * (slope * input - std::pow(slope * input, 3.0) / 3.0) / 2.0 + (lowerLim / 2.0);
            }
            if(std::isnan(out)) {
                output[i] = input < 0? -1.0 : 1.0;
            }
            else
                output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    DspFloatType crush (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = pre_gain*(gain + 1.0);
        auto out = 0.0;
        
        gain /= 100.0;
        
        DspFloatType dry = input;
        DspFloatType wet = std::round(input * std::pow(2, gain));
        out = (wet + dry)  * std::asin(gain) + dry;
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);        
    }
    inline void crush_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {               
            DspFloatType gain = (G + 1.0);
            auto out = 0.0;
            
            gain /= 100.0;
            
            DspFloatType dry = input[i];
            DspFloatType wet = std::round(input[i] * std::pow(2, gain));
            out = (wet + dry)  * std::asin(gain) + dry;
            if(std::isnan(out)) 
                output[i] = input < 0? -1.0 : 1.0;            
            else
                output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }
    DspFloatType tuboid (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = pre_gain*(gain + 1.0);
        auto ktp = 1.0;
        auto ktn = 3.0;
        auto sfn = 0.0;
        
        auto threshPos = 0.3;
        auto threshNeg = -0.7;
        
        auto out = 0.0;
        
        gain /= 10.0;
        
        auto so = input * gain;
        
        if (so >= threshPos)
            sfn = ktp * std::pow(so - threshPos, 3.0);
        else if (so <= threshNeg)
            sfn = -ktn * abs(std::pow(so - threshNeg, 3.0));
        else
            sfn = 0.0;
        
        so = (input - sfn) * gain;
        out = so;    
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);        
    }

    inline void tuboid_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType gain = (G + 1.0);           
            auto ktp = 1.0;
            auto ktn = 3.0;
            auto sfn = 0.0;
            
            auto threshPos = 0.3;
            auto threshNeg = -0.7;
            
            auto out = 0.0;
            
            gain /= 10.0;
            
            auto so = input[i] * gain;
            
            if (so >= threshPos)
                sfn = ktp * std::pow(so - threshPos, 3.0);
            else if (so <= threshNeg)
                sfn = -ktn * abs(std::pow(so - threshNeg, 3.0));
            else
                sfn = 0.0;
            
            so = (input[i] - sfn) * gain;
            out = so;    
            if(std::isnan(out)) 
                output[i] = input < 0? -1.0 : 1.0;
            else
                output[i] = amp_clamp(post_gain*out,-1.0,1.0);        
        }
    }

    DspFloatType pakarinen_yeh (DspFloatType input, DspFloatType gain, DspFloatType X =-1, DspFloatType Y=1)
    {
        gain = pre_gain*(gain + 1.0);
        auto out = 0.0;
        
        gain /= 100.0;
        
        auto x = input * gain;
        
        if ((x >= 0.320018) && (x <= 1.0))
            out = 0.630035;
        else if ((x >= -0.08905) && (x < 0.320018))
            out = (-6.153 * std::pow(x, 2.0)) + (3.9375 * x);
        else if ((x >= -1.0) && (x < -0.08905))
            out = (-0.75 * (1.0 - std::pow(1.0 - (abs(x) - 0.029), 12.0) + (0.333 * (abs(x) - 0.029)))) + 0.01;
        else
            out = -0.9818;
        
        out *= 1.5;
        if(std::isnan(out)) {
            return input < 0? -1.0 : 1.0;
        }
        out = out < X? X : out > Y? Y : out;
        return amp_clamp(post_gain*out,-1.0,1.0);
    }
    
    inline void pakarinen_yeh_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto gain = pre_gain*(G + 1.0);
            auto out = 0.0;            
            gain /= 100.0;            
            auto x = input[i] * gain;            
            if ((x >= 0.320018) && (x <= 1.0))
                out = 0.630035;
            else if ((x >= -0.08905) && (x < 0.320018))
                out = (-6.153 * std::pow(x, 2.0)) + (3.9375 * x);
            else if ((x >= -1.0) && (x < -0.08905))
                out = (-0.75 * (1.0 - std::pow(1.0 - (abs(x) - 0.029), 12.0) + (0.333 * (abs(x) - 0.029)))) + 0.01;
            else
                out = -0.9818;            
            out *= 1.5;
            if(std::isnan(out)) {
                output[i] = input < 0? -1.0 : 1.0;
            }            
            else output[i] = amp_clamp(post_gain*out,-1.0,1.0);
        }
    }

    template<typename T>
    T logiclip (T x) noexcept
    {
        return 2.0 / (1.0 + std::exp (-2.0 * x)) - 1.0;
    }
    inline void logiclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            output [i] = amp_clamp(post_gain*(2.0 / (1.0 + std::exp (-2.0 * pre_gain*G*input[i])) - 1.0),-1.0,1.0);
        }
    }
    template<typename T>
    T hardclip (T x) noexcept
    {
        return sgn (x) * std::fminf (std::fabs(x), 1.0);
    }
    inline void hardclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain;
            auto s = x < 0? -1.0 : 1.0;
            output[i] = amp_clamp(post_gain* s * std::fminf(std::fabs(x), 1.0),-1.0,1.0);
        }
    }
    template<typename T>
    T tanhclip (T x) noexcept
    {
        DspFloatType soft = 0.0;
        return std::tanh ((1.0 - 0.5 * soft) * x);
    }
    inline void tanhclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain;
            output[i] = amp_clamp(post_gain*std::tanh ((1.0 - 0.5 * G) * x),-1.0,1.0);
        }
    }

    template<typename T>
    T quintic (T x) noexcept
    {
        if (std::fabs (x) < 1.25)
        {
            return x - (256.0 / 3125.0) * std::pow (x, 5.0);
        } else
        {
            return sgn (x) * 1.0;
        }
    }
    inline void quintic_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain;
            auto s = x < 0? -1.0:1.0;
            if (std::fabs (x) < 1.25)
            {
                output[i] = amp_clamp(post_gain * x - (256.0 / 3125.0) * std::pow (G*x, 5.0),-1.0,1.0);
            } else
            {
                output[i] = amp_clamp(post_gain * s * 1.0,-1.0,1.0);
            }
        }
    }

    template<typename T>
    T cubicBasic (T x) noexcept
    {
        if (std::fabs (x) < 1.5)
        {
            return x - (4.0 / 27.0) * std::pow (x, 3.0);
        } else
        {
            return sgn (x) * 1.0;
        }
    }
    inline void basic_cubic_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain*G;
            auto s = x < 0? -1.0:1.0;
            if (std::fabs(x) < 1.5)
            {
                x = x - (4.0 / 27.0) * std::pow (x, 3.0);
            } else
            {
                x= s * 1.0;
            }
            output[i] = amp_clamp(post_gain * x, -1.0, 1.0);
        }
    }
    
    template<typename T>
    T algClip (T x) noexcept
    {
        DspFloatType soft = 0.0;
        return x / std::sqrt ((1.0 + 2.0 * soft + std::pow (x, 2.0)));
    }
    
    inline void algclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain;
            DspFloatType soft = G;
            output[i] =  amp_clamp( post_gain * x / std::sqrt ((1.0 + 2.0 * soft + std::pow (x, 2.0))), -1.0, 1.0);
        }
        }

    template<typename T>
    T arcClip (T x) noexcept
    {
        DspFloatType soft = 0.0;
        return (2.0 / M_PI) * std::atan ((1.6 - soft * 0.6) * x);
    }
    inline void arcclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain;
            DspFloatType soft = G;
            output[i] = amp_clamp(post_gain * (2.0 / M_PI) * std::atan ((1.6 - soft * 0.6) * x),-1.0,1.0);
        }
    }

    template<typename T>
    T sinclip (T x) noexcept
    {
        if (std::fabs (x) < M_PI)
        {
            return std::sin (x);
        }
        else
        {
            return sgn (x) * 1.0;
        }
    }
    inline void sinclip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            DspFloatType x = input[i]*pre_gain*G;
            auto s = x < 0? -1.0:1.0;
            if (std::fabs (x) < M_PI)
            {
                x = std::sin (x);
            }
            else
            {
                x = s * 1.0;
            }
            output[i] = amp_clamp(post_gain * x, -1.0, 1.0);
        }
    }
    inline DspFloatType FuzzCtrTable(const DspFloatType x)
    {
        static auto gen = std::minstd_rand(2112);
        static const DspFloatType b = 20;
        static auto dist = std::uniform_real_distribution<DspFloatType>(-1.0, 1.0);
        auto g = exp(-x * x * b);
        auto xadj = x + g * dist(gen);
        return xadj;
    }
    inline void fuzzctrtable_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        static auto gen = std::minstd_rand(2112);
        static const DspFloatType b = 20;
        static auto dist = std::uniform_real_distribution<DspFloatType>(-1.0, 1.0);
        #pragma omp simd                
        for(size_t i = 0; i < n; i++) {    
            auto x = input[i];
            auto g = exp(-x * x * b);
            auto xadj = x + g * dist(gen);
            output[i] = xadj;
        }
    }
    template <int scale> DspFloatType FuzzTable(const DspFloatType x)
    {
        static auto gen = std::minstd_rand(2112);        
        static const DspFloatType range = 0.1 * scale;
        static auto dist = std::uniform_real_distribution<DspFloatType>(-range, range);
        auto xadj = x * (1 - range) + dist(gen);
        return xadj;
    }
    inline void fuzztable_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType G = 1.0) {
        static auto gen = std::minstd_rand(2112);
        static const DspFloatType range = 0.1 * 10;
        static auto dist = std::uniform_real_distribution<DspFloatType>(-range, range);
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            output[i] = input[i] * (1 - range) + dist(gen);
        }
    }

    inline  DspFloatType absfunc(DspFloatType x) {
        return (x >= 0.0 ? x : -x);
    }
    // cube function
    inline  DspFloatType cubefunc(DspFloatType x) {
        return (x * x * x);
    }
    // use this to process audio via the rectification algorithm
    inline  DspFloatType Rectify(DspFloatType sample, DspFloatType RectifierThreshold=0.9) {
        return ((1 - RectifierThreshold) * sample) + (absfunc(sample) * RectifierThreshold);
    };
    inline void rectify_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType RectifierThreshold=0.9) {        
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            output[i] = ((1 - RectifierThreshold) * sample) + (absfunc(sample) * RectifierThreshold);
        }
    }
    // hard clip of input signal
    inline  DspFloatType HardClip(DspFloatType sample, DspFloatType thresh) {
        return 0.5 * (absfunc(sample + thresh) - absfunc(sample - thresh));
    };
    inline void hardclip2_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType thresh = 0.5) {        
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            output[i] = 0.5 * (absfunc(sample + thresh) - absfunc(sample - thresh));
        }
    }

    // cubic soft clip function
    inline DspFloatType SoftCubicClip(DspFloatType sample, DspFloatType thresh=0.5) {
        DspFloatType threshInv = 1 / thresh;
        return threshInv * ((thresh * 1.5 * HardClip(sample, thresh)) -
            (0.5 * cubefunc(HardClip(sample, thresh)) * threshInv));
    };
    inline void soft_cubic_clip_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType thresh = 0.5) {        
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            DspFloatType threshInv = 1 / thresh;
            output[i] = threshInv * ((thresh * 1.5 * HardClip(sample, thresh)) -
                        (0.5 * cubefunc(HardClip(sample, thresh)) * threshInv));
        }
    }

    // use this to process audio via the SoftCubicClip algorithm
    inline DspFloatType SoftCubic(DspFloatType sample, DspFloatType CubicSoftClipThreshold=0.9, DspFloatType CubicHarmonicBalance=0.5) {
        DspFloatType invsqrt2 = 1.0/sqrt(2);
        return (invsqrt2 / 3) * (SoftCubicClip(sample, CubicSoftClipThreshold) +
            (CubicHarmonicBalance * SoftCubicClip(absfunc(sample), CubicSoftClipThreshold)));
    };
    
    inline void soft_cubic_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType CubicSoftClipThreshold=0.9, DspFloatType CubicHarmonicBalance=0.5) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            DspFloatType invsqrt2 = 1.0/sqrt(2);
            output[i] = (invsqrt2 / 3) * (SoftCubicClip(sample, CubicSoftClipThreshold) +
            (CubicHarmonicBalance * SoftCubicClip(absfunc(sample), CubicSoftClipThreshold)));
        }
    }

    // soft clip function with adjustable knee
    inline  DspFloatType SKClip(DspFloatType sample, DspFloatType knee) {
        return sample / (knee * absfunc(sample) + 1.0);
    };

    
    // use this to process audio via the SKClip algorithm
    inline  DspFloatType SoftKnee(DspFloatType sample, DspFloatType SoftClipKnee=0.9) {
        return 0.5 * (SKClip(sample, SoftClipKnee) + ((SoftClipKnee / 2.0) * SKClip(absfunc(sample), SoftClipKnee)));
    };
    inline void soft_knee_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType SoftClipKnee=0.9) {
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            output[i] = 0.5 * (SKClip(sample, SoftClipKnee) + ((SoftClipKnee / 2.0) * SKClip(absfunc(sample), SoftClipKnee)));
        }
    }

    // use this to process audio via the leaky integrator algorithm
    inline  DspFloatType LeakyInt(DspFloatType sample, DspFloatType previous_sample, DspFloatType TcRise = 0.5, DspFloatType TcFall=0.5) {
        DspFloatType invsqrt2 = 1.0/sqrt(2);
        if (sample > previous_sample) 
        {
        return invsqrt2 * (((1.0 - TcRise) * sample) + (TcRise * previous_sample));
        }
        else {
            return invsqrt2 * (((1.0 - TcFall) * sample) + (TcFall * previous_sample));
        } 
    }
    
    inline DspFloatType linearScale( DspFloatType in, DspFloatType min, DspFloatType max )
    {
        DspFloatType ret;
        if ( min == 0.0 && max == 0.0 )
        {
            ret = 0.0;
        }
        else if ( min > max )
        {
            ret = min - ( in * ( min - max ) );
        }
        else
        {
            ret = min + ( in * ( max - min ) );
        }
        return ret;
    }
    
    inline DspFloatType linearDescale( DspFloatType in, DspFloatType min, DspFloatType max )
    {
        DspFloatType ret;
        if ( min == 0.0 && max == 0.0 )
        {
            ret = 0.0;
        }
        else if ( min > max )
        {
            ret = ( min - in ) / ( min - max );
        }
        else
        {
            ret = ( in - min ) / ( max - min );
        }
        return ret;
    }
    inline DspFloatType expoScale( DspFloatType in, DspFloatType min, DspFloatType max )
    {
        // negative log makes no sense...
        if ( min < 0.0 || max < 0.0 )
        {
            return 0.0;
        }
        // not handling min > max (inverse) case yet
        // figure out how many "octaves" (doublings) it takes to get from min to
        // max
        // we only have log & log10 so we have to do change of base
        // note this uses + instead of * so we can handle min == 0
        DspFloatType octaves = log( max - min + 1 ) / log( 2.0 );
        return ( min - 1 ) + std::pow( 2.0, in * octaves );
    }
    inline DspFloatType expoDescale( DspFloatType in, DspFloatType min, DspFloatType max )
    {
        // see above
        if ( min < 0.0 || max < 0.0 )
        {
            return 0.0;
        }
        // again, not handling min > max (inverse) case yet
        
        // note this was derived by simply inverting the previous function
        DspFloatType log2 = log( 2.0 );
        return ( log( in - min + 1 ) / log2 ) / ( log( max - min + 1 ) / log2 );
    }
    inline DspFloatType floorScale( DspFloatType in, DspFloatType min, DspFloatType max )
    {
        if ( min > max )
        {
            return ceil( linearScale( in, min, max ) );
        }
        else
        {
            return floor( linearScale( in, min, max ) );
        }
    }
    inline DspFloatType expoShape( DspFloatType in, DspFloatType amount )
    {
        if ( in == 0.0 )
            return in;
        DspFloatType flip = in < 0.0 ? -1.0 : 1.0;
        return std::pow( in * flip, amount ) * flip;
    }
    inline void exposhape_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType amount) {    
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto in = input[i];            
            if ( in == 0.0 )
                output[i] = 0.0;
            else {
                DspFloatType flip = in < 0.0 ? -1.0 : 1.0;
                output[i] = std::pow( in * flip, amount ) * flip;
            }
        }
    }

    inline DspFloatType softClipShape( DspFloatType in, DspFloatType amount )
    {
        return in / ( 1 + fabs( in ) );
    }
    inline void softclipshape_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType amount) {    
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto in = input[i];        
            output[i] = in / ( 1 + fabs( in ) );
        }
    }    

    inline DspFloatType sineShape( DspFloatType in, DspFloatType amount )
    {
        return in * cos( in * amount );
    }
    inline void sineshape_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType amount) {    
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto in = input[i];        
            output[i] = in * cos( in * amount );
        }
    }

    inline DspFloatType chebyshevRec( DspFloatType sample, int depth, DspFloatType softClipThreshold=0.9 )
    {
        if ( depth == 0 )
        {
            return 1.0;
        }
        // lastval represents C(k-1)DspFloatType Distortion::softClip(DspFloatType sample)
        {
            if (sample < -1.f) {
                return -softClipThreshold; //2/3
            }
            else if (sample > 1.f) {
                return softClipThreshold;
            }
            else {
                return sample - ((sample * sample * sample) / 3.f);
            }
        }
    }
    
    inline DspFloatType chebyshevShape( DspFloatType in, DspFloatType amount )
    {
        return chebyshevRec( in, (int)amount );
    }

    // Arctangent nonlinearity
    inline DspFloatType arctangent(DspFloatType sample, DspFloatType alpha)
    {
        // f(x) = (2 / M_PI) * arctan(alpha * x[n]), where alpha >> 1 (drive param)
        return (2.f / M_PI)* atan(alpha * sample);
    }
    // Hard-clipping nonlinearity
    inline DspFloatType hardClip(DspFloatType sample)
    {
        if (sample < -1.f) {
            return -1.f;
        }
        else if (sample > 1.f) {
            return 1.f;
        }
        else {
            return sample;
        }
    }
    // Square law series expansion
    inline DspFloatType squareLaw(DspFloatType sample, DspFloatType alpha)
    {
        return sample + alpha * sample * sample;
    }
    inline void squarelaw_vector(size_t n, DspFloatType * input, DspFloatType *output, DspFloatType alpha) {    
        #pragma omp simd        
        for(size_t i = 0; i < n; i++) {    
            auto sample = input[i];
            output[i] = sample + alpha * sample * sample;
        }
    }

    inline DspFloatType cubicWaveShaper(DspFloatType sample)
    {
        return 1.5 * sample - 0.5 * sample * sample * sample;
    }
    

    // Foldback nonlinearity, input range: (-inf, inf)
    inline DspFloatType foldback(DspFloatType sample, DspFloatType threshold=0.96)
    {
        // Threshold should be > 0.f
        if (sample > threshold || sample < -threshold) {
            sample = fabs(fabs(fmod(sample - threshold,
                                    threshold * 4))
                        - threshold * 2) - threshold;
        }
        return sample;
    }
    
    // A nonlinearity by Partice Tarrabia and Bram de Jong
    inline DspFloatType waveShaper1(DspFloatType sample, DspFloatType alpha)
    {
        const DspFloatType k = 2.f * alpha / (1.f - alpha);
        return (1.f + k) * sample / (1.f + k * fabs(sample));
    }
    // A nonlinearity by Jon Watte
    inline DspFloatType waveShaper2(DspFloatType sample, DspFloatType alpha)
    {
        const DspFloatType z = M_PI * alpha;
        const DspFloatType s = 1.f / sin(z);
        const DspFloatType b = 1.f / alpha;
        
        if (sample > b) {
            return 1.f;
        }
        else {
            return sin(z * sample) * s;
        }
    }
    // A nonlinearity by Bram de Jong, input range: [-1, 1]
    inline DspFloatType waveShaper3(DspFloatType sample, DspFloatType alpha)
    {
        // original design requires sample be positive
        // alpha: [0, 1]
        bool isNegative = false;
        DspFloatType output = sample;
        if (sample < 0.f) {
            isNegative = true;
            output = -output;
        }
        
        if (output > alpha) {
            output = alpha + (output - alpha)
                / (1.f + std::pow(((output - alpha)) / (1.f - alpha), 2.f));
        }
        if (output > 1.f) {
            output = (alpha + 1.f) / 2.f;
        }
        
        if (isNegative) {
            output = -output;
        }
        
        return output;
    }
    
    
    inline DspFloatType gloubiBoulga(DspFloatType sample)
    {
        const DspFloatType x = sample * 0.686306;
        const DspFloatType a = 1 + exp(sqrt(fabs(x)) * -0.75);
        return (exp(x) - exp(-x * a)) / (exp(x) + exp(-x));
    }
    // Approximation based on description in gloubiBoulga
    inline DspFloatType gloubiApprox(DspFloatType sample)
    {
        return sample - (0.15 * sample * sample) - (0.15 * sample * sample * sample);
    }
    inline DspFloatType FuzzEdgeTable(const DspFloatType x)
    {
        static auto gen = std::minstd_rand(2112);
        static auto dist = std::uniform_real_distribution<DspFloatType>(-1.0, 1.0);
        auto g = x * x * x * x;
        auto xadj = 0.85 * x + 0.15 * g * dist(gen);
        return xadj;
    }
    template<typename T>
    T limitclip (T x) noexcept
    {
        return clamp(x,-0.1, 0.1);
    }
           
    inline DspFloatType hardClipping(const DspFloatType& _in)
    {
        DspFloatType out;
        DspFloatType threshold = 0.5;
        
        if (_in > threshold)
            out = threshold;
        else if (_in < -threshold)
            out = -threshold;
        else
            out = _in;
        
        return out;
    }
    inline DspFloatType softClipping(const DspFloatType& _in)
    {
        DspFloatType out;
        DspFloatType threshold1 = 1.0 / 3.0;
        DspFloatType threshold2 = 2.0 / 3.0;
        
        if (_in > threshold2)
            out = 1.0;
        else if (_in > threshold1)
            out = 1.0 - std::pow (2.0 - 3.0 * _in, 2.0) / 3.0;
        else if (_in < -threshold2)
            out = -1.0;
        else if (_in < -threshold1)
            out = -1.0 + std::pow (2.0 + 3.0 * _in, 2.0) / 3.0;
        else
            out = 2.0 * _in;
            out *= 0.5;
        
        return out;
    }
    inline DspFloatType exponential(const DspFloatType& _in)
    {
        DspFloatType out;
        
        if (_in > 0.0)
            out = 1.0 - expf (-_in);
        else
            out = -1.0 + expf (_in);
        
        return out;
    }
    inline DspFloatType fullWaveRectifier(const DspFloatType& _in)
    {
        DspFloatType out;
        
        out = fabs (_in);
        
        return out;
    }
    inline DspFloatType halfWaveRectifier(const DspFloatType& _in)
    {
        DspFloatType out;
        
        if (_in > 0.0)
            out = _in;
        else
            out = 0.0;
        
        return out;
        
    }
    inline DspFloatType ArayaAndSuyama(const DspFloatType &_in)
    {
        DspFloatType out;
        
        out = (3/2) * (_in) * (1 - std::pow(_in, 2)/3);
        out = (3/2) * (out) * (1 - std::pow(out, 2)/3);
        out = (3/2) * (out) * (1 - std::pow(out, 2)/3);
        return out;
        
    }
    inline DspFloatType doidicSymmetric(const DspFloatType& _in)
    {
        DspFloatType out;
        
        out =  ( (2*fabs(_in))  - std::pow(_in, 2)) * copysignf(1, _in);
        
        return out;
    }
    inline DspFloatType doidicAssymetric(const DspFloatType& _in)
    {
        DspFloatType out;
        
        if (_in >= -1 && _in < -0.08905) {
            out = -(0.75)*( 1 - (1 - std::pow(std::fabs(_in) - 0.032847, 12)) + 1/3*(std::fabs(_in) - 0.032847)) + 0.01;
        }
        else if (_in >= -0.08905 && _in < 0.320018)
        {
            out = -6.153 * std::pow(_in,2) + 3.9375 * _in;
        }
        else if (_in >= 0.320018 && _in <= 1)
        {
            out = 0.630035;
        }
        
        return out;        
    }
    inline void InitRandom() {
        std::srand(std::time(NULL));
    }
}

    

