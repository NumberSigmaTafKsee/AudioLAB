#pragma once

/* FxDSP Processor
Biquad
CircularBuffer
Decimator
DiodeRectifier
DiodeSaturator
FFT
FIRFilter
LadderFilter
LinkwitzRileyFilter
Metering
MidiUtils
MultibandBank
OnePole
Optocoupler
PanLaw
PolySaturation
Polyphase
RBJFilter
RMSEstimator
SpectrumAnalyzer
Stereo
Tape
Upsampler
Utilities
WindowFunctions
bs1770
*/

#include "FxDSP/BiquadFilter.h"
#include "FxDSP/Decimator.h"
#include "FxDSP/DiodeRectifier.h"
#include "FxDSP/DiodeSaturator.h"
#include "FxDSP/FFT.h"
#include "FxDSP/FIRFilter.h"

namespace DSPFX
{
    struct FxBiquadFilter
    {
        BiquadFilter * filter;

        FxBiquadFilter(const float * Bc, const float * Ac)
        {
            filter = BiquadFilterInit(Bc,Ac);
            assert(filter != NULL);
        }    
        ~FxBiquadFilter() {
            if(filter) BiquadFilterFree(filter);
        }
        void flush() {
            BiquadFilterFlush(filter);
        }
        void Process(size_t n, const float * in, float * out)
        {
            BiquadFilterProcess(filter,out,in,n);
        }
        float Tick(float In) {
            return BiquadFilterTick(filter,In);
        }
        void updateCoefficients(const float * B, const float * A) {
            BiquadFilterUpdateKernel(B,A);
        }
    };

    struct FxRBJFilter
    {
        RBJFilter * filter;
        float Cutoff,Q;
        
        enum Type
        {
        /** Lowpass */
        LOWPASS,

        /** Highpass */
        HIGHPASS,

        /** Bandpass */
        BANDPASS,

        /** Allpass */
        ALLPASS,

        /** Notch */
        NOTCH,

        /** Peaking */
        PEAK,

        /** Low Shelf */
        LOW_SHELF,

        /** High Shelf */
        HIGH_SHELF,

        /** Number of Filter types */
        N_FILTER_TYPES
        };

        FxRBJFilter(float cutoff, float sampleRate)
        {
            Cutoff = cutoff;
            Q      = 0.5;
            filter = RBJFilterInit(cutoff,sampleRate);
            assert(filter != NULL);
        }    
        ~FxRBJFilter() {
            if(filter) RBJFilterFree(filter);
        }
        void setType(Type type) {
            RBJFilterSetType(filter,(Filter_t)type);
        }
        void setCutoff(float cut) {
            Cutoff = cut;
            RBJFilterSetCutoff(filter,cut);
        }
        void setQ(float q) {
            Q = q;
            RBJFilterSetQ(filter,Q);
        }
        void flush() {
            RBJFilterFlush(filter);
        }
        void Process(size_t n, const float * in, float * out)
        {
            RBJFilterProcess(filter,out,in,n);
        }
        
        float Tick(float I) 
        {
            float out = 0;
            Process(1,&I,&out);
            return out;
        }
    };

    struct FxDecimator
    {
        Decimator * deci;

        /** Resampling Factor constants */
        enum factor
        {
            /** 2x resampling */
            X2 = 0,

            /** 4x resampling */
            X4,

            /** 8x resampling */
            X8,

            /** 16x resampling */
            /*X16,*/

            /** number of resampling factors */
            N_FACTORS
        };
        FxDecimator(factor fac) {
            deci = DecimatorInit(fac);
            assert(deci != NULL);
        }
        ~FxDecimator() {
            if(deci) DecimatorFree(deci);
        }
        void flush() {
            DecimatorFlush(deci);
        }
        void Process(size_t n, const float * in, float * out) {
            Error_t err = DecimatorProcess(deci,out,in,n);
        }
    };


    struct FxDiodeRectifier
    {
        enum _bias_t
        {
        /** Pass positive signals, clamp netagive signals to 0 */
        FORWARD_BIAS,

        /** Pass negative signals, clamp positive signals to 0 */
        REVERSE_BIAS,

        /** Full-wave rectification. */
        FULL_WAVE
        };
        DiodeRectifier * rect;

        FxDiodeRectifier( _bias_t bias, float thresh)
        {
            rect = DiodeRectifierInit(bias,thresh);
            assert(rect != NULL);
        }
        ~FxDiodeRectifier() {
            if(rect) DiodeRectifierFree(rect);
        }
        void setThreshold(float t) {
            Error_t err = DiodeRectifierSetThreshold(rect,t);
        }
        void ProcessBlock(size_t n, const float * in, float * out) {
            Error_t err = DiodeRectifierProcessBlock(rect,out,in,n);
        }
        float Tick(float in) {
            return DiodeRectifierTick(rect,in);
        }
    };

    struct FxDiodeSaturator : public AmplifierProcessor
    {
        enum _bias_t
        {
        /** Pass positive signals, clamp netagive signals to 0 */
        FORWARD_BIAS,

        /** Pass negative signals, clamp positive signals to 0 */
        REVERSE_BIAS,

        /** Full-wave rectification. */
        FULL_WAVE
        };
        DiodeSaturator * rect;

        FxDiodeSaturator( _bias_t bias, float amt) : AmplifierProcessor()
        {
            rect = DiodeSaturatorInit(bias,amt);
            assert(rect != NULL);
        }
        ~FxDiodeSaturator() {
            if(rect) DiodeSaturatorFree(rect);
        }
        void setThreshold(float t) {
            Error_t err = DiodeSaturatorSetAmount(rect,t);
        }
        void ProcessBlock(size_t n, const float * in, float * out) {
            Error_t err = DiodeSaturatorProcessBlock(rect,out,in,n);
        }
        double Tick(double in, double A=1, double X=1, double Y=1) {
            return DiodeSaturatorTick(rect,in);
        }        
    };

    struct FxFFT : public SpectrumProcessor
    {
        FFTConfig * config;

        FxFFT(size_t n) : SpectrumProcessor() {
            config = FFTInit(n);
            assert(config != NULL);
        }
        ~FxFFT() {
            if(config) FFTFree(config);
        }
        void R2C(const float * in, float * real, float * imag) {
            Error_t err = FFT_R2C(config, in, real, imag);
        }
        void C2R(const float * real, const float * imag, float * out) {
            Error_t err = FFT_R2C(config, real, imag, out);
        }
        void convolve(const float * in1, size_t len1, const float *in2, size_t len2, float * out) {
            FFTConvolve(config,in1,len1,in2,len2,out);
        }
        void filterConvolve(const float * in, size_t len, float * real, float * imag, float * out)
        {
            FFTSplitComplex split;
            split.realp = real;
            split.imagp = imag;
            FFTFilterConvolve(config,in,len,split,out);
        }        
        void print(float * buffer) {
            FFTdemo(config,buffer);
        }        
    };


    struct FxFIRFilter : public MonoFXProcessor
    {
        /** Convolution Algorithm to use */
        enum 
        {
            /** Choose the best algorithm based on filter size */
            BEST    = 0,

            /** Use direct convolution */
            DIRECT  = 1,

            /** Use FFT Convolution (Better for longer filter kernels */
            FFT     = 2

        };

        FIRFilter * fir;

        FxFIRFilter(float * kernel, size_t len, int mode) : MonoFXProcessor() {
            fir = FIRFilterInit(kernel,len,mode);
            assert(fir != NULL);
        }
        ~FxFIRFilter() {
            if(fir) FIRFilterFree(fir);
        }
        void flush() { FIRFilterFlush(fir); }
        void ProcessBlock(size_t n, float * in, float * out) {
            FIRFilterProcess(fir,out,in,n);
        }
        void updateKernel(float * kernel) {
            FIRFilterUpdateKernel(fir,kernel);
        }
    };

    struct FxLadderFilter : public MonoFXProcessor
    {
        LadderFilter * filter;

        FxLadderFilter(float sr) {
            filter = LadderFilterInit(sr);
            assert(filter != NULL);
        }
        ~FxLadderFilter() {
            if(filter) LadderFilterFree(filter);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            LadderFilterProcess(filter,out,in,n);
        }
        void setCutoff(float c) {
            LadderFilterSetCutoff(filter,c);
        }
        void setResonance(float q) {
            LadderFilterSetResonance(filter,q);
        }
        void setTemperature(float t) {
            LadderFilterSetTemperature(filter,t);
        }
    };

    struct FxLRFilter
    {
        LRFilter * filter;
        enum Type
        {
        /** Lowpass */
        LOWPASS,

        /** Highpass */
        HIGHPASS,

        /** Bandpass */
        BANDPASS,

        /** Allpass */
        ALLPASS,

        /** Notch */
        NOTCH,

        /** Peaking */
        PEAK,

        /** Low Shelf */
        LOW_SHELF,

        /** High Shelf */
        HIGH_SHELF,

        /** Number of Filter types */
        N_FILTER_TYPES
        };
        float cutoff, Q;
        int type = LOWPASS;
        FxLRFilter(int type, float cut, float q, float sr = sampleRate) {
            filter = LRFilterInit(type,cut,q,sr);
            assert(filter != NULL);
            cutoff = cut;
            Q = q;
            this->type = type;
        }
        ~FxLRFilter() {
            if(filter) LRFilterFree(filter);
        }
        void flush() { LRFilterFlush(filter); }
        void setParams(int type, float cut, float q) {
            this->type = type;
            cutoff = cut;
            Q = q;
            LRFilterSetParams(filter,type,cutoff,Q);
        }
        void setCutoff(float cut) {
            cutoff = cut;
            LRFilterSetParams(filter,type,cutoff,Q);
        }
        void setQ(float q) {
            Q = q;
            LRFilterSetParams(filter,type,cutoff,Q);
        }
        void setType(int type) {
            this->type = type;
            LRFilterSetParams(filter,type,cutoff,Q);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            LRFilterProcessBlock(filter,out,in,n);
        }
    };

    struct FxMetering
    {

        enum
        {
            FULL_SCALE,
            K_12,
            K_14,
            K_20
        };

        static float phase_correlation(float * left, float * right, size_t n) {
            return ::phase_correlation(left,right,n);
        }
        static float balance(float * left, float * right, size_t n) {
            return ::balance(left,right,n);
        }
        static float vu_peak(float * signal, size_t n, int scale);
            return ::vu_peak(left,right,n);
        }

    };

    struct FxMidi
    {
        static float midiNoteToFrequency(unsigned note) {
            return ::midiNoteToFrequency(note);
        }
        static unsigned frequencyToMidiNote(float f) {
            return ::frequencyToMidiNote(f);
        }
    };

    struct FxMultibandFilter
    {
        MultibandFilter * filter;

        FxMultiBandFilter(float low, float high, float sr = sampleRate)
        {
            filter = MultibandFilterInit(low,high,sr);
            assert(filter != NULL);
        }
        ~FxMultiBandFilter() {
            if(filter) MultibandFilterFree(filter);
        }
        void flush() { MultibandFilterFlush(filter); }
        void setHighCutoff(float c) {
            MultibandFilterSetHighCutoff(filter,c);
        }
        void setLowCutoff(float c) {
            MultibandFilterSetLowCutoff(filter,c);
        }
        void update(float low, float high) {
            MultibandFilterUpdate(filter,low,high);
        }
        void ProcessBlock(size_t n, float * in, float * low, float * mid, float * high) {
            MultibandFilterProcess(filter,low,mid,high,in,n);
        }
    };

    struct FxOnePoleFilter : FilterProcessor
    {
        OnePole * filter;

        enum Type
        {
        /** Lowpass */
        LOWPASS,

        /** Highpass */
        HIGHPASS,

        /** Bandpass */
        BANDPASS,

        /** Allpass */
        ALLPASS,

        /** Notch */
        NOTCH,

        /** Peaking */
        PEAK,

        /** Low Shelf */
        LOW_SHELF,

        /** High Shelf */
        HIGH_SHELF,

        /** Number of Filter types */
        N_FILTER_TYPES
        };
        float cutoff;
        // only low and highpass are valie
        FxOnePoleFilter(float cut, int type=LOWPASS, float sr=sampleRate)
        : FilterProcessor()
        {
            cutoff = cut;
            filter = OnePoleInit(cut,sr,type);
            assert(filter != NULL);
        }
        ~FxOnePoleFilter() {
            if(filter) OnePoleFree(filter);
        }
        void flush() { OnePoleFlush(filter); }
        void setType(int type) {
            OnePoleSetType(filter,type);
        }
        void setCutoff(float c) {
            cutoff = c;
            OnePoleSetCutoff(filer,c);
        }
        void setCoefficients(float beta, float alpha) {
            OnePoleSetCoefficients(&beta, &alpha);
        }
        void ProcessBlock(size_t n, float * in, float * out)
        {
            OnePoleProcess(filter,out,in,n);
        }
        double Tick(double I, double A=1, double X=1, double Y=1)
        {
            double c = cutoff;            
            setCutoff(cutoff * X * Y);
            double out = OnePoleTick(I);
            setCutoff(c);
            return out * A;
        }
    };

    struct FxOptoCoupler
    {
        enum 
        {
            /** Light-Dependent-Resistor output. Based
             on Vactrol VTL series. datasheet:
            http://pdf.datasheetcatalog.com/datasheet/perkinelmer/VT500.pdf
            Midpoint Delay values:
            Turn-on delay:   ~10ms
            Turn-off delay:  ~133ms
            */
            OPTO_LDR,

            /** TODO: Add Phototransistor/voltage output opto model*/
            OPTO_PHOTOTRANSISTOR
        };

        OptoCoupler * oc;
        double delay;

        FxOptoCoupler(int type, float delay, float sr=sampleRate)
        {
            this->delay = delay;
            oc = OptoInit(type,delay,sr);
            assert(oc != NULL);
        }
        ~FxOptoCoupler() {
            OptoFree(oc);
        }
        void setDelay(float d) {
            delay = d;
            OptoSetDelay(oc,d);
        }
        void ProcessBlock(size_t n, float * in, float * out)
        {
            OptoProcess(oc,out,in,n);
        }
        double Tick(double I, double A=1, double X=1, double Y=1)
        {
            double d = delay;
            setDelay(d*X*Y);
            double r = OptoTick(I);
            setDelay(d);
            return A*r;
        }
    };

    struct FxPanLaw
    {
        void linear_pan(float control, float * l_gain, float * r_gain) {
            ::linear_pan(control,l_gain,r_gain);
        }
        void equal_power_3db_pan(float control, float * l_gain, float * r_gain) 
        {
            ::equal_power_3db_pan(control,l_gain,r_gain);
        }
        void equal_power_6db_pan(float control, float * l_gain, float * r_gain) 
        {
            ::equal_power_6db_pan(control,l_gain,r_gain);
        }
    };

    struct FxPolySaturator
    {
        PolySaturator * sat;
        double N;

        FxPolySatuator(float n) {
            N = n;
            sat = PolySaturatorInit(n);
            assert(sat != NULL);
        }
        ~FxPolySaturator() {
            if(sat) PolySaturatorFree(sat);
        }
        void setN(float n) {
            N = n;
            PolySaturatorSetN(sat,n);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            PolySaturatorProcess(sat,out,in,n);
        }
        double Tick(double I, double A=1, double X=1, double Y=1)
        {
            double t = N;
            setN(N*X*Y);
            double r = PolySaturatorTick(I);
            setN(t);
            return A*r;
        }
    };

    struct FxRMSEstimator
    {
        RMSEstimator *rms;

        RMSEstimator(float avg, float sr=sampleRate)
        {
            rms = RMSEstimatorInit(avg,sr);
            assert(rms != NULL);
        }
        ~RMSEstimator()
        {
            if(rms) RMSEstimatorFree(rms);
        }
        void flush() { RMSEstimatorFlush(rms); }
        void ProcessBlock(size_t n, float * in, float * out) {
            RMSEstimatorProcess(rms,out,in,n);
        }
        // this is never modulated
        double Tick(double I, double A=1, double X=1, double Y=1) {
            return RMSEstimatorTick(I);
        }
    };

    struct FxSpectrumAnalyzer
    {
        SpectrumAnalyzer *spc;

        FxSpectrumAnalyzer(size_t len, float sr=sampleRate) {
            spc = SpectrumAnalyzerInit(len,sr);
            assert(spc != NULL);
        }
        ~FxSpectrumAnalyzer() {
            if(spc) SpectrumAnalyzerFree(spc);
        }
        float getCentroid() {
            return SpectrumAnalyzerCentroid(spc);
        }
        float getSpread() {
            return SpectrumAnalyzerSpread(spc);
        }
        float getSkewness() {
            return SpectrumAnalyzerSkewnes(spc);
        }
        float getKurtosis() {
            return SpectrumAnalyzerKurtosis(spc);
        }
        // todo : should be able to get magnitude and phase
        // should be able to set the bins and inverse it to samples too    
        void Analyze(float * signal) {
            SpectrumAnalyzerAnalyze(spc,signal);
        }
    };


}
