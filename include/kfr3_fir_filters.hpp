#pragma once

namespace kfr3
{
    template<typename T>
    struct FIRFilter {
    private:
        kfr::filter_fir<SampleType> * filter;
        kfr::univector<T> taps;
        
    public:
        
        // need to be able to input coefficients
        FIRFilter(size_t num_taps) { 
            taps.resize(num_taps); 
            filter = nullptr;
        }
        FIRFilter(const kfr::univector<T> & taps) {
            filter = new kfr::filter_fir<T>(taps);
        }
        ~FIRFilter() {
            if(filter != NULL) delete filter;
        }
        void set_taps(const kfr::univector<T> & taps) {
            filter = new kfr::filter_fir<T>(taps);
        }
        void bandpass(T x, T y, kfr::expression_pointer<T> & window, bool normalize=true ) {        
            kfr::fir_bandpass(taps,x,y,window,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void bandstop(T x, T y, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_bandstop(taps,x,y, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void highpass(T cutoff, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_highpass(taps,cutoff, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void lowpass(T cutoff, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_lowpass(taps,cutoff, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void apply(T * data, size_t s) {
            filter->apply(data,s);
        }
        void apply(size_t n, T * in, T * out) {
            filter->apply(out,in,n);
        }
        void reset() { filter->reset(); }

    };
    
    template<typename T>
    struct FIRBandpassFilter : public FIRFilter<T>
    {        
        FIRBandpassFilter(size_t num_taps, T x, T y, kfr::expression_pointer<T> & window, bool normalize = true) : FIRFilter<T>() {            
            bandpass(x,y,window,normalize);
        }
        ~FIRBandpassFilter() {
            
        }        
    };
    
    template<typename T>
    struct FIRBandstopFilter : public FIRFilter<T>
    {       
        FIRBandstopFilter(size_t num_taps, T x, T y, kfr::expression_pointer<T> & window, bool normalize = true) : FIRFilter<T>() {
            bandstop(x,y,window,normalize);
        }
        ~FIRBandstopFilter() {
            
        }
        
    };
    template<typename T>
    struct FIRHighpassFilter : public FIRFilter<T>
    {        
        FIRHighpassFilter(size_t num_taps, T x, kfr::expression_pointer<T> & window, bool normalize = true) : FIRFilter<T>() {
            highpass(x,window,normalize);
        }
        ~FIRHighpassFilter() {
            
        }
    };
    template<typename T>
    struct FIRLowpassFilter : public FIRFilter<T>
    {        
        FIRLowpassFilter(size_t num_taps, T x, kfr::expression_pointer<T> & window, bool normalize = true) : FIRFilter<T>() {
            lowpass(x,window,normalize);
        }
        ~FIRLowpassFilter() {
            
        }        
    };
    template<typename T>
    void filter(FIRFilter<T> * filter, size_t n, T * in, T * out)
    {
		filter->apply(n,in,out);
    }
    
    // parallel and serial fir are easier than iir as it is just a single op
    // the iir must manipulate the numerator and denominator to add and multiply in the frequency response
    template<typename T>
    // H = H1 + H2 is low shelf = G*lp + hp
    kfr::univector<T> parallel_fir(const kfr::univector<T> & a, const kfr::univector<T> & b) {		
		return a + b;
	}
	template<typename T>
    // H = H1 * H2
    // ie a cascade of 2 biquads = Hb1 * Hb2 = polynomial multiply in frequency domain
    kfr::univector<T> serial_fir(const kfr::univector<T> & a, const kfr::univector<T> & b) {		
		return a * b;
	}
}    
