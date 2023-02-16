#pragma once

namespace KfrDSP1
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
        void reset() { filter->reset(); }

    };
}