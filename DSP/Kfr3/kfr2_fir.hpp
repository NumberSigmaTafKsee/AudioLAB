#pragma once
#include "KfrWindows.hpp"
#include "KfrFIRFilter.hpp"
#include "KfrFIRBandpassFilter.hpp"
#include "KfrFIRBandstopFilter.hpp"
#include "KfrFIRLowpassFilter.hpp"
#include "KfrFIRHighpassFilter.hpp"
   
namespace KfrDSP1
{
    ///////////////////////////////////////////////////////////////
    // FIR
    ///////////////////////////////////////////////////////////////
    template<typename T>
    kfr::univector<T> fir_lowpass(kfr::univector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);        
        return DSP::fir_lowpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    kfr::univector<T> fir_highpass(kfr::univector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);        
        return DSP::fir_highpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    kfr::univector<T> fir_bandpass(kfr::univector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);
        return DSP::fir_bandpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    kfr::univector<T> fir_bandstop(kfr::univector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);
        return DSP::fir_bandstop(in,num_taps,window,normalize);        
    }
}