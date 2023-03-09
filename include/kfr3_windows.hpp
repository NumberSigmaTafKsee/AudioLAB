#pragma once

namespace kfr3
{
    ///////////////////////////////////////////////////////////////
    // Windows
    ///////////////////////////////////////////////////////////////

    enum WindowType 
    {
        rectangular     = 1,
        triangular      = 2,
        bartlett        = 3,
        cosine          = 4,
        hann            = 5,
        bartlett_hann   = 6,
        hamming         = 7,
        bohman          = 8,
        blackman        = 9,
        blackman_harris = 10,
        kaiser          = 11,
        flattop         = 12,
        gaussian        = 13,
        lanczos         = 14,
    };

    template<typename T>
    kfr::expression_pointer<T> make_window(WindowType type, size_t s, const T x = -1)
    {
        switch(type)
        {
        case hann: return DSP::make_window_hann_ptr<T>(s);
        case hamming: return DSP::make_window_hamming_ptr<T>(s);
        case blackman: return DSP::make_window_blackman_ptr<T>(s);
        case blackman_harris: return DSP::make_window_blackman_ptr<T>(s);
        case gaussian: return DSP::make_window_gaussian_ptr<T>(s);
        case triangular: return DSP::make_window_triangular_ptr<T>(s);
        case bartlett: return DSP::make_window_bartlett_ptr<T>(s);
        case cosine: return DSP::make_window_cosine_ptr<T>(s);
        case bartlett_hann: return DSP::make_window_bartlett_hann_ptr<T>(s);
        case bohman: return DSP::make_window_bohman_ptr<T>(s);
        case lanczos: return DSP::make_window_lanczos_ptr<T>(s);
        case flattop: return DSP::make_window_flattop_ptr<T>(s);
        case rectangular: return DSP::make_window_rectangular_ptr<T>(s);
        case kaiser: return DSP::make_window_kaiser_ptr<T>(s);
        }
    }
    template<typename T>
    kfr::univector<T> window(WindowType type,size_t s, const T x = -1) {
        return kfr::univector<T>(make_window<T>(type,s));
    }    
}        
