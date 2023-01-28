#pragma once

namespace KfrDSP1
{
    /*
    template<typename T> using CDFTPlan = DSP::DFTPlan<T>;
    template<typename T> using RDFT = DSP::DFTRealPlan<T>;
    template<typename T> using DCT = DSP::DCTPlan<T>;
    */
    ///////////////////////////////////////////////////////////////
    // DFT/DCT
    ///////////////////////////////////////////////////////////////

    template<typename T>
    kfr::univector<kfr::complex<T>> dft_forward(kfr::univector<kfr::complex<T>> input)
    {        
        return DSP::run_dft(input);        
    }
    template<typename T>
    kfr::univector<kfr::complex<T>> dft_inverse(kfr::univector<kfr::complex<T>> input)
    {        
        return  DSP::run_idft(input);     
    }
    template<typename T>
    kfr::univector<kfr::complex<T>> real_dft_forward(kfr::univector<T> input)
    {        
        return run_realdft(input);     
    }
    template<typename T>
    kfr::univector<T> real_dft_inverse(kfr::univector<kfr::complex<T>> in)
    {        
        return run_irealdft(in);        
    }

    template<typename T>
    kfr::univector<T> dct_forward(kfr::univector<T> in)
    {        
        DSP::DCTPlan<T> dct(in.size());
        kfr::univector<T> out(in.size());
        dct.execute(out,in,false);
        return out;
    }
    template<typename T>
    kfr::univector<T> dct_inverse(kfr::univector<T> in)
    {        
        DSP::DCTPlan<T> dct(in.size());
        kfr::univector<T> out(in.size());
        dct.execute(out,in,true);
        return out;
    }
}        