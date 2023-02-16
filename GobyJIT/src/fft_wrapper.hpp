#ifndef FFT_WRAPPER_
#define FFT_WRAPPER_

#include <fftw3.h>

/// Basic wrapper for FFT operations form real data
/*
 * Holds two buffers of the needed size. To transform data, copy in into the right buffer
 * and call the corresponding function. Beware that c2r changes input data.
 * 
 * The wrapper remembers the current "state" of the data,
 * so consecutive calls to the same transform transforms the data only once.
 * This also ensures correct rescaling of the data.
 * Rescaling occurs only on c2r transform.
 */
class FFTwrapper{
public:
    ///Default constructor that wraps transforms of size n
    /*
     * Expects the first transform to be from real to complex
     */
    FFTwrapper(int n);
    ///Constructor for both directions of first transform
    /*
     * The is_complex_state parameter indicates wheter the initial data will be complex or real.
     * Beware that the data are always rescaled after the c2r transform, so starting from complex
     * may introduce some scaling error
     */
    FFTwrapper(int n, bool is_complex_state);
    
    //Delete copying
    FFTwrapper(const FFTwrapper &fft) = delete;
    FFTwrapper &operator=(const FFTwrapper &fft) = delete;

    //Move constructors and desctructors
    FFTwrapper(FFTwrapper &&fft);
    FFTwrapper &operator=(FFTwrapper &&fft);

    ~FFTwrapper();
    ///Returns the pointer to the real time-donain data
    double *get_real_ptr();
    ///Returns the length of the real buffer
    size_t get_real_c() const;
    ///Returns the pointer to the complex frequency-domain data
    fftw_complex *get_complex_ptr();
    ///Returns the length of the complex buffer
    size_t get_complex_c() const;
    ///Transforms data from real buffer into complex buffer
    void r2c_transform();
    ///Transforms data from complex buffer into real buffer
    void c2r_transform();

    //Can be changed manually to avoid unnecesarry transforms
    bool is_complex_state_ = false;

private:
    int n_r_, n_c_;
    double *real_;
    fftw_complex *complex_;
    fftw_plan r2c_, c2r_;
};

#endif

