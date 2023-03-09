#pragma once
#include <iostream>
namespace AudioDSP
{
    #define SQR(x) ((x)*(x))

    /** Window types */
    typedef enum
    {
    win_ones,
    win_rectangle,
    win_hamming,
    win_hanning,
    win_hanningz,
    win_blackman,
    win_blackman_harris,
    win_gaussian,
    win_welch,
    win_parzen,
    win_default = win_hanningz,
    } window_type;

    template<typename T>
    T unwrap2pi (T phase)
    {
        /* mod(phase+pi,-2pi)+pi */
        return phase + 2.0*M_PI * (1. + std::floor(-(phase + M_PI) / 2.0*M_PI));
    }


    template<typename T>
    T mean (const size_t n, T * s)
    {
        T tmp = 0.0;
        #if defined(HAVE_INTEL_IPP)
            ippsMean(s, (int)n, &tmp);
            return tmp;
        #elif defined(HAVE_ACCELERATE)
            vDSP_meanv(s, 1, &tmp, n);
            return tmp;
        #else
            size_t j;
            #pragma omp simd
            for (j = 0; j < n; j++) {
                tmp += s[j];
            }
            return tmp / (T)(n);
        #endif
    }

    template<typename T>
    T sum (const size_t n, T * s)
    {
        T tmp = 0.0;
        #if defined(HAVE_INTEL_IPP)
            ippsSum(s, (int)n, &tmp);
        #elif defined(HAVE_ACCELERATE)
            vDSP_sve(s, 1, &tmp, n);
        #else
            size_t j;
            #pragma omp simd
            for (j = 0; j < n; j++) {
                tmp += s[j];
            }
        #endif
        return tmp;
    }

    template<typename T>
    T max (const size_t n, T * s)
    {
        #if defined(HAVE_INTEL_IPP)
            T tmp = 0.;
            ippsMax( s, (int)n, &tmp);
        #elif defined(HAVE_ACCELERATE)
            T tmp = 0.;
            vDSP_maxv( s, 1, &tmp, n );
        #else
            size_t j;
            T tmp = s[0];
            for (j = 1; j < n; j++) {
                tmp = (tmp > s[j]) ? tmp : s[j];
            }
        #endif
        return tmp;
    }

    template<typename T>
    T min (const size_t n, T * s)
    {
        #if defined(HAVE_INTEL_IPP)
            T tmp = 0.;
            ippsMin(s, (int)n, &tmp);
        #elif defined(HAVE_ACCELERATE)
            T tmp = 0.;
            vDSP_minv(s, 1, &tmp, n);
        #else
            size_t j;
            T tmp = s[0];
            #pragma omp simd
            for (j = 1; j < n; j++) {
                tmp = (tmp < s[j]) ? tmp : s[j];
            }
        #endif
        return tmp;
    }
   
    template<typename T>
    void swap(T & a, T & b) {
        T x = a;
        a = b;
        b = x;
    }
    template<typename T>
    void shift (size_t n, T * s)
    {
        size_t half = n / 2, start = half, j;
        // if length is odd, middle element is moved to the end
        if (2 * half < n) start ++;
        #ifndef HAVE_BLAS
            for (j = 0; j < half; j++) {
                swap(s[j], s[j + start]);
            }
        #else
            cblas_swap(half, s, 1, s + start, 1);
        #endif
        if (start != half) {
            for (j = 0; j < half; j++) {
            swap(s[j + start - 1], s[j + start]);
            }
        }
    }

    template<typename T>
    void ishift (const size_t n, T * s)
    {
        size_t half = n / 2, start = half, j;
        // if length is odd, middle element is moved to the beginning
        if (2 * half < n) start ++;
        #ifndef HAVE_BLAS
            #pragma omp simd
            for (j = 0; j < half; j++) {
                swap(s[j], s[j + start]);
            }
        #else
            cblas_swap(half, s, 1, s + start, 1);
        #endif            
            if (start != half) {
                #pragma omp simd
                for (j = 0; j < half; j++) {
                    swap(s[half], s[j]);
                }
            }
        }

    
    template<typename T>
    void clamp(const size_t n, T * in, T absmax) {
        size_t i;
        absmax = fabs(absmax);
        #pragma omp simd  
        for (i = 0; i < n; i++) in[i] = std::clamp(in[i],-absmax,absmax);  
    }


    template<typename T>
    T level_lin (const size_t n, T * f)
    {
        T energy = 0.;
        #ifndef HAVE_BLAS
            size_t j;
            #pragma omp simd
            for (j = 0; j < n; j++) {
                energy += SQR (f[j]);
            }
        #else
            energy = cblas_dot(n, f, 1, f, 1);
        #endif
        return energy / n;
    }

    template<typename T>
    T local_hfc (const size_t n, T * v)
    {
        T hfc = 0.;
        size_t j;
        #pragma omp simd
        for (j = 0; j < n; j++) {
            hfc += (j + 1) * v[j];
        }
        return hfc;
        }

    template<typename T>
    void min_removal (const size_t n, T * v)
    {
        T v_min = min(v);
        v += -v_min;  
    }

    template<typename T>
    T alpha_norm (const size_t n, T * o, T alpha)
    {
        size_t j;
        T tmp = 0.;
        #pragma omp simd
        for (j = 0; j < n; j++) {
            tmp += std::pow(std::fabs(o[j]), alpha);
        }
        return std::pow(tmp / n, 1. / alpha);
    }

    template<typename T>
    void alpha_normalise (const size_t n, T * o, T alpha)
    {
        size_t j;
        T norm = alpha_norm (o, alpha);
        o /= norm;  
    }


    template<typename T>
    T median (const size_t n, T * input) {
        T * arr = input;
        size_t low, high ;
        size_t median;
        size_t middle, ll, hh;

        low = 0 ; high = n-1 ; median = (low + high) / 2;
        
        for (;;) {
            if (high <= low) /* One element only */
            return arr[median] ;

            if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                std::swap(arr[low], arr[high]) ;
            return arr[median] ;
            }

            /* Find median of low, middle and high items; swap into position low */
            middle = (low + high) / 2;
            if (arr[middle] > arr[high])    std::swap(arr[middle], arr[high]);
            if (arr[low]    > arr[high])    std::swap(arr[low],    arr[high]);
            if (arr[middle] > arr[low])     std::swap(arr[middle], arr[low]) ;

            /* Swap low item (now in position middle) into position (low+1) */
            std::swap(arr[middle], arr[low+1]) ;

            /* Nibble from each end towards middle, swapping items when stuck */
            ll = low + 1;
            hh = high;
            for (;;) {
            do ll++; while (arr[low] > arr[ll]) ;
            do hh--; while (arr[hh]  > arr[low]) ;

            if (hh < ll)
                break;

            std::swap(arr[ll], arr[hh]) ;
            }

            /* Swap middle item (in position low) back into correct position */
            std::swap(arr[low], arr[hh]) ;

            /* Re-set active partition */
            if (hh <= median)
            low = ll;
            if (hh >= median)
            high = hh - 1;
        }
    }


    template<typename T>
    T moving_thres (const size_t n, T * vec, T * tmpvec,
        size_t post, size_t pre, size_t pos)
    {
        size_t k;
        T *medar = (T *) tmpvec;
        size_t win_length = post + pre + 1;
        size_t length = n;
        /* post part of the buffer does not exist */
        if (pos < post + 1) {
            for (k = 0; k < post + 1 - pos; k++)
            medar[k] = 0.;            /* 0-padding at the beginning */
            #pragma omp simd
            for (k = post + 1 - pos; k < win_length; k++)
            medar[k] = vec[k + pos - post];
            /* the buffer is fully defined */
        } else if (pos + pre < length) {
            #pragma omp simd
            for (k = 0; k < win_length; k++)
            medar[k] = vec[k + pos - post];
            /* pre part of the buffer does not exist */
        } else {
            #pragma omp simd
            for (k = 0; k < length - pos + post; k++)
            medar[k] = vec[k + pos - post];
            #pragma omp simd
            for (k = length - pos + post; k < win_length; k++)
            medar[k] = 0.;            /* 0-padding at the end */
        }
        return median (tmpvec);
    }

    template<typename T>
    void sample_vector_adapt_thres(const size_t n, T * vec, T * tmp,
        size_t post, size_t pre) {
        size_t length = n, j;  
        for (j=0;j<length;j++) {
            vec[j] -= moving_thres(vec, tmp, post, pre, j);
        }
    }

    template<typename T>
    T quadratic_peak_pos (const size_t n, T * x, size_t pos) {
        T s0, s1, s2; 
        size_t x0, x2;
        T half = .5, two = 2.;
        if (pos == 0 || pos == n - 1) return pos;
        x0 = (pos < 1) ? pos : pos - 1;
        x2 = (pos + 1 < n) ? pos + 1 : pos;
        if (x0 == pos) return (x[pos] <= x[x2]) ? pos : x2;
        if (x2 == pos) return (x[pos] <= x[x0]) ? pos : x0;
        s0 = x[x0];
        s1 = x[pos];
        s2 = x[x2];
        return pos + half * (s0 - s2 ) / (s0 - two * s1 + s2);
    }


    template<typename T>
    T quadratic_peak_mag (const size_t n, T * x, T pos) {
        T x0, x1, x2;
        size_t index = (size_t)(pos - .5) + 1;
        if (pos >= n || pos < 0.) return 0.;
        if ((T)index == pos) return x[index];
        x0 = x[index - 1];
        x1 = x[index];
        x2 = x[index + 1];
        return x1 - .25 * (x0 - x2) * (pos - index);
    }

    template<typename T>
    size_t peakpick(const size_t n, T * onset, size_t pos) {
        size_t tmp=0;
        tmp = (onset[pos] > onset[pos-1]
            &&  onset[pos] > onset[pos+1]
            &&  onset[pos] > 0.);
        return tmp;
    }

    template<typename T>
    T quadfrac (T s0, T s1, T s2, T pf)
    {
        T tmp =
            s0 + (pf / 2.) * (pf * (s0 - 2. * s1 + s2) - 3. * s0 + 4. * s1 - s2);
        return tmp;
    }

    template<typename T>
    T freqtomidi (T freq)
    {
        T midi;
        if (freq < 2. || freq > 100000.) return 0.; // avoid nans and infs
        /* log(freq/A-2)/log(2) */
        midi = freq / 6.875;
        midi = std::log (midi) / 0.6931471805599453;
        midi *= 12;
        midi -= 3;
        return midi;
    }

    template<typename T>
    T miditofreq (T midi)
    {
        T freq;
        if (midi > 140.) return 0.; // avoid infs
        freq = (midi + 3.) / 12.;
        freq = EXP (freq * 0.6931471805599453);
        freq *= 6.875;
        return freq;
    }

    template<typename T>
    T bintofreq (T bin, T samplerate, T fftsize)
    {
        T freq = samplerate / fftsize;
        return freq * MAX(bin, 0);
    }

    template<typename T>
    T bintomidi (T bin, T samplerate, T fftsize)
    {
        T midi = bintofreq (bin, samplerate, fftsize);
        return freqtomidi (midi);
    }

    template<typename T>
    T freqtobin (T freq, T samplerate, T fftsize)
    {
        T bin = fftsize / samplerate;
        return MAX(freq, 0) * bin;
    }

    template<typename T>
    T miditobin (T midi, T samplerate, T fftsize)
    {
        T freq = miditofreq (midi);
        return freqtobin (freq, samplerate, fftsize);
    }


    size_t is_power_of_two (size_t a)
    {
        if ((a & (a - 1)) == 0) {
            return 1;
        } else {
            return 0;
        }
    }


    size_t next_power_of_two (size_t a)
    {
        size_t i = 1;
        while (i < a) i <<= 1;
            return i;
        }

    size_t power_of_two_order (size_t a)
    {
        int order = 0;
        int temp = next_power_of_two(a);
        while (temp >>= 1) {
            ++order;
        }
        return order;
    }

    template<typename T>
    T db_spl (const size_t n, T * o)
    {
        return 10. * LOG10 (level_lin (o));
    }

    template<typename T>
    size_t silence_detection (const size_t n, T * o, T threshold)
    {
        return (db_spl (o) < threshold);
    }

    template<typename T>
    T level_detection (const size_t n, T * o, T threshold)
    {
        T db_spl = db_spl (o);
        if (db_spl < threshold) {
            return 1.;
        } else {
            return db_spl;
        }
    }

    template<typename T>
    T zero_crossing_rate (size_t n, T * input)
    {
        size_t j;
        size_t zcr = 0;
        for (j = 1; j < n; j++) {
            // previous was strictly negative
            if (input[j - 1] < 0.) {
            // current is positive or null
            if (input[j] >= 0.) {
                zcr += 1;
            }
            // previous was positive or null
            } else {
            // current is strictly negative
            if (input[j] < 0.) {
                zcr += 1;
            }
            }
        }
        return zcr / (T) n;
    }

    template<typename T>
    void autocorr (const size_t n, const T * input, T * output)
    {
        size_t i, j, length = n;
        T *data, *acf;
        T tmp = 0;
        data = input;
        acf = output->data;
        #pragma omp simd
        for (i = 0; i < length; i++) {
            tmp = 0.;
            for (j = i; j < length; j++) {
                tmp += data[j - i] * data[j];
            }
            acf[i] = tmp / (T) (length - i);
        }
    }

    template<typename T>
    void create_window (size_t n, T * win, window_type wintype) {        
        T * w = win;
        size_t i, size =n;
    
        switch(wintype) {
        case win_ones:
            win.ones();      
            break;
        case win_rectangle:
            win.fill(.5);
            break;
        case win_hamming:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.54 - 0.46 * std::cos(2.0*M_PI * i / (size));
            break;
        case win_hanning:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.5 - (0.5 * std::cos(2.0*M_PI * i / (size)));
            break;
        case win_hanningz:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.5 * (1.0 - std::cos(2.0*M_PI * i / (size)));
            break;
        case win_blackman:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.42
                - 0.50 * std::cos(    2.0*M_PI*i/(size-1.0))
                + 0.08 * std::cos(2.0*2.0*M_PI*i/(size-1.0));
            break;
        case win_blackman_harris:
        #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.35875
                - 0.48829 * std::cos(    2.0*M_PI*i/(size-1.0))
                + 0.14128 * std::cos(2.0*2.0*M_PI*i/(size-1.0))
                - 0.01168 * std::cos(3.0*2.0*M_PI*i/(size-1.0));
            break;
        case win_gaussian:
            {
            T a, b, c = 0.5;
            size_t n;
            #pragma omp simd
            for (n = 0; n < size; n++)
            {
                a = (n-c*(size-1))/(SQR(c)*(size-1));
                b = -c*SQR(a);
                w[n] = std::exp(b);
            }
            }
            break;
        case win_welch:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 1.0 - SQR((2.*i-size)/(size+1.0));
            break;
        case win_parzen:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 1.0 - std::fabs((2.f*i-size)/(size+1.0f));
            break;
        default:
            break;
        }  
    }

    template<typename T>
    T hztomel (T freq)
    {
        const T lin_space = 3./200.;
        const T split_hz = 1000.;
        const T split_mel = split_hz * lin_space;
        const T log_space = 27./std::log(6400/1000.);
        if (freq < 0) {
            std::cerr << ("hztomel: input frequency should be >= 0\n");
            return 0;
        }
        if (freq < split_hz)
        {
            return freq * lin_space;
        } else {
            return split_mel + log_space * std::log (freq / split_hz);
        }

    }

    template<typename T>
    T meltohz (T mel)
    {
        const T lin_space = 200./3.;
        const T split_hz = 1000.;
        const T split_mel = split_hz / lin_space;
        const T logSpacing = std::pow(6400/1000., 1/27.);
        if (mel < 0) {
            std::cerr << "meltohz: input mel should be >= 0\n";
            return 0;
        }
        if (mel < split_mel) {
            return lin_space * mel;
        } else {
            return split_hz * std::pow(logSpacing, mel - split_mel);
        }
    }

    template<typename T>
    T hztomel_htk (T freq)
    {
        const T split_hz = 700.;
        const T log_space = 1127.;
        if (freq < 0) {
            std::cerr << "hztomel_htk: input frequency should be >= 0\n";
            return 0;
        }
        return log_space * std::log(1 + freq / split_hz);
    }

    template<typename T>
    T meltohz_htk (T mel)
    {
        const T split_hz = 700.;
        const T log_space = 1./1127.;
        if (mel < 0) {
            std::cerr << "meltohz_htk: input frequency should be >= 0\n";
            return 0;
        }
        return split_hz * (std::exp( mel * log_space) - 1.);
    }
}