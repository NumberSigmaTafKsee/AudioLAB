#pragma once 

#define AUBIO_UNSTABLE 1
#include "aubio.h"
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <string>

namespace Aubio {

    typedef float Sample;

    template<typename T> 
    struct BufferBase {
        T *      buffer;
        size_t   len;

        BufferBase() { buffer = NULL; len = 0; }
        BufferBase(T * p, size_t n) : buffer(p),len(n) {} 
        BufferBase(size_t n) {
            buffer = (T*)calloc(n, sizeof(T));
            assert(buffer != NULL);
            len = n;
        }
        virtual ~BufferBase() {
            
        }

        T& operator[](size_t i) { return buffer[i]; }    
        T  __getitem(size_t i) { return buffer[i]; }
        void    __setitem(size_t i, const T s) { buffer[i] = s; }

        BufferBase<T>& operator = (const BufferBase<T> & s) {
            copy(s);
            return *this;
        }
        void copy(const BufferBase<T> & s) {
            len = s.b;
            if(buffer) free(buffer);
            buffer = (T*)calloc(len,sizeof(T));
            memcpy(buffer,s.buffer,len*sizeof(T));
        }
        void copy(const T * p, size_t len) {
            if(buffer) free(buffer);
            this->len = len;
            buffer = (T*)calloc(len,sizeof(T));
            memcpy(buffer,p,len*sizeof(T));
        }

    };

    template<typename T>
    struct SampleBuffer : public BufferBase<T> {
        
        SampleBuffer() : BufferBase<T>() {}
        SampleBuffer(size_t len) : BufferBase<T>(len) {}
        ~SampleBuffer() {
            if(this->buffer) free(this->buffer);
        }

        SampleBuffer<T>& operator = (const SampleBuffer<T> & s) {
            copy(s);
            return *this;
        }
    };

    template<typename T>
    struct Buffer : public BufferBase<T> 
    {
        Buffer() : BufferBase<T>() {}
        Buffer(T * p, size_t n) : BufferBase<T>(p,n) {}
        Buffer(size_t size) : BufferBase<T>(size) {} 

        Buffer<T>& operator = (const Buffer<T> & s) {
            copy(s);
            return *this;
        }
    };

    struct FVec 
    {
        fvec_t * pvec;

        FVec(const FVec & c) { 
            pvec = new_fvec(c.pvec->length);
            fvec_copy(c.pvec,pvec);
        }
        FVec(fvec_t * p) : pvec(p) {}
        FVec(size_t len) { pvec = new_fvec(len); assert(pvec != NULL); }
        ~FVec() { if(pvec) del_fvec(pvec); }

        FVec& operator = (const FVec & v) { 
            copy(v);
            return *this; 
        }
        void copy(const FVec & v) { 
            if(pvec) del_fvec(pvec);
            pvec = new_fvec(v.pvec->length);
            fvec_copy(v.pvec,pvec); 
        }

        void weighted_copy( const FVec & in, const FVec & w) { 
            copy(in);
            fvec_weighted_copy(in.pvec, w.pvec, pvec); 
        }
        
        Sample get_sample(uint_t position) { return fvec_get_sample(pvec, position); }
        void   set_sample(uint_t pos, Sample v) { fvec_set_sample(pvec,v, pos); }
        Sample& operator[](size_t index) { return pvec->data[index]; }
        Sample __getitem(size_t index) { return pvec->data[index]; }
        void   __setitem(size_t index, Sample value) { pvec->data[index] = value; }

        size_t size() const { return pvec->length; }

        void resize(size_t n) {
            if( n == pvec->length) return;
            del_fvec(pvec);
            pvec = new_fvec(n);
        }
        Buffer<Sample> get_data() { 
            Buffer<Sample> r(pvec->data,pvec->length);
            return r;
        }
        void set_data(const Buffer<Sample> & v) { 
            memcpy(pvec->data,v.buffer,v.len*sizeof(float));
        }
        void set_data(const SampleBuffer<Sample> & v) { 
            memcpy(pvec->data,v.buffer,v.len*sizeof(float));
        }

        void print() { fvec_print(pvec); }
        void setall(Sample v) { fvec_set_all(pvec, v); }
        void zeros() { fvec_zeros(pvec); }
        void ones() { fvec_ones(pvec); }
        void reverse() { fvec_rev(pvec); }
        void weight(FVec & v) { fvec_weight(pvec,v.pvec); }
        
        Sample zero_crossing_rate() { return aubio_zero_crossing_rate(pvec); }
        Sample level_lin() { return aubio_level_lin(pvec); }
        Sample db_spl() { return aubio_db_spl(pvec); }
        size_t silence_detection(Sample threshold) { return aubio_silence_detection(pvec, threshold); }
        Sample level_detection(Sample threshold) { return aubio_level_detection(pvec, threshold); }
        void clamp(Sample absmax) { fvec_clamp(pvec, absmax); }

        // can add more functions like tanh, sinh, cosh, etc
        void exp() { fvec_exp(pvec); }
        void cos() { fvec_cos(pvec); }
        void sin() { fvec_sin(pvec); }
        void abs() { fvec_abs(pvec); }
        void sqrt() { fvec_sqrt(pvec); }
        void log10() { fvec_log10(pvec); }
        void log() { fvec_log(pvec); }
        void floor() { fvec_floor(pvec); }
        void ceil() { fvec_ceil(pvec); }
        void round() { fvec_round(pvec); }
        void pow(Sample pow) { fvec_pow(pvec,pow); }
    };


    struct CVec
    {
        cvec_t * cvec;

        CVec(size_t len) { cvec = new_cvec(len); assert(cvec != NULL); }
        CVec(cvec_t * c) : cvec(c) {}
        ~CVec() { del_cvec(cvec); }

        CVec& operator = (const CVec & a) { copy(a); return *this; }
        
        void copy(const CVec & a) { 
            if(cvec) del_cvec(cvec);
            cvec = new_cvec(a.cvec->length);
            cvec_copy(a.cvec, cvec); 
        }

        void norm_zeros() { cvec_norm_zeros(cvec); }
        void norm_ones() { cvec_norm_ones(cvec); }
        void phas_set_all(Sample v) { cvec_phas_set_all(cvec,v); }
        void phas_zeros() { cvec_phas_zeros(cvec); }
        void phas_ones() { cvec_phas_ones(cvec); }
        void zeros() { cvec_zeros(cvec); }
        void logmag(Sample lambda) { cvec_logmag(cvec,lambda); }

        size_t size() const { return cvec->length; }
        
        void norm_set_sample( Sample v, size_t p) { cvec_norm_set_sample(cvec, v, p); }        
        void phas_set_sample( Sample v, size_t p) { cvec_phas_set_sample(cvec, v, p); }
        
        Sample norm_get_sample(size_t p) { return cvec_norm_get_sample(cvec, p); }
        Sample phas_get_sample(size_t p) { return cvec_phas_get_sample(cvec, p); }
        
        Buffer<Sample> norm_get_data()
        {
            Buffer<Sample> r(cvec->length);
            memcpy(r.buffer,cvec->norm, cvec->length*sizeof(Sample));
            return r;
        }
        Buffer<Sample> phas_get_data()
        {
            Buffer<Sample> r(cvec->length);
            memcpy(r.buffer, cvec->phas, cvec->length*sizeof(Sample));
            return r;
        }
    };

    struct LVec 
    {
        lvec_t * pvec;

        LVec(size_t len) { pvec = new_lvec(len); assert(pvec != NULL); }
        LVec(lvec_t * p) { pvec = p; }
        
        ~LVec() { del_lvec(pvec); }
        
        double get_sample(uint_t position) { return lvec_get_sample(pvec, position); }
        void   set_sample(uint_t pos, double v) { lvec_set_sample(pvec,v, pos); }
        
        double& operator[](size_t index) { return pvec->data[index]; }
        double __getitem(size_t index) { return pvec->data[index]; }
        void  __setitem(size_t index, double value) { pvec->data[index] = value; }

        Buffer<double> get_data() { 
            Buffer<double> r(pvec->length);
            memcpy(r.buffer,pvec->data,pvec->length*sizeof(double));
            return r;
        }

        void set_data(const Buffer<double> & v) { 
            if(pvec->length != v.len) {
                del_lvec(pvec);
                pvec = new_lvec(v.len);            
            }
            memcpy(pvec->data, v.buffer, sizeof(double) * v.len);
        }
        
        void set_data(const SampleBuffer<double> & v) { 
            if(pvec->length != v.len) {
                del_lvec(pvec);
                pvec = new_lvec(v.len);            
            }
            memcpy(pvec->data, v.buffer, sizeof(double) * v.len);
        }
            
        void print() { lvec_print(pvec); }
        void setall(double v) { lvec_set_all(pvec, v); }
        void zeros() { lvec_zeros(pvec); }
        void ones() { lvec_ones(pvec); }
    };


    struct FMat
    {
        fmat_t * m;

        FMat(size_t height, size_t width) { m = new_fmat(height,width); assert(m != NULL); }
        FMat(fmat_t * q) : m(q) {}
        ~FMat() { del_fmat(m); }

        Sample get_sample(size_t channel, size_t pos) { return fmat_get_sample(m, channel, pos); }           

        void set_sample(Sample data, size_t channel, size_t pos) { fmat_set_sample(m, data, channel, pos); }
        
        void get_channel(size_t channel, FVec & output) { fmat_get_channel(m,channel,output.pvec); }
        
        void set_channel_data(size_t channel, const Buffer<Sample> & buffer) {
            Sample * data = fmat_get_channel_data(m, channel);
            memcpy(data,buffer.buffer,m->length*sizeof(Sample));
        }
        void set_channel_data(size_t channel, const SampleBuffer<Sample> & buffer) {
            Sample * data = fmat_get_channel_data(m, channel);
            memcpy(data,buffer.buffer,m->length*sizeof(Sample));
        }
        Buffer<Sample> get_channel_data(size_t channel) { 
            Sample * p = fmat_get_channel_data(m,channel);        
            return Buffer<Sample>(p,m->length);
        }

        size_t size() const { return m->height*m->length; }
        size_t rows() const { return m->height; }
        size_t cols() const { return m->length; }
        void resize(size_t h, size_t l) {
            if(h == m->height && l == m->length) return;
            del_fmat(m);
            m = new_fmat(h,l);
        }
        void print() { fmat_print(m); }
        void set(Sample v) { fmat_set(m,v); }
        void zeros() { fmat_zeros(m); }
        void ones() { fmat_ones(m); }
        void reverse() { fmat_rev(m); }
        void weight(const FMat & w) { fmat_weight(m, w.m); }
        
        Sample operator()(size_t h, size_t w) { return m->data[h][w]; }

        FMat& operator = (const FMat & a) { copy(a); return *this; }
        
        void  copy(const FMat & a) { 
            if(m) del_fmat(m);
            m = new_fmat(a.m->height,a.m->length);
            fmat_copy(a.m, m); 
        }

        void  vecmul(const FVec & scale, const FVec & output) { fmat_vecmul(m, scale.pvec, output.pvec); }

        FVec operator * (const FVec & input)
        {
            FVec v(input);
            vecmul(input,v);
            return v;
        }
    };


    inline FVec*  new_window(char* name, size_t size) { FVec *v = new FVec(size); v->pvec = new_aubio_window(name, size); return v; }
    inline size_t set_window(FVec & window, char * window_type) { return fvec_set_window(window.pvec, window_type); }
        
    enum ResamplerType
    {
        SRC_SINC_BEST_QUALITY		= 0,
        SRC_SINC_MEDIUM_QUALITY		= 1,
        SRC_SINC_FASTEST			= 2,
        SRC_ZERO_ORDER_HOLD			= 3,
        SRC_LINEAR					= 4,
    } ;

    struct Resampler
    {
        aubio_resampler_t * resampler;

        Resampler(Sample ratio, ResamplerType type) {
            resampler = new_aubio_resampler(ratio, type);
            assert(resampler != NULL);
        }
        
        ~Resampler()
        {
            del_aubio_resampler(resampler);
        }

        void process(const FVec & input, FVec & output) { 
            output.resize(input.size());
            aubio_resampler_do(resampler, input.pvec, output.pvec); 
        }
        
    };




    struct Filter 
    {
        aubio_filter_t  *filter;

        Filter() { filter = NULL; }
        Filter(uint32_t order) { 
            filter = new_aubio_filter(order);
            assert(filter != NULL);
        }
        ~Filter()
        {
            if(filter) del_aubio_filter(filter);
        }

        void process(FVec & input)
        {
            aubio_filter_do(filter, input.pvec);
        }
        void do_outplace(const FVec & in, FVec & out)
        {
            out.resize(in.size());
            aubio_filter_do_outplace(filter, in.pvec, out.pvec);
        }
        void do_filtfilt(FVec & input, FVec & temp)
        {
            temp.resize(input.size());
            aubio_filter_do_filtfilt(filter, input.pvec, temp.pvec);
        }
        LVec& get_feedback(LVec & vec)
        {
            vec.pvec = aubio_filter_get_feedback(filter);
            return vec;
        }
        LVec&  get_feedforward(LVec & vec) { 
            vec.pvec = aubio_filter_get_feedforward(filter);            
            return vec;
        }
        uint32_t get_order() { return aubio_filter_get_order(filter); }
        uint32_t get_samplerate() { return aubio_filter_get_samplerate(filter); }
        void set_samplerate(uint32_t samplerate) { aubio_filter_set_samplerate(filter,samplerate); }
        void do_reset() { aubio_filter_do_reset(filter); }
    };



    struct BiQuad : public Filter
    {
        BiQuad(double b0, double b1, double b2, double a1, double a2)
        {
            filter = new_aubio_filter_biquad(b0,b1,b2,a1,a2);
            assert(filter != NULL);            
        }
        
        void set_biquad(double b0, double b1, double b2, double a1, double a2)
        {
            aubio_filter_set_biquad(filter,b0,b1,b2,a1,a2);
        }
        
    };

    struct AWeighting : public Filter
    {
        AWeighting(uint32_t samplerate)
        {
            filter = new_aubio_filter_a_weighting(samplerate);
        }

        size_t set_a_weighting(uint32_t samplerate)
        {
            return aubio_filter_set_a_weighting(filter, samplerate);
        }
    };

    struct CWeighting : public Filter
    {
        CWeighting(uint32_t samplerate)
        {
            filter = new_aubio_filter_c_weighting(samplerate);
        }

        size_t set_c_weighting(uint32_t samplerate)
        {
            return aubio_filter_set_c_weighting(filter, samplerate);
        }
    };

    struct FFT
    {
        aubio_fft_t  *fft;

        FFT(uint32_t size)
        {
            fft = new_aubio_fft(size);
            assert(fft != NULL);
        }
        ~FFT()
        {
            if(fft) del_aubio_fft(fft);
        }

        void forward(const FVec & input, CVec & spectrum)
        {
            aubio_fft_do(fft,input.pvec, spectrum.cvec);
        }

        void reverse(const CVec & spectrum, FVec & output)
        {
            aubio_fft_rdo(fft,spectrum.cvec, output.pvec);
        }

        void forward_complex(const FVec & input, FVec & compspec)
        {
            aubio_fft_do_complex(fft,input.pvec, compspec.pvec);
        }
        void inverse_complex(const FVec & compspec, FVec & output)
        {
            aubio_fft_rdo_complex(fft, compspec.pvec, output.pvec);
        }
        static void get_spectrum( const FVec & compspec, CVec & spectrum)
        {
            aubio_fft_get_spectrum(compspec.pvec,spectrum.cvec);
        }
        static void get_realimag( const CVec & spectrum, FVec & compspec)
        {
            aubio_fft_get_realimag(spectrum.cvec, compspec.pvec);
        }
        static void get_phase( const FVec & compspec, CVec & spectrum)
        {
            aubio_fft_get_phas( compspec.pvec, spectrum.cvec);
        }
        static void get_imaginary( const CVec & spectrum, FVec & compspec)
        {
            aubio_fft_get_imag(spectrum.cvec, compspec.pvec);
        }
        static void get_norm(const FVec & compspec, CVec & spectrum)
        {
            aubio_fft_get_norm(compspec.pvec,spectrum.cvec);
        }
        static void get_real(const CVec & spectrum, FVec & compspec)
        {
            aubio_fft_get_real(spectrum.cvec, compspec.pvec);
        }
    };

    struct DCT
    {
        aubio_dct_t * dct; 

        DCT(uint32_t size) { 
            dct = new_aubio_dct(size);
            assert(dct != NULL);
        }

        ~DCT() { if(dct != NULL) del_aubio_dct(dct); }

        void forward( const FVec & input, const FVec &output)
        {
            aubio_dct_do(dct,input.pvec, output.pvec);
        }
        void reverse( const FVec & input, const FVec &output)
        {
            aubio_dct_rdo(dct,input.pvec, output.pvec);
        }
    };

    struct PhaseVocoder 
    {
        aubio_pvoc_t * pvoc; 

        PhaseVocoder(uint32_t win_s, uint32_t hop_s)
        {
            pvoc = new_aubio_pvoc(win_s,hop_s);
            assert(pvoc != NULL);
        }
        ~PhaseVocoder()
        {
            if(pvoc) del_aubio_pvoc(pvoc);
        }

        void forward(const FVec & in, CVec& fftgrain)
        {
            aubio_pvoc_do(pvoc,in.pvec, fftgrain.cvec);
        }
        void reverse(const CVec & fftgrain, FVec& out)
        {
            aubio_pvoc_rdo(pvoc,fftgrain.cvec,out.pvec);
        }
        size_t get_win() { return aubio_pvoc_get_win(pvoc); }
        size_t get_hop() { return aubio_pvoc_get_hop(pvoc); }

        size_t set_window(const char * window_type)
        {
            return aubio_pvoc_set_window(pvoc, window_type);
        }
    };

    struct FilterBank
    {
        aubio_filterbank_t * bank;

        FilterBank(uint32_t nfilters, uint32_t win_s)
        {
            bank = new_aubio_filterbank(nfilters, win_s);
            assert(bank != NULL);
        }
        ~FilterBank()
        {
            if(bank != NULL) del_aubio_filterbank(bank);
        }
        void process(const CVec & in, FVec & out)
        {
            out.resize(in.size());
            aubio_filterbank_do(bank,in.cvec, out.pvec);
        }
        FMat get_coeffs() 
        {
            return FMat(aubio_filterbank_get_coeffs(bank));
        }
        uint32_t set_coeffs(FMat & m)
        {
            return aubio_filterbank_set_coeffs(bank, m.m);
        }
        uint32_t set_norm(double norm) { return aubio_filterbank_set_norm(bank,norm); }
        uint32_t set_power(double power) { return aubio_filterbank_set_power(bank,power); }
        double get_power() { return aubio_filterbank_get_power(bank); }        
    };

    struct MelFilterBank
    {
        aubio_filterbank_t * fb;

        MelFilterBank(uint32_t nfilters, uint32_t win_s)
        {
            fb = new_aubio_filterbank(nfilters, win_s);
            assert(fb != NULL);
        }
        ~MelFilterBank() 
        {
            if(fb != NULL) del_aubio_filterbank(fb);
        }
        size_t set_triangle_bands(const FVec & freqs, double samplerate)
        {
            return aubio_filterbank_set_triangle_bands(fb, freqs.pvec, samplerate);
        }
        size_t set_mel_coeffs_slaney( double samplerate)
        {
            return aubio_filterbank_set_mel_coeffs_slaney(fb,samplerate);
        }
        size_t set_mel_coeffs( double samplerate, double fmin, double fmax)
        {
            return aubio_filterbank_set_mel_coeffs(fb,samplerate,fmin,fmax);
        }
        size_t set_mel_coeffs_htk( double samplerate, double fmin, double fmax)
        {
            return aubio_filterbank_set_mel_coeffs_htk(fb,samplerate,fmin,fmax);
        }
    };


    struct MFCC
    {
        aubio_mfcc_t * mfcc; 

        MFCC(size_t buf_size, uint32_t nfilters, size_t size_coeffs, uint32_t samplerate)
        {
            mfcc = new_aubio_mfcc(buf_size,nfilters,size_coeffs, samplerate);
            assert(mfcc != NULL);
        }
        ~MFCC()
        {
            if(mfcc != NULL) del_aubio_mfcc(mfcc);
        }
        void process(const CVec &in, FVec & out)
        {
            out.resize(in.size());
            aubio_mfcc_do(mfcc, in.cvec, out.pvec);
        }
        size_t set_power(double power) { return aubio_mfcc_set_power(mfcc, power); }
        double get_power() { return aubio_mfcc_get_power(mfcc); }
        size_t set_scale(double power) { return aubio_mfcc_set_scale(mfcc, power); }
        double get_scale() { return aubio_mfcc_get_scale(mfcc); }

        size_t set_mel_coeffs(double fmin, double fmax)
        {
            return aubio_mfcc_set_mel_coeffs(mfcc,fmin, fmax);
        }
        size_t set_mel_coeffs_htk(double fmin, double fmax)
        {
            return aubio_mfcc_set_mel_coeffs_htk(mfcc,fmin, fmax);
        }
        size_t set_mel_coeffs_slaney()
        {
            return aubio_mfcc_set_mel_coeffs_slaney(mfcc);
        }
    };


    struct AWhitening
    {

        aubio_spectral_whitening_t * w;

        AWhitening(size_t buf_size, uint32_t hop_size, uint32_t samplerate)
        {
            w = new_aubio_spectral_whitening(buf_size,hop_size, samplerate);
            assert(w != NULL);
        }
        ~AWhitening() 
        {
            if(w != NULL) del_aubio_spectral_whitening(w);
        }

        void reset() { aubio_spectral_whitening_reset(w); }
        void set_relax_time(double relax_time )
        {
            aubio_spectral_whitening_set_relax_time(w,relax_time);
        }
        double get_relax_time( )
        {
            return aubio_spectral_whitening_get_relax_time(w);
        }
        void set_floor(double floor )
        {
            aubio_spectral_whitening_set_relax_time(w,floor);
        }
        double get_floor( )
        {
            return aubio_spectral_whitening_get_floor(w);
        }
    };



    struct TSS 
    {
        aubio_tss_t *tss;

        TSS(size_t buf_size, uint32_t hop_size) { 
            tss = new_aubio_tss(buf_size,hop_size);
            assert(tss != NULL);
        }
        ~TSS() { if(tss != NULL) del_aubio_tss(tss); }

        void process(const CVec & in, CVec & out, CVec & stead)
        {
            aubio_tss_do(tss,in.cvec,out.cvec,stead.cvec);
        }
        size_t set_threshold(double thrs)
        {
            return aubio_tss_set_threshold(tss,thrs);
        }
        size_t set_alpha(double alpha)
        {
            return aubio_tss_set_alpha(tss,alpha);
        }
        size_t set_beta(double beta)
        {
            return aubio_tss_set_beta(tss,beta);
        }        
    };




    struct SpecDesc
    {
        aubio_specdesc_t * sd;

        SpecDesc(const char * method, size_t buf_size)
        {
            sd = new_aubio_specdesc(method, buf_size);
            assert(sd != NULL);
        }
        ~SpecDesc() 
        {
            if(sd != NULL) del_aubio_specdesc(sd);
        }
        void process(const CVec& fftgrain, FVec & desc)
        {
            aubio_specdesc_do(sd,fftgrain.cvec, desc.pvec);
        }
    };


    struct Pitch
    {
        aubio_pitch_t * p;

        Pitch(const char * method, size_t buf_size, uint32_t hop_size, uint32_t samplerate)
        {
            p = new_aubio_pitch(method, buf_size, hop_size, samplerate);
            assert(p != NULL);
        }
        ~Pitch()
        {
            if(p != NULL) del_aubio_pitch(p);
        }

        void process(const FVec  &in, FVec &out)
        {
            out.resize(in.size());
            aubio_pitch_do(p,in.pvec,out.pvec);
        }
        size_t set_tolerance(double tol)
        {
            return aubio_pitch_set_tolerance(p,tol);
        }
        double get_tolerance()
        {
            return aubio_pitch_get_tolerance(p);
        }
        size_t set_unit(const char * mode)
        {
            return aubio_pitch_set_unit(p,mode);
        }
        size_t set_silence(const Pitch& silence, double s)
        {
            return aubio_pitch_set_silence(silence.p,s);
        }
        double get_silence()
        {
            return aubio_pitch_get_silence(p);
        }
        double get_confidence()
        {
            return aubio_pitch_get_confidence(p);
        }
    };

    struct PitchFComb
    {
        aubio_pitchfcomb_t * pitch;

        PitchFComb(uint32_t buf_size, uint32_t hop_size) {
            pitch = new_aubio_pitchfcomb(buf_size, hop_size);
            assert(pitch != NULL);
        }
        ~PitchFComb() {
            if(pitch) {
                del_aubio_pitchfcomb(pitch);
            }
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchfcomb_do(pitch, input.pvec,output.pvec);
        }
    };

    struct PitchMComb
    {
        aubio_pitchmcomb_t * pitch;

        PitchMComb(uint32_t buf_size, uint32_t hop_size) {
            pitch = new_aubio_pitchmcomb(buf_size, hop_size);
            assert(pitch != NULL);
        }
        ~PitchMComb() {
            if(pitch) {
                del_aubio_pitchmcomb(pitch);
            }
        }

        void process(const CVec & input, FVec & output) {
            aubio_pitchmcomb_do(pitch,input.cvec,output.pvec);
        }
    };

    struct PitchSchmitt
    {
        aubio_pitchschmitt_t * pitch;

        PitchSchmitt(uint32_t buf_size) {
            pitch = new_aubio_pitchschmitt(buf_size);
            assert(pitch != NULL);
        }
        ~PitchSchmitt() {
            if(pitch) {
                del_aubio_pitchschmitt(pitch);
            }
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchschmitt_do(pitch,input.pvec,output.pvec);
        }
    };

    struct PitchSpecACF
    {
        aubio_pitchspecacf_t * pitch;

        PitchSpecACF(uint32_t buf_size) {
            pitch = new_aubio_pitchspecacf(buf_size);
            assert(pitch != NULL);
        }
        ~PitchSpecACF() {
            if(pitch) {
                del_aubio_pitchspecacf(pitch);
            }
        }

        Sample get_tolerance() {
            return aubio_pitchspecacf_get_tolerance(pitch);
        }
        uint32_t set_tolerance(Sample tol) {
            return aubio_pitchspecacf_set_tolerance(pitch,tol);
        }

        Sample get_confidence() {
            return aubio_pitchspecacf_get_confidence(pitch);
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchspecacf_do(pitch,input.pvec,output.pvec);
        }
    };

    struct PitchYin
    {
        aubio_pitchyin_t * pitch;

        PitchYin(uint32_t buf_size) {
            pitch = new_aubio_pitchyin(buf_size);
            assert(pitch != NULL);
        }
        ~PitchYin() {
            if(pitch) {
                del_aubio_pitchyin(pitch);
            }
        }

        Sample get_tolerance() {
            return aubio_pitchyin_get_tolerance(pitch);
        }
        uint32_t set_tolerance(Sample tol) {
            return aubio_pitchyin_set_tolerance(pitch,tol);
        }

        Sample get_confidence() {
            return aubio_pitchyin_get_confidence(pitch);
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchyin_do(pitch,input.pvec,output.pvec);
        }
    };

    struct PitchYinFast
    {
        aubio_pitchyinfast_t * pitch;

        PitchYinFast(uint32_t buf_size) {
            pitch = new_aubio_pitchyinfast(buf_size);
            assert(pitch != NULL);
        }
        ~PitchYinFast() {
            if(pitch) {
                del_aubio_pitchyinfast(pitch);
            }
        }

        Sample get_tolerance() {
            return aubio_pitchyinfast_get_tolerance(pitch);
        }
        uint32_t set_tolerance(Sample tol) {
            return aubio_pitchyinfast_set_tolerance(pitch,tol);
        }

        Sample get_confidence() {
            return aubio_pitchyinfast_get_confidence(pitch);
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchyinfast_do(pitch,input.pvec,output.pvec);
        }
    };

    struct PitchYinFFT
    {
        aubio_pitchyinfft_t * pitch;

        PitchYinFFT(uint32_t sr, uint32_t buf_size) {
            pitch = new_aubio_pitchyinfft(sr,buf_size);
            assert(pitch != NULL);
        }
        ~PitchYinFFT() {
            if(pitch) {
                del_aubio_pitchyinfft(pitch);
            }
        }

        Sample get_tolerance() {
            return aubio_pitchyinfft_get_tolerance(pitch);
        }
        uint32_t set_tolerance(Sample tol) {
            return aubio_pitchyinfft_set_tolerance(pitch,tol);
        }

        Sample get_confidence() {
            return aubio_pitchyinfft_get_confidence(pitch);
        }

        void process(const FVec & input, FVec & output) {
            aubio_pitchyinfft_do(pitch,input.pvec,output.pvec);
        }        
    };


    struct Tempo
    {
        aubio_tempo_t * tempo;

        Tempo(const char * method, size_t buf_size, uint32_t hop_size, uint32_t samplerate)
        {
            tempo = new_aubio_tempo(method,buf_size,hop_size,samplerate);
            assert(tempo != NULL);
        }
        ~Tempo()
        {
            if(tempo != NULL) del_aubio_tempo(tempo);
        }
        void process(const FVec &in, FVec & out)
        {
            out.resize(in.size());
            aubio_tempo_do(tempo, in.pvec,out.pvec);
        }
        size_t get_last()
        {
            return aubio_tempo_get_last(tempo);
        }
        double get_last_s() 
        {
            return aubio_tempo_get_last_s(tempo);
        }
        double get_last_ms() 
        {
            return aubio_tempo_get_last_ms(tempo);
        }
        size_t set_silence(double silence)
        {
            return aubio_tempo_set_silence(tempo,silence);
        }
        double get_silence()
        {
            return aubio_tempo_get_silence(tempo);
        }
        size_t set_threshold(double threshold)
        {
            return aubio_tempo_set_threshold(tempo,threshold);
        }
        double get_threshold() 
        {
            return aubio_tempo_get_threshold(tempo);
        }
        double get_period() 
        {
            return aubio_tempo_get_period(tempo);
        }
        double get_period_s() 
        {
            return aubio_tempo_get_period_s(tempo);
        }
        double get_bpm() 
        {
            return aubio_tempo_get_bpm(tempo);
        }
        double get_confidence() 
        {
            return aubio_tempo_get_confidence(tempo);        
        }
        size_t get_set_tatum_signature(uint_t x) 
        {
            return aubio_tempo_set_tatum_signature(tempo,x);
        }
        size_t was_tatum()
        {
            return aubio_tempo_was_tatum(tempo); 
        }
        double get_last_tatum() 
        {
            return aubio_tempo_get_last_tatum(tempo);
        }
        double get_delay()
        {
            return aubio_tempo_get_delay(tempo);
        }
        double get_delay_s()
        {
            return aubio_tempo_get_delay_s(tempo);
        }
        double get_delay_ms()
        {
            return aubio_tempo_get_delay_ms(tempo);
        }
        size_t set_delay(int delay)
        {
            return aubio_tempo_set_delay(tempo,delay);
        }
        size_t set_delay_s(int delay)
        {
            return aubio_tempo_set_delay_s(tempo,delay);
        }
        size_t set_delay_ms(int delay)
        {
            return aubio_tempo_set_delay_ms(tempo,delay);
        }
    };

    struct BeatTrack 
    {
        aubio_beattracking_t * bt;

        BeatTrack(size_t winlen, uint32_t hop_size, uint32_t samplerate)
        {
            bt = new_aubio_beattracking(winlen,hop_size,samplerate);
            assert(bt != NULL);
        }
        ~BeatTrack()
        {
            if(bt != NULL) del_aubio_beattracking(bt);
        }
        void process(const FVec & dfframes, FVec & out)
        {
            out.resize(dfframes.size());
            aubio_beattracking_do(bt,dfframes.pvec,out.pvec);
        }  
        double get_period()
        {
            return aubio_beattracking_get_period(bt);
        }
        double get_period_s()
        {
            return aubio_beattracking_get_period_s(bt);
        }
        double get_bpm()
        {
            return aubio_beattracking_get_bpm(bt);
        }
        double get_confidence()
        {
            return aubio_beattracking_get_confidence(bt);
        }
    };

    struct Sink
    {
        aubio_sink_t * sink;

        Sink(const char * uri, uint32_t samplerate)
        {
            sink = new_aubio_sink(uri,samplerate);
            assert(sink != NULL);
        }
        ~Sink()
        {
            if(sink != NULL) del_aubio_sink(sink);
        }
        size_t preset_samplerate(size_t samplerate)
        {
            return aubio_sink_preset_samplerate(sink,samplerate);
        }
        size_t get_samplerate() 
        {
            return aubio_sink_get_samplerate(sink);
        }
        size_t get_channels() 
        {
            return aubio_sink_get_channels(sink);
        }
        void process(FVec & write_data, size_t write)
        {
            aubio_sink_do(sink,write_data.pvec, write);
        }
        void process_multi(FMat & write_data,size_t write)
        {
            aubio_sink_do_multi(sink,write_data.m,write);
        }
        size_t close() 
        {
            return aubio_sink_close(sink);
        }
    };

    enum FormatType {
        FORMAT_WAV,
        FORMAT_AIFF,
        FORMAT_FLAC,
        FORMAT_OGG
    };

    struct SinkSoundFile 
    {
        aubio_sink_sndfile_t * file;

        SinkSoundFile(const char * uri, uint32_t samplerate) {
            file = new_aubio_sink_sndfile(uri,samplerate);
            assert(file != NULL);
        }
        ~SinkSoundFile() {
            if(file) aubio_sink_sndfile_close(file);
            if(file) del_aubio_sink_sndfile(file);
        }

        uint32_t preset_samplerate(uint32_t samplerate) {
            return aubio_sink_sndfile_preset_samplerate(file,samplerate);
        }
        uint32_t preset_channels(uint32_t channels) {
            return aubio_sink_sndfile_preset_channels(file,channels);
        }
        
        uint32_t get_samplerate() {
            return aubio_sink_sndfile_get_samplerate(file);
        }
        uint32_t get_channels() {
            return aubio_sink_sndfile_get_channels(file);
        }
        void process(FVec & write_data, uint32_t write) {
            aubio_sink_sndfile_do(file,write_data.pvec,write);
        }
        void do_multi(FMat & write_data, uint32_t write) {
            aubio_sink_sndfile_do_multi(file,write_data.m, write);
        }
    };

    struct SinkWavWrite {
        aubio_sink_wavwrite_t * wav;

        SinkWavWrite(const char* uri, uint32_t samplerate) {
            wav = new_aubio_sink_wavwrite(uri,samplerate);
            assert(wav != NULL);
        }
        ~SinkWavWrite() {
            if(wav) {
                aubio_sink_wavwrite_close(wav);
                del_aubio_sink_wavwrite(wav);
            }
        }

        uint32_t preset_samplerate(uint32_t samplerate) {
            return aubio_sink_wavwrite_preset_samplerate(wav,samplerate);
        }
        uint32_t preset_channels(uint32_t channels) {
            return aubio_sink_wavwrite_preset_channels(wav,channels);
        }

        uint32_t get_samplerate() {
            return aubio_sink_wavwrite_get_samplerate(wav);
        }
        uint32_t get_channels() {
            return aubio_sink_wavwrite_get_channels(wav);
        }
        void process(FVec & write_data, uint32_t write) {
            aubio_sink_wavwrite_do(wav,write_data.pvec,write);
        }
        void do_multi(FMat & write_data, uint32_t write) {
            aubio_sink_wavwrite_do_multi(wav,write_data.m, write);
        }
    };

    struct Sampler
    {
        aubio_sampler_t * sampler;

        Sampler(uint32_t samplerate, size_t hop_size)
        {
            sampler = new_aubio_sampler(samplerate, hop_size);
            assert(sampler != NULL);
        }
        ~Sampler() {
            if(sampler) del_aubio_sampler(sampler);
        }

        uint32_t load(const char * uri) {
            return aubio_sampler_load(sampler,uri);
        }
        void process(const FVec & input, FVec & output) {
            output.resize(input.size());
            aubio_sampler_do(sampler,input.pvec,output.pvec);
        }
        void process_multi(const FMat & input, FMat & output) {
            output.resize(input.rows(),input.cols());
            aubio_sampler_do_multi(sampler,input.m,output.m);
        }

    };

    struct Source {
        aubio_source_t * source;

        Source(const char * uri, uint32_t samplerate, uint32_t hop_size) {
            source = new_aubio_source(uri,samplerate,hop_size);
            assert(source != NULL);
        }
        ~Source() {
            if(source) {
                aubio_source_close(source);
                del_aubio_source(source);
            }
        }

        uint32_t process(FVec & read_to) {
            uint32_t r = 0;
            aubio_source_do(source, read_to.pvec, &r );
            return r;
        }
        uint32_t process_multi(FMat & read_to) {
            uint32_t r = 0;
            aubio_source_do_multi(source, read_to.m, &r );
            return r;
        }

        uint32_t get_samplerate() {
            return aubio_source_get_samplerate(source);
        }
        uint32_t get_channels() {
            return aubio_source_get_channels(source);
        }

        uint32_t seek(uint32_t pos) {
            return aubio_source_seek(source,pos);            
        }

        uint32_t get_duration() {
            return aubio_source_get_duration(source);
        }
    };

    struct SourceSoundFile 
    {        
        aubio_source_sndfile_t * file; 

        SourceSoundFile(const char * uri, uint32_t samplerate, size_t hop_size)
        {
            file = new_aubio_source_sndfile(uri, samplerate, hop_size);
            assert(file != NULL);
        }
        ~SourceSoundFile()
        {
            if(file != NULL) del_aubio_source_sndfile(file);
        }

        size_t process(FVec & read_to)
        {
            uint_t read=0;
            aubio_source_sndfile_do(file, read_to.pvec, &read);
            return read;
        }
        size_t multi_process(FMat & read_to)
        {
            uint_t read=0;
            aubio_source_sndfile_do_multi(file, read_to.m, &read);
            return read;
        }
        uint32_t get_samplerate() { return aubio_source_sndfile_get_samplerate(file); }
        uint32_t get_channels() { return aubio_source_sndfile_get_channels(file); }
        uint64_t seek(size_t pos) { return aubio_source_sndfile_seek(file,pos); }
        uint32_t get_duration() { return aubio_source_sndfile_get_duration(file); }
        uint32_t close() { return aubio_source_sndfile_close(file); }

    };

    struct SourceWavFile 
    {        
        aubio_source_wavread_t * file; 

        SourceWavFile(const char * uri, uint32_t samplerate, size_t hop_size)
        {
            file = new_aubio_source_wavread(uri, samplerate, hop_size);
            assert(file != NULL);
        }
        ~SourceWavFile()
        {
            if(file != NULL) del_aubio_source_wavread(file);
        }

        size_t process(FVec & read_to)
        {
            uint_t read=0;
            aubio_source_wavread_do(file, read_to.pvec, &read);
            return read;
        }
        size_t multi_process(FMat & read_to)
        {
            uint_t read=0;
            aubio_source_wavread_do_multi(file, read_to.m, &read);
            return read;
        }
        uint32_t get_samplerate() { return aubio_source_wavread_get_samplerate(file); }
        uint32_t get_channels() { return aubio_source_wavread_get_channels(file); }
        uint64_t seek(size_t pos) { return aubio_source_wavread_seek(file,pos); }
        uint32_t get_duration() { return aubio_source_wavread_get_duration(file); }
        uint32_t close() { return aubio_source_wavread_close(file); }

    };

    struct Notes {
        aubio_notes_t * notes;

        Notes(uint32_t buf_size, uint32_t hop_size, uint32_t sample_rate) {
            notes = new_aubio_notes("default",buf_size,hop_size,sample_rate);
            assert(notes != NULL);
        }
        ~Notes() {
            if(notes) {
                del_aubio_notes(notes);
            }
        }

        void process(const FVec & input, FVec & output) {
            aubio_notes_do(notes,input.pvec,output.pvec);
        }
        uint32_t set_silence(Sample silence) {
            return aubio_notes_set_silence(notes,silence);
        }
        Sample get_silence() {
            return aubio_notes_get_silence(notes);
        }
        Sample get_minioi_ms() {
            return aubio_notes_get_minioi_ms(notes);
        }
        uint32_t set_minitoi_ms(Sample minioi_ms) {
            return aubio_notes_set_minioi_ms(notes,minioi_ms);
        }

        uint32_t release_drop(Sample release_drop) {
            return aubio_notes_set_release_drop(notes,release_drop);
        }
    };

    struct Onset {
        aubio_onset_t * onset;

        Onset(uint32_t buf_size, uint32_t hop_size, uint32_t samplerate) {
            onset = new_aubio_onset("default",buf_size,hop_size,samplerate);
        }
        ~Onset() {
            if(onset) {
                del_aubio_onset(onset);
            }
        }

        void process(const FVec & input, FVec & output) {
            aubio_onset_do(onset,input.pvec,output.pvec);
        }

        uint32_t get_last() {
            return aubio_onset_get_last(onset);
        }
        Sample get_last_s() {
            return aubio_onset_get_last_s(onset);
        }
        Sample get_last_ms() {
            return aubio_onset_get_last_ms(onset);            
        }
        uint32_t set_awhitening(uint32_t enable) {
            return aubio_onset_set_awhitening(onset,enable);
        }
        Sample get_awhitening() {
            return aubio_onset_get_awhitening(onset);
        }
        uint32_t set_compression(Sample lambda) {
            return aubio_onset_set_compression(onset,lambda);
        }
        Sample get_compression() {
            return aubio_onset_get_compression(onset);
        }
        uint32_t set_silence(Sample silence) {
            return aubio_onset_set_silence(onset,silence);
        }
        Sample get_silence() {
            return aubio_onset_get_silence(onset);
        }
        Sample get_thresholded_descriptor() {
            return aubio_onset_get_thresholded_descriptor(onset);
        }
        uint32_t set_threshold(Sample thresh) {
            return aubio_onset_set_threshold(onset,thresh);
        }
        uint32_t set_minioi(uint32_t minioi) {
            return aubio_onset_set_minioi(onset,minioi);
        }
        uint32_t set_minioi_s(Sample mini) {
            return aubio_onset_set_minioi_s(onset,mini);
        }
        uint32_t set_minioi_ms(Sample mini) {
            return aubio_onset_set_minioi_ms(onset,mini);
        }
        uint32_t set_delay(uint32_t delay) {
            return aubio_onset_set_delay(onset,delay);
        }
        uint32_t set_delay_s(Sample delay) {
            return aubio_onset_set_delay_s(onset,delay);
        }
        uint32_t set_delay_ms(Sample ms) {
            return aubio_onset_set_delay_ms(onset,ms);
        }
        uint32_t get_minioi() {
            return aubio_onset_get_minioi(onset);
        }
        Sample get_minioi_ms() {
            return aubio_onset_get_minioi_ms(onset);
        }
        uint32_t get_delay() {
            return aubio_onset_get_delay(onset);
        }
        Sample get_delay_s() {
            return aubio_onset_get_delay_s(onset);
        }
        Sample get_delay_ms() {
            return aubio_onset_get_delay_ms(onset);
        }
        Sample get_threshold() {
            return aubio_onset_get_threshold(onset);
        }        
        uint32_t set_default_parameters(const char* onset_mode) {
            return aubio_onset_set_default_parameters(onset,onset_mode);
        }
        void reset() {
            aubio_onset_reset(onset);
        }
    };

    struct PeakPicker {
        aubio_peakpicker_t * pick;

        PeakPicker() {
            pick = new_aubio_peakpicker();
            assert(pick != NULL);
        }
        ~PeakPicker() {
            if(pick) del_aubio_peakpicker(pick);
        }

        void process(FVec & in, FVec & out) {
            aubio_peakpicker_do(pick,in.pvec,out.pvec);
        }
        FVec get_thresholded_input() {
            return FVec(aubio_peakpicker_get_thresholded_input(pick));
        }
        uint32_t set_threshold(Sample thresh) {
            return aubio_peakpicker_set_threshold(pick,thresh);
        }
        Sample get_threshold() {
            return aubio_peakpicker_get_threshold(pick);
        }
    };

    struct Histogram
    {
        aubio_hist_t * hist;

        Histogram(Sample flow, Sample fhig, uint32_t nelems) {
            hist = new_aubio_hist(flow,fhig,nelems);
            assert(hist != NULL);
        }
        ~Histogram() {
            if(hist != NULL) del_aubio_hist(hist);
        }

        void process(FVec & input) {
            aubio_hist_do(hist,input.pvec);
        }
        void process_notnull(FVec & input) {
            aubio_hist_do_notnull(hist,input.pvec);
        }
        Sample hist_mean() {
            return aubio_hist_mean(hist);
        }
        void hist_weight() {
            aubio_hist_weight(hist);
        }
        void dyn_notnull(FVec & input) {
            aubio_hist_dyn_notnull(hist,input.pvec);
        }
    };

    inline Sample mean(FVec & s) { return fvec_mean(s.pvec); }
    inline Sample max(FVec & s) { return fvec_max(s.pvec); }
    inline Sample min(FVec & s) { return fvec_min(s.pvec); }
    inline uint32_t min_elem(FVec & s) { return fvec_min_elem(s.pvec); }
    inline uint32_t max_elem(FVec & s) { return fvec_max_elem(s.pvec); }
    inline void shift(FVec & s) { fvec_shift(s.pvec); }
    inline void ishift(FVec & s) { fvec_ishift(s.pvec); }
    inline void push(FVec & s, Sample new_elem) { fvec_push(s.pvec, new_elem); }
    inline Sample sum(FVec & s) { return fvec_sum(s.pvec); }
    inline Sample local_hfc(FVec & s) { return fvec_local_hfc(s.pvec); }
    inline Sample alpha_norm(FVec & s, Sample p) { return fvec_alpha_norm(s.pvec, p); }
    inline void alpha_normalise(FVec & s, Sample p) { fvec_alpha_normalise(s.pvec, p); }
    inline void add(FVec & v, Sample c) { fvec_add(v.pvec,c); }
    inline void mul(FVec & v, Sample c) { fvec_mul(v.pvec,c); }
    inline void remove_min(FVec & v) { fvec_min_removal(v.pvec); }
    inline Sample moving_threshold(FVec & v, FVec & tmp, uint32_t post, uint32_t pre, uint32_t pos) { return fvec_moving_thres(v.pvec, tmp.pvec, post, pre, pos); }
    inline void adapt_threshold(FVec & v, FVec & tmp, uint32_t post, uint32_t pre) { fvec_adapt_thres(v.pvec,tmp.pvec,post,pre); }
    inline Sample median(FVec & v) { return fvec_median(v.pvec); }
    inline Sample quadratic_peak_pos(FVec & x, uint32_t p) { return fvec_quadratic_peak_pos(x.pvec,p); }
    inline Sample quadratic_peak_mag(FVec &x, Sample p) { return fvec_quadratic_peak_mag(x.pvec, p); }
    inline Sample quadfrac(Sample s0, Sample s1, Sample s2, Sample s3) { return aubio_quadfrac(s0,s1,s2,s3); }
    inline Sample peakpick(FVec & v, uint32_t p) { return fvec_peakpick(v.pvec,p); }
    inline bool is_power_of_two(uint32_t x) { return aubio_is_power_of_two(x); }
    inline uint32_t next_power_of_two(uint32_t x) { return aubio_next_power_of_two(x); }
    inline uint32_t power_of_two_order(uint32_t x) { return aubio_power_of_two_order(x); }
    inline void autocorr(FVec & input, FVec & output) { aubio_autocorr(input.pvec,output.pvec); }


}