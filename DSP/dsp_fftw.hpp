#pragma once
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <fftw3.h>
#include <cstring>

namespace Casino
{

    template<typename T>
    struct Window : public sample_vector<T>
    {   
        sample_vector<T> window;
        enum WindowType {
            RECTANGLE,
            HANNING,
            HAMMING,
            BLACKMAN,
            BLACKMANHARRIS,
            GAUSSIAN,
            WELCH,
            PARZEN,
            TUKEY,
        } type = HANNING;

        Window() = default;
        Window(WindowType type, size_t i) { window.resize(i); make_window(i,type); }
        virtual ~Window() = default;

        void make_window(size_t n, WindowType type)
        {
            this->type = type;
            this->resize(n);
            switch(type)
            {
                case RECTANGLE: rectangle(n); break;
                case HANNING: hanning(n); break;
                case HAMMING: hamming(n); break;
                case BLACKMAN: blackman(n); break;
                case BLACKMANHARRIS: blackmanharris(n); break;
                case GAUSSIAN: gaussian(n); break;
                case WELCH: welch(n); break;
                case PARZEN: parzen(n); break;
                case TUKEY: throw std::runtime_error("Can't init tukey window with make_window");
            }            
        }
        void rectangle(size_t n) {
            
            std::fill(this->begin(),this->end(),(T)1.0);
        }
        void hamming(size_t n) {
            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)[i] = 0.54 - (0.46 * std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
        void hanning(size_t n)
        {
            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)[i] = 0.5*(1 - std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
        void blackman(size_t n)
        {        
            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)                    
                (*this)[i] = 0.42 - (0.5* std::cos(2*M_PI*i/(n)) + (0.08*std::cos(4*M_PI*i/n)));        
        }
        void blackmanharris(size_t n)
        {
            
            #pragma omp simd        
            for(size_t i = 0; i < this->size(); i++)            
            {   
                double ci = (double) i / (double) n;
                (*this)[i] = 0.35875 
                        - 0.48829*std::cos(2*M_PI*(ci))
                        + 0.14128*std::cos(4.0*M_PI*(ci)) 
                        - 0.01168*std::cos(6.0*M_PI*(ci));
            }
        }
        void gaussian(size_t n)
        {
            
            T a,b,c=0.5;
            #pragma omp simd        
            for(size_t i = 0; i < this->size(); i++)
            {
                a = ((double)i - c*(this->size()-1)/(std::sqrt(c)*this->size()-1));
                b = -c * std::sqrt(a);
                (*this)[i] = std::exp(b);
            }    
        }
        void welch(size_t n)
        {
            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
                (*this)[i] = 1.0 - std::sqrt((2.0*(double)i-(double)this->size()-1)/((double)this->size()));        
        }
        void parzen(size_t n)
        {
            
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
                (*this)[i] = 1.0 - std::abs((2.0*(double)i-this->size()-1)/(this->size()));        
        }
        void tukey(size_t num_samples, T alpha)
        {
            
            T value = (-1*(num_samples/2)) + 1;
            double n2 = (double)num_samples / 2.0;
            #pragma omp simd
            for(size_t i = 0; i < this->size(); i++)
            {    
                if(value >= 0 && value <= (alpha * (n2))) 
                    (*this)[i] = 1.0; 
                else if(value <= 0 && (value >= (-1*alpha*(n2)))) 
                    (*this)[i] = 1.0;
                else 
                    (*this)[i] = 0.5 * (1 + std::cos(M_PI *(((2.0*value)/(alpha*(double)num_samples))-1)))        ;
                value = value + 1;
            }     
        }
    };

    struct FFTPlanComplexDouble
    {
        fftw_complex *  x=nullptr;    
        fftw_complex *  y=nullptr;
        size_t          size;
        fftw_plan       pf,pb;

        FFTPlanComplexDouble() = default;
        FFTPlanComplexDouble(size_t n) 
        {
            init(n);
        }
        ~FFTPlanComplexDouble()
        {
            deinit();
        }
        void deinit() {
            if(x) fftw_free(x);
            if(y) fftw_free(y);
            if(pf) fftw_destroy_plan(pf);
            if(pb) fftw_destroy_plan(pb);
        }
        void init(size_t n)
        {
            if(x != nullptr) deinit();
            x = fftw_alloc_complex(n);
            y = fftw_alloc_complex(n);        
            size = n;
            pf = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            pb = fftw_plan_dft_1d(n, x, y, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        void setInput(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void getOutput(std::complex<double> * r) {
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<double>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (double)size;    
                y[i][1] /= (double)size;
            }
        }
    };

    struct FFTPlanComplexDouble2D
    {
        fftw_complex *  x=nullptr;    
        fftw_complex *  y;
        size_t          size,M,N;
        fftw_plan       pf,pb;

        FFTPlanComplexDouble2D() = default;
        FFTPlanComplexDouble2D(size_t m,size_t n) 
        {
            init(m,n);
        }
        ~FFTPlanComplexDouble2D()
        {
            deinit();
        }
        void deinit() {
            if(x) fftw_free(x);
            if(y) fftw_free(y);
            if(pf) fftw_destroy_plan(pf);
            if(pb) fftw_destroy_plan(pb);
        }

        void init(size_t m,size_t n) 
        {
            if(x != nullptr) deinit();
            size = m*n;
            M = m;
            N = n;
            x = fftw_alloc_complex(size);
            y = fftw_alloc_complex(size);                    
            pf = fftw_plan_dft_2d(m,n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            pb = fftw_plan_dft_2d(m,n, x, y, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        void setInput(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void getOutput(std::complex<double> * r) {
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<double>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (double)size;    
                y[i][1] /= (double)size;
            }
        }
    };

    struct FFTPlanComplexFloat
    {
        fftwf_complex * x=nullptr;    
        fftwf_complex * y;
        size_t size;
        fftwf_plan pf,pb;

        FFTPlanComplexFloat() = default;
        FFTPlanComplexFloat(size_t n)     
        {
            init(n);
        }
        ~FFTPlanComplexFloat()
        {
            deinit();
        }
        void deinit() {
            if(x) fftwf_free(x);
            if(y) fftwf_free(y);
            if(pf) fftwf_destroy_plan(pf);
            if(pb) fftwf_destroy_plan(pb);
        }        
        void init(size_t n)
        {
            if(x != nullptr) deinit();
            x = fftwf_alloc_complex(n);
            y = fftwf_alloc_complex(n);        
            size = n;
            pf = fftwf_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            pb = fftwf_plan_dft_1d(n, x, y, FFTW_BACKWARD, FFTW_ESTIMATE);
        }

        void setInput(const std::complex<float> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void getOutput(std::complex<float> * r) {
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<float>(y[i][0],y[i][1]);
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (float)size;    
                y[i][1] /= (float)size;
            }
        }
    };

    struct FFTPlanComplexFloat2D
    {
        fftwf_complex * x=nullptr;    
        fftwf_complex * y;
        size_t size,M,N;
        fftwf_plan pf,pb;

        FFTPlanComplexFloat2D() = default;
        FFTPlanComplexFloat2D(size_t m, size_t n)     
        {
            init(m,n);
        }
        ~FFTPlanComplexFloat2D()
        {
            deinit();
        }
               
        void deinit()
        {
            if(x) fftwf_free(x);
            if(y) fftwf_free(y);
            if(pf) fftwf_destroy_plan(pf);
            if(pb) fftwf_destroy_plan(pb);
        }


        void init(size_t m, size_t n)     
        {
            if(x != nullptr) deinit();
            size = n;
            M = m;
            N = n;
            x = fftwf_alloc_complex(size);
            y = fftwf_alloc_complex(size);        
            
            pf = fftwf_plan_dft_2d(m,n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            pb = fftwf_plan_dft_2d(m,n, x, y, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        void setInput(const std::complex<float> * input)
        {            
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void getOutput(std::complex<float> * r) {        
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<float>(y[i][0],y[i][1]);
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (float)size;    
                y[i][1] /= (float)size;
            }
        }
    };

    struct FFTPlanRealDouble
    {
        double * x = nullptr;
        fftw_complex * y;
        size_t size;
        fftw_plan pf,pb;

        FFTPlanRealDouble() = default;
        FFTPlanRealDouble(size_t n)
        {
            init(n);
        }
        ~FFTPlanRealDouble() {
            deinit();
        }


        void deinit()
        {
            if(x) fftw_free(x);
            if(y) fftw_free(y);
            if(pf) fftw_destroy_plan(pf);
            if(pb) fftw_destroy_plan(pb);
        }
        void init(size_t n)     
        {
            if(x != nullptr) deinit();
            x = fftw_alloc_real(n);
            y = fftw_alloc_complex(n);        
            size = n;
            pf = fftw_plan_dft_r2c_1d(n, x, y, FFTW_ESTIMATE);
            pb = fftw_plan_dft_c2r_1d(n, y, x, FFTW_ESTIMATE);
        }
        void setReal(const double * input)
        {
            memcpy(x,input,size*sizeof(double));
        }
        void setComplex(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void getReal(double * r) {        
            memcpy(r,x,size*sizeof(double));        
        }
        void getComplex(std::complex<double> * r) {        
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<double>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (double)size;    
                y[i][1] /= (double)size;
            }
        }
    };

    struct FFTPlanRealDouble2D
    {
        double * x = nullptr;
        fftw_complex * y;
        size_t size,M,N;
        fftw_plan pf,pb;

        FFTPlanRealDouble2D() = default;
        FFTPlanRealDouble2D(size_t m, size_t n)     
        {
            init(m,n);
        }
        ~FFTPlanRealDouble2D() {
            deinit();
        }
        
        void deinit()
        {
            if(x) fftw_free(x);
            if(y) fftw_free(y);
            if(pf) fftw_destroy_plan(pf);
            if(pb) fftw_destroy_plan(pb);
        }
        void init(size_t m, size_t n)     
        {
            if(x != nullptr) deinit();
            size = m*n;
            M = m;
            N = n;
            x = fftw_alloc_real(size);
            y = fftw_alloc_complex(size);        
            
            pf = fftw_plan_dft_r2c_2d(m,n, x, y, FFTW_ESTIMATE);
            pb = fftw_plan_dft_c2r_2d(m,n, y, x, FFTW_ESTIMATE);
        }
        void setReal(const double * input)
        {
            memcpy(x,input,size*sizeof(double));
        }
        void setComplex(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void getReal(double * r) {        
            memcpy(r,x,size*sizeof(double));        
        }
        void getComplex(std::complex<double> * r) {        
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<double>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (double)size;    
                y[i][1] /= (double)size;
            }
        }
    };

    struct FFTPlanRealFloat
    {
        float * x = nullptr;
        fftwf_complex * y;
        size_t size;
        fftwf_plan pf,pb;

        FFTPlanRealFloat() = default;
        FFTPlanRealFloat(size_t n)     
        {
            init(n);
        }
        ~FFTPlanRealFloat()
        {
            deinit();
        }
        void deinit()
        {
            if(x) fftwf_free(x);
            if(y) fftwf_free(y);
            if(pf) fftwf_destroy_plan(pf);
            if(pb) fftwf_destroy_plan(pb);
        }
        void init(size_t n)
        {
            if(x != nullptr) deinit();
            x = fftwf_alloc_real(n);
            y = fftwf_alloc_complex(n);        
            size = n;
            pf = fftwf_plan_dft_r2c_1d(n, x, y, FFTW_ESTIMATE);
            pb = fftwf_plan_dft_c2r_1d(n, y, x, FFTW_ESTIMATE);
        }

        void setReal(const float * input)
        {
            memcpy(x,input,size*sizeof(float));
        }
        void setComplex(const std::complex<float>* input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void getReal(float * r) {        
            memcpy(r,x,size*sizeof(float));        
        }
        void getComplex(std::complex<float> * r) {        
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<float>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (float)size;    
                y[i][1] /= (float)size;
            }
        }
    };

    struct FFTPlanRealFloat2D
    {
        float * x = nullptr;
        fftwf_complex * y;
        size_t size,M,N;
        fftwf_plan pf,pb;

        FFTPlanRealFloat2D() = default;

        FFTPlanRealFloat2D(size_t m, size_t n)     
        {
            init(m,n);
        }
        
        ~FFTPlanRealFloat2D() {
            deinit();
        }
        
        void deinit()
        {
            if(x) fftwf_free(x);
            if(y) fftwf_free(y);
            if(pf) fftwf_destroy_plan(pf);
            if(pb) fftwf_destroy_plan(pb);
        }
        void init(size_t m, size_t n)     
        {
            size = m*n;
            M = m;
            N = n;
            x = fftwf_alloc_real(size);
            y = fftwf_alloc_complex(size);                    
            pf = fftwf_plan_dft_r2c_2d(m,n, x, y, FFTW_ESTIMATE);
            pb = fftwf_plan_dft_c2r_2d(m,n, y, x, FFTW_ESTIMATE);
        }
        void setReal(const float * input)
        {
            memcpy(x,input,size*sizeof(float));
        }
        void setComplex(const std::complex<float>* input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void getReal(float * r) {        
            memcpy(r,x,size*sizeof(float));        
        }
        void getComplex(std::complex<float> * r) {        
            for(size_t i = 0; i < size; i++)        
                r[i] = std::complex<float>(y[i][0],y[i][1]);        
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] /= (float)size;    
                y[i][1] /= (float)size;
            }
        }
    };



    void fft(FFTPlanComplexDouble & plan, const std::complex<double> * in, std::complex<double> * out, bool norm=true)
    {
        plan.setInput(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.getOutput(out);
    }

    void ifft(FFTPlanComplexDouble & plan, std::complex<double> * in, std::complex<double> * out, bool norm=false)
    {
        plan.setInput(in);
        fftw_execute(plan.pb);    
        plan.getOutput(out);
        if(norm) {
            double max = std::abs(out[0]);
            for(size_t i = 1; i < plan.size; i++)
            {
                double temp = std::abs(out[i]);
                if(temp > max) max = temp;
            }
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }


    void fft2(FFTPlanComplexDouble2D & plan, const std::complex<double> * in, std::complex<double> * out, bool norm=true)
    {
        plan.setInput(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.getOutput(out);
    }

    void ifft2(FFTPlanComplexDouble2D & plan, std::complex<double> * in, std::complex<double> * out, bool norm=false)
    {
        plan.setInput(in);
        fftw_execute(plan.pb);    
        plan.getOutput(out);
        if(norm) {
            double max = std::abs(out[0]);
            for(size_t i = 1; i < plan.size; i++)
            {
                double temp = std::abs(out[i]);
                if(temp > max) max = temp;
            }
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft(FFTPlanComplexFloat & plan, const std::complex<float> * in, std::complex<float> * out, bool norm = true)
    {
        plan.setInput(in);
        fftwf_execute(plan.pf);
        if(norm) plan.normalize();
        plan.getOutput(out);
    }

    void ifft(FFTPlanComplexFloat & plan, const std::complex<float> * in, std::complex<float> * out, bool norm=false)
    {
        plan.setInput(in);
        fftwf_execute(plan.pb);    
        plan.getOutput(out);
        if(norm) {
            double max = std::abs(out[0]);
            for(size_t i = 1; i < plan.size; i++)
            {
                double temp = std::abs(out[i]);
                if(temp > max) max = temp;
            }
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft2(FFTPlanComplexFloat2D & plan, const std::complex<float> * in, std::complex<float> * out, bool norm=true)
    {
        plan.setInput(in);
        fftwf_execute(plan.pf);
        if(norm) plan.normalize();
        plan.getOutput(out);
    }

    void ifft2(FFTPlanComplexFloat2D & plan, const std::complex<float> * in, std::complex<float> * out, bool norm=false)
    {
        plan.setInput(in);
        fftwf_execute(plan.pb);    
        plan.getOutput(out);
        if(norm) {
            double max = std::abs(out[0]);
            for(size_t i = 1; i < plan.size; i++)
            {
                double temp = std::abs(out[i]);
                if(temp > max) max = temp;
            }
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft(FFTPlanRealDouble & plan, const double * in, std::complex<double> * out, bool norm=true)
    {
        plan.setReal(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.getComplex(out);
    }

    void ifft(FFTPlanRealDouble & plan, std::complex<double> * in, double * out, bool norm=false)
    {
        plan.setComplex(in);
        fftw_execute(plan.pb);    
        plan.getReal(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }


    void fft2(FFTPlanRealDouble2D & plan, const double * in, std::complex<double> * out, bool norm=true)
    {
        plan.setReal(in);
        fftw_execute(plan.pf);
        plan.normalize();
        plan.getComplex(out);
    }

    void ifft2(FFTPlanRealDouble2D & plan, std::complex<double> * in, double * out,bool norm=false)
    {
        plan.setComplex(in);
        fftw_execute(plan.pb);    
        plan.getReal(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft(FFTPlanRealFloat & plan, const float * in, std::complex<float> * out,bool norm=true)
    {
        plan.setReal(in);
        fftwf_execute(plan.pf);
        plan.normalize();
        plan.getComplex(out);
    }

    void ifft(FFTPlanRealFloat & plan, const std::complex<float> * in, float * out,bool norm=false)
    {
        plan.setComplex(in);
        fftwf_execute(plan.pb);    
        plan.getReal(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft2(FFTPlanRealFloat2D & plan, const float * in, std::complex<float> * out,bool norm=true)
    {
        plan.setReal(in);
        fftwf_execute(plan.pf);
        plan.normalize();
        plan.getComplex(out);
    }

    void ifft2(FFTPlanRealFloat2D & plan, const std::complex<float> * in, float * out,bool norm=false)
    {
        plan.setComplex(in);
        fftwf_execute(plan.pb);    
        plan.getReal(out);
        if(norm) {
            for(size_t i = 0; i < plan.size; i++)
                out[i] /= (double)plan.size;
        }
    }


    // Code adapted from gsl/fft/factorize.c
    void factorize (const int n,
                    int *n_factors,
                    int factors[],
                    int * implemented_factors)
    {
        int nf = 0;
        int ntest = n;
        int factor;
        int i = 0;

        if (n == 0)
        {
            printf("Length n must be positive integer\n");
            return ;
        }

        if (n == 1)
        {
            factors[0] = 1;
            *n_factors = 1;
            return ;
        }

        /* deal with the implemented factors */

        while (implemented_factors[i] && ntest != 1)
        {
            factor = implemented_factors[i];
            while ((ntest % factor) == 0)
            {
                ntest = ntest / factor;
                factors[nf] = factor;
                nf++;
            }
            i++;
        }

        // Ok that's it
        if(ntest != 1)
        {
            factors[nf] = ntest;
            nf++;
        }

        /* check that the factorization is correct */
        {
            int product = 1;

            for (i = 0; i < nf; i++)
            {
                product *= factors[i];
            }

            if (product != n)
            {
                printf("factorization failed");
            }
        }

        *n_factors = nf;
    }



    bool is_optimal(int n, int * implemented_factors)
    {
        // We check that n is not a multiple of 4*4*4*2
        if(n % 4*4*4*2 == 0)
            return false;

        int nf=0;
        int factors[64];
        int i = 0;
        factorize(n, &nf, factors,implemented_factors);

        // We just have to check if the last factor belongs to GSL_FACTORS
        while(implemented_factors[i])
        {
            if(factors[nf-1] == implemented_factors[i])
                return true;
            ++i;
        }
        return false;
    }

    int find_closest_factor(int n)
    {
        int j;
        int FFTW_FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array
        if(is_optimal(n,FFTW_FACTORS))
            return n;
        else
        {
            j = n+1;
            while(!is_optimal(j,FFTW_FACTORS))
                ++j;
            return j;
        }
    }

    
    
    struct FFTConvolutionDouble
    {
        FFTPlanRealDouble fftPlan;
        size_t  length,fftSize,blockSize;
 
        std::vector<std::complex<double>> t1,t2,tempC;
        std::vector<double> H,i1,ola,temp;
        
        FFTConvolutionDouble(size_t len, double * h, size_t blocks)
        {           
            length  = len;            
            blockSize = blocks;
            fftSize = find_closest_factor(len + blocks -1);
            
            fftPlan.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            tempC.resize(fftSize);

            H.resize(fftSize);
            memset(H.data(),0,H.size()*sizeof(double));
            memcpy(H.data(),h,len*sizeof(double));     

            i1.resize(fftSize);
            memset(i1.data(),0,i1.size()*sizeof(double));

            ola.resize(fftSize);
            memset(ola.data(),0,ola.size()*sizeof(double));
            
            temp.resize(fftSize);            
        }

        void ProcessBlock(size_t n, double * in, double * out)
        {
            memcpy(i1.data(),in,n*sizeof(double));

            
            fft(fftPlan,H.data(),t1.data(),false);
            fft(fftPlan,i1.data(),t2.data(),false);
                        
            tempC[0] = std::complex<double>(0,0);
            tempC[t2.size()/2-1] = std::complex<double>(0,0);

            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                tempC[i] = t1[i] * t2[i];
            }
              
            ifft(fftPlan,tempC.data(),temp.data(),true);

            for(size_t i = 0; i < n; i++)
            {
                out[i] = (temp[i] + ola[i]);                
                ola[i] = ola[i+n];                
            }            
            for(size_t i = n; i < temp.size(); i++)
                ola[i] = temp[i];            
        }
    };

    struct FFTConvolutionFloat
    {
        FFTPlanRealFloat fftPlan;
        size_t  length,fftSize,blockSize;
 
        std::vector<std::complex<float>> t1,t2,tempC;
        std::vector<float> H,i1,ola,temp;
        
        FFTConvolutionFloat(size_t len, float * h, size_t blocks)
        {           
            length  = len;            
            blockSize = blocks;
            fftSize = find_closest_factor(len + blocks -1);
            
            fftPlan.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            tempC.resize(fftSize);

            H.resize(fftSize);
            memset(H.data(),0,H.size()*sizeof(float));
            memcpy(H.data(),h,len*sizeof(float));     

            i1.resize(fftSize);
            memset(i1.data(),0,i1.size()*sizeof(float));

            ola.resize(fftSize);
            memset(ola.data(),0,ola.size()*sizeof(float));
            
            temp.resize(fftSize);            
        }

        void ProcessBlock(size_t n, float * in, float * out)
        {
            memcpy(i1.data(),in,n*sizeof(float));

            
            fft(fftPlan,H.data(),t1.data(),false);
            fft(fftPlan,i1.data(),t2.data(),false);
                        
            tempC[0] = std::complex<double>(0,0);
            tempC[t2.size()/2-1] = std::complex<float>(0,0);

            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                tempC[i] = t1[i] * t2[i];
            }
              
            ifft(fftPlan,tempC.data(),temp.data(),true);

            for(size_t i = 0; i < n; i++)
            {
                out[i] = (temp[i] + ola[i]);                
                ola[i] = ola[i+n];                
            }            
            for(size_t i = n; i < temp.size(); i++)
                ola[i] = temp[i];            
        }
    };
}