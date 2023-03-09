#pragma once
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <fftw3.h>
#include <cstring>

#include "audiodsp_fft.hpp"

namespace AudioDSP
{
   
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
        void set_complex_input(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void get_complex_output(std::complex<double> * r) {
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
        void forward() {
            fftw_execute(pf);
        }
        void backward() {
            fftw_execute(pb);
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
        void set_complex_input(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void get_complex_output(std::complex<double> * r) {
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
        void forward() {
            fftw_execute(pf);
        }
        void backward() {
            fftw_execute(pb);
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

        void set_complex_input(const std::complex<float> * input)
        {
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void get_complex_output(std::complex<float> * r) {
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
        void forward() {
            fftwf_execute(pf);
        }
        void backward() {
            fftwf_execute(pb);
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
        void set_complex_input(const std::complex<float> * input)
        {            
            for(size_t i = 0; i < size; i++) {
                x[i][0] = input[i].real();
                x[i][1] = input[i].imag();
            }
        }
        void get_complex_output(std::complex<float> * r) {        
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
        void forward() {
            fftwf_execute(pf);
        }
        void backward() {
            fftwf_execute(pb);
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
        void set_input(const double * input)
        {
            memcpy(x,input,size*sizeof(double));
        }
        void set_complex_input(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void get_output(double * r) {        
            memcpy(r,x,size*sizeof(double));        
        }
        void get_complex_output(std::complex<double> * r) {        
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
        void forward() {
            fftw_execute(pf);
        }
        void backward() {
            fftw_execute(pb);
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
        void set_input(const double * input)
        {
            memcpy(x,input,size*sizeof(double));
        }
        void set_complex_input(const std::complex<double> * input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void get_output(double * r) {        
            memcpy(r,x,size*sizeof(double));        
        }
        void get_complex_output(std::complex<double> * r) {        
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
        void forward() {
            fftw_execute(pf);
        }
        void backward() {
            fftw_execute(pb);
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

        void set_input(const float * input)
        {
            memcpy(x,input,size*sizeof(float));
        }
        void set_complex_input(const std::complex<float>* input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void get_output(float * r) {        
            memcpy(r,x,size*sizeof(float));        
        }
        void get_complex_output(std::complex<float> * r) {        
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
        void forward() {
            fftwf_execute(pf);
        }
        void backward() {
            fftwf_execute(pb);
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
        void set_input(const float * input)
        {
            memcpy(x,input,size*sizeof(float));
        }
        void set_complex_input(const std::complex<float>* input)
        {
            for(size_t i = 0; i < size; i++) {
                y[i][0] = input[i].real();
                y[i][1] = input[i].imag();
            }
        }
        void get_output(float * r) {        
            memcpy(r,x,size*sizeof(float));        
        }
        void get_complex_output(std::complex<float> * r) {        
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
        void forward() {
            fftwf_execute(pf);
        }
        void backward() {
            fftwf_execute(pb);
        }
    };



    void fft(FFTPlanComplexDouble & plan, const std::complex<double> * in, std::complex<double> * out, bool norm=true)
    {
        plan.set_complex_input(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlanComplexDouble & plan, std::complex<double> * in, std::complex<double> * out, bool norm=false)
    {
        plan.set_complex_input(in);
        fftw_execute(plan.pb);    
        plan.get_complex_output(out);
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
        plan.set_complex_input(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft2(FFTPlanComplexDouble2D & plan, std::complex<double> * in, std::complex<double> * out, bool norm=false)
    {
        plan.set_complex_input(in);
        fftw_execute(plan.pb);    
        plan.get_complex_output(out);
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
        plan.set_complex_input(in);
        fftwf_execute(plan.pf);
        if(norm) plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlanComplexFloat & plan, const std::complex<float> * in, std::complex<float> * out, bool norm=false)
    {
        plan.set_complex_input(in);
        fftwf_execute(plan.pb);    
        plan.get_complex_output(out);
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
        plan.set_complex_input(in);
        fftwf_execute(plan.pf);
        if(norm) plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft2(FFTPlanComplexFloat2D & plan, const std::complex<float> * in, std::complex<float> * out, bool norm=false)
    {
        plan.set_complex_input(in);
        fftwf_execute(plan.pb);    
        plan.get_complex_output(out);
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
        plan.set_input(in);
        fftw_execute(plan.pf);
        if(norm) plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlanRealDouble & plan, std::complex<double> * in, double * out, bool norm=false)
    {
        plan.set_complex_input(in);
        fftw_execute(plan.pb);    
        plan.get_output(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }


    void fft2(FFTPlanRealDouble2D & plan, const double * in, std::complex<double> * out, bool norm=true)
    {
        plan.set_input(in);
        fftw_execute(plan.pf);
        plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft2(FFTPlanRealDouble2D & plan, std::complex<double> * in, double * out,bool norm=false)
    {
        plan.set_complex_input(in);
        fftw_execute(plan.pb);    
        plan.get_output(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft(FFTPlanRealFloat & plan, const float * in, std::complex<float> * out,bool norm=true)
    {
        plan.set_input(in);
        fftwf_execute(plan.pf);
        plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlanRealFloat & plan, const std::complex<float> * in, float * out,bool norm=false)
    {
        plan.set_complex_input(in);
        fftwf_execute(plan.pb);    
        plan.get_output(out);
        if(norm) {
            double max = *std::max_element(out,out+plan.size);
            if(max != 0.0) for(size_t i = 0; i < plan.size; i++) out[i] /= max;
        }
    }

    void fft2(FFTPlanRealFloat2D & plan, const float * in, std::complex<float> * out,bool norm=true)
    {
        plan.set_input(in);
        fftwf_execute(plan.pf);
        plan.normalize();
        plan.get_complex_output(out);
    }

    void ifft2(FFTPlanRealFloat2D & plan, const std::complex<float> * in, float * out,bool norm=false)
    {
        plan.set_complex_input(in);
        fftwf_execute(plan.pb);    
        plan.get_output(out);
        if(norm) {
            for(size_t i = 0; i < plan.size; i++)
                out[i] /= (double)plan.size;
        }
    }
}