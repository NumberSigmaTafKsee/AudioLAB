#pragma once

namespace AudioDSP
{
    ////////////////////////////////////////////////////////////////
    // FFTW Real 2 Complex
    ////////////////////////////////////////////////////////////////
    struct R2CD
    {
        double       * in;    
        fftw_complex * out;
        size_t         size;
        fftw_plan      p;

        R2CD() {
            in = NULL;
            out = NULL;
            size= 0;
        }
        R2CD(size_t n) {
            in = NULL;
            out = NULL;
            size= 0;
            init(n);            
        }
        ~R2CD() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void init(size_t n) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) fftw_free(in);
            if(out != NULL) fftw_free(out);
            in = NULL;
            out = NULL;
            size = n;
            in = fftw_alloc_real(n);
            out= fftw_alloc_complex(n);                    
            p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(sample_vector<double> & input) {
            memcpy(in,input.data(),size*sizeof(double));
        }        
        complex_vector<double> get_output() {
            complex_vector<double> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i].real(out[i][0]);
                r[i].imag(out[i][1]);
            }
            return r;
        }        
        void get_output(complex_vector<double> & output) {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++)
            {
                output[i].real(out[i][0]);
                output[i].imag(out[i][1]);
            }            
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                out[i][0] /= (double)size;    
                out[i][1] /= (double)size;
            }
        }
        void Execute() {
            fftw_execute(p);            
        }
    };
