#pragma once

namespace AudioDSP
{
    struct R2CF
    {
        float * in;    
        fftwf_complex * out;
        size_t size;
        fftwf_plan p;

        R2CF() {
            in = NULL;
            out = NULL;
            size = 0;
        }
        R2CF(size_t n) {
            in = NULL;
            out = NULL;
            size = 0;
            init(n);            
        }
        ~R2CF() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = NULL;
            out = NULL;
            size = n;
            in = fftwf_alloc_real(n);
            out= fftwf_alloc_complex(n);                    
            p = fftwf_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(sample_vector<float> & input) {
            memcpy(in,input.data(),size*sizeof(float));
        }
        void set_input(float * buffer) {
            memcpy(in,buffer,size*sizeof(float));
        }
        complex_vector<float> get_output() {
            complex_vector<float> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i].real(out[i][0]);
                r[i].imag(out[i][1]);
            }
                
            return r;
        }    
        void get_output(float * output)
        {
            for(size_t i = 0; i < size; i++ )
            {
                output[2*i]   = out[i][0];
                output[2*i+1] = out[i][1];
            }
        }
        void get_output( complex_vector<float> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(out[i][0]);
                output[i].imag(out[i][1]);
            }
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                out[i][0] /= (float)size;    
                out[i][1] /= (float)size;
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
