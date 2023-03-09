#pragma once

namespace AudioDSP
{
    struct R2CF2D
    {
        float * in;    
        fftwf_complex * out;
        size_t M,N;
        fftwf_plan p;

        R2CF2D() {
            in = NULL;
            out = NULL;
            M = N = 0;
        }
        R2CF2D(size_t m, size_t n) {
            in = NULL;
            out = NULL;
            M = N = 0;
            init(m,n);            
        }
        ~R2CF2D() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t m, size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = NULL;
            out = NULL;
            M = m;
            N = n;
            in = fftwf_alloc_real(n);
            out= fftwf_alloc_complex(n);                    
            p = fftwf_plan_dft_r2c_2d(m,n, in, out, FFTW_ESTIMATE);
        }
        size_t size() { return M*N; }

        void set_input(sample_matrix<float> & input) {
            for(size_t i = 0; i < input.rows(); i++)
                for(size_t j = 0; j < input.cols(); j++) {
                    in[i*N+j] = input(i,j);                    
                }
        }
        complex_matrix<float> get_output() {
            complex_matrix<float> r(M,N);
            for(size_t i = 0; i < r.rows(); i++)
                for(size_t j = 0; j < r.cols(); j++) 
                {
                    r(i,j).real(out[i*N+j][0]);
                    r(i,j).imag(out[i*N+j][1]);
                }
            return r;
        }
        
        void get_output(complex_matrix<float>&  output)
        {
            if(output.size() != size()) output.resize(M,N);
            for(size_t i = 0; i < output.rows(); i++)
                for(size_t j = 0; j < output.cols(); j++) 
                {
                    output(i,j).real(out[i*N+j][0]);
                    output(i,j).imag(out[i*N+j][1]);
                }
        }
        void normalize() {
            for(size_t i = 0; i < M*N; i++)                
            {
                out[i][0] /= (float)size();    
                out[i][1] /= (float)size();
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
