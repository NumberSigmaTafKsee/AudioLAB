#pragma once

namespace AudioDSP
{
    struct C2RF2D
    {
        fftwf_complex * in;    
        float * out;
        size_t M,N;
        fftwf_plan p;

        C2RF2D() {
            in = NULL;
            out = NULL;            
            M = N = 0;
        }
        C2RF2D(size_t m,size_t n) {
            in = NULL;
            out = NULL;            
            init(m,n);
        }
        ~C2RF2D() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t m,size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = NULL;
            out = NULL;
            M = m;
            N = n;
            in = fftwf_alloc_complex(m*n);
            out= fftwf_alloc_real(m*n);                    
            p = fftwf_plan_dft_c2r_2d(m,n, in, out, FFTW_ESTIMATE);
        }
        size_t size() { return M*N; }

        void set_input(complex_matrix<float> & input) {
            for(size_t i = 0; i < input.rows(); i++)
                for(size_t j = 0; j < input.cols(); j++) {
                    in[i*N+j][0] = input(i,j).real();
                    in[i*N+j][1] = input(i,j).imag();
                }
        }
        sample_matrix<float> get_output() {
            sample_matrix<float> r(M,N);
            for(size_t i = 0; i < r.rows(); i++)
                for(size_t j = 0; j < r.cols(); j++) 
                {
                    r(i,j) = (out[i*N+j]);                    
                }
            return r;
        }
        
        void get_output(sample_matrix<float>&  output)
        {
            if(output.size() != size()) output.resize(M,N);
            for(size_t i = 0; i < output.rows(); i++)
                for(size_t j = 0; j < output.cols(); j++) 
                {
                    output(i,j) = (out[i*N+j]);                    
                }
        }
        void normalize() {
            for(size_t i = 0; i < M*N; i++)                
            {
                out[i] /= (float)size();    
                out[i] /= (float)size();
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
}