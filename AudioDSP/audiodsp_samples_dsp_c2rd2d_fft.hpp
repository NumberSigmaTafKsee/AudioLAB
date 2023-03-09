#pragma once

namespace AudioDSP
{
    struct C2RD2D
    {
        fftw_complex * in;    
        double * out;
        size_t M,N;
        fftw_plan p;

        C2RD2D() {
            in = NULL;
            out = NULL;
            M = N = 0;
        }
        C2RD2D(size_t m, size_t n) {
            init(m,n);
        }
        ~C2RD2D() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void init(size_t m,size_t n) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) fftw_free(in);
            if(out!= NULL) fftw_free(out);
            in = NULL;
            out = NULL;
            M = m;
            N = n;
            in = fftw_alloc_complex(m*n);
            out= fftw_alloc_real(m*n);                    
            p = fftw_plan_dft_c2r_2d(m,n, in, out, FFTW_ESTIMATE);
        }
        size_t size() { return M*N; }

        void set_input(complex_matrix<double> & input) {
            for(size_t i = 0; i < input.rows(); i++)
                for(size_t j = 0; j < input.cols(); j++) {
                    in[i*N+j][0] = input(i,j).real();
                    in[i*N+j][1] = input(i,j).imag();
                }
        }
        sample_matrix<double> get_output() {
            sample_matrix<double> r(M,N);
            for(size_t i = 0; i < r.rows(); i++)
                for(size_t j = 0; j < r.cols(); j++) 
                {
                    r(i,j) = (out[i*N+j]);                    
                }
            return r;
        }
        
        void get_output(sample_matrix<double>&  output)
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
                out[i] /= (double)size();    
                out[i] /= (double)size();
            }
        }
        void Execute() {
            fftw_execute(p);
        }
    };
}