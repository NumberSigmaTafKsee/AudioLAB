#pragma once

namespace AudioDSP
{
    struct C2CF2D
    {
        fftwf_complex * in;    
        fftwf_complex * out;        
        fftwf_plan p;
        size_t M,N;

        C2CF2D(size_t m,size_t n, Direction dir = FORWARD) {
            in = fftwf_alloc_complex(m*n);
            out= fftwf_alloc_complex(m*n);        
            M = m;
            N = n;
            p = fftwf_plan_dft_2d(m,n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CF2D() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        size_t size() { return M*N; }

        void set_input(complex_matrix<float> & input) {
            for(size_t i = 0; i < input.rows(); i++)
                for(size_t j = 0; j < input.cols(); j++) {
                    in[i*N+j][0] = input(i,j).real();
                    in[i*N+j][1] = input(i,j).imag();
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
}