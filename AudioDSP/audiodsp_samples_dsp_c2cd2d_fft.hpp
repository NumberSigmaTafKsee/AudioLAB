#pragma once

namespace AudioDSP
{
    struct C2CD2D
    {
        fftw_complex * in;    
        fftw_complex * out;
        size_t M,N;
        fftw_plan p;

        C2CD2D(size_t m, size_t n, Direction dir = FORWARD) {
            in = fftw_alloc_complex(m*n);
            out= fftw_alloc_complex(m*n);        
            M = m;
            N = n;
            p = fftw_plan_dft_2d(m,n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CD2D() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        size_t size() { return M*N; }
        void set_input(complex_matrix<double> & input) {
            for(size_t i = 0; i < input.rows(); i++)
                for(size_t j = 0; j < input.cols(); j++) {
                    in[i*N+j][0] = input(i,j).real();
                    in[i*N+j][1] = input(i,j).imag();
                }
        }
        complex_vector<double> get_output() {
            complex_vector<double> r(M,N);
            for(size_t i = 0; i < r.rows(); i++)
                for(size_t j = 0; j < r.cols(); j++) 
                {
                    r(i,j).real(out[i*N+j][0]);
                    r(i,j).imag(out[i*N+j][1]);
                }
            return r;
        }
        
        void get_output(complex_matrix<double>&  output)
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
                out[i][0] /= (double)size();    
                out[i][1] /= (double)size();
            }
        }
        void Execute() {
            fftw_execute(p);
        }
    };
}