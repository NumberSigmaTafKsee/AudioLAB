#pragma once

namespace AudioDSP
{
    ////////////////////////////////////////////////////////////////
    // FFTW Complex 2 Complex
    ////////////////////////////////////////////////////////////////
    struct C2CD
    {
        fftw_complex * in;    
        fftw_complex * out;
        size_t size;
        fftw_plan p;

        C2CD(size_t n, Direction dir = FORWARD) {
            in = fftw_alloc_complex(n);
            out= fftw_alloc_complex(n);        
            size = n;
            p = fftw_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CD() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void set_input(complex_vector<double> & input) {
            for(size_t i = 0; i < size; i++) {
                in[i][0] = input[i].real();
                in[i][1] = input[i].imag();
            }
        }
        complex_vector<double> get_output() {
            complex_vector<double> r(size);
            for(size_t i = 0; i < size; i++ )
            {
                r[i].real(out[i][0]);
                r[i].imag(out[i][1]);
            }
            return r;
        }
        void get_output(complex_vector<double>&  output)
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
                out[i][0] /= (double)size;    
                out[i][1] /= (double)size;
            }
        }
        void Execute() {
            fftw_execute(p);
        }
    };
}