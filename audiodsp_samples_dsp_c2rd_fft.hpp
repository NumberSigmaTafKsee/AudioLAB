#pragma once

namespace AudioDSP
{
    ////////////////////////////////////////////////////////////////
    // FFTW Complex 2 Real
    ////////////////////////////////////////////////////////////////
    struct C2RD
    {
        fftw_complex * in;    
        double * out;
        size_t size;
        fftw_plan p;

        C2RD() {
            in = NULL;
            out = NULL;
            size = 0;
        }
        C2RD(size_t n) {
            init(n);
        }
        ~C2RD() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void init(size_t n) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) fftw_free(in);
            if(out!= NULL) fftw_free(out);
            in = NULL;
            out = NULL;
            size = n;
            in = fftw_alloc_complex(n);
            out= fftw_alloc_real(n);                    
            p = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(complex_vector<double> & input) {
            for(size_t i = 0; i < size; i++) {
                in[i][0] = input[i].real();
                in[i][1] = input[i].imag();
            }
        }        
        sample_vector<double> get_output() {
            sample_vector<double> r(size);
            memcpy(r.data(),out, size * sizeof(double));
            return r;
        }        
        void get_output( sample_vector<double> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = out[i];                
            }
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                out[i] /= (double)size;                    
            }
        }
        void Execute() {
            fftw_execute(p);
        }
    };
}