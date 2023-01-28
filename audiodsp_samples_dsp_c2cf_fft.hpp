#pragma once

namespace AudioDSP
{
    struct C2CF
    {
        fftwf_complex * in;    
        fftwf_complex * out;
        size_t size;
        fftwf_plan p;

        C2CF(size_t n, Direction dir = FORWARD) {
            in = fftwf_alloc_complex(n);
            out= fftwf_alloc_complex(n);        
            size = n;
            p = fftwf_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CF() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void set_input(complex_vector<float> & input) {
            for(size_t i = 0; i < size; i++) {
                in[i][0] = input[i].real();
                in[i][1] = input[i].imag();
            }
        }
        void set_input(float * buffer) {
            memcpy(in,buffer,2*size*sizeof(float));
        }
        complex_vector<float> get_output() {
            complex_vector<float> r(size);
            for(size_t i = 0; i < size; i++ )
            {
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
        void get_output(complex_vector<float>& output)
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
}