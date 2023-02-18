#pragma once

namespace AudioDSP
{
    struct C2RF
    {
        fftwf_complex * in;    
        float * out;
        size_t size;
        fftwf_plan p;

        C2RF() {
            in = NULL;
            out = NULL;
            size = 0;
        }
        C2RF(size_t n) {
            in = NULL;
            out = NULL;
            size = 0;
            init(n);
        }
        ~C2RF() {
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
            in = fftwf_alloc_complex(n);
            out= fftwf_alloc_real(n);                    
            p = fftwf_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(complex_vector<float> & input) {
            for(size_t i = 0; i < size; i++) {
                in[i][0] = input[i].real();
                in[i][1] = input[i].imag();
            }
        }        
        sample_vector<float> get_output() {
            sample_vector<float> r(size);
            memcpy(r.data(),out, size*sizeof(float));
            return r;
        }        
        void get_output( sample_vector<float> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = out[i];
            }
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) 
                out[i] /= (float)size;                
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
}