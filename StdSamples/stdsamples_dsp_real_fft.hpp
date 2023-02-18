#pragma once

namespace AudioDSP
{
   ////////////////////////////////////////////////////////////////
    // FFTW R2R DCT/DST
    ////////////////////////////////////////////////////////////////
    enum R2RKind
    {
        DCTI = FFTW_REDFT00,
        DCTII= FFTW_REDFT01,
        DCTIII=FFTW_REDFT10,
        DCTIV=FFTW_REDFT11,
        DSTI=FFTW_REDFT11,
        DSTII=FFTW_REDFT00,
        DSTIII=FFTW_REDFT10,
        DSTIV=FFTW_REDFT11,
    };

    struct R2RD
    {
        double       * in;    
        double       * out;
        size_t size;
        fftw_plan p;

        R2RD() {
            in = out = NULL;
            size = 0;
        }
        R2RD(size_t n, R2RKind type = DCTI) {
            init(n,type);
        }
        ~R2RD() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void init(size_t n, R2RKind type) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) fftw_free(in);
            if(out != NULL) fftw_free(out);
            in = out = NULL;            
            size = n;            
            in = fftw_alloc_real(n);
            out= fftw_alloc_real(n);                    
            p = fftw_plan_r2r_1d(n, in, out, (fftw_r2r_kind)type, FFTW_ESTIMATE);
        }
        void set_input(sample_vector<double> & input) {
            memcpy(in,input.data(),size*sizeof(double));
        }
        sample_vector<double> get_output() {
            sample_vector<double> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i] = out[i];                
            }
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
            for(size_t i = 0; i < size; i++)
                out[i] /= (double)size;    
        }
        void Execute() {
            fftw_execute(p);            
        }
    };

    struct R2RD2D
    {
        double       * in;    
        double       * out;
        size_t M,N;
        fftw_plan p;

        R2RD2D() {
            in = out = NULL;
            M  = N   = 0;
        }
        R2RD2D(size_t m, size_t n, R2RKind type1 = DCTI, R2RKind type2 = DCTI) {
            init(m,n,type1,type2);
        }
        ~R2RD2D() {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);    
        }
        void init(size_t m, size_t n, R2RKind type1, R2RKind type2,unsigned flag = FFTW_ESTIMATE) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) fftw_free(in);
            if(out != NULL) fftw_free(out);
            in = out = NULL;            
            M = m;
            N = n;     
            in = fftw_alloc_real(m*n);
            out= fftw_alloc_real(m*n);                    
            p = fftw_plan_r2r_2d(m,n, in, out, (fftw_r2r_kind)type1, (fftw_r2r_kind)type2,flag);
        }
        
        size_t size() const { return M*N; }

        void set_input(sample_matrix<double> & input) {
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                    in[i*N+j] = input(i,j);
        }
        
        sample_matrix<double> get_output() {
            sample_matrix<double> r(M,N);
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) = out[i*N+j];
                }
            return r;
        }        
        void get_output( sample_matrix<double> & output)
        {
            if(output.size() != size()) output.resize(M,N);
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                {
                    output(i,j) = out[i*N+j];
                }
        }
        void normalize() {
            for(size_t i = 0; i < size(); i++)
                out[i] /= (double)size();    
        }
        void Execute() {
            fftw_execute(p);            
        }
    };

    struct R2RF
    {
        float * in;    
        float * out;
        size_t size;
        fftwf_plan p;

        R2RF() {
            in = out = NULL;
            size = 0;
        }
        R2RF(size_t n, R2RKind type = DCTI) {
            in = out = NULL;
            size = 0;
            init(n,type);            
        }
        ~R2RF() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t n, R2RKind type) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = out = NULL;            
            size = n;            
            in = fftwf_alloc_real(n);
            out= fftwf_alloc_real(n);                    
            p = fftwf_plan_r2r_1d(n, in, out, (fftw_r2r_kind)type, FFTW_ESTIMATE);
        }
        void set_input(sample_vector<float> & input) {
            memcpy(in,input.data(),size*sizeof(float));
        }
        sample_vector<float> get_output() {
            sample_vector<float> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i] = out[i];                
            }                
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
        void normalize() {
            for(size_t i = 0; i < size; i++)
                out[i] /= (float)size;    
        }    
        void Execute() {
            fftwf_execute(p);            
        }
    };


    struct R2RF2D
    {
        float * in;    
        float * out;
        size_t M,N;
        fftwf_plan p;

        R2RF2D() {
            in = out = NULL;
            M=N=0;
        }
        R2RF2D(size_t m,size_t n, R2RKind type1 = DCTI, R2RKind type2 = DCTI) {
            in = out = NULL;
            M=N=0;
            init(m,n,type1,type2);     
        }
        ~R2RF2D() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t m,size_t n, R2RKind type1,R2RKind type2,int flags=FFTW_ESTIMATE) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = out = NULL;            
            M=m;
            N=n;           
            in = fftwf_alloc_real(m*n);
            out= fftwf_alloc_real(m*n);                    
            p = fftwf_plan_r2r_2d(m,n, in, out, (fftw_r2r_kind)type1,(fftw_r2r_kind)type2,flags);
        }

        size_t size() const { return M*N; }

        void set_input(sample_matrix<float> & input) {
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                    in[i*N+j] = input(i,j);
        }
        
        sample_matrix<float> get_output() {
            sample_matrix<float> r(M,N);
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) = out[i*N+j];
                }
            return r;
        }        
        void get_output( sample_matrix<float> & output)
        {
            if(output.size() != size()) output.resize(M,N);
            for(size_t i = 0; i < M; i++)
                for(size_t j = 0; j < N; j++)
                {
                    output(i,j) = out[i*N+j];
                }
        }
        void normalize() {
            for(size_t i = 0; i < size(); i++)
                out[i] /= (float)size();    
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
}