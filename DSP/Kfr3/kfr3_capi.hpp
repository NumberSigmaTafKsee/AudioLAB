#pragma once

namespace kfr3
{
        enum KfrType {
        DFT,
        DFT_REAL,
        DCT,
        FIR,
        IIR,
        CONVOLUTION,
        RESAMPLER,
    };

    
    struct KFRBase
    {
        virtual KfrType getType() const = 0;
        virtual void*   getPlan() = 0;
    };

    
    struct DFTPlanFloat : public KFRBase
    {
        KFR_DFT_PLAN_F32 * plan;
        sample_vector<uint8_t> temp;
        
        DFTPlanFloat(size_t n) {
            plan = kfr_dft_create_plan_f32(n);        
            temp.resize(kfr_dft_get_temp_size_f32(plan));
        }
        ~DFTPlanFloat() {
            if(plan) kfr_dft_delete_plan_f32(plan);
        }        

        KfrType getType() { return DFT; }
        void*   getPlan() { return (void*)plan; }
    };

    
    struct DFTRealPlanFloat : public KFRBase
    {
        KFR_DFT_REAL_PLAN_F32 * plan;
        sample_vector<uint8_t> temp;
        
        KfrType getType() { return DFT_REAL; }
        void*   getPlan() { return (void*)plan; }

        DFTRealPlanFloat(size_t n, KFR_DFT_PACK_FORMAT flags=Perm) {
            plan = kfr_dft_real_create_plan_f32(n,flags);        
            temp.resize(kfr_dft_real_get_temp_size_f32(plan));
        }
        ~DFTRealPlanFloat() {
            if(plan) kfr_dft_real_delete_plan_f32(plan);
        }        
    };

    
    struct DCTPlanFloat : public KFRBase
    {
        KFR_DCT_PLAN_F32 * plan;
        sample_vector<uint8_t> temp;    

        KfrType getType() { return DCT; }
        void*   getPlan() { return (void*)plan; }

        DCTPlanFloat(size_t n) {  
            plan = kfr_dct_create_plan_f32(n);
            temp.resize(kfr_dct_get_temp_size_f32(plan));
        }
        ~DCTPlanFloat() {
            if(plan) kfr_dct_delete_plan_f32(plan);
        }        
    };

    // kfr filter not good for music they are not time-variant
    // they are good if the parameters never change.    
    struct FIRFilterFloat : public KFRBase
    {
        KFR_FILTER_F32 * filter;

        KfrType getType() { return FIR; }
        void*   getPlan() { return (void*)filter; }

        FIRFilterFloat(const float* taps, size_t size)
        {
            filter = kfr_filter_create_fir_plan_f32(taps,size);
        }
        ~FIRFilterFloat() {
            if(filter) kfr_filter_delete_plan_f32(filter);
        }
    };  

    struct IIRFilterFloat : public KFRBase
    {
        KFR_FILTER_F32 * filter;

        KfrType getType() { return IIR; }
        void*   getPlan() { return (void*)filter; }

        IIRFilterFloat(const float* sos, size_t size)
        {
            filter = kfr_filter_create_iir_plan_f32(sos,size);
        }
        ~IIRFilterFloat() {
            if(filter) kfr_filter_delete_plan_f32(filter);
        }
    };    

    struct ConvolutionFilterFloat
    {
        KFR_FILTER_F32* filter;
        
        KfrType getType() { return CONVOLUTION; }
        void*   getPlan() { return (void*)filter; }
        
        ConvolutionFilterFloat(float * taps, size_t size, size_t block_size)
        {
            kfr_filter_create_convolution_plan_f32(taps,size,block_size);
        }
        ~ConvolutionFilterFloat() {
            if(filter) kfr_filter_delete_plan_f32(filter);
        }
    };
    

    struct DFTPlanDouble : public KFRBase
    {
        KFR_DFT_PLAN_F64 * plan;
        sample_vector<uint8_t> temp;
        
        DFTPlanDouble(size_t n) {
            plan = kfr_dft_create_plan_f64(n);        
            temp.resize(kfr_dft_get_temp_size_f64(plan));
        }
        ~DFTPlanDouble() {
            if(plan) kfr_dft_delete_plan_f64(plan);
        }        

        KfrType getType() { return DFT; }
        void*   getPlan() { return (void*)plan; }
    };

    
    struct DFTRealPlanDouble : public KFRBase
    {
        KFR_DFT_REAL_PLAN_F64 * plan;
        sample_vector<uint8_t> temp;
        
        KfrType getType() { return DFT_REAL; }
        void*   getPlan() { return (void*)plan; }

        DFTRealPlanDouble(size_t n, KFR_DFT_PACK_FORMAT flags=Perm) {
            plan = kfr_dft_real_create_plan_f64(n,flags);        
            temp.resize(kfr_dft_real_get_temp_size_f64(plan));
        }
        ~DFTRealPlanDouble() {
            if(plan) kfr_dft_real_delete_plan_f64(plan);
        }        
    };

    
    struct DCTPlanDouble : public KFRBase
    {
        KFR_DCT_PLAN_F64 * plan;
        sample_vector<uint8_t> temp;    

        KfrType getType() { return DCT; }
        void*   getPlan() { return (void*)plan; }

        DCTPlanDouble(size_t n) {  
            plan = kfr_dct_create_plan_f64(n);
            temp.resize(kfr_dct_get_temp_size_f64(plan));
        }
        ~DCTPlanDouble() {
            if(plan) kfr_dct_delete_plan_f64(plan);
        }        
    };

    // kfr filter not good for music they are not time-variant
    // they are good if the parameters never change.    
    struct FIRFilterDouble : public KFRBase
    {
        KFR_FILTER_F64 * filter;

        KfrType getType() { return FIR; }
        void*   getPlan() { return (void*)filter; }

        FIRFilterDouble(const double* taps, size_t size)
        {
            filter = kfr_filter_create_fir_plan_f64(taps,size);
        }
        ~FIRFilterDouble() {
            if(filter) kfr_filter_delete_plan_f64(filter);
        }
    };  

    struct IIRFilterDouble : public KFRBase
    {
        KFR_FILTER_F64 * filter;

        KfrType getType() { return IIR; }
        void*   getPlan() { return (void*)filter; }

        IIRFilterDouble(const double* sos, size_t size)
        {
            filter = kfr_filter_create_iir_plan_f64(sos,size);
        }
        ~IIRFilterDouble() {
            if(filter) kfr_filter_delete_plan_f64(filter);
        }
    };    

    struct ConvolutionFilterDouble
    {
        KFR_FILTER_F64* filter;
        
        KfrType getType() { return CONVOLUTION; }
        void*   getPlan() { return (void*)filter; }
        
        ConvolutionFilterDouble(double * taps, size_t size, size_t block_size)
        {
            kfr_filter_create_convolution_plan_f64(taps,size,block_size);
        }
        ~ConvolutionFilterDouble() {
            if(filter) kfr_filter_delete_plan_f64(filter);
        }
    };

    void dump(DFTPlanFloat & plan)
    {
        kfr_dft_dump_f32(plan.plan);
    }
    size_t get_size(DFTPlanFloat & plan)
    {
        return kfr_dft_get_size_f32(plan.plan);
    }
    size_t get_temp_size(DFTPlanFloat & plan)
    {
        return kfr_dft_get_temp_size_f32(plan.plan);
    }
    void forward(DFTPlanFloat & plan, std::complex<float> * out, std::complex<float> * in)
    {        
        kfr_dft_execute_f32(plan.plan,(kfr_c32*)out,(kfr_c32*)in,plan.temp.data());
    }
    void inverse(DFTPlanFloat & plan, std::complex<float> * out, std::complex<float> * in)
    {        
        kfr_dft_execute_inverse_f32(plan.plan,(kfr_c32*)out,(kfr_c32*)in,plan.temp.data());
    }
    void dft(DFTPlanFloat & plan, std::complex<float> * out, std::complex<float> * in)
    {        
        kfr_dft_execute_f32(plan.plan,(kfr_c32*)out,(kfr_c32*)in,plan.temp.data());
    }
    void idft(DFTPlanFloat & plan, std::complex<float> * out, std::complex<float> * in)
    {        
        kfr_dft_execute_inverse_f32(plan.plan,(kfr_c32*)out,(kfr_c32*)in,plan.temp.data());
    }



    void dump(DFTPlanDouble & plan)
    {
        kfr_dft_dump_f64(plan.plan);
    }
    size_t get_size(DFTPlanDouble & plan)
    {
        return kfr_dft_get_size_f64(plan.plan);
    }
    size_t get_temp_size(DFTPlanDouble & plan)
    {
        return kfr_dft_get_temp_size_f64(plan.plan);
    }
    void forward(DFTPlanDouble & plan, std::complex<double> * out, std::complex<double> * in)
    {        
        kfr_dft_execute_f64(plan.plan,(kfr_c64*)out,(kfr_c64*)in,plan.temp.data());
    }
    void inverse(DFTPlanDouble & plan, std::complex<double> * out, std::complex<double> * in)
    {        
        kfr_dft_execute_inverse_f64(plan.plan,(kfr_c64*)out,(kfr_c64*)in,plan.temp.data());
    }
    void dft(DFTPlanDouble & plan, std::complex<double> * out, std::complex<double> * in)
    {        
        kfr_dft_execute_f64(plan.plan,(kfr_c64*)out,(kfr_c64*)in,plan.temp.data());
    }
    void idft(DFTPlanDouble & plan, std::complex<double> * out, std::complex<double> * in)
    {        
        kfr_dft_execute_inverse_f64(plan.plan,(kfr_c64*)out,(kfr_c64*)in,plan.temp.data());
    }

    void dump(DFTRealPlanFloat & plan)
    {
        kfr_dft_real_dump_f32(plan.plan);
    }
    size_t get_size(DFTRealPlanFloat & plan)
    {
        return kfr_dft_real_get_size_f32(plan.plan);
    }
    size_t get_temp_size(DFTRealPlanFloat & plan)
    {
        return kfr_dft_real_get_temp_size_f32(plan.plan);
    }
    void forward(DFTRealPlanFloat & plan, float * out, std::complex<float> * in)
    {        
        kfr_dft_real_execute_f32(plan.plan,out,(kfr_c32*)in,plan.temp.data());
    }
    void inverse(DFTRealPlanFloat & plan, std::complex<float> * out, float * in)
    {        
        kfr_dft_real_execute_inverse_f32(plan.plan,(kfr_c32*)out,in,plan.temp.data());
    }
    void dft(DFTRealPlanFloat & plan, std::complex<float> * out, float * in)
    {        
        kfr_dft_real_execute_f32(plan.plan,(kfr_c32*)out,in,plan.temp.data());
    }
    void idft(DFTRealPlanFloat & plan, float * out, std::complex<float> * in)
    {     
        kfr_dft_real_execute_inverse_f32(plan.plan,out,(kfr_c32*)in,plan.temp.data());
    }

    void dump(DFTRealPlanDouble & plan)
    {
        kfr_dft_real_dump_f64(plan.plan);
    }
    size_t get_size(DFTRealPlanDouble & plan)
    {
        return kfr_dft_real_get_size_f64(plan.plan);
    }
    size_t get_temp_size(DFTRealPlanDouble & plan)
    {
        return kfr_dft_real_get_temp_size_f64(plan.plan);
    }
    void forward(DFTRealPlanDouble & plan, double * out, std::complex<double> * in)
    {        
        kfr_dft_real_execute_f64(plan.plan,out,(kfr_c64*)in,plan.temp.data());
    }
    void inverse(DFTRealPlanDouble & plan, std::complex<double> * out, double * in)
    {        
        kfr_dft_real_execute_inverse_f64(plan.plan,(kfr_c64*)out,in,plan.temp.data());
    }
    void dft(DFTRealPlanDouble & plan, std::complex<double> * out, double * in)
    {        
        kfr_dft_real_execute_f64(plan.plan,(kfr_c64*)out,in,plan.temp.data());
    }
    void idft(DFTRealPlanDouble & plan, double * out, std::complex<double> * in)
    {        
        kfr_dft_real_execute_inverse_f64(plan.plan,out,(kfr_c64*)in,plan.temp.data());
    }

    void dump(DCTPlanFloat & plan)
    {
        kfr_dct_dump_f32(plan.plan);
    }
    size_t get_size(DCTPlanFloat & plan)
    {
        return kfr_dct_get_size_f32(plan.plan);
    }
    size_t get_temp_size(DCTPlanFloat & plan)
    {
        return kfr_dct_get_temp_size_f32(plan.plan);
    }
    void forward(DCTPlanFloat & plan, float * out, float * in)
    {        
        kfr_dct_execute_f32(plan.plan,out,in,plan.temp.data());
    }
    void inverse(DCTPlanFloat & plan, float * out, float * in)
    {        
        kfr_dct_execute_inverse_f32(plan.plan,out,in,plan.temp.data());
    }


    void filter(FIRFilterFloat& filter, float * output, const float * input, size_t size)
    {
        kfr_filter_process_f32(filter.filter,output,input,size);
    }
    void filter(IIRFilterFloat& filter, float * output, const float * input, size_t size)
    {
        kfr_filter_process_f32(filter.filter,output,input,size);
    }
    void filter(ConvolutionFilterFloat& filter, float * output, const float * input, size_t size)
    {
        kfr_filter_process_f32(filter.filter,output,input,size);
    }
    void reset(FIRFilterFloat& filter)
    {
        kfr_filter_reset_f32(filter.filter);
    }
    void reset(IIRFilterFloat& filter)
    {
        kfr_filter_reset_f32(filter.filter);
    }
    void reset(ConvolutionFilterFloat& filter)
    {
        kfr_filter_reset_f32(filter.filter);
    }

    void dump(DCTPlanDouble & plan)
    {
        kfr_dct_dump_f64(plan.plan);
    }
    size_t get_size(DCTPlanDouble & plan)
    {
        return kfr_dct_get_size_f64(plan.plan);
    }
    size_t get_temp_size(DCTPlanDouble & plan)
    {
        return kfr_dct_get_temp_size_f64(plan.plan);
    }
    void forward(DCTPlanDouble & plan, double * out, double * in)
    {        
        kfr_dct_execute_f64(plan.plan,out,in,plan.temp.data());
    }
    void inverse(DCTPlanDouble & plan, double * out, double * in)
    {        
        kfr_dct_execute_inverse_f64(plan.plan,out,in,plan.temp.data());
    }


    void filter(FIRFilterDouble& filter, double * output, const double * input, size_t size)
    {
        kfr_filter_process_f64(filter.filter,output,input,size);
    }
    void filter(IIRFilterDouble& filter, double * output, const double * input, size_t size)
    {
        kfr_filter_process_f64(filter.filter,output,input,size);
    }
    void filter(ConvolutionFilterDouble& filter, double * output, const double * input, size_t size)
    {
        kfr_filter_process_f64(filter.filter,output,input,size);
    }
    void reset(FIRFilterDouble& filter)
    {
        kfr_filter_reset_f64(filter.filter);
    }
    void reset(IIRFilterDouble& filter)
    {
        kfr_filter_reset_f64(filter.filter);
    }
    void reset(ConvolutionFilterDouble& filter)
    {
        kfr_filter_reset_f64(filter.filter);
    }


    template<typename T>
    sample_vector<T> convolve(size_t n, const T *a, const T * b) {
        kfr::univector<T> v1(n);
        kfr::univector<T> v2(n);
        memcpy(v1.data(),a,sizeof(T)*n);
        memcpy(v2.data(),b,sizeof(T)*n);
        kfr::univector<T> r = (kfr::convolve(v1,v2));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }

    template<typename T>
    sample_vector<T> convolve(const sample_vector<T> &a, const sample_vector<T> &b) {        
        assert(a.size() == b.size());
        kfr::univector<T> r = (kfr::convolve(a,b));        
        return sample_vector<T>(r);
    }
    
    // there isn't lag in kfr
    template<typename T>
    sample_vector<T> xcorr(size_t n, const T *src1, const T * src2) {
        kfr::univector<T> v1(n);
        kfr::univector<T> v2(n);
        memcpy(v1.data(),src1.data(),sizeof(T)*n);
        memcpy(v2.data(),src2.data(),sizeof(T)*n);
        kfr::univector<T> r = (kfr::correlate(v1,v2));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }
    
    template<typename T>
    sample_vector<T> xcorr(const sample_vector<T> &a, const sample_vector<T> &b) {        
        assert(a.size() == b.size());
        kfr::univector<T> r = (kfr::correlate(a,b));        
        return sample_vector<T>(r);
    }
    
    // there isn't lag in kfr
    template<typename T>
    sample_vector<T> acorr(size_t n, const T *src) {
        kfr::univector<T> v1(n);        
        memcpy(v1.data(),src.data(),sizeof(T)*n);        
        kfr::univector<T> r = (kfr::autocorrelate(v1));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }
    
    template<typename T>
    sample_vector<T> acorr(const sample_vector<T> &a) {        
        kfr::univector<T> r = (kfr::autocorrelate(a));        
        return sample_vector<T>(r);
    }

}