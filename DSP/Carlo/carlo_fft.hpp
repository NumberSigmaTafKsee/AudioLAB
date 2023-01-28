#pragma once

namespace Casino::IPP
{

    
    template<typename T>
    void FFTInitR(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTInitR<float>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_R_32f((IppsFFTSpec_R_32f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitR<double>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_R_64f((IppsFFTSpec_R_64f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }

    template<typename T>
    void FFTInitC(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    
    template<>
    void FFTInitC<float>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_32f((IppsFFTSpec_C_32f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitC<double>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_64f((IppsFFTSpec_C_64f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    
    template<>
    void FFTInitC<std::complex<float>>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_32fc((IppsFFTSpec_C_32fc**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitC<std::complex<double>>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_64fc((IppsFFTSpec_C_64fc**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }

    template<typename T>
    void FFTFwd_RToPack(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTFwd_RToPack<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPack_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToPack<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPack_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }


    template<typename T>
    void FFTFwd_RToPerm(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTFwd_RToPerm<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPerm_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToPerm<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPerm_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void FFTFwd_RToCCS(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTFwd_RToCCS<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToCCS_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToCCS<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToCCS_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    
    
    template<typename T>
    void FFTInv_PackToR(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    
    template<>
    void FFTInv_PackToR<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PackToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_PackToR<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PackToR_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void FFTInv_PermToR(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTInv_PermToR<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PermToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_PermToR<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PermToR_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void FFTInv_CCSToR(const T * pSrc, T * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTInv_CCSToR<float>(const float * pSrc, float * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_CCSToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_CCSToR<double>(const double * pSrc, double * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_CCSToR_64f(pSrc,pDst,(IppsFFTSpec_R_64f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<typename T>
    void FFTGetSizeR(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void FFTGetSizeR<float>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_R_32f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeR<double>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_R_64f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }

    template<typename T>
    void FFTGetSizeC(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void FFTGetSizeC<float>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_32f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<double>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_64f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<std::complex<float>>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_32fc(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<std::complex<double>>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_64fc(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }

    template<typename T>
    void FFTFwd_C2C(const T * pSrcRe, const T * pSrcIm,  T* pDstRe, T* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        // sometime later will try to put some more useful exceptions here
        // it should never be called but I dont think there is a virtual function without a class
        assert(1==0);
    }
    
    template<>
    void FFTFwd_C2C<float>(const float * pSrcRe, const float * pSrcIm, float* pDstRe, float* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_32f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTFwd_C2C<double>(const double * pSrcRe, const double * pSrcIm, double* pDstRe, double* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_64f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_64f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void FFTFwd_C2C(const T * pSrc, T* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        assert(1==0);
    }

    template<>
    void FFTFwd_C2C<std::complex<float>>(const std::complex<float> * pSrc,std::complex<float>* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_32fc((Ipp32fc*)pSrc,(Ipp32fc*)pDst,(IppsFFTSpec_C_32fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTFwd_C2C<std::complex<double>>(const std::complex<double> * pSrc,std::complex<double>* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_64fc((Ipp64fc*)pSrc,(Ipp64fc*)pDst,(IppsFFTSpec_C_64fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void FFTInv_C2C(const T * pSrcRe, const T * pSrcIm,  T* pDstRe, T* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        // sometime later will try to put some more useful exceptions here
        // it should never be called but I dont think there is a virtual function without a class
        assert(1==0);
    }
    template<>
    void FFTInv_C2C<float>(const float * pSrcRe, const float * pSrcIm, float* pDstRe, float* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_32f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<>
    void FFTInv_C2C<double>(const double * pSrcRe, const double * pSrcIm, double* pDstRe, double* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_64f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_64f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void FFTInv_C2C(const T * pSrc, T* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        assert(1==0);
    }

    template<>
    void FFTInv_C2C<std::complex<float>>(const std::complex<float> * pSrc,std::complex<float>* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_32fc((Ipp32fc*)pSrc,(Ipp32fc*)pDst,(IppsFFTSpec_C_32fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTInv_C2C<std::complex<double>>(const std::complex<double> * pSrc,std::complex<double>* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_64fc((Ipp64fc*)pSrc,(Ipp64fc*)pDst,(IppsFFTSpec_C_64fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    
    template<typename T>
    struct CFFT
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        void * fft;
        
        CFFT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            FFTGetSizeC<T>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = std::log2(n);
            FFTInitC<T>(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);                        
        }
        ~CFFT() {
            if(pSpec) Free(pSpec);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
        void Forward(const T* pSrc, T * pDst)
        {                               
            FFTFwd_C2C<T>(pSrc,pDst, fft, pBuffer);                                    
        }
        void Inverse(const T* pSrc, T * pDst)
        {                               
            FFTInv_C2C<T>(pSrc,pDst, fft, pBuffer);            
        }        
    };

    
    struct CFFT32 
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        void * fft;
        size_t blocks;
        
        CFFT32(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            FFTGetSizeC<std::complex<float>>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = std::log2(n);
            blocks = n;
            FFTInitC<std::complex<float>>(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);                        
        }
        ~CFFT32() {
            if(pSpec) Free(pSpec);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
        void Forward(const std::complex<float>* pSrc, std::complex<float> * pDst)
        {                               
            FFTFwd_C2C<std::complex<float>>(pSrc,pDst, fft, pBuffer);                                    
        }
        void Inverse(const std::complex<float>* pSrc, std::complex<float> * pDst)
        {                               
            FFTInv_C2C<std::complex<float>>(pSrc,pDst, fft, pBuffer);            
        }        
    };

    struct CFFT64
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        void * fft;
        size_t blocks;
        
        CFFT64(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            FFTGetSizeC<std::complex<double>>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = std::log2(n);
            blocks = n;
            FFTInitC<std::complex<double>>(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);                        
        }
        ~CFFT64() {
            if(pSpec) Free(pSpec);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
        void Forward(const std::complex<double>* pSrc, std::complex<double> * pDst)
        {                               
            FFTFwd_C2C<std::complex<double>>(pSrc,pDst, fft, pBuffer);                                    
        }
        void Inverse(const std::complex<double>* pSrc, std::complex<double> * pDst)
        {                               
            FFTInv_C2C<std::complex<double>>(pSrc,pDst, fft, pBuffer);            
        }        
    };

    
    void fft(size_t n, const std::complex<float> * pSrc, std::complex<float> * pDst)
    {
        CFFT32 f(n);
        f.Forward(pSrc,pDst);
    }
    void fft(CFFT32 &f,const std::complex<float> * pSrc, std::complex<float> * pDst)
    {                
        f.Forward(pSrc,pDst);
    }
    void ifft(size_t n, const std::complex<float> * pSrc, std::complex<float> * pDst)
    {
        CFFT32 f(n);
        f.Inverse(pSrc,pDst);
    }
    void ifft(CFFT32 &f,const std::complex<float> * pSrc, std::complex<float> * pDst)
    {        
        f.Inverse(pSrc,pDst);
    }
    
        void fft(size_t n, const float * pSrc, std::complex<float> * pDst)
    {
        CFFT32 d(n);
        std::vector<std::complex<float>> in(n);
        for(size_t i = 0; i < n; i++) {
            in[i].real(pSrc[i]);
            in[i].imag(0);
        }
        d.Forward(in.data(),pDst);        
    }
    void fft(CFFT32 &d, const float * pSrc, std::complex<float> * pDst)
    {        
        std::vector<std::complex<float>> in(d.blocks);
        for(size_t i = 0; i < d.blocks; i++) {
            in[i].real(pSrc[i]);
            in[i].imag(0);
        }
        d.Forward(in.data(),pDst);        
    }
    void ifft(size_t n, const std::complex<float> * pSrc, float * pDst)
    {
        CFFT32 d(n);
        std::vector<std::complex<float>> out(n);
        d.Inverse(pSrc,out.data());
        for(size_t i = 0; i < n; i++) {
            pDst[i] = (out[i].real());            
        }        
    }
    void ifft(CFFT32 &d, const std::complex<float> * pSrc, float * pDst)
    {        
        std::vector<std::complex<float>> out(d.blocks);
        d.Inverse(pSrc,out.data());
        for(size_t i = 0; i < d.blocks; i++) {
            pDst[i] = (out[i].real());
        }
    }   

    void fft(size_t n, const std::complex<double> * pSrc, std::complex<double> * pDst)
    {
        CFFT64 f(n);
        f.Forward(pSrc,pDst);
    }
    void fft(CFFT64 &f,const std::complex<double> * pSrc, std::complex<double> * pDst)
    {        
        f.Forward(pSrc,pDst);
    }
    void ifft(size_t n, const std::complex<double> * pSrc, std::complex<double> * pDst)
    {
        CFFT64 f(n);
        f.Inverse(pSrc,pDst);
    }
    void ifft(CFFT64 &f, std::complex<double>* pSrc, std::complex<double> * pDst)
    {        
        f.Inverse(pSrc,pDst);
    }

    void fft(size_t n, const double * pSrc, std::complex<double> * pDst)
    {
        CFFT64 d(n);
        std::vector<std::complex<double>> in(n);
        for(size_t i = 0; i < n; i++) {
            in[i].real(pSrc[i]);
            in[i].imag(0);
        }
        d.Forward(in.data(),pDst);        
    }
    void fft(CFFT64 &d, const double * pSrc, std::complex<double> * pDst)
    {        
        std::vector<std::complex<double>> in(d.blocks);
        for(size_t i = 0; i < d.blocks; i++) {
            in[i].real(pSrc[i]);
            in[i].imag(0);
        }
        d.Forward(in.data(),pDst);        
    }
    void ifft(size_t n, const std::complex<double> * pSrc, double * pDst)
    {
        CFFT64 d(n);
        std::vector<std::complex<double>> out(n);
        d.Inverse(pSrc,out.data());
        for(size_t i = 0; i < n; i++) {
            pDst[i] = (out[i].real());            
        }        
    }
    void ifft(CFFT64 &d, const std::complex<double> * pSrc, double * pDst)
    {        
        std::vector<std::complex<double>> out(d.blocks);
        d.Inverse(pSrc,out.data());
        for(size_t i = 0; i < d.blocks; i++) {
            pDst[i] = (out[i].real());
        }
    }   
}