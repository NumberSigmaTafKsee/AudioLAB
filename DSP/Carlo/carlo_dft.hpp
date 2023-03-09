#pragma once

namespace Casino::IPP
{

    template<typename T>
    void DFTInitR(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        assert(1==0);
    }
    template<>
    void DFTInitR<float>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_R_32f(length,flag,hint,(IppsDFTSpec_R_32f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitR<double>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_R_64f(length,flag,hint,(IppsDFTSpec_R_64f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<typename T>
    void DFTInitC(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        assert(1==0);
    }
    template<>
    void DFTInitC<float>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_32f(length,flag,hint,(IppsDFTSpec_C_32f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<double>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_64f(length,flag,hint,(IppsDFTSpec_C_64f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<std::complex<float>>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_32fc(length,flag,hint,(IppsDFTSpec_C_32fc*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<std::complex<double>>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_64fc(length,flag,hint,(IppsDFTSpec_C_64fc*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }

    template<typename T>
    void DFTGetSizeR(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        assert(1==0);
    }
    template<>
    void DFTGetSizeR<float>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_R_32f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeR<double>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_R_64f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<typename T>
    void DFTGetSizeC(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        assert(1==0);
    }
    template<>
    void DFTGetSizeC<float>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_32f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<double>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_64f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<std::complex<float>>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_32fc(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<std::complex<double>>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_64fc(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    
    template<typename T>
    void DFTFwd_C2C(const T * pSrcRe, const T * pSrcIm,  T* pDstRe, T* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        // sometime later will try to put some more useful exceptions here
        // it should never be called but I dont think there is a virtual function without a class
        assert(1==0);
    }
    
    template<>
    void DFTFwd_C2C<float>(const float * pSrcRe, const float * pSrcIm, float* pDstRe, float* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_32f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTFwd_C2C<double>(const double * pSrcRe, const double * pSrcIm, double* pDstRe, double* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_64f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_64f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void DFTFwd_C2C(const T * pSrc, T* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        assert(1==0);
    }

    template<>
    void DFTFwd_C2C<std::complex<float>>(const std::complex<float> * pSrc,std::complex<float>* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_32fc((Ipp32fc*)pSrc,(Ipp32fc*)pDst,(IppsDFTSpec_C_32fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTFwd_C2C<std::complex<double>>(const std::complex<double> * pSrc,std::complex<double>* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_64fc((Ipp64fc*)pSrc,(Ipp64fc*)pDst,(IppsDFTSpec_C_64fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void DFTInv_C2C(const T * pSrcRe, const T * pSrcIm,  T* pDstRe, T* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        // sometime later will try to put some more useful exceptions here
        // it should never be called but I dont think there is a virtual function without a class
        assert(1==0);
    }
    template<>
    void DFTInv_C2C<float>(const float * pSrcRe, const float * pSrcIm, float* pDstRe, float* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_32f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<>
    void DFTInv_C2C<double>(const double * pSrcRe, const double * pSrcIm, double* pDstRe, double* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_64f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_64f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void DFTInv_C2C(const T * pSrc, T* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        assert(1==0);
    }

    template<>
    void DFTInv_C2C<std::complex<float>>(const std::complex<float> * pSrc,std::complex<float>* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_32fc((Ipp32fc*)pSrc,(Ipp32fc*)pDst,(IppsDFTSpec_C_32fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTInv_C2C<std::complex<double>>(const std::complex<double> * pSrc,std::complex<double>* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_64fc((Ipp64fc*)pSrc,(Ipp64fc*)pDst,(IppsDFTSpec_C_64fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    // this just sacks massive dong just unpack it please
    template<typename T>
    void DFTFwd_RToPack(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTFwd_RToPack<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPack_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTFwd_RToPack<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPack_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void DFTFwd_RToPerm(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTFwd_RToPerm<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPerm_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTFwd_RToPerm<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPerm_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    
    template<typename T>
    void DFTFwd_RToCCS(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTFwd_RToCCS<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToCCS_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTFwd_RToCCS<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToCCS_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    
    template<typename T>
    void DFTInv_PackToR(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTInv_PackToR<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PackToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_PackToR<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PackToR_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    
    template<typename T>
    void DFTInv_PermToR(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTInv_PermToR<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PermToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_PermToR<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PermToR_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void DFTInv_CCSToR(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTInv_CCSToR<float>(const float * pSrc, float * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_CCSToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_CCSToR<double>(const double * pSrc, double * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_CCSToR_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    // doesn't work with swig at all
    // it's complex that is a complete mess to wrap
    template<typename T>
    struct CDFT
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;        
        void  * fft;
                              
        CDFT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeC<T>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            DFTInitC<T>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~CDFT() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
                        
        virtual void Forward(const T* pSrc, T* pDst)
        {                   
            DFTFwd_C2C<T>(pSrc,pDst,fft,pBuffer);
        }
        virtual void Inverse(const T* pSrc, T* pDst)
        {                   
            DFTInv_C2C<T>(pSrc,pDst,fft,pBuffer);
        }
    };


    
    struct CDFT32 
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;        
        void  * fft;
        size_t blocks;            
        CDFT32(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeC<std::complex<float>>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            blocks = n;
            DFTInitC<std::complex<float>>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~CDFT32() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
                        
        void Forward(const std::complex<float>* pSrc, std::complex<float>* pDst)
        {                   
            DFTFwd_RToCCS<std::complex<float>>((std::complex<float>*)pSrc,(std::complex<float>*)pDst,fft,pBuffer);
        }
        void Inverse(const std::complex<float>* pSrc, std::complex<float>* pDst)
        {                   
            DFTInv_C2C<std::complex<float>>((std::complex<float>*)pSrc,(std::complex<float>*)pDst,fft,pBuffer);
        }
    };
	
	struct RDFT32 
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;        
        void  * fft;
        size_t blocks;            
        RDFT32(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeR<float>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            blocks = n;
            DFTInitR<float>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~RDFT32() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
                        
        void Forward(const float* pSrc, std::complex<float>* pDst)
        {                
			std::vector<float> temp(blocks*2);
            DFTFwd_RToCCS<float>(pSrc,temp.data(),fft,pBuffer);
            for(size_t i = 0; i < blocks; i++)
            {
				pDst[i].real(temp[i*2]);
				pDst[i].imag(temp[i*2+1]);
			}
        }
        void Inverse(const std::complex<float>* pSrc, float* pDst)
        {   
			std::vector<float> temp(blocks*2);
            temp[0] = (pSrc[0].real());
            temp[1] = 0;
            for(size_t i = 1; i < blocks; i++) {
				temp[i*2] =(pSrc[i].real());
				temp[i*2+1] = (pSrc[i].imag());
			}
            DFTInv_CCSToR<float>(temp.data(),pDst,fft,pBuffer);
        }
    };
    
    struct CDFT64 
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;        
        void  * fft;
        size_t blocks;
                        
        CDFT64(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeC<std::complex<double>>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            blocks = n;
            DFTInitC<std::complex<double>>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~CDFT64() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
                        
        void Forward(const std::complex<double>* pSrc, std::complex<double>* pDst)
        {                   
            DFTFwd_C2C<std::complex<double>>((std::complex<double>*)pSrc,(std::complex<double>*)pDst,fft,pBuffer);
        }
        void Inverse(const std::complex<double>* pSrc, std::complex<double>* pDst)
        {                   
            DFTInv_C2C<std::complex<double>>((std::complex<double>*)pSrc,(std::complex<double>*)pDst,fft,pBuffer);
        }
    };
    
    struct RDFT64
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;        
        void  * fft;
        size_t blocks;            
        RDFT64(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeR<double>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            blocks = n;
            DFTInitR<double>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~RDFT64() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);
        }
                        
        void Forward(const double* pSrc, std::complex<double>* pDst)
        {                
			std::vector<double> temp(blocks*2);
            DFTFwd_RToCCS<double>(pSrc,temp.data(),fft,pBuffer);
            for(size_t i = 0; i < blocks; i++)
            {
				pDst[i].real(temp[i*2]);
				pDst[i].imag(temp[i*2+1]);
			}
        }
        void Inverse(const std::complex<double>* pSrc, double* pDst)
        {   
			std::vector<double> temp(blocks*2);
            temp[0]=(pSrc[0].real());
            temp[0]=0;
            for(size_t i = 1; i < blocks; i++) {
				temp[i*2]=(pSrc[i].real());
				temp[i*2+1]=(pSrc[i].imag());
			}
            DFTInv_CCSToR<double>(temp.data(),pDst,fft,pBuffer);
        }
    };
    
	void dft(size_t n, const std::complex<float> * pSrc, std::complex<float> * pDst)
    {
        CDFT32 d(n);
        d.Forward(pSrc,pDst);
    }
    void dft(CDFT32 &d, const std::complex<float> * pSrc, std::complex<float> * pDst)
    {
        d.Forward(pSrc,pDst);
    }
    void idft(size_t n, const std::complex<float> * pSrc, std::complex<float> * pDst)
    {
        CDFT32 d(n);
        d.Inverse(pSrc,pDst);
    }
    void idft(CDFT32 &d,const std::complex<float> * pSrc, std::complex<float> * pDst)
    {        
        d.Inverse(pSrc,pDst);
    }

    void dft(size_t n, const float * pSrc, std::complex<float> * pDst)
    {
        RDFT32 d(n);        
        d.Forward(pSrc,pDst);        
    }
    void dft(RDFT32 &d, const float * pSrc, std::complex<float> * pDst)
    {        
        d.Forward(pSrc,pDst);        
    }
    void idft(size_t n, const std::complex<float> * pSrc, float * pDst)
    {
        RDFT32 d(n);
        d.Inverse(pSrc,pDst);        
    }
    void idft(RDFT32 &d, const std::complex<float> * pSrc, float * pDst)
    {        
        d.Inverse(pSrc,pDst);
    }

    void dft(size_t n, const std::complex<double> * pSrc, std::complex<double> * pDst)
    {
        CDFT64 d(n);
        d.Forward(pSrc,pDst);
    }
    void dft(CDFT64 &d, const std::complex<double> * pSrc, std::complex<double> * pDst)
    {
        d.Forward(pSrc,pDst);
    }
    void idft(size_t n, const std::complex<double> * pSrc, std::complex<double> * pDst)
    {
        CDFT64 d(n);
        d.Inverse(pSrc,pDst);
    }
    void idft(CDFT64 &d,const std::complex<double> * pSrc, std::complex<double> * pDst)
    {        
        d.Inverse(pSrc,pDst);
    }    


    void dft(size_t n, const double * pSrc, std::complex<double> * pDst)
    {
        RDFT64 d(n);        
        d.Forward(pSrc,pDst);        
    }
    void dft(RDFT64 &d, const double * pSrc, std::complex<double> * pDst)
    {                
        d.Forward(pSrc,pDst);        
    }
    void idft(size_t n, const std::complex<double> * pSrc, double * pDst)
    {
        RDFT64 d(n);        
        d.Inverse(pSrc,pDst);        
    }
    void idft(RDFT64 &d, const std::complex<double> * pSrc, double * pDst)
    {                
        d.Inverse(pSrc,pDst);        
    }   
}
