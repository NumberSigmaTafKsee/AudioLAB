#pragma once


namespace Casino::IPP
{
    template<typename T>
    void AutoCorrNormGetBufferSize(size_t srcLen, size_t dstLen, IppEnum algType, int * pBufferSize)
    {
        IppDataType dataType = GetDataType<T>();        
        IppStatus status = ippsAutoCorrNormGetBufferSize(srcLen,dstLen,dataType,algType,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void AutoCorrNorm(const T* pSrc, int srcLen, T * pDst, int dstLen,IppEnum algType,Ipp8u * pBuffer)
    {
        std::runtime_error("Called the abstract AutoCorrNorm");                                
    }
    template<>
    void AutoCorrNorm<Ipp32f>(const Ipp32f* pSrc, int srcLen, Ipp32f * pDst, int dstLen,IppEnum algType,Ipp8u * pBuffer)
    {
        IppStatus status;
        status = ippsAutoCorrNorm_32f(pSrc,srcLen,pDst,dstLen,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void AutoCorrNorm<Ipp64f>(const Ipp64f* pSrc, int srcLen, Ipp64f * pDst, int dstLen,IppEnum algType,Ipp8u * pBuffer)
    {
        IppStatus status;
        status = ippsAutoCorrNorm_64f(pSrc,srcLen,pDst,dstLen,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void AutoCorrNorm<Ipp32fc>(const Ipp32fc* pSrc, int srcLen, Ipp32fc * pDst, int dstLen,IppEnum algType,Ipp8u * pBuffer)
    {
        IppStatus status;
        status = ippsAutoCorrNorm_32fc(pSrc,srcLen,pDst,dstLen,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void AutoCorrNorm<Ipp64fc>(const Ipp64fc* pSrc, int srcLen, Ipp64fc * pDst, int dstLen,IppEnum algType,Ipp8u * pBuffer)
    {
        IppStatus status;
        status = ippsAutoCorrNorm_64fc(pSrc,srcLen,pDst,dstLen,algType,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    struct AutoCorrelation
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrelation(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            AutoCorrNormGetBufferSize<T>(srcLen,dstLen,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = Malloc<Ipp8u>(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrelation() {
            if(buffer) Free(buffer);
        }
        void Process(T * src, T * dst) {
            AutoCorrNorm<T>(src,srcLen,dst,dstLen,type,buffer);
        }
    };

    void acorr(size_t srcLen, Ipp32f * src, size_t dstLen, Ipp32f * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrelation<Ipp32f> a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    void acorr(AutoCorrelation<float> &a, Ipp32f * src, Ipp32f * dst)
    {        
        a.Process(src,dst);
    }    
    void acorr(size_t srcLen, Ipp32fc * src, size_t dstLen, Ipp32fc * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrelation<Ipp32fc> a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    void acorr(AutoCorrelation<Ipp32fc> &a,Ipp32fc * src, Ipp32fc * dst)
    {        
        a.Process(src,dst);
    }    
    void acorr(AutoCorrelation<Ipp32fc> &a,std::complex<float> * src, std::complex<float> * dst)
    {        
        a.Process((Ipp32fc*)src,(Ipp32fc*)dst);
    }    
    void acorr(size_t srcLen, Ipp64f * src, size_t dstLen, Ipp64f * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrelation<double> a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    void acorr(AutoCorrelation<double> &a, Ipp64f * src, Ipp64f * dst)
    {        
        a.Process(src,dst);
    }
    
    void acorr(size_t srcLen, Ipp64fc * src, size_t dstLen, Ipp64fc * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrelation<Ipp64fc> a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    void acorr(AutoCorrelation<Ipp64fc> &a,Ipp64fc * src, Ipp64fc * dst)
    {        
        a.Process(src,dst);
    }
    void acorr(AutoCorrelation<Ipp64fc> &a, std::complex<double> * src, std::complex<double> * dst)
    {        
        a.Process((Ipp64fc*)src,(Ipp64fc*)dst);
    }
    
}