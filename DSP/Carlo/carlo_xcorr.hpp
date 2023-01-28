#pragma once

namespace Casino::IPP
{
    template<typename T>
    void CrossCorrNormGetBufferSize(int src1Len, int src2Len, int dstLen, int lowLag, IppEnum algType, int * pBufferSize)
    {
        IppDataType dataType = GetDataType<T>();        
        ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,dataType,algType,pBufferSize);    
    }

    template<typename T>
    void CrossCorrNorm(const T * pSrc1, int src1Len, const T * pSrc2, int src2Len, T * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        throw std::runtime_error("Called abstract CrossCorrNorm");
    }
    
    template<>
    void CrossCorrNorm<float>(const float * pSrc1, int src1Len, const float * pSrc2, int src2Len, float * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_32f(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<double>(const double * pSrc1, int src1Len, const double * pSrc2, int src2Len, double * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_64f(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<std::complex<float>>(const std::complex<float> * pSrc1, int src1Len, const std::complex<float> * pSrc2, int src2Len, std::complex<float> * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_32fc((Ipp32fc*)pSrc1,src1Len,(Ipp32fc*)pSrc2,src2Len,(Ipp32fc*)pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<std::complex<double>>(const std::complex<double> * pSrc1, int src1Len, const std::complex<double> * pSrc2, int src2Len, std::complex<double> * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_64fc((Ipp64fc*)pSrc1,src1Len,(Ipp64fc*)pSrc2,src2Len,(Ipp64fc*)pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    
    template<typename T>
    struct CrossCorrelation
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrelation(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            CrossCorrNormGetBufferSize<T>(src1Len,src2Len,dstLen,lowLag,algorithm,&bufferLen);        
            if(bufferLen > 0) buffer = Malloc<Ipp8u>(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->lowLag = lowLag;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
        }
        ~CrossCorrelation() {
            if(buffer) Free(buffer);
        }
        void Process(T * src1, T * src2, T * dst) {
            CrossCorrNorm<T>(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);            
        }
    };

    void xcorr(size_t srcLen, float * src1, size_t srcLen2, float* src2, size_t dstLen, float * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrelation<float> c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    void xcorr(CrossCorrelation<float> &c, float * src1, float* src2, size_t dstLen, float * dst)
    {        
        c.Process(src1,src2,dst);
    }        
    void xcorr(size_t srcLen, std::complex<float> * src1, size_t srcLen2, std::complex<float>* src2, size_t dstLen, std::complex<float> * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrelation<std::complex<float>> c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }            
    void xcorr(CrossCorrelation<std::complex<float>> &c, std::complex<float> * src1, std::complex<float>* src2, std::complex<float> * dst)
    {    
        c.Process(src1,src2,dst);
    }            
    void xcorr(size_t srcLen, double * src1, size_t srcLen2, double* src2, size_t dstLen, double * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrelation<double> c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    void xcorr(CrossCorrelation<double> &c, double * src1, double* src2, double * dst)
    {        
        c.Process(src1,src2,dst);
    }        
    void xcorr(size_t srcLen, std::complex<double> * src1, size_t srcLen2, std::complex<double>* src2, size_t dstLen, std::complex<double> * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrelation<std::complex<double>> c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    void xcorr(CrossCorrelation<std::complex<double>> &c, std::complex<double> * src1, std::complex<double>* src2, std::complex<double> * dst)
    {    
        c.Process(src1,src2,dst);
    }        


}