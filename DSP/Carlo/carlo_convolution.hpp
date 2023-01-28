#pragma once

namespace Casino::IPP
{
    template<typename T>
    void ConvolveGetBufferSize(int src1Len, int src2Len, IppEnum algType, int * pBufferSize)
    {
        IppDataType dataType = GetDataType<T>();        
        ippsConvolveGetBufferSize(src1Len,src2Len,dataType,algType,pBufferSize);
    }
    template<typename T>
    void Convolve(const T* pSrc1, int src1Len, const T* pSrc2, int src2Len, T* pDst, IppEnum algType, Ipp8u* pBuffer)
    {
        throw std::runtime_error("Called abstract CrossCorrNorm");
    }
    template<>
    void Convolve<Ipp32f>(const Ipp32f* pSrc1, int src1Len, const Ipp32f* pSrc2, int src2Len, Ipp32f* pDst, IppEnum algType, Ipp8u* pBuffer)
    {
        IppStatus status = ippsConvolve_32f(pSrc1,src1Len,pSrc2,src2Len,pDst,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void Convolve<Ipp64f>(const Ipp64f* pSrc1, int src1Len, const Ipp64f* pSrc2, int src2Len, Ipp64f* pDst, IppEnum algType, Ipp8u* pBuffer)
    {
        IppStatus status = ippsConvolve_64f(pSrc1,src1Len,pSrc2,src2Len,pDst,algType,pBuffer);
        checkStatus(status);
    }

    template<typename T>
    void ConvBiased(const T* pSrc1, int src1Len, const T* pSrc2, int src2Len, T* pDst, int dstLen, int bias)
    {
        throw std::runtime_error("Called abstract CrossCorrNorm");
    }
    template<>
    void ConvBiased<Ipp32f>(const Ipp32f* pSrc1, int src1Len, const Ipp32f* pSrc2, int src2Len, Ipp32f* pDst, int dstLen, int bias)
    {
        IppStatus status = ippsConvBiased_32f(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,bias);
        checkStatus(status);
    }
    
    template<typename T>
    struct Convolver
    {
        Ipp8u * buffer;        
        int bufferLen;
        int      type;
        int src1Len,dstLen,src2Len;

        Convolver(size_t src1Len,size_t src2Len, int algorithm = (int)Algorithm::ALG_AUTO) {
            ConvolveGetBufferSize<T>(src1Len,src2Len,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = Malloc<Ipp8u>(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;        
        }
        ~Convolver() {
            if(buffer) ippsFree(buffer);
        }
        void Process(T * src1, T * src2, T * dst) {
            Convolve<T>(src1,src1Len,src2,src2Len,dst,type,buffer);            
        }
    };

    template<typename T>
    struct ConvolutionFilter
    {
        Convolver<T> *filter;
        std::shared_ptr<T> h;
        size_t            len,block;
        std::vector<T> ola,temp;

        ConvolutionFilter(size_t n, T * impulse, size_t block_size)
        {
            filter = new Convolver<T>(block_size,block_size);
            block = block_size;
            T *hi = Malloc<T>(len);
            h = std::shared_ptr<T>(hi,[](T * p) { Free(p); }); 
            len = n;
            assert(filter != NULL);
            assert(hi != nullptr);
            ola.resize(block_size);
            temp.resize(block_size + n -1);
        }        
        void ProcessBlock(T * signal, T * dest) {
            filter->Process(h.get(),signal,dest);
            for(size_t i = 0; i < block; i++) {
                dest[i] = temp[i] + ola[i];
                ola[i]  = temp[i+block];
            }
        }
    };

    void conv(size_t src1Len, Ipp32f * src1, size_t src2Len, Ipp32f * src2, Ipp32f * dst, int algorithm = (int)Algorithm::ALG_AUTO)
    {
        Convolver<Ipp32f> c(src1Len,src2Len,algorithm);
        c.Process(src1,src2,dst);
    }
    void conv(Convolver<Ipp32f> &c, Ipp32f * src1, Ipp32f * src2, Ipp32f * dst)
    {        
        c.Process(src1,src2,dst);
    }
    void conv(size_t src1Len, Ipp64f * src1, size_t src2Len, Ipp64f * src2, Ipp64f * dst, int algorithm = (int)Algorithm::ALG_AUTO)
    {
        Convolver<Ipp64f> c(src1Len,src2Len,algorithm);
        c.Process(src1,src2,dst);
    }
    void conv(Convolver<Ipp64f> &c, Ipp64f * src1, Ipp64f * src2, Ipp64f * dst)
    {    
        c.Process(src1,src2,dst);
    }

}