    
    
    


    
    
    
    

    
    
    
    

    
        struct AutoCorrDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrDouble(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src, Ipp64f * dst) {
            IppStatus status = ippsAutoCorrNorm_64f(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    struct AutoCorrDoubleComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrDoubleComplex(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp64fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrDoubleComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64fc * src, Ipp64fc * dst) {
            IppStatus status = ippsAutoCorrNorm_64fc(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    
    struct CrossCorrDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrDouble(size_t src1Len,  size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
            this->lowLag = lowLag;
        }
        ~CrossCorrDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src1, Ipp64f * src2, Ipp64f * dst) {
            IppStatus status = ippsCrossCorrNorm_64f(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    struct CrossCorrDoubleComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrDoubleComplex(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp64fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
            this->lowLag = lowLag;
        }
        ~CrossCorrDoubleComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64fc * src1, Ipp64fc * src2, Ipp64fc * dst) {
            IppStatus status = ippsCrossCorrNorm_64fc(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    
    struct ConvolveDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,dstLen,src2Len;

        ConvolveDouble(size_t src1Len,  size_t src2Len, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsConvolveGetBufferSize(src1Len,src2Len,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;            
        }
        ~ConvolveDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src1, Ipp64f * src2, Ipp64f * dst) {
            IppStatus status = ippsConvolve_64f(src1,src1Len,src2,src2Len,dst,type,buffer);
            checkStatus(status);
        }
    };


    struct ConvolutionFilterDouble
    {
        ConvolveDouble *filter;
        std::shared_ptr<Ipp64f> h;
        size_t            len,block;
        std::vector<Ipp64f> ola,temp;

        ConvolutionFilterDouble(size_t n, Ipp64f * impulse, size_t block_size)
        {
            filter = new ConvolveDouble(block_size,block_size);
            block = block_size;
            Ipp64f *hi = ippsMalloc_64f(len);
            h = std::shared_ptr<Ipp64f>(hi,[](Ipp64f * p) { ippsFree(p); }); 
            len = n;
            assert(filter != NULL);
            assert(hi != nullptr);
            ola.resize(block_size);
            temp.resize(block_size + n -1);
        }
        ~ConvolutionFilterDouble() {

        }
        void ProcessBlock(Ipp64f * signal, Ipp64f * dest) {
            filter->Process(h.get(),signal,dest);
            for(size_t i = 0; i < block; i++) {
                dest[i] = temp[i] + ola[i];
                ola[i]  = temp[i+block];
            }
        }
    };
    
    
    void Goertzal(int len, const Ipp64fc * src, Ipp64fc * dst, Ipp64f freq) {
        IppStatus status = ippsGoertz_64fc(src,len,dst,freq);
    }
    void Goertzal(int len, const Ipp64f * src, Ipp64fc * dst, Ipp64f freq) {
        IppStatus status = ippsGoertz_64f(src,len,dst,freq);
    }

    
    void UpSample(const Ipp16s * pSrc, int srcLen, Ipp16s* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_16s(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    
    void UpSample(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    
    void UpSample(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp16s * pSrc, int srcLen, Ipp16s* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_16s(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    
    
