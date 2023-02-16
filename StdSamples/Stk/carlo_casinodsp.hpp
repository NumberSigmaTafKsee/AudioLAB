// going to specialize everything
#pragma once


namespace Casino
{
    template<typename T1, typename T2>
    void RealToComplex(const T1 * pSrcR, const T1 * pSrcD, T2* pDst, int len) {
        throw std::runtime_error("called the abstract realtocomplex");
    }
    template<>
    void RealToComplex<Ipp32f,Ipp32fc>(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32fc* pDst, int len) {
        IppStatus status = ippsRealToCplx_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    
    template<typename T1, typename T2>
    void ComplexToReal(const T1 * src, T2 * real, T2 * imag, int len) {
        throw std::runtime_error("called the abstract complextoreal");
    }

    template<>
    void ComplexToReal<Ipp32fc,Ipp32f>(const Ipp32fc * src, Ipp32f * real, Ipp32f * imag, int len) {
        IppStatus status = ippsCplxToReal_32fc(src,real,imag,len);
        checkStatus(status);
    }

    template<typename T>
    void Magnitude(const T * pSrcR, const T * pSrcD, T* pDst, int len)
    {
        throw std::runtime_error("called the abstract magnitude");
    }

    template<>
    void Magnitude<Ipp32f>(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32f* pDst, int len) {
        IppStatus status = ippsMagnitude_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }

    template<typename T>
    void Phase(const T * pSrcR, const T * pSrcD, T* pDst, int len) {
        throw std::runtime_error("called the abstract phase");
    }

    template<>
    void Phase<Ipp32f>(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32f* pDst, int len) {
        IppStatus status = ippsPhase_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }

    template<typename T>
    void CartToPolar(const T* pSrcR, const T* pSrcI, T* pDstMag, T* pDstPhase, int len)
    {
        throw std::runtime_error("called the abstract CartToPolar");
    }

    template<>
    void CartToPolar<Ipp32f>(const Ipp32f * pSrcR, const Ipp32f * pSrcI, Ipp32f* pDstMag, Ipp32f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_32f(pSrcR,pSrcI,pDstMag,pDstPhase,len);
        checkStatus(status);
    }

    template<typename T>
    void PolarToCart(const T* pMag, const T* pPhase, T* real, T* imag, int len) {
        throw std::runtime_error("called the abstract PolarToCart");
    }

    template<>
    void PolarToCart<Ipp32f>(const Ipp32f* pMag, const Ipp32f* pPhase, Ipp32f* real, Ipp32f* imag, int len) {
        IppStatus status = ippsPolarToCart_32f(pMag,pPhase,real,imag,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void Magnitude(const T1 * pSrc, T2* pDst, int len) {
        throw std::runtime_error("called the abstract Magnitude");
    }

    template<>
    void Magnitude<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsMagnitude_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    
    template<typename T1, typename T2>
    void Phase(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        throw std::runtime_error("called the abstract Phase");
    }

    template<>
    void Phase<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsPhase_32fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void CartToPolar(const T1 * src, T2* pDstMag, T2* pDstPhase, int len) {
        throw std::runtime_error("called the abstract CartToPolar");
    }

    template<>
    void CartToPolar<Ipp32fc,Ipp32f>(const Ipp32fc * src, Ipp32f* pDstMag, Ipp32f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_32fc(src,pDstMag,pDstPhase,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void PolarToCart(const T1* pMag, const T1* pPhase, T2 * dst, int len) {
        throw std::runtime_error("called the abstract PolarToCar");
    }

    template<>
    void PolarToCart<Ipp32f,Ipp32fc>(const Ipp32f* pMag, const Ipp32f* pPhase, Ipp32fc * dst, int len) {
        IppStatus status = ippsPolarToCart_32fc(pMag,pPhase,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Conj(const T * pSrc, T* pDst, int len)
    {
        throw std::runtime_error("called the abstract Conj");
    }

    template<>
    void Conj<Ipp32fc>(const Ipp32fc * pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsConj_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    void Conj<Ipp64fc>(const Ipp64fc * pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsConj_64fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T>
    void ConjFlip(const T * pSrc, T* pDst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void ConjFlip<Ipp32fc>(const Ipp32fc * pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsConjFlip_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    void ConjFlip<Ipp64fc>(const Ipp64fc * pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsConjFlip_64fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void PowerSpectrum(const T1 * pSrc, T2* pDst, int len) {
        throw std::runtime_error("foobar");
    }

    template<>
    void PowerSpectrum<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsPowerSpectr_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    void PowerSpectrum<Ipp64fc,Ipp64f>(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsPowerSpectr_64fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void Real(const T1 * pSrc, T2* pDst, int len) {
        throw std::runtime_error("foobar");
    }

    template<>
    void Real<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsReal_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    void Real<Ipp64fc,Ipp64f>(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsReal_64fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void Imag(const T1 * pSrc, T2* pDst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Imag<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsImag_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    void Imag<Ipp64fc,Ipp64f>(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsImag_64fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T, typename T2>
    void Threshold(const T * pSrc, T * pDst, size_t len, T2 level, IppCmpOp op = ippCmpGreater) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Threshold<Ipp32f,Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, size_t len, Ipp32f level, IppCmpOp op) {
        IppStatus status = ippsThreshold_32f(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    template<>
    void Threshold<Ipp32fc,Ipp32f>(const Ipp32fc * pSrc, Ipp32fc * pDst, size_t len, Ipp32f level, IppCmpOp op) {
        IppStatus status = ippsThreshold_32fc(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    template<>
    void Threshold<Ipp64f,Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, size_t len, Ipp64f level, IppCmpOp op) {
        IppStatus status = ippsThreshold_64f(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    template<>
    void Threshold<Ipp64fc,Ipp64f>(const Ipp64fc * pSrc, Ipp64fc * pDst, size_t len, Ipp64f level, IppCmpOp op) {
        IppStatus status = ippsThreshold_64fc(pSrc,pDst,level,len, op);
        checkStatus(status);
    }

    template<typename T>
    void WinBartlett(const T* src, T * dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void WinBartlett<Ipp32f>(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp32fc>(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32fc(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp64f>(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp64fc>(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64fc(src,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void WinBartlett(T* src, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void WinBartlett<Ipp32f>(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp32fc>(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32fc_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp64f>(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinBartlett<Ipp64fc>(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64fc_I(dst,len);
        checkStatus(status);
    }

    template<typename T>
    void WinBlackman(const T* src, T* dst, int len, Ipp32f alpha) {
        throw std::runtime_error("foobar");
    }
    template<>
    void WinBlackman<Ipp32f>(const Ipp32f* src, Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32f(src,dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinBlackman<Ipp32fc>(const Ipp32fc* src, Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32fc(src,dst,len,alpha);
        checkStatus(status);
    }

    template<typename T1,typename T2>
    void WinBlackman(T1 * dst, int len, T2 alpha)
    {
        throw std::runtime_error("foobar");
    }

    template<>
    void WinBlackman<Ipp32f,Ipp32f>(Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32f_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinBlackman<Ipp32fc,Ipp32f>(Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32fc_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinBlackman<Ipp64f,Ipp64f>(Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64f_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinBlackman<Ipp64fc,Ipp64f>(Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64fc_I(dst,len,alpha);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void WinKaiser(const T1* src, T1 * dst, int len, T2 alpha) {
           throw std::runtime_error("foobar");    
    }

    template<>
    void WinKaiser<Ipp32f,Ipp32f>(const Ipp32f* src, Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32f(src,dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp32fc,Ipp32f>(const Ipp32fc* src, Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32fc(src,dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp64f,Ipp64f>(const Ipp64f* src, Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64f(src,dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp64fc,Ipp64f>(const Ipp64fc* src, Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64fc(src,dst,len,alpha);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void WinKaiser(T1 * dst, int len, T2 alpha) {
        throw std::runtime_error("foobar");    
    }

    template<>
    void WinKaiser<Ipp32f,Ipp32f>(Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32f_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp32fc,Ipp32f>(Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32fc_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp64f,Ipp64f>(Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64f_I(dst,len,alpha);
        checkStatus(status);
    }
    template<>
    void WinKaiser<Ipp64fc,Ipp64f>(Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64fc_I(dst,len,alpha);
        checkStatus(status);
    }

    template<typename T>
    void WinHamming(const T* src, T* dst, int len) {
        throw std::runtime_error("foobar");    
    }

    template<>
    void WinHamming<Ipp32f>(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHamming_32f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp32fc>(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_32fc(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp64f>(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHamming_64f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp64fc>(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_64fc(src,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void WinHamming(T * dst, int len)
    {
        throw std::runtime_error("foobar");    
    }
    template<>
    void WinHamming<Ipp32f>(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHamming_32f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp32fc>(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_32fc_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp64f>(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHamming_64f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHamming<Ipp64fc>(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_64fc_I(dst,len);
        checkStatus(status);
    }

    template<typename T>
    void WinHann(const T* src, T* dst, int len)
    {
        throw std::runtime_error("foobar");    
    }

    template<>
    void WinHann<Ipp32f>(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHann_32f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp32fc>(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHann_32fc(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp64f>(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHann_64f(src,dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp64fc>(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHann_64fc(src,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void WinHann(T * dst, int len)
    {
        throw std::runtime_error("foobar");    
    }

    template<>
    void WinHann<Ipp32f>(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHann_32f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp32fc>(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHann_32fc_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp64f>(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHann_64f_I(dst,len);
        checkStatus(status);
    }
    template<>
    void WinHann<Ipp64fc>(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHann_64fc_I(dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Sum(const T * pSrc, int len, T* pDst, int hint = 0) {
        throw std::runtime_error("Called abstract Sum");
    }

    template<>
    void Sum<Ipp32f>(const Ipp32f * pSrc, int len, Ipp32f* pDst, int hint) {
        IppStatus status = ippsSum_32f(pSrc,len,pDst,(IppHintAlgorithm)hint);
        checkStatus(status);
    }
    template<>
    void Sum<Ipp64f>(const Ipp64f * pSrc, int len, Ipp64f* pDst, int hint) {
        IppStatus status = ippsSum_64f(pSrc,len,pDst);
        checkStatus(status);
    }

    template<typename T>
    void AddC(const T* src, T val, T* dst, int len)
    {
        throw std::runtime_error("foobar");    
    }

    template<>
    void AddC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsAddC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsAddC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsAddC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsAddC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void AddC(T val, T* dst, int len)
    {
        throw std::runtime_error("foobar");    
    }
    template<>
    void AddC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsAddC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsAddC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsAddC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsAddC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Add(const T * a, T *b, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Add<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsAdd_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsAdd_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsAdd_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsAdd_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Add(const T * pSrc, T * pDst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Add<Ipp32f>(const Ipp32f * src, Ipp32f* dst, int len) {
        IppStatus status = ippsAdd_32f_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp32fc>(const Ipp32fc * src, Ipp32fc* dst, int len) {
        IppStatus status = ippsAdd_32fc_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64f>(const Ipp64f * src, Ipp64f* dst, int len) {
        IppStatus status = ippsAdd_64f_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64fc>(const Ipp64fc * src, Ipp64fc* dst, int len) {
        IppStatus status = ippsAdd_64fc_I(src,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubC(const T* src, T val, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void SubC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubC(T val, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void SubC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Sub(const T * a, T *b, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void Sub<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsSub_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsSub_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsSub_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsSub_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Sub(const T * a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sub<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsSub_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsSub_32fc_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsSub_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsSub_64fc_I(a,dst,len);
        checkStatus(status);
    }
    
    template<typename T>
    void SubCRev(const T * src, T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void SubCRev<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubCRev_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubCRev_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubCRev_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubRev(const T * a, T b, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void SubRev<Ipp32f>(const Ipp32f * a, Ipp32f b, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp32fc>(const Ipp32fc * a, Ipp32fc b, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubCRev_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp64f>(const Ipp64f * a, Ipp64f b, Ipp64f* dst, int len) {
        IppStatus status = ippsSubCRev_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp64fc>(const Ipp64fc * a, Ipp64fc b, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubCRev_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void MulC(const T* src, T val, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void MulC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsMulC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsMulC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsMulC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsMulC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void MulC(T val, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void MulC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsMulC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsMulC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsMulC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsMulC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Mul(const T* a, T*b, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Mul<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsMul_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsMul_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsMul_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsMul_64fc(a,b,dst,len);
        checkStatus(status);
    }
    
    template<typename T>
    void Mul(const T* a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Mul<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsMul_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsMul_32fc_I(a,dst,len);
        checkStatus(status);
    }    
    template<>
    void Mul<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsMul_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsMul_64fc_I(a,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivC(const T* src, T val, T* dst, int len)
    {
        assert(1==0);
    }

    template<>
    void DivC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsDivC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsDivC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsDivC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivC(T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsDivC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsDivC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsDivC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Div(const T* a, T *b, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Div<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsDiv_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsDiv_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsDiv_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsDiv_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Div(const T* a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Div<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsDiv_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsDiv_32fc_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsDiv_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsDiv_64fc_I(a,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivCRev(const T* src, T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivCRev<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivCRev_32f(src,val,dst,len);
        checkStatus(status);
    }    
     
    template<typename T>
    void DivCRev(T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivCRev<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivCRev_32f_I(val,dst,len);
        checkStatus(status);
    }    
    

    
    void AddProductC(const Ipp32f * array, int len, const Ipp32f val, Ipp32f * dst) {            
        IppStatus status = ippsAddProductC_32f(array,val,dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void AddProduct(const T* array1, int len, const T* array2, T* dst) {
        assert(1==0);
    }
    template<>
    void AddProduct<Ipp32f>(const Ipp32f * array1, int len, const Ipp32f *array2, Ipp32f * dst) {            
        IppStatus status = ippsAddProduct_32f(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp64f>(const Ipp64f * array1, int len, const Ipp64f *array2, Ipp64f * dst) {            
        IppStatus status = ippsAddProduct_64f(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp32fc>(const Ipp32fc * array1, int len, const Ipp32fc *array2, Ipp32fc * dst) {            
        IppStatus status = ippsAddProduct_32fc(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp64fc>(const Ipp64fc * array1, int len, const Ipp64fc *array2, Ipp64fc * dst) {            
        IppStatus status = ippsAddProduct_64fc(array1,array2,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Abs(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Abs<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsAbs_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Abs<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsAbs_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Abs(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Abs<Ipp32f>( Ipp32f * dst, int len) {
        IppStatus status = ippsAbs_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Abs<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsAbs_64f_I(dst,len);
        checkStatus(status);        
    }


    template<typename T>
    void Sqr(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqr<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqr_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp32fc>(const Ipp32fc * src, Ipp32fc * dst, int len) {
        IppStatus status = ippsSqr_32fc(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsSqr_64f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64fc>(const Ipp64fc * src, Ipp64fc * dst, int len) {
        IppStatus status = ippsSqr_64fc(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqr(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqr<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsSqr_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp32fc>(Ipp32fc * dst, int len) {
        IppStatus status = ippsSqr_32fc_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsSqr_64f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64fc>(Ipp64fc * dst, int len) {
        IppStatus status = ippsSqr_64fc_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqrt(const T * src, T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqrt<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqrt_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp32fc>(const Ipp32fc * src, Ipp32fc * dst, int len) {
        IppStatus status = ippsSqrt_32fc(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsSqrt_64f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64fc>(const Ipp64fc * src, Ipp64fc * dst, int len) {
        IppStatus status = ippsSqrt_64fc(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqrt(T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqrt<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsSqrt_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp32fc>(Ipp32fc * dst, int len) {
        IppStatus status = ippsSqrt_32fc_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsSqrt_64f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64fc>(Ipp64fc * dst, int len) {
        IppStatus status = ippsSqrt_64fc_I(dst,len);
        checkStatus(status);        
    }

    void Cubrt(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsCubrt_32f(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Exp(const T * src, T * dst, int len) {
        assert(1==0);
    }

    template<>
    void Exp<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsExp_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Exp<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsExp_64f(src,dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void Exp( T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Exp<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsExp_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Exp<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsExp_64f_I(dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void Ln(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Ln<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Ln<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Ln(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Ln<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Ln<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void SumLn(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void SumLn<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void SumLn<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Arctan(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Arctan<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsArctan_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Arctan<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsArctan_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Arctan(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Arctan<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsArctan_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Arctan<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsArctan_64f_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    T Max(const T* src,int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Max<Ipp32f>(const Ipp32f * src,int len) {
        Ipp32f max=0;
        IppStatus status = ippsMax_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    template<>
    Ipp64f Max<Ipp64f>(const Ipp64f * src,int len) {
        Ipp64f max=0;
        IppStatus status = ippsMax_64f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    
    template<typename T>
    T Min(const T* src,int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Min<Ipp32f>(const Ipp32f * src,int len) {
        Ipp32f min=0;
        IppStatus status = ippsMin_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    template<>
    Ipp64f Min<Ipp64f>(const Ipp64f * src,int len) {
        Ipp64f min=0;
        IppStatus status = ippsMin_64f(src,len,&min);
        checkStatus(status);        
        return min;
    }

    template<typename T>
    int MaxIndex(const T* src, T * max, int len) {
        assert(1==0);
    }
    template<>
    int MaxIndex<Ipp32f>(const Ipp32f * src,Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxIndx_32f(src,len,max,&index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MaxIndex<Ipp64f>(const Ipp64f * src,Ipp64f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxIndx_64f(src,len,max,&index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    int MinIndex(const T * src, int len, T * min) {
        assert(1==0);
    }

    template<>
    int MinIndex<Ipp32f>(const Ipp32f * src, int len, Ipp32f * min) {
        int index=-1;
        IppStatus status = ippsMinIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MinIndex<Ipp64f>(const Ipp64f * src, int len, Ipp64f * min) {
        int index=-1;
        IppStatus status = ippsMinIndx_64f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    T MaxAbs(const T* src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f MaxAbs<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f max=0;
        IppStatus status = ippsMaxAbs_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    template<>
    Ipp64f MaxAbs<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f max=0;
        IppStatus status = ippsMaxAbs_64f(src,len,&max);
        checkStatus(status);        
        return max;
    }

    template<typename T>
    T MinAbs(const T * src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f MinAbs<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f min=0;
        IppStatus status = ippsMinAbs_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    template<>
    Ipp64f MinAbs<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f min=0;
        IppStatus status = ippsMinAbs_64f(src,len,&min);
        checkStatus(status);        
        return min;
    }

    template<typename T>
    int MaxAbsIndex(const T* src, T* max, int len) {
        assert(1==0);
    }

    template<>
    int MaxAbsIndex<Ipp32f>(const Ipp32f * src, Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxAbsIndx_32f(src,len,max, &index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MaxAbsIndex<Ipp64f>(const Ipp64f * src, Ipp64f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxAbsIndx_64f(src,len,max, &index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    int MinAbsIndex(const T* src, T* min, int len) {
        assert(1==0);
    }
    template<>
    int MinAbsIndex<Ipp32f>(const Ipp32f * src, Ipp32f * min, int len) {
        int index=-1;
        IppStatus status = ippsMinAbsIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    
    template<typename T>
    void MinMax(const T* src, int len, T* min, T* max) {    
        assert(1==0);
    }
    template<>
    void MinMax<Ipp32f>(const Ipp32f * src, int len, Ipp32f * min, Ipp32f * max) {    
        IppStatus status = ippsMinMax_32f(src,len,min,max);
        checkStatus(status);            
    }

    template<typename T>
    void MinMaxIndex(const T* src, int len, T *min, int * min_index, T *max, int * max_index) {
        assert(1==0);
    }
        
    template<>
    void MinMaxIndex<Ipp32f>(const Ipp32f * src, int len, Ipp32f *min, int * min_index, Ipp32f *max, int * max_index) {
        IppStatus status = ippsMinMaxIndx_32f(src,len,min,min_index,max,max_index);
        checkStatus(status);            
    }


    void ReplaceNAN(Ipp32f * p, Ipp32f v, int len) {
        IppStatus status = ippsReplaceNAN_32f_I(p,len,v);
        checkStatus(status);            
    }
    
    template<typename T>
    T Mean(const T * src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Mean<Ipp32f>(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsMean_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    template<>
    Ipp64f Mean<Ipp64f>(const Ipp64f * src, int len)
    {
        Ipp64f mean;
        IppStatus status = ippsMean_64f(src,len,&mean);
        checkStatus(status);            
        return mean;
    }

    template<typename T>
    T StdDev(const T* src, int len)
    {
        assert(1==0);
    }

    template<>
    Ipp32f StdDev<Ipp32f>(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsStdDev_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    template<>
    Ipp64f StdDev<Ipp64f>(const Ipp64f * src, int len)
    {
        Ipp64f mean;
        IppStatus status = ippsStdDev_64f(src,len,&mean);
        checkStatus(status);            
        return mean;
    }

    template<typename T>
    void MeanStdDev(const T* src, int len, T* mean, T* dev)
    {
        assert(1==0);
    }
    template<>
    void MeanStdDev<Ipp32f>(const Ipp32f * src, int len, Ipp32f * mean, Ipp32f * dev)
    {        
        IppStatus status = ippsMeanStdDev_32f(src,len,mean,dev,ippAlgHintNone);
        checkStatus(status);                    
    }
    template<>
    void MeanStdDev<Ipp64f>(const Ipp64f * src, int len, Ipp64f * mean, Ipp64f * dev)
    {        
        IppStatus status = ippsMeanStdDev_64f(src,len,mean,dev);
        checkStatus(status);                    
    }

    template<typename T>
    T NormInf(const T* src, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormInf<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_Inf_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormInf<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_Inf_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormL1(const T* src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f NormL1<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L1_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormL1<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_L1_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormL2(const T* src, int len) {
        assert(1==0);
    }

    
    template<>
    Ipp32f NormL2<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L2_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormL2<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_L2_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffInf(const T* src1,const T* src2, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f NormDiffInf<Ipp32f>(const Ipp32f * src1,const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_Inf_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffInf<Ipp64f>(const Ipp64f * src1,const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_Inf_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffL1(const T* src1, const T* src2, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormDiffL1<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L1_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffL1<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_L1_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffL2(const T* src1, const T* src2, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormDiffL2<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L2_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffL2<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_L2_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T DotProduct(const T * src1, const T * src2, int len) 
    {
        assert(1==0);
    }

    template<>
    Ipp32f DotProduct<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f dp;
        IppStatus status = ippsDotProd_32f(src1,src2,len,&dp);
        checkStatus(status);     
        return dp;
    }
    template<>
    Ipp64f DotProduct<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f dp;
        IppStatus status = ippsDotProd_64f(src1,src2,len,&dp);
        checkStatus(status);     
        return dp;
    }

    void MaxEvery(const Ipp32f * src1, const Ipp32f * src2, Ipp32f * dst, int len) {        
        IppStatus status = ippsMaxEvery_32f(src1,src2,dst,len);
        checkStatus(status);             
    }
    void MinEvery(const Ipp32f * src1, const Ipp32f * src2, Ipp32f * dst, int len) {        
        IppStatus status = ippsMinEvery_32f(src1,src2,dst,len);
        checkStatus(status);             
    }
    Ipp32f ZeroCrossing(const Ipp32f * src1, int len, IppsZCType zcType = ippZCR) {
        Ipp32f zc;
        IppStatus status = ippsZeroCrossing_32f(src1,len,&zc,zcType);
        checkStatus(status);     
        return zc;
    }
    Ipp32f MSE(const Ipp32f * thetaAbs, const Ipp32f * thetaModel, int len) {
        Ipp32f r = 0;
        for(size_t i = 0; i <  len; i++)
            r += pow(thetaAbs[i] - thetaModel[i],2.0);
        return (1.0f / (Ipp32f)len) * r;
    }
    Ipp32f RMSE(const Ipp32f * thetaAbs, const Ipp32f * thetaModel, int len) {
        Ipp32f r = MSE(thetaAbs,thetaModel,len);
        return sqrt(r);
    }
    Ipp32f MeanSquare(const Ipp32f * x, int len) {
        Ipp32f r = 0;
        for(size_t i = 0; i < len; i++) r += pow(x[i],2.0);
        return (1.0f / (Ipp32f)len) * r;
    }
    Ipp32f AmpToDb(Ipp32f a) {
        return pow(10.0,a/20.0);
    }
    Ipp32f DbToAmp(Ipp32f a) {
        return 20.0*log10(a);
    }
    void Tone(Ipp32f * array, int len, Ipp32f mag, Ipp32f freq, Ipp32f * phase) {
            IppStatus status = ippsTone_32f(array,len,mag,freq,phase,ippAlgHintNone);
            checkStatus(status);
    }
    void Triangle(Ipp32f * array, int len, Ipp32f m, Ipp32f f, Ipp32f a, Ipp32f * p) {
        IppStatus status = ippsTriangle_32f(array,len,m,f,a,p);
        checkStatus(status);
    }
    void Sort(Ipp32f * array, int len, int dir=1) {
        IppStatus status;
        if(dir >= 0)
            status = ippsSortAscend_32f_I(array,len);
        else
            status = ippsSortDescend_32f_I(array,len);
        checkStatus(status);
    }
    
    template<typename T>
    void UpSample(const T * pSrc, int srcLen, T* pDst, int *dstLen, int fact, int * pPhase)
    {
        throw std::runtime_error("");
    }

    template<>
    void UpSample<Ipp32f>(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp32fc>(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp64f>(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp64fc>(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }

    template<typename T>
    void DownSample(const T * pSrc, int srcLen, T* pDst, int *dstLen, int fact, int * pPhase) {
        throw std::runtime_error("");
    }

    template<>
    void DownSample<Ipp32f>(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp32fc>(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp64f>(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp64fc>(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }

    template<typename T>
    T kaiserBeta(T As)
    {
        if(As > (T)50.0)
            return (T)0.1102 * (As - (T)8.7);
        else if(As >= (T)21.)
            return (T)0.5842 * std::pow(As - (T)21.0, (T)0.4) + (T)0.07886 * (As - (T)21.0);
        else
            return (T)0.;
    }

    
    void Goertzal(int len, const Ipp32fc * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32fc(src,len,dst,freq);
    }
    
    void Goertzal(int len, const Ipp32f * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32f(src,len,dst,freq);
    }

    void DoubleToFloat(size_t len, double * src, Ipp32f * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (Ipp32f)src[i];
    }
    void FloatToDouble(size_t len, Ipp32f * src, double * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (double)src[i];
    }

    template<typename T>
    IppDataType GetDataType()
    {
        IppDataType dataType;
        if(typeid(T) == typeid(Ipp32f)) dataType = ipp32f;
        else if(typeid(T) == typeid(Ipp64f)) dataType = ipp64f;
        else if(typeid(T) == typeid(Ipp32fc)) dataType = ipp32fc;
        else if(typeid(T) == typeid(Ipp64fc)) dataType = ipp64fc;
        else throw std::runtime_error("Type not supported yet");
        return dataType;
    }

    template<typename T>
    void RandUniformInit(void * state, T low, T high, int seed)
    {
        assert(1==0);
    }
    template<>
    void RandUniformInit<Ipp32f>(void * state, Ipp32f low, Ipp32f high, int seed)
    {
        IppStatus status = ippsRandUniformInit_32f((IppsRandUniState_32f*)state,low,high,seed);
        checkStatus(status);
    }
    template<typename T>
    void RandUniform(T * array, int len, void * state)
    {
        assert(1==0);
    }
    template<>
    void RandUniform<Ipp32f>(Ipp32f * array, int len, void * state)
    {
        IppStatus status = ippsRandUniform_32f(array,len,(IppsRandUniState_32f*)state);
        checkStatus(status);
    }
    template<typename T>
    void RandUniformGetSize(int * p)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status;
        if(dType ==ipp32f) status = ippsRandUniformGetSize_32f(p);
        else if(dType == ipp64f) status = ippsRandUniformGetSize_64f(p);        
        checkStatus(status);
    }
    template<typename T>
    struct RandomUniform
    {
        Ipp8u * state;
        RandomUniform(T high, T low, int seed=-1)
        {
            if(seed == -1) seed = time(NULL);            
            int size=0;
            RandUniformGetSize<T>(&size);
            if(size==0) throw std::runtime_error("RandomUniform size is 0");
            state = Malloc<Ipp8u>(size);
            RandUniformInit<T>(state,high,low,seed);    
        }
        ~RandomUniform() {
            if(state) Free(state);
        }
        void fill(T* array, int len) { 
            RandUniform<T>(array,len,state);
        }
        
    };
    

    template<typename T>
    struct IPPArray
    {
        std::shared_ptr<T> ptr;
        T * array;
        size_t   len;
        int      r,w;

        IPPArray(const IPPArray<T> & a) {
            *this = a;
        }
        virtual ~IPPArray() = default;

        IPPArray(size_t n) {
            array = Malloc<T>(len = n);
            ptr = std::shared_ptr<T>(array,[](T* p) { Free<T>(p); });
            assert(array != NULL);
            Zero<T>(array,n);
            r = 0;
            w = 0;
        }
                
        
        void resize(size_t n) {
            T * p = Malloc<T>(n);
            Move<T>(array,p,n);
            Free<T>(array);
            array  = p;
            len    = n;
        }
        void fill(T value) {
            if(array == NULL) return;
            Set<T>(value,array,len);
        }
        T sum() {
            T r = 0;
            Sum<T>(array,len,&r);
            return r;
        }
        
        T& operator[] (size_t i) { return array[i]; }

        T      __getitem__(size_t i) { return array[i]; }
        void   __setitem__(size_t i, T v) { array[i] = v; }

        IPPArray<T>& operator = (const IPPArray & x) {
            ptr.reset();
            ptr = x.ptr;
            array = x.array;
            len = x.len;
            return *this;
        }
        
        void ring_push(const T& value) {
            array[w++] = value;
            w = w % len;
        }
        T ring_pop() {
            T v = array[r++];
            r = r % len;
            return v;
        }
        T ring_linear_pop() {
            T v1 = array[r];
            r = (r+1) % len;
            T v2 = array[r];            
            T frac = v1 - std::floor(v1);
            return v1 + frac*(v2-v1);
        }

        IPPArray<T> operator + (const T& value) {
            IPPArray<T> r(*this);
            AddC<T>(array,value,r.array,len);
            return r;
        }
        IPPArray<T> operator + (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Add<T>(array,b.array,r.array,len);
            return r;
        }
        IPPArray<T> operator - (const T& value) {
            IPPArray<T> r(*this);
            SubC<T>(array,value,r.array,len);
            return r;
        }
        IPPArray<T> operator - (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Sub<T>(array,b.array,r.array,len);
            return r;
        }        
        IPPArray<T> operator * (const T& value) {
            IPPArray<T> r(*this);
            MulC<T>(array,value,r.array,len);
            return r;
        }        
        IPPArray<T> operator * (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Mul<T>(array,b.array,r.array,len);
            return r;
        }
        IPPArray<T> operator / (const T& value) {
            IPPArray<T> r(*this);
            DivC<T>(array,value,r.array,len);
            return r;
        }        
        IPPArray<T> operator / (const IPPArray<T>& b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Div<T>(array,b.array,r.array,len);
            return r;
        }        

        void print() {
            std::cout << "Array[" << len << "]=";
            for(size_t i = 0; i < len-1; i++) std::cout << array[i] << ",";
            std::cout << array[len-1] << std::endl;
        }
        IPPArray<T>& copy(const IPPArray<T> & a) {
            ptr.reset();
            array = Malloc<T>(a.len);
            memcpy(array,a.array,a.len*sizeof(T));
            return *this; 
       }        
    };
    /*
    sample_vector<Ipp32f> ConvertToSampleVector(const Ipp32fArray &a) {
        sample_vector<Ipp32f> r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(Ipp32f));
        return r;
    };
    */
    /*
    Vector::Vector4f ConvertToVector4f(const Ipp32fArray &a) {
        Vector::Vector4f r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(Ipp32f));
        return r;
    }
    Vector::Vector8f ConvertToVector8f(const Ipp32fArray &a) {
        Vector::Vector8f r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(Ipp32f));
        return r;
    }
    sample_vector<double> ConvertToSampleVectorDouble(const Ipp32fArray &a) {
        sample_vector<double> r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    Vector::Vector4f ConvertToVector2d(const Ipp32fArray &a) {
        Vector::Vector4f r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    Vector::Vector8f ConvertToVector4d(const Ipp32fArray &a) {
        Vector::Vector8f r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    */

    
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
    struct AutoCorr
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorr(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            AutoCorrNormGetBufferSize<T>(srcLen,dstLen,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = Malloc<Ipp8u>(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorr() {
            if(buffer) Free(buffer);
        }
        void Process(T * src, T * dst) {
            AutoCorrNorm<T>(src,srcLen,dst,dstLen,type,buffer);
        }
    };

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
    void CrossCorrNorm<Ipp32f>(const Ipp32f * pSrc1, int src1Len, const Ipp32f * pSrc2, int src2Len, Ipp32f * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_32f(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<Ipp64f>(const Ipp64f * pSrc1, int src1Len, const Ipp64f * pSrc2, int src2Len, Ipp64f * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_64f(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<Ipp32fc>(const Ipp32fc * pSrc1, int src1Len, const Ipp32fc * pSrc2, int src2Len, Ipp32fc * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_32fc(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    template<>
    void CrossCorrNorm<Ipp64fc>(const Ipp64fc * pSrc1, int src1Len, const Ipp64fc * pSrc2, int src2Len, Ipp64fc * pDst, int dstLen, int lowLag, IppEnum algType, Ipp8u* pBuffer)
    {   
        IppStatus status;
        status = ippsCrossCorrNorm_64fc(pSrc1,src1Len,pSrc2,src2Len,pDst,dstLen,lowLag,algType,pBuffer);
        checkStatus(status);
    }
    
    template<typename T>
    struct CrossCorr
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorr(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            CrossCorrNormGetBufferSize<T>(src1Len,src2Len,dstLen,lowLag,algorithm,&bufferLen);        
            if(bufferLen > 0) buffer = Malloc<Ipp8u>(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->lowLag = lowLag;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
        }
        ~CrossCorr() {
            if(buffer) Free(buffer);
        }
        void Process(T * src1, T * src2, T * dst) {
            CrossCorrNorm<T>(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);            
        }
    };

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

    template<typename T1, typename T2>
    void ConjPack(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjPack<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPack_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjPack<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPack_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void ConjPerm(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjPerm<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPerm_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjPerm<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPerm_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void ConjCcs(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjCcs<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjCcs_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjCcs<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjCcs_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }


    template<typename T>
    void MulPack(const T* pSrc1, const T * pSrc2, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPack<Ipp32f>(const Ipp32f* pSrc1, const Ipp32f * pSrc2, Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPack_32f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPack<Ipp64f>(const Ipp64f* pSrc1, const Ipp64f * pSrc2, Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPack_64f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T>
    void MulPerm(const T* pSrc1, const T* pSrc2, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPerm<Ipp32f>(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPerm_32f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPerm<Ipp64f>(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPerm_64f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T>
    void MulPackConj(const T* pSrc1, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPackConj<Ipp32f>(const Ipp32f* pSrc1,Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPackConj_32f_I(pSrc1,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPackConj<Ipp64f>(const Ipp64f* pSrc1,Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPackConj_64f_I(pSrc1,pDst,dstLen);
        checkStatus(status);
    }
    
    template<typename T>
    void FFTInitR(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>
    void FFTInitR<Ipp32f>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_R_32f((IppsFFTSpec_R_32f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitR<Ipp64f>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
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
    void FFTInitC<Ipp32f>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_32f((IppsFFTSpec_C_32f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitC<Ipp64f>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_64f((IppsFFTSpec_C_64f**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    
    template<>
    void FFTInitC<Ipp32fc>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsFFTInit_C_32fc((IppsFFTSpec_C_32fc**)ppFFTSpec,order,flag,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInitC<Ipp64fc>(void ** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u * pSpec, Ipp8u* pSpecBuffer)
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
    void FFTFwd_RToPack<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPack_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToPack<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTFwd_RToPerm<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToPerm_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToPerm<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTFwd_RToCCS<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTFwd_RToCCS_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTFwd_RToCCS<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTInv_PackToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PackToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_PackToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTInv_PermToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_PermToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_PermToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTInv_CCSToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFFTInv_CCSToR_32f(pSrc,pDst,(IppsFFTSpec_R_32f*)pFFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void FFTInv_CCSToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pFFTSpec, Ipp8u * pBuffer)
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
    void FFTGetSizeR<Ipp32f>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_R_32f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeR<Ipp64f>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
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
    void FFTGetSizeC<Ipp32f>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_32f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<Ipp64f>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_64f(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<Ipp32fc>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsFFTGetSize_C_32fc(order,flag,hint,pSpecSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void FFTGetSizeC<Ipp64fc>(int order, int flag, IppHintAlgorithm hint, int * pSpecSize, int * pSpecBufferSize, int * pBufferSize)
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
    void FFTFwd_C2C<Ipp32f>(const Ipp32f * pSrcRe, const Ipp32f * pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_32f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTFwd_C2C<Ipp64f>(const Ipp64f * pSrcRe, const Ipp64f * pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
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
    void FFTFwd_C2C<Ipp32fc>(const Ipp32fc * pSrc,Ipp32fc* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_32fc(pSrc,pDst,(IppsFFTSpec_C_32fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTFwd_C2C<Ipp64fc>(const Ipp64fc * pSrc,Ipp64fc* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTFwd_CToC_64fc(pSrc,pDst,(IppsFFTSpec_C_64fc*)pFFTSpec,pBuffer);        
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
    void FFTInv_C2C<Ipp32f>(const Ipp32f * pSrcRe, const Ipp32f * pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsFFTSpec_C_32f*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<>
    void FFTInv_C2C<Ipp64f>(const Ipp64f * pSrcRe, const Ipp64f * pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, void * pFFTSpec, Ipp8u* pBuffer)
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
    void FFTInv_C2C<Ipp32fc>(const Ipp32fc * pSrc,Ipp32fc* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_32fc(pSrc,pDst,(IppsFFTSpec_C_32fc*)pFFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void FFTInv_C2C<Ipp64fc>(const Ipp64fc * pSrc,Ipp64fc* pDst, void * pFFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsFFTInv_CToC_64fc(pSrc,pDst,(IppsFFTSpec_C_64fc*)pFFTSpec,pBuffer);        
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


    template<typename T>
    struct RFFT
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        void  * fft;
                        
        RFFT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            FFTGetSizeR<T>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = std::log2(n);
            FFTInitR<T>(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);                        
        }
        ~RFFT() {
            if(pSpec) Free(pSpec);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);            
        }
        
        void Forward(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            FFTFwd_RToPack<T>(pSrc,pDst,fft,pBuffer);            
        }
        void Inverse(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            FFTInv_PackToR<T>(pSrc,pDst,fft,pBuffer);            
        }
        void Unpack(int len, Ipp32f * pSrc, Ipp32fc * pDst) {
            ConjPack<Ipp32f,Ipp32fc>(pSrc,pDst,len);
        }
    };

    template<typename T>
    void DFTInitR(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        assert(1==0);
    }
    template<>
    void DFTInitR<Ipp32f>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_R_32f(length,flag,hint,(IppsDFTSpec_R_32f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitR<Ipp64f>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
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
    void DFTInitC<Ipp32f>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_32f(length,flag,hint,(IppsDFTSpec_C_32f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<Ipp64f>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_64f(length,flag,hint,(IppsDFTSpec_C_64f*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<Ipp32fc>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
    {
        IppStatus status = ippsDFTInit_C_32fc(length,flag,hint,(IppsDFTSpec_C_32fc*)pDFTSpec,pMemInit);        
        checkStatus(status);            
    }
    template<>
    void DFTInitC<Ipp64fc>(int length, int flag, IppHintAlgorithm hint, void * pDFTSpec, Ipp8u * pMemInit)
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
    void DFTGetSizeR<Ipp32f>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_R_32f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeR<Ipp64f>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
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
    void DFTGetSizeC<Ipp32f>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_32f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<Ipp64f>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_64f(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<Ipp32fc>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
    {
        IppStatus status = ippsDFTGetSize_C_32fc(length,flag,hint,pSizeSpec,pSizeInt,pSizeBuf);
        checkStatus(status);            
    }
    template<>
    void DFTGetSizeC<Ipp64fc>(int length, int flag, IppHintAlgorithm hint, int * pSizeSpec, int * pSizeInt, int * pSizeBuf)
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
    void DFTFwd_C2C<Ipp32f>(const Ipp32f * pSrcRe, const Ipp32f * pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_32f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTFwd_C2C<Ipp64f>(const Ipp64f * pSrcRe, const Ipp64f * pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
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
    void DFTFwd_C2C<Ipp32fc>(const Ipp32fc * pSrc,Ipp32fc* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_32fc(pSrc,pDst,(IppsDFTSpec_C_32fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTFwd_C2C<Ipp64fc>(const Ipp64fc * pSrc,Ipp64fc* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTFwd_CToC_64fc(pSrc,pDst,(IppsDFTSpec_C_64fc*)pDFTSpec,pBuffer);        
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
    void DFTInv_C2C<Ipp32f>(const Ipp32f * pSrcRe, const Ipp32f * pSrcIm, Ipp32f* pDstRe, Ipp32f* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,(IppsDFTSpec_C_32f*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<>
    void DFTInv_C2C<Ipp64f>(const Ipp64f * pSrcRe, const Ipp64f * pSrcIm, Ipp64f* pDstRe, Ipp64f* pDstIm, void * pDFTSpec, Ipp8u* pBuffer)
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
    void DFTInv_C2C<Ipp32fc>(const Ipp32fc * pSrc,Ipp32fc* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_32fc(pSrc,pDst,(IppsDFTSpec_C_32fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }
    template<>
    void DFTInv_C2C<Ipp64fc>(const Ipp64fc * pSrc,Ipp64fc* pDst, void * pDFTSpec, Ipp8u* pBuffer)
    {
        IppStatus status = ippsDFTInv_CToC_64fc(pSrc,pDst,(IppsDFTSpec_C_64fc*)pDFTSpec,pBuffer);        
        checkStatus(status);
    }

    template<typename T>
    void DFTFwd_RToPack(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTFwd_RToPack<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPack_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTFwd_RToPack<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
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
    void DFTFwd_RToPerm<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPerm_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTFwd_RToPerm<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTFwd_RToPerm_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    
    template<typename T>
    void DFTInv_PackToR(const T * pSrc, T * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DFTInv_PackToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PackToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_PackToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
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
    void DFTInv_PermToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_PermToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_PermToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
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
    void DFTInv_CCSToR<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_CCSToR_32f(pSrc,pDst,(IppsDFTSpec_R_32f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DFTInv_CCSToR<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, void * pDFTSpec, Ipp8u * pBuffer)
    {
        IppStatus status = ippsDFTInv_CCSToR_64f(pSrc,pDst,(IppsDFTSpec_R_64f*)pDFTSpec,pBuffer);
        checkStatus(status);
    }

    
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
        
        void Forward(const T* pSrc, T* pDst)
        {                   
            DFTFwd_C2C<T>(pSrc,pDst,fft,pBuffer);
        }
        void Inverse(const T* pSrc, T* pDst)
        {                   
            DFTInv_C2C<T>(pSrc,pDst,fft,pBuffer);
        }
    };

    
    template<typename T>
    struct RDFT
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;
        void  * fft;
                        
        RDFT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            DFTGetSizeR<T>(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = Malloc<Ipp8u>(spec);
            pSpecBuffer = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer = Malloc<Ipp8u>(size);
            int order = n;
            DFTInitR<T>(order,flag,ippAlgHintNone,fft,pSpecBuffer);                        
        }
        ~RDFT() {
            if(fft) Free(fft);
            if(pSpecBuffer) Free(pSpecBuffer);
            if(pBuffer) Free(pBuffer);            
        }
        void Forward(const T* pSrc, T* pDst)
        {                   
            DFTFwd_RToPack<T>(pSrc,pDst,fft,pBuffer);            
        }
        void Inverse(const T* pSrc, T* pDst)
        {                   
            DFTInv_PackToR<T>(pSrc,pDst,fft,pBuffer);            
        }
        void Unpack(int len, Ipp32f * pSrc, Ipp32fc * pDst) {
            ConjPack<Ipp32f,Ipp32fc>(pSrc,pDst,len);
        }
        void Unpack(int len, Ipp64f * pSrc, Ipp64fc * pDst) {
            ConjPack<Ipp64f,Ipp64fc>(pSrc,pDst,len);
        }
    };
    
    template<typename T>
    void DCTFwdInit(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>        
    void DCTFwdInit<Ipp32f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTFwdInit_32f((IppsDCTFwdSpec_32f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>        
    void DCTFwdInit<Ipp64f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTFwdInit_64f((IppsDCTFwdSpec_64f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTInvInit(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>        
    void DCTInvInit<Ipp32f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTInvInit_32f((IppsDCTInvSpec_32f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>        
    void DCTInvInit<Ipp64f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTInvInit_64f((IppsDCTInvSpec_64f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTFwdGetSize(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBuffersize)
    {
        assert(1==0);
    }
    template<>
    void DCTFwdGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTFwdGetSize_32f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void DCTFwdGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTFwdGetSize_64f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void DCTInvGetSize(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void DCTInvGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTInvGetSize_32f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void DCTInvGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTInvGetSize_64f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void DCTFwd(const T* pSrc, T * pDst, void *, Ipp8u* pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DCTFwd<Ipp32f>(const Ipp32f* pSrc, Ipp32f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTFwd_32f(pSrc,pDst,(const IppsDCTFwdSpec_32f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DCTFwd<Ipp64f>(const Ipp64f* pSrc, Ipp64f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTFwd_64f(pSrc,pDst,(const IppsDCTFwdSpec_64f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTInv(const T* pSrc, T * pDst, void *, Ipp8u* pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DCTInv<Ipp32f>(const Ipp32f* pSrc, Ipp32f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTInv_32f(pSrc,pDst,(const IppsDCTInvSpec_32f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DCTInv<Ipp64f>(const Ipp64f* pSrc, Ipp64f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTInv_64f(pSrc,pDst,(const IppsDCTInvSpec_64f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }




    template<typename T>
    struct DCT
    {
        Ipp8u * pBuffer[2];        
        Ipp8u * pSpec[2];
        Ipp8u * pSpecBuffer[2];
        void * forward;
        void *inverse;
                        
        DCT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;

            DCTFwdGetSize<T>(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[0]       = Malloc<Ipp8u>(spec);
            pSpecBuffer[0] = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer[0]     = Malloc<Ipp8u>(size);

            int order = n;
            DCTFwdInit<T>(&forward,n,ippAlgHintNone,pSpec[0],pSpecBuffer[0]);            
            

            DCTInvGetSize<T>(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[1]       = Malloc<Ipp8u>(spec);
            pSpecBuffer[1] = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer[1]     = Malloc<Ipp8u>(size);

            DCTInvInit<T>(&inverse,n,ippAlgHintNone,pSpec[1],pSpecBuffer[1]);                        
        }
        ~DCT() {
            if(forward) Free(forward);
            if(inverse) Free(inverse);
            if(pSpecBuffer[0]) Free(pSpecBuffer[0]);
            if(pBuffer[0]) Free(pBuffer[0]);
            if(pSpecBuffer[1]) Free(pSpecBuffer[1]);
            if(pBuffer[1]) Free(pBuffer[1]);        
        }        
        void Forward(const T* pSrc, T* pDst)
        {                   
            DCTFwd<T>(pSrc,pDst,forward,pSpecBuffer[0]);            
        }
        void Inverse(const T* pSrc, T* pDst)
        {                   
            DCTInv<T>(pSrc,pDst,inverse,pSpecBuffer[1]);            
        }
    };
    
        
    struct HilbertTransform32
    {   
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;        

        HilbertTransform32(size_t n) {
            int specSize,bufferSize;
            size = n;
            IppStatus status = ippsHilbertGetSize_32f32fc(n,ippAlgHintNone,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            status = ippsHilbertInit_32f32fc(n,ippAlgHintNone,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);
        }
        ~HilbertTransform32() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
        }
        void Execute(const Ipp32f * pSrc, Ipp32fc * pDst) {
            IppStatus status = ippsHilbert_32f32fc(pSrc,pDst,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);            
        }
    };

    struct HilbertTransform64 
    {    
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;        
        HilbertTransform64(size_t n) {
            int specSize,bufferSize;
            size = n;
            IppStatus status = ippsHilbertGetSize_64f64fc(n,ippAlgHintNone,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            status = ippsHilbertInit_64f64fc(n,ippAlgHintNone,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);
        }
        ~HilbertTransform64() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
        }
        void Execute(const Ipp64f * pSrc, Ipp64fc * pDst) {
            IppStatus status = ippsHilbert_64f64fc(pSrc,pDst,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);            
        }
    };


    
    template<typename T>
    void FIRMRGetSize(int tapsLen, int upFactor, int downFactor, int * pSpecSize, int *pBufSize)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status = ippsFIRMRGetSize(tapsLen,upFactor,downFactor,dType,pSpecSize,pBufSize);
        checkStatus(status);
    }
    template<typename T1>
    void FIRMRInit(const T1* pTaps, int tapsLen, int upFactor, int upPhase, int downFact, int downPhase, void* pSpec)
    {
        assert(1==0);
    }
    template<>
    void FIRMRInit<Ipp32f>(const Ipp32f* pTaps, int tapsLen, int upFactor, int upPhase,int downFactor, int downPhase, void * pSpec)
    {        
        IppStatus status = ippsFIRMRInit_32f(pTaps,tapsLen,upFactor,upPhase,downFactor,downPhase,(IppsFIRSpec_32f*)pSpec);
        checkStatus(status);
    }
    template<>
    void FIRMRInit<Ipp64f>(const Ipp64f* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, void * pSpec)
    {        
        IppStatus status = ippsFIRMRInit_64f(pTaps,tapsLen,upFactor,upPhase,downFactor,downPhase,(IppsFIRSpec_64f*)pSpec);
        checkStatus(status);
    }

    template<typename T>
    void FIRMR_(const T * pSrc, T * pDst, int numIters, void* pSpec, const T* pDlySrc, T* pDlyDst, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FIRMR_<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, int numIters, void * pSpec, const Ipp32f* pDlySrc, Ipp32f* pDlyDst, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFIRMR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpec,pDlySrc,pDlyDst,pBuffer);
        checkStatus(status);
    }
    template<>
    void FIRMR_<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, int numIters, void * pSpec, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFIRMR_64f(pSrc,pDst,numIters,(IppsFIRSpec_64f*)pSpec,pDlySrc,pDlyDst,pBuffer);
        checkStatus(status);
    }


    template<typename T>
    struct FIRMR
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;
        T     * pTaps;
        FIRMR(size_t n, int up, int down, int upPhase=0,int downPhase=0) {
            int specSize,bufferSize;
            FIRMRGetSize<T>(n,up,down,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = Malloc<T>(n);
            FIRMRInit<T>(pTaps,n,up,upPhase,down,downPhase,(void*)pSpec);            
        }
        ~FIRMR() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) Free(pTaps);
        }
        void Execute(const T* pSrc, T* pDst, int numIters) {
            FIRMR_<T>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }
    };    

    
    template<typename T>
    void FIRSRGetSize(int tapsLen, int * pSpecSize, int *pBufSize)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status = ippsFIRSRGetSize(tapsLen,dType,pSpecSize,pBufSize);
        checkStatus(status);
    }
    template<typename T1>
    void FIRSRInit(const T1* pTaps, int tapsLen, IppAlgType algType, void* pSpec)
    {
        assert(1==0);
    }
    template<>
    void FIRSRInit<Ipp32f>(const Ipp32f* pTaps, int tapsLen, IppAlgType algType, void * pSpec)
    {
        IppStatus status = ippsFIRSRInit_32f(pTaps,tapsLen, algType,(IppsFIRSpec_32f*)pSpec);
        checkStatus(status);
    }
    template<>
    void FIRSRInit<Ipp64f>(const Ipp64f* pTaps, int tapsLen, IppAlgType algType, void * pSpec)
    {
        IppStatus status = ippsFIRSRInit_64f(pTaps,tapsLen, algType,(IppsFIRSpec_64f*)pSpec);
        checkStatus(status);
    }

    template<typename T>
    void FIRSR_(const T * pSrc, T * pDst, int numIters, void * pSpect, const T * pDlySrc, T * pDlyDst, Ipp8u* pBuf)
    {
        assert(1==0);
    }
    template<>
    void FIRSR_<Ipp32f>(const Ipp32f * pSrc, Ipp32f * pDst, int numIters, void * pSpect, const Ipp32f * pDlySrc, Ipp32f * pDlyDst, Ipp8u* pBuf)
    {
        IppStatus status = ippsFIRSR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpect,pDlySrc,pDlyDst,pBuf);
        checkStatus(status);
    }
    template<>
    void FIRSR_<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, int numIters, void * pSpect, const Ipp64f * pDlySrc, Ipp64f * pDlyDst, Ipp8u* pBuf)
    {
        IppStatus status = ippsFIRSR_64f(pSrc,pDst,numIters,(IppsFIRSpec_64f*)pSpect,pDlySrc,pDlyDst,pBuf);
        checkStatus(status);
    }

    template<typename T>
    struct FIRSR
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        T *pTaps;

        enum {
            LP,
            HP,
            BP,
            BS
        };
        int type = LP;

        FIRSR(size_t n, T* taps, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            FIRSRGetSize<T>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            nTaps   = n;
            pTaps   = Malloc<T>(n);
            memcpy(pTaps,taps,n*sizeof(T));
            FIRSRInit<T>(pTaps,n,algType,(void*)pSpec);            
        }
        FIRSR(int type, size_t n, T fc1, T fc2=0, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            FIRSRGetSize<T>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = Malloc<T>(n);
            makeFilter(type,fc1,fc2);
            FIRSRInit<T>(pTaps,n,algType,(void*)pSpec);            
        }
        ~FIRSR() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) Free(pTaps);
        }

        // bandpass/bandstop - fc1 = lower, fc2 = upper
        void makeFilter(int filterType, T fc1, T fc2=0) {
            type = filterType;
            switch(type) {
                case LP: makeLowPass(fc1);  break;
                case HP: makeHighPass(fc1);  break;
                case BP: makeBandPass(fc1,fc2);  break;
                case BS: makeBandStop(fc1,fc2);  break;
            }
        }
        void makeLowPass(T fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                Ipp32f x= fc*(T)(i-nTaps/2);            
                pTaps[i] = fc * sin(M_PI*x)/(M_PI*x);
            }
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeHighPass(T fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                Ipp32f x= fc*(Ipp32f)(i-nTaps/2);            
                pTaps[i] = (sin(M_PI*x)/(M_PI*x)) - fc * sin(M_PI*x)/(M_PI*x);
            }            
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeBandPass(T fc1, T fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                Ipp32f xl= fc1*(T)(i-nTaps/2);            
                Ipp32f xu= fc2*(T)(i-nTaps/2);            
                pTaps[i] = (xu * sin(M_PI*xu)/(M_PI*xu)) - (xl * sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }        
        void makeBandStop(T fc1, T fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                Ipp32f xl= fc1*(T)(i-nTaps/2);            
                Ipp32f xu= fc2*(T)(i-nTaps/2);            
                pTaps[i] = (xu * sin(M_PI*xu)/(M_PI*xu)) - (xl * sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }        
        void Execute(const T* pSrc, T* pDst, int numIters) {
            FIRSR_<T>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }
    };



    
    
    // warning - do not know if this works yet
    struct FIRLMS32F
    {
        Ipp32f * pDlyLine;
        Ipp8u * pBuffer;
        int  size,nTaps;
        Ipp32f *pTaps;
        IppsFIRLMSState_32f * pState;
    
        FIRLMS32F(size_t n, Ipp32f * taps) {            
            IppStatus status = ippsFIRLMSGetStateSize_32f(n,0,&size);
            checkStatus(status);
            nTaps = n;
            pTaps = ippsMalloc_32f(n);
            pDlyLine = ippsMalloc_32f(2*n);
            pBuffer = ippsMalloc_8u(size);            
            status = ippsFIRLMSInit_32f(&pState,pTaps,nTaps,pDlyLine,0,pBuffer);
            checkStatus(status);
        }
        ~FIRLMS32F() {            
            if(pBuffer) ippsFree(pBuffer);
            if(pDlyLine) ippFree(pDlyLine);            
            if(pTaps) ippsFree(pTaps);
        }        
        void Execute(const Ipp32f * pSrc, Ipp32f * pRef, Ipp32f* pDst, int len, Ipp32f mu) {
            IppStatus status = ippsFIRLMS_32f(pSrc,pRef,pDst,len,mu,pState);
            checkStatus(status);
        }
    };
    
    inline int kaiserTapsEstimate(double delta, double As)
    {
        return int((As - 8.0) / (2.285 * delta) + 0.5);
    }
    struct Resample32
    {
        IppsResamplingPolyphaseFixed_32f * pSpec;
        size_t N;
        Ipp64f Time;
        int Outlen;
        Resample32(int inRate, int outRate, Ipp32f rollf=0.9f, Ipp32f as = 80.0f) {
            int size,height,length;
            Ipp32f alpha = kaiserBeta(as);
            Ipp32f delta = (1.0f - rollf) * M_PI;
            int n = kaiserTapsEstimate(delta,as);
            IppStatus status = ippsResamplePolyphaseFixedGetSize_32f(inRate,outRate,n,&size,&length,&height,ippAlgHintFast);
            checkStatus(status);
            pSpec = (IppsResamplingPolyphaseFixed_32f*)ippsMalloc_8u(size);
            status = ippsResamplePolyphaseFixedInit_32f(inRate,outRate,n,rollf,alpha,(IppsResamplingPolyphaseFixed_32f*)pSpec,ippAlgHintFast);
            checkStatus(status);            
        }
        ~Resample32() {
            if(pSpec) ippsFree(pSpec);
        }
        void Execute(const Ipp32f * pSrc, int len, Ipp32f* pDst, Ipp64f norm=0.98) {
            IppStatus status = ippsResamplePolyphaseFixed_32f(pSrc,len,pDst,norm,&Time,&Outlen,(IppsResamplingPolyphaseFixed_32f*)pSpec);
            checkStatus(status);
        }
    };

    template<typename T>
    void IIRInit(void ** ppState, const T* pTaps, int order, const T* pDlyLine, Ipp8u * pBuf)
    {
        assert(1==0);
    }
    template<>
    void IIRInit<Ipp32f>(void ** ppState, const Ipp32f* pTaps, int order, const Ipp32f* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_32f((IppsIIRState_32f**)ppState,pTaps,order,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<>
    void IIRInit<Ipp64f>(void ** ppState, const Ipp64f* pTaps, int order, const Ipp64f* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_64f((IppsIIRState_64f**)ppState,pTaps,order,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<typename T>
    void IIRInit_BiQuad(void ** ppState, const T* pTaps, int numBq, const T* pDlyLine, Ipp8u * pBuf)
    {
        assert(1==0);
    }
    template<>
    void IIRInit_BiQuad<Ipp32f>(void ** ppState, const Ipp32f* pTaps, int numBq, const Ipp32f* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_BiQuad_32f((IppsIIRState_32f**)ppState,pTaps,numBq,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<>
    void IIRInit_BiQuad<Ipp64f>(void ** ppState, const Ipp64f* pTaps, int numBq, const Ipp64f* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_BiQuad_64f((IppsIIRState_64f**)ppState,pTaps,numBq,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetStateSize(int order, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void IIRGetStateSize<Ipp32f>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_32f(order,pBufferSize);
        checkStatus(status);
    }
    template<>
    void IIRGetStateSize<Ipp64f>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_64f(order,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetStateSize_BiQuad(int order, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void IIRGetStateSize_BiQuad<Ipp32f>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_BiQuad_32f(order,pBufferSize);
        checkStatus(status);
    }
    template<>
    void IIRGetStateSize_BiQuad<Ipp64f>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_BiQuad_64f(order,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetDlyLine(void * pState, T * pDlyLine)
    {
        assert(1==0);
    }
    template<>
    void IIRGetDlyLine<Ipp32f>(void * pState, Ipp32f * pDlyLine)
    {
        IppStatus status = ippsIIRGetDlyLine_32f((IppsIIRState_32f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<>
    void IIRGetDlyLine<Ipp64f>(void * pState, Ipp64f * pDlyLine)
    {
        IppStatus status = ippsIIRGetDlyLine_64f((IppsIIRState_64f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<typename T>
    void IIRSetDlyLine(void * pState, T * pDlyLine)
    {
        assert(1==0);
    }
    template<>
    void IIRSetDlyLine<Ipp32f>(void * pState, Ipp32f * pDlyLine)
    {
        IppStatus status = ippsIIRSetDlyLine_32f((IppsIIRState_32f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<>
    void IIRSetDlyLine<Ipp64f>(void * pState, Ipp64f * pDlyLine)
    {
        IppStatus status = ippsIIRSetDlyLine_64f((IppsIIRState_64f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<typename T>
    void IIR_(const T* pSrc, T * pDst, int len, void * pState)
    {
        assert(1==0);
    }
    template<>
    void IIR_<Ipp32f>(const Ipp32f* pSrc, Ipp32f * pDst, int len, void * pState)
    {
        IppStatus status = ippsIIR_32f(pSrc,pDst,len,(IppsIIRState_32f*)pState);
        checkStatus(status);
    }
    template<>
    void IIR_<Ipp64f>(const Ipp64f* pSrc, Ipp64f * pDst, int len, void * pState)
    {
        IppStatus status = ippsIIR_64f(pSrc,pDst,len,(IppsIIRState_64f*)pState);
        checkStatus(status);
    }    
    void IIRGenGetBufferSize(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGenGetBufferSize(order,pBufferSize);
        checkStatus(status);
    }
    
    void IIRGenLowpassButterworth(Ipp64f freq, Ipp64f ripple, int order, Ipp64f * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenLowpass_64f(freq,ripple,order,pTaps,ippButterworth,pBuffer);
        checkStatus(status);
    }
    void IIRGenHighpassButterworth(Ipp64f freq, Ipp64f ripple, int order, Ipp64f * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenHighpass_64f(freq,ripple,order,pTaps,ippButterworth,pBuffer);
        checkStatus(status);
    }
    void IIRGenLowpassChebyshev1(Ipp64f freq, Ipp64f ripple, int order, Ipp64f * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenLowpass_64f(freq,ripple,order,pTaps,ippChebyshev1,pBuffer);
        checkStatus(status);
    }
    void IIRGenHighpassChebyshev1(Ipp64f freq, Ipp64f ripple, int order, Ipp64f * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenHighpass_64f(freq,ripple,order,pTaps,ippChebyshev1,pBuffer);
        checkStatus(status);
    }


    
    template<typename T>
    struct IIR
    {
        Ipp8u * pBuffer;
        size_t  len;                
        void  * pState;
        T*      dlyLine;
        
        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIR(size_t n, int order, const T * taps) {
            int bufferSize;
            IIRGetStateSize<T>(n,&bufferSize);            
            len = n;
            pBuffer = Malloc<Ipp8u>(bufferSize);            
            IIRInit<T>(&pState,taps,order,NULL,pBuffer);            
            dlyLine = Malloc<T>(2*(order+1));
        }
        ~IIR() {
            if(pBuffer) Free(pBuffer);      
            if(dlyLine) Free(dlyLine);  
        }
        void setCoefficients(int n, int order, const T* taps) {
            if(pBuffer) Free(pBuffer);
            int bufferSize;
            len = n;
            // have to save it so it doesn't pop
            IIRGetDlyLine<T>(pState,dlyLine);            
            IIRGetStateSize<T>(n,&bufferSize);            
            pBuffer = Malloc<Ipp8u>(bufferSize);            
            IIRInit<T>(&pState,taps,order,NULL,pBuffer);            
            IIRSetDlyLine<T>(pState,dlyLine);            
        }
        void Execute(const T* pSrc, T* pDst)
        {
            IIR_<T>(pSrc,pDst,len,pState);            
        }
    };

    template<typename T>
    struct IIRBiquad
    {        
        Ipp8u * pBuffer;
        size_t  len;                
        void  * pState;
        T*      dlyLine;
        

        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIRBiquad(size_t n, int numBiquads, const T* taps) {
            int bufferSize;
            IIRGetStateSize<T>(n,&bufferSize);            
            len = n;
            pBuffer = Malloc<Ipp8u>(bufferSize);            
            IIRInit_BiQuad<T>(&pState,taps,numBiquads,NULL,pBuffer);            
            dlyLine = Malloc<T>(2*numBiquads);
        }
        ~IIRBiquad() {
            if(pBuffer) Free(pBuffer);        
            if(dlyLine) Free(dlyLine);
        }
        void setCoefficients(int order, int numBiquads, const T* taps) {
            if(pBuffer) Free(pBuffer);
            int bufferSize;
            IIRGetDlyLine<T>(pState,dlyLine);            
            IIRGetStateSize<T>(numBiquads,&bufferSize);            
            pBuffer = Malloc<Ipp8u>(bufferSize);            
            IIRInit_BiQuad<T>(&pState,taps,numBiquads,NULL,pBuffer);            
            IIRSetDlyLine<T>(pState,dlyLine);            
        }
        void Execute(const T* pSrc, T* pDst)
        {
            IIR_<T>(pSrc,pDst,len,pState);            
        }
    };

    /*
    struct STFTF32
    {
        RFFTF fft;    
        size_t blockSize,windowSize,hopSize;

        STFTF32(size_t window_size, size_t hop_size) : fft(window_size) {            
            windowSize= window_size;
            hopSize = hop_size;
        }
        std::vector<std::vector<std::complex<Ipp32f>>> 
        stft(std::vector<float> data) {
            std::vector<std::vector<std::complex<float>>> output;
            std::vector<float> _data(windowSize);
            int data_size = data.size();    
            size_t result_size = (data_size / hopSize);                       
            size_t idx = 0;
            for (size_t pos = 0; pos < data_size; pos += hopSize) {
                for (size_t i = 0; i < window_size; ++i) {
                    if (pos + i < data_size)
                        _data[i] = data[pos + i];
                    else
                        _data[i] = 0;
                }
                WindowHann(_data,_data,window_size);
                fft.set_input(_data);
                fft.Execute();
                fft.Normalize();
                temp = fft.get_output();
                output.push_back(temp);
                idx++;
            }            
            return output;
        }

        std::vector<std::vector<float>> 
        istft(complex_vector<float> data, size_t window_size, size_t hop_size) {
            int data_size = data.size();
            size_t result_size = data_size * hop_size + (window_size - hop_size);
            std::vector<std::vector<float>> output;
            std::vector<float> temp;
            std::vector<float> frame(window_size);
            std::vector<float> result(result_size);
            std::vector<std::complex<float>> slice(window_size);
                                    
            for (size_t i = 0; i < data_size; ++i) {                
                memcpy(slice.data(), data.data() + i * window_size, sizeof(fftw_complex) * window_size);
                fft.set_input(slice);
                fft.Execute();
                fft.Normalize();
                result = fft.get_output();
                output.push_back(result);
            }                
            return output;
        }
    };
    */

}