/*

    Templatized DSP vector operations
    Copyright (C) 2001-2005 Jussi Laako

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#include <dsp/DSPOp.hh>
#include <dsp/ReBufferT.hh>


#ifndef DSPVECTOR_HH
    #define DSPVECTOR_HH

    #define DSPV_SHORTER(x)     ((Size() <= x.Size()) ? Size() : x.Size())
    #define DSPV_SHORTER2(x, y) ((x.Size() <= y.Size()) ? x.Size() : y.Size())


    template <class TDSPVector_t> class clDSPVector : 
        public clReBufferT<TDSPVector_t>
    {
            long lFFTSize;
            clDSPOp DSP;
        protected:
            TDSPVector_t * GetPtr ()
                { return clReBufferT<TDSPVector_t>::GetPtr(); }
            TDSPVector_t * GetRawPtr ()
                { return clReBufferT<TDSPVector_t>::GetRawPtr(); }
        public:
            long Size ()
                { return clReBufferT<TDSPVector_t>::Size(); }
            void SetSize (long lNewSize)
                { clReBufferT<TDSPVector_t>::SetSize(lNewSize); }

            clDSPVector ()
                {
                    lFFTSize = 0;
                }
            clDSPVector (const clDSPVector &CopySrc)
                {
                    lFFTSize = 0;
                    *this = CopySrc;
                }
            clDSPVector (long lNewSize) : 
                clReBufferT<TDSPVector_t>(lNewSize)
                {
                    lFFTSize = 0;
                }
            clDSPVector (const TDSPVector_t *fpSrcData, long lSrcCount) :
                clReBufferT<TDSPVector_t>(fpSrcData, lSrcCount)
                {
                    lFFTSize = 0;
                }
            virtual ~clDSPVector () {}
            TDSPVector_t * Ptr()
                {
                    return GetPtr();
                }

            clDSPVector operator+ (const TDSPVector_t &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Add(Temp.Ptr(), Src, Temp.Size());
                    return Temp;
                }
            clDSPVector operator+ (clDSPVector &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Add(Temp.Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return Temp;
                }
            clDSPVector operator- (const TDSPVector_t &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Sub(Temp.Ptr(), Src, Temp.Size());
                    return Temp;
                }
            clDSPVector operator- (clDSPVector &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Sub(Temp.Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return Temp;
                }
            clDSPVector operator* (const TDSPVector_t &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Mul(Temp.Ptr(), Src, Temp.Size());
                    return Temp;
                }
            clDSPVector operator* (clDSPVector &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Mul(Temp.Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return Temp;
                }
            clDSPVector operator/ (const TDSPVector_t &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Div(Temp.Ptr(), Src, Temp.Size());
                    return Temp;
                }
            clDSPVector operator/ (clDSPVector &Src)
                {
                    clDSPVector<TDSPVector_t> Temp(*this);
                    DSP.Div(Temp.Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return Temp;
                }
            clDSPVector & operator+= (const TDSPVector_t &Src)
                {
                    DSP.Add(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & operator+= (clDSPVector &Src)
                {
                    DSP.Add(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & operator-= (const TDSPVector_t &Src)
                {
                    DSP.Sub(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & operator-= (clDSPVector &Src)
                {
                    DSP.Sub(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & operator*= (const TDSPVector_t &Src)
                {
                    DSP.Mul(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & operator*= (clDSPVector &Src)
                {
                    DSP.Mul(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & operator/= (const TDSPVector_t &Src)
                {
                    DSP.Div(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & operator/= (clDSPVector &Src)
                {
                    DSP.Div(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }

            clDSPVector & Zero ()
                {
                    DSP.Zero(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Set (const TDSPVector_t &Src)
                {
                    DSP.Set(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & Set (const TDSPVector_t *Src, long lSrcCount)
                {
                    SetSize(lSrcCount);
                    DSP.Set(Ptr(), Src, lSrcCount);
                    return (*this);
                }
            clDSPVector & Set (const TDSPVector_t &Src, long lStartIdx,
                long lSetCount)
                {
                    DSP.Set(Ptr(), Src, lStartIdx, lSetCount, Size());
                    return (*this);
                }
            clDSPVector & Clip (const TDSPVector_t &Src)
                {
                    DSP.Clip(Ptr(), Src, Size());
                    return (*this);
                }
            clDSPVector & Clip (const TDSPVector_t &Src1, 
                const TDSPVector_t &Src2)
                {
                    DSP.Clip(Ptr(), Src1, Src2, Size());
                    return (*this);
                }
            clDSPVector & ClipZero ()
                {
                    DSP.ClipZero(Ptr(), Size());
                    return (*this);
                }

            clDSPVector & MulC (clDSPVector &Src)
                {
                    DSP.MulC(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & MulC (clDSPVector &Src1, clDSPVector &Src2)
                {
                    *this = Src1;
                    DSP.MulC(Ptr(), Src2.Ptr(), DSPV_SHORTER(Src2));
                    return (*this);
                }
            clDSPVector & Div1x ()
                {
                    DSP.Div1x(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & MulAdd (TDSPVector_t fMul, TDSPVector_t fAdd)
                {
                    DSP.MulAdd(Ptr(), fMul, fAdd, Size());
                }
            clDSPVector & MulAddC (TDSPVector_t spSrc1, TDSPVector_t spSrc2)
                {
                    DSP.MulAddC(Ptr(), spSrc1, spSrc2, Size());
                }
            clDSPVector & Square ()
                {
                    DSP.Square(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Abs ()
                {
                    DSP.Abs(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Sqrt ()
                {
                    DSP.Sqrt(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Negate ()
                {
                    DSP.Negate(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Normalize ()
                {
                    DSP.Normalize(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Reverse ()
                {
                    DSP.Reverse(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Reverse (clDSPVector &Src)
                {
                    SetSize(Src.Size());
                    DSP.Reverse(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Scale ()
                {
                    DSP.Scale(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Scale (clDSPVector &Src)
                {
                    SetSize(Src.Size());
                    DSP.Scale(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Scale01 ()
                {
                    DSP.Scale01(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Scale01 (clDSPVector &Src)
                {
                    SetSize(Src.Size());
                    DSP.Scale01(Ptr(), Src.Ptr(), DSPV_SHOTER(Src));
                    return (*this);
                }
            clDSPVector & Sort ()
                {
                    DSP.Sort(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & Sort (clDSPVector &Src)
                {
                    *this = Src;
                    DSP.Sort(Ptr(), Size());
                    return (*this);
                }

            TDSPVector_t Sum ()
                {
                    return DSP.Sum(Ptr(), Size());
                }
            TDSPVector_t Conv (clDSPVector &Src)
                {
                    return DSP.Convolve(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                }
            clDSPVector & Conv (clDSPVector &Src1, clDSPVector &Src2)
                {
                    SetSize(DSPV_SHORTER2(Src1, Src2));
                    DSP.Convolve(Ptr(), Src1.Ptr(), Src2.Ptr(),
                        DSPV_SHORTER2(Src1, Src2));
                    return (*this);
                }
            TDSPVector_t Corr (clDSPVector &Src)
                {
                    return DSP.Correlate(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                }
            clDSPVector & Corr (clDSPVector &Src1, clDSPVector &Src2)
                {
                    SetSize(DSPV_SHORTER2(Src1, Src2));
                    DSP.Correlate(Ptr(), Src1.Ptr(), Src2.Ptr(),
                        DSPV_SHORTER2(Src1, Src2));
                    return (*this);
                }
            TDSPVector_t AutoCorr ()
                {
                    return DSP.AutoCorrelate(Ptr(), Size());
                }
            clDSPVector & AutoCorr (clDSPVector &Src)
                {
                    SetSize(Src.Size());
                    DSP.AutoCorrelate(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                }
            TDSPVector_t CrossCorr (clDSPVector &Src)
                {
                    return DSP.CrossCorr(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                }
            TDSPVector_t CrossCorr (clDSPVector &Src, long lDelay)
                {
                    return DSP.DelCrossCorr(Ptr(), Src.Ptr(), lDelay,
                        DSPV_SHORTER(Src));
                }
            clDSPVector & CrossCorr (clDSPVector &Src1, clDSPVector &Src2,
                long *lpDelays, long lDelayCount)
                {
                    SetSize(lDelayCount);
                    DSP.DelCrossCorr(Ptr(), Src1.Ptr(), Src2.Ptr(),
                        DSPV_SHORTER2(Src1, Src2), lpDelays, lDelayCount);
                    return (*this);
                }
            TDSPVector_t DotProduct (clDSPVector &Src)
                {
                    return DSP.DotProduct(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                }
            TDSPVector_t Mean ()
                {
                    return DSP.Mean(Ptr(), Size());
                }
            TDSPVector_t Median ()
                {
                    return DSP.Median(Ptr(), Size());
                }
            TDSPVector_t Energy ()
                {
                    return DSP.Energy(Ptr(), Size());
                }
            TDSPVector_t RMS ()
                {
                    return DSP.RMS(Ptr(), Size());
                }
            TDSPVector_t PeakLevel ()
                {
                    return DSP.PeakLevel(Ptr(), Size());
                }
            void MinMax (TDSPVector_t &fMin, TDSPVector_t &fMax)
                {
                    DSP.MinMax(&fMin, &fMax, Ptr(), Size());
                }
            void StdDev (TDSPVector_t &fStdDev, TDSPVector_t &fMean)
                {
                    DSP.StdDev(&fStdDev, &fMean, Ptr(), Size());
                }
            void Variance (TDSPVector_t &fVariance, TDSPVector_t &fMean)
                {
                    DSP.Variance(&fVariance, &fMean, Ptr(), Size());
                }

            clDSPVector & Convert (const unsigned char *ucpSrc, 
                long lSrcCount)
                {
                    SetSize(lSrcCount);
                    DSP.Convert(Ptr(), ucpSrc, lSrcCount);
                    return (*this);
                }
            clDSPVector & Convert (const signed short *ipSrc, long lSrcCount,
                bool b12bit = false)
                {
                    SetSize(lSrcCount);
                    DSP.Convert(Ptr(), ipSrc, lSrcCount, b12bit);
                    return (*this);
                }
            clDSPVector & Convert (const signed int *ipSrc, long lSrcCount,
                bool b24bit = false)
                {
                    SetSize(lSrcCount);
                    DSP.Convert(Ptr(), ipSrc, lSrcCount, b24bit);
                    return (*this);
                }
            clDSPVector & Convert (const float *fpSrc, long lSrcCount)
                {
                    SetSize(lSrcCount);
                    DSP.Convert(Ptr(), fpSrc, lSrcCount);
                    return (*this);
                }
            clDSPVector & Convert (const double *dpSrc, long lSrcCount)
                {
                    SetSize(lSrcCount);
                    DSP.Convert(Ptr(), dpSrc, lSrcCount);
                    return (*this);
                }
            void Convert (unsigned char *ucpDest)
                {
                    DSP.Convert(ucpDest, Ptr(), Size());
                }
            void Convert (signed short *ipDest, bool b12bit = false)
                {
                    DSP.Convert(ipDest, Ptr(), Size(), b12bit);
                }
            void Convert (signed int *ipDest, bool b24bit = false)
                {
                    DSP.Convert(ipDest, Ptr(), Size(), b24bit);
                }
            void Convert (float *fpDest)
                {
                    DSP.Convert(fpDest, Ptr(), Size());
                }
            void Convert (double *dpDest)
                {
                    DSP.Convert(dpDest, Ptr(), Size());
                }
            clDSPVector & CartToPolar ()
                {
                    DSP.CartToPolar(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & CartToPolar (clDSPVector &Cart)
                {
                    SetSize(Cart.Size());
                    DSP.CartToPolar(Ptr(), Cart.Ptr(), Size());
                    return (*this);
                }
            void CartToPolar (clDSPVector<float> &Magn, 
                clDSPVector<float> &Phase)
                {
                    Magn.SetSize(Size());
                    Phase.SetSize(Size());
                    DSP.CartToPolar(Magn.Ptr(), Phase.Ptr(), Ptr(), Size());
                }
            void CartToPolar (clDSPVector<double> &Magn, 
                clDSPVector<double> &Phase)
                {
                    Magn.SetSize(Size());
                    Phase.SetSize(Size());
                    DSP.CartToPolar(Magn.Ptr(), Phase.Ptr(), Ptr(), Size());
                }
            clDSPVector & PolarToCart ()
                {
                    DSP.PolarToCart(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & PolarToCart (clDSPVector &Polar)
                {
                    SetSize(Polar.Size());
                    DSP.PolarToCart(Ptr(), Polar.Ptr(), Size());
                    return (*this);
                }
            void PolarToCart (clDSPVector<float> &Real, 
                clDSPVector<float> &Imag)
                {
                    Real.SetSize(Size());
                    Imag.SetSize(Size());
                    DSP.PolarToCart(Real.Ptr(), Imag.Ptr(), Ptr(), Size());
                }
            void PolarToCart (clDSPVector<double> &Real, 
                clDSPVector<double> &Imag)
                {
                    Real.SetSize(Size());
                    Imag.SetSize(Size());
                    DSP.PolarToCart(Real.Ptr(), Imag.Ptr(), Ptr(), Size());
                }
            clDSPVector & Magnitude (clDSPVector<stSCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Magnitude(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Magnitude (clDSPVector<stDCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Magnitude(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Power (clDSPVector<stSCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Power(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Power (clDSPVector<stDCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Power(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Phase (clDSPVector<stSCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Phase(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            clDSPVector & Phase (clDSPVector<stDCplx> &Src)
                {
                    SetSize(Src.Size());
                    DSP.Phase(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            void PowerPhase (clDSPVector<float> &Power, 
                clDSPVector<float> &Phase)
                {
                    Power.SetSize(Size());
                    Phase.SetSize(Size());
                    DSP.PowerPhase(Power.Ptr(), Phase.Ptr(), Ptr(), Size());
                }
            void PowerPhase (clDSPVector<double> &Power, 
                clDSPVector<double> &Phase)
                {
                    Power.SetSize(Size());
                    Phase.SetSize(Size());
                    DSP.PowerPhase(Power.Ptr(), Phase.Ptr(), Ptr(), Size());
                }
            clDSPVector & Decimate (long lFactor)
                {
                    DSP.Decimate(Ptr(), Ptr(), lFactor, Size());
                    return (*this);
                }
            clDSPVector & Decimate (clDSPVector &Src, long lFactor)
                {
                    SetSize(Src.Size() / lFactor);
                    DSP.Decimate(Ptr(), Src.Ptr(), lFactor, Src.Size());
                    return (*this);
                }
            clDSPVector & DecimateAvg (clDSPVector &Src, long lFactor)
                {
                    SetSize(Src.Size() / lFactor);
                    DSP.DecimateAvg(Ptr(), Src.Ptr(), lFactor, Src.Size());
                    return (*this);
                }
            clDSPVector & Interpolate (clDSPVector &Src, long lFactor)
                {
                    SetSize(Src.Size() * lFactor);
                    DSP.Interpolate(Ptr(), Src.Ptr(), lFactor, Src.Size());
                    return (*this);
                }
            clDSPVector & InterpolateAvg (clDSPVector &Src, long lFactor)
                {
                    SetSize(Src.Size() * lFactor);
                    DSP.InterpolateAvg(Ptr(), Src.Ptr(), lFactor, Src.Size());
                    return (*this);
                }
            clDSPVector & Resample (clDSPVector &Src)
                {
                    DSP.Resample(Ptr(), Size(), Src.Ptr(), Src.Size());
                    return (*this);
                }
            clDSPVector & ResampleAvg (clDSPVector &Src)
                {
                    DSP.ResampleAvg(Ptr(), Size(), Src.Ptr(), Src.Size());
                    return (*this);
                }

            clDSPVector & WinBartlett (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinBartlett(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinBlackman (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinBlackman(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinBlackman (long lWinSize, TDSPVector_t fAlpha)
                {
                    SetSize(lWinSize);
                    DSP.WinBlackman(Ptr(), Size(), fAlpha);
                    return (*this);
                }
            clDSPVector & WinBlackmanHarris (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinBlackmanHarris(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinCosTapered (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinCosTapered(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinExactBlackman (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinExactBlackman(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinExp (TDSPVector_t fZ, long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinExp(Ptr(), fZ, Size());
                    return (*this);
                }
            clDSPVector & WinFlatTop (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinFlatTop(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinGenericCos (clDSPVector &Coeffs, long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinGenericCos(Ptr(), Size(), 
                        Coeffs.Ptr(), Coeffs.Size());
                    return (*this);
                }
            clDSPVector & WinHamming (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinHamming(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinHanning (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinHanning(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinKaiser (TDSPVector_t fBeta, long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinKaiser(Ptr(), fBeta, Size());
                    return (*this);
                }
            clDSPVector & WinKaiserBessel (TDSPVector_t fAlpha, long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinKaiserBessel(Ptr(), fAlpha, Size());
                    return (*this);
                }
            clDSPVector & WinTukey (long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinTukey(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & WinDolphChebyshev (TDSPVector_t fGamma, long lWinSize)
                {
                    SetSize(lWinSize);
                    DSP.WinDolphChebyshev(Ptr(), fGamma, Size());
                    return (*this);
                }

            clDSPVector & Mix (clDSPVector &Src, long lChCount)
                {
                    SetSize(Src.Size() / lChCount);
                    DSP.Mix(Ptr(), Src.Ptr(), lChCount, Size());
                    return (*this);
                }
            clDSPVector & Extract (clDSPVector &Src, long lCh, long lChCount)
                {
                    SetSize(Src.Size() / lChCount);
                    DSP.Extract(Ptr(), Src.Ptr(), lCh, lChCount, Src.Size());
                    return (*this);
                }
            clDSPVector & Pack (clDSPVector &Src, long lCh, long lChCount)
                {
                    Resize(Src.Size() * lChCount);
                    DSP.Pack(Ptr(), Src.Ptr(), lCh, lChCount, Src.Size());
                    return (*this);
                }
            /// \note Vector must be sized before calling this!
            clDSPVector & FFTWConvert (clDSPVector &Src)
                {
                    DSP.FFTWConvert(Ptr(), Src.Ptr(), Src.Size());
                    return (*this);
                }

            void FIRAllocate (clDSPVector &Src)
                {
                    DSP.FIRAllocate(Src.Ptr(), Src.Size());
                }
            clDSPVector & FIRFilter ()
                {
                    DSP.FIRFilter(Ptr(), Size());
                    return (*this);
                }
            clDSPVector & FIRFilter (clDSPVector &Src)
                {
                    SetSize(Src.Size());
                    DSP.FIRFilter(Ptr(), Src.Ptr(), DSPV_SHORTER(Src));
                    return (*this);
                }
            void FIRFree ()
                {
                    DSP.FIRFree();
                }

            void FFTInitialize (long lWinSize, bool bReal)
                {
                    lFFTSize = lWinSize;
                    DSP.FFTInitialize(lFFTSize, bReal);
                }
            clDSPVector & FFT (clDSPVector<float> &Src)
                {
                    SetSize(lFFTSize / 2 + 1);
                    DSP.FFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & FFT (clDSPVector<double> &Src)
                {
                    SetSize(lFFTSize / 2 + 1);
                    DSP.FFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & FFT (clDSPVector<stSCplx> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.FFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & FFT (clDSPVector<stDCplx> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.FFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & IFFT (clDSPVector<stSCplx> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.IFFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & IFFT (clDSPVector<stDCplx> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.IFFTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            void FFTUninitialize ()
                {
                    DSP.FFTUninitialize();
                }
                
            void DCTInitialize (long lWinSize)
                {
                    lFFTSize = lWinSize;
                    DSP.DCTInitialize(lWinSize);
                }
            void DCTUninitialize ()
                {
                    DSP.DCTUninitialize();
                }
            clDSPVector & DCT ()
                {
                    DSP.DCTi(Ptr());
                    return (*this);
                }
            clDSPVector & DCT (clDSPVector<float> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.DCTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & DCT (clDSPVector<double> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.DCTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & IDCT ()
                {
                    DSP.IDCTi(Ptr());
                    return (*this);
                }
            clDSPVector & IDCT (clDSPVector<float> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.IDCTo(Ptr(), Src.Ptr());
                    return (*this);
                }
            clDSPVector & IDCT (clDSPVector<double> &Src)
                {
                    SetSize(lFFTSize);
                    DSP.IDCTo(Ptr(), Src.Ptr());
                    return (*this);
                }
    };


#endif
