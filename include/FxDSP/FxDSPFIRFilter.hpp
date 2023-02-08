#pragma once

namespace DSPFX
{

    struct FxFIRFilter : public MonoFXProcessor
    {
        /** Convolution Algorithm to use */
        enum 
        {
            /** Choose the best algorithm based on filter size */
            BEST    = 0,

            /** Use direct convolution */
            DIRECT  = 1,

            /** Use FFT Convolution (Better for longer filter kernels */
            FFT     = 2

        };

        #ifdef DSPFLOATDOUBLE
        FIRFilterD * fir;
        #else
        FIRFilter * fir;
        #endif

        FxFIRFilter(DspFloatType * kernel, size_t len, int mode) : MonoFXProcessor() {
            #ifdef DSPFLOATDOUBLE
            fir = FIRFilterInitD(kernel,len,(ConvolutionMode_t)mode);
            #else
            fir = FIRFilterInit(kernel,len,(ConvolutionMode_t)mode);
            #endif
            assert(fir != NULL);
        }
        ~FxFIRFilter() {
            #ifdef DSPFLOATDOUBLE
            if(fir) FIRFilterFreeD(fir);
            #else
            if(fir) FIRFilterFree(fir);
            #endif
        }
        void flush() { 
            #ifdef DSPFLOATDOUBLE
            FIRFilterFlushD(fir); 
            #else
            FIRFilterFlush(fir); 
            #endif
        }

        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            FIRFilterProcessD(fir,out,in,n);
            #else
            FIRFilterProcess(fir,out,in,n);
            #endif
        }
        void updateKernel(const DspFloatType * kernel) {
            #ifdef DSPFLOATDOUBLE
            FIRFilterUpdateKernelD(fir,kernel);
            #else
            FIRFilterUpdateKernel(fir,kernel);
            #endif
        }
        enum {
            PORT_KERNEL,
        };
        void setPortV(int port, const std::vector<DspFloatType> & k) {
            if(port == PORT_KERNEL) updateKernel(k.data());
            else printf("No port %d",port);
        }
    };
}
