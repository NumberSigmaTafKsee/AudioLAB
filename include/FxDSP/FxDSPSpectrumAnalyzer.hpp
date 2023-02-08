#pragma once

namespace DSPFX
{

    struct FxSpectrumAnalyzer
    {
        #ifdef DSPFLOATDOUBLE
        SpectrumAnalyzerD *spc;
        #else
        SpectrumAnalyzer *spc;
        #endif

        FxSpectrumAnalyzer(size_t len, DspFloatType sr=44100) {
            #ifdef DSPFLOATDOUBLE
            spc = SpectrumAnalyzerInitD(len,sr);
            #else
            spc = SpectrumAnalyzerInit(len,sr);
            #endif
            assert(spc != NULL);
        }
        ~FxSpectrumAnalyzer() {
            #ifdef DSPFLOATDOUBLE
            if(spc) SpectrumAnalyzerFreeD(spc);
            #else
            if(spc) SpectrumAnalyzerFree(spc);
            #endif
        }
        DspFloatType getCentroid() {
            #ifdef DSPFLOATDOUBLE
            return SpectralCentroidD(spc);
            #else
            return SpectralCentroid(spc);
            #endif
        }
        DspFloatType getSpread() {
            #ifdef DSPFLOATDOUBLE
            return SpectralSpreadD(spc);
            #else
            return SpectralSpread(spc);
            #endif
        }
        DspFloatType getSkewness() {
            #ifdef DSPFLOATDOUBLE
            return SpectralSkewnessD(spc);
            #else
            return SpectralSkewness(spc);
            #endif
        }
        DspFloatType getKurtosis() {
            #ifdef DSPFLOATDOUBLE
            return SpectralKurtosisD(spc);
            #else
            return SpectralKurtosis(spc);
            #endif
        }
        // todo : should be able to get magnitude and phase
        // should be able to set the bins and inverse it to samples too    
        void Analyze(DspFloatType * signal) {
            #ifdef DSPFLOATDOUBLE
            SpectrumAnalyzerAnalyzeD(spc,signal);
            #else
            SpectrumAnalyzerAnalyze(spc,signal);
            #endif
        }
    };
}
