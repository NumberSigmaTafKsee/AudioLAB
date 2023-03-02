// https://github.com/jpcima/rezonateur
#pragma once

#include <complex>
#include <memory>
#include <cstring>
#include <cassert>

#include "caps.hpp"

namespace Filters::Formant::Rezonateur
{
    //==============================================================================

    /** The type of filter that the State Variable Filter will output. */
    enum SVFType {
        SVFLowpass = 0,
        SVFBandpass,
        SVFHighpass,
        SVFUnitGainBandpass,
        SVFBandShelving,
        SVFNotch,
        SVFAllpass,
        SVFPeak
    };

    //==============================================================================
    class VAStateVariableFilter {
    public:
        /** Create and initialize the filter with default values defined in constructor. */
        VAStateVariableFilter();

        //------------------------------------------------------------------------------

        ~VAStateVariableFilter();

        //------------------------------------------------------------------------------

        /**    Sets the type of the filter that processAudioSample() or processAudioBlock() will
            output. This filter can choose between 8 different types using the enums listed
            below or the int given to each.
            0: SVFLowpass
            1: SVFBandpass
            2: SVFHighpass
            3: SVFUnitGainBandpass
            4: SVFBandShelving
            5: SVFNotch
            6: SVFAllpass
            7: SVFPeak
        */
        void setFilterType(int newType);

        //------------------------------------------------------------------------------
        /**    Used for changing the filter's cutoff parameter linearly by frequency (Hz) */
        void setCutoffFreq(double newCutoffFreq);

        //------------------------------------------------------------------------------
        /** Used for setting the resonance amount. This is then converted to a Q
            value, which is used by the filter.
            Range: (0-1)
        */
        void setResonance(double newResonance);

        //------------------------------------------------------------------------------
        /** Used for setting the filter's Q amount. This is then converted to a
            damping parameter called R, which is used in the original filter.
        */
        void setQ(double newQ);

        //------------------------------------------------------------------------------
        /**    Sets the gain of the shelf for the BandShelving filter only. */
        void setShelfGain(double newGain);

        //------------------------------------------------------------------------------
        /**    Statically set the filters parameters. */
        void setFilter(int newType, double newCutoff,
                    double newResonance, double newShelfGain);

        //------------------------------------------------------------------------------
        /**    Set the sample rate used by the host. Needs to be used to accurately
            calculate the coefficients of the filter from the cutoff.
            Note: This is often used in AudioProcessor::prepareToPlay
        */
        void setSampleRate(double newSampleRate);

        //------------------------------------------------------------------------------
        /**    Performs the actual processing.
        */
        void process(float gain, const float *input, float *output, unsigned count);

        //------------------------------------------------------------------------------
        /**    Reset the state variables.
        */
        void clear();

        //------------------------------------------------------------------------------
        /**    Compute the transfer function at given frequency.
        */
        std::complex<double> calcTransfer(double freq) const;

        //------------------------------------------------------------------------------


        double getCutoffFreq() const { return cutoffFreq; }

        double getFilterType() const { return filterType; }

        double getQ() const { return Q; }

        double getShelfGain() const { return shelfGain; }

    private:
        //==============================================================================
        //    Calculate the coefficients for the filter based on parameters.
        void calcFilter();

        //
        template <int FilterType>
        void processInternally(float gain, const float *input, float *output, unsigned count);

        //    Parameters:
        int filterType;
        double cutoffFreq;
        double Q;
        double shelfGain;

        double sampleRate;

        //    Coefficients:
        double gCoeff;        // gain element
        double RCoeff;        // feedback damping element
        double KCoeff;        // shelf gain element

        double z1_A, z2_A;        // state variables (z^-1)

    private:
        static std::complex<double> calcTransferLowpass(double w, double wc, double r);
        static std::complex<double> calcTransferBandpass(double w, double wc, double r);
        static std::complex<double> calcTransferHighpass(double w, double wc, double r);
        static std::complex<double> calcTransferUnitGainBandpass(double w, double wc, double r);
        static std::complex<double> calcTransferBandShelving(double w, double wc, double r, double k);
        static std::complex<double> calcTransferNotch(double w, double wc, double r);
        static std::complex<double> calcTransferAllpass(double w, double wc, double r);
        static std::complex<double> calcTransferPeak(double w, double wc, double r);
    };

    class Rezonateur {
    public:
        void init(double samplerate);

        void setFilterMode(int mode);
        void setFilterGain(unsigned nth, float gain);
        void setFilterCutoff(unsigned nth, float cutoff);
        void setFilterEmph(unsigned nth, float emph);
        int getFilterMode() const;
        float getFilterGain(unsigned nth) const;
        float getFilterCutoff(unsigned nth) const;
        float getFilterEmph(unsigned nth) const;

        unsigned getOversampling() const;
        void setOversampling(unsigned oversampling);

        void process(const float *input, float *output, unsigned count);

        double getResponseGain(double f) const;

        enum Mode {
            LowpassMode,
            BandpassMode,
            HighpassMode,
            BandpassNotchMode,
        };

    private:
        template <class Oversampler> void processOversampled(Oversampler &oversampler, const float *input, float *output, unsigned count);
        template <class Oversampler> void processWithinBufferLimit(Oversampler &oversampler, const float *input, float *output, unsigned count);
        void getEffectiveFilterGains(float gains[3]) const;
        static int getFilterTypeForMode(int mode);

    private:
        int fMode;
        float fFilterGains[3];
        float fFilterCutoffFreqs[3];
        float fFilterQ[3];
        VAStateVariableFilter fFilters[3];

        unsigned fOversampling;

        DSP::Oversampler<2, 32> fOversampler2x;
        DSP::Oversampler<4, 64> fOversampler4x;
        DSP::Oversampler<8, 64> fOversampler8x;
        enum { MaximumOversampling = 8 };

    private:
        void allocateWorkBuffers(unsigned count);
        float *getWorkBuffer(unsigned index);

    private:
        unsigned fNumWorkBuffers = 0;
        std::unique_ptr<float[]> fWorkBuffers;
        static constexpr unsigned sBufferLimit = 256;
    };

    #if __cplusplus >= 201703L
    # define if_constexpr if constexpr
    #else
    # define if_constexpr if
    #endif

    //==============================================================================

    static double resonanceToQ(double resonance)
    {
        return 1.0 / (2.0 * (1.0 - resonance));
    }

    //==============================================================================

    VAStateVariableFilter::VAStateVariableFilter()
    {
        sampleRate = 44100.0;                // default sample rate when constructed
        filterType = SVFLowpass;            // lowpass filter by default

        gCoeff = 1.0;
        RCoeff = 1.0;
        KCoeff = 0.0;

        cutoffFreq = 1000.0;
        Q = resonanceToQ(0.5);
        shelfGain = 1.0;

        z1_A = 0.0;
        z2_A = 0.0;
    }

    VAStateVariableFilter::~VAStateVariableFilter()
    {
    }

    // Member functions for setting the filter's parameters (and sample rate).
    //==============================================================================
    void VAStateVariableFilter::setFilterType(int newType)
    {
        filterType = newType;
    }

    void VAStateVariableFilter::setCutoffFreq(double newCutoffFreq)
    {
        if (cutoffFreq == newCutoffFreq)
            return;

        cutoffFreq = newCutoffFreq;
        calcFilter();
    }

    void VAStateVariableFilter::setResonance(double newResonance)
    {
        setQ(resonanceToQ(newResonance));
    }

    void VAStateVariableFilter::setQ(double newQ)
    {
        if (Q == newQ)
            return;

        Q = newQ;
        calcFilter();
    }

    void VAStateVariableFilter::setShelfGain(double newGain)
    {
        if (shelfGain == newGain)
            return;

        shelfGain = newGain;
        calcFilter();
    }

    void VAStateVariableFilter::setFilter(int newType, double newCutoffFreq,
                                        double newResonance, double newShelfGain)
    {
        double newQ = resonanceToQ(newResonance);

        if (filterType == newType && cutoffFreq == newCutoffFreq && Q == newQ && shelfGain == newShelfGain)
            return;

        filterType = newType;
        cutoffFreq = newCutoffFreq;
        Q = newQ;
        shelfGain = newShelfGain;
        calcFilter();
    }

    void VAStateVariableFilter::setSampleRate(double newSampleRate)
    {
        if (sampleRate == newSampleRate)
            return;

        sampleRate = newSampleRate;
        calcFilter();
    }

    //==============================================================================
    void VAStateVariableFilter::calcFilter()
    {
        // prewarp the cutoff (for bilinear-transform filters)
        double wd = cutoffFreq * (2.0 * M_PI);
        double T = 1.0 / sampleRate;
        double wa = (2.0 / T) * std::tan(wd * T / 2.0);

        // Calculate g (gain element of integrator)
        gCoeff = wa * T / 2.0;            // Calculate g (gain element of integrator)

        // Calculate Zavalishin's R from Q (referred to as damping parameter)
        RCoeff = 1.0 / (2.0 * Q);

        // Gain for BandShelving filter
        KCoeff = shelfGain;
    }

    static double analogSaturate(double x)
    {
        // simple filter analog saturation

        if (x > +1)
            x = 2. / 3.;
        else if (x < -1)
            x = -2. / 3.;
        else
            x = x - (x * x * x) * (1.0 / 3.0);

        return x;
    }

    template <int FilterType>
    void VAStateVariableFilter::processInternally(float gain, const float *input, float *output, unsigned count)
    {
        const double gCoeff = this->gCoeff;
        const double RCoeff = this->RCoeff;
        const double KCoeff = this->KCoeff;

        double z1_A = this->z1_A;
        double z2_A = this->z2_A;

        for (unsigned i = 0; i < count; ++i) {
            double in = gain * input[i];

            double HP = (in - ((2.0 * RCoeff + gCoeff) * z1_A) - z2_A)
                * (1.0 / (1.0 + (2.0 * RCoeff * gCoeff) + gCoeff * gCoeff));
            double BP = HP * gCoeff + z1_A;
            double LP = BP * gCoeff + z2_A;

            z1_A = analogSaturate(gCoeff * HP + BP);        // unit delay (state variable)
            z2_A = analogSaturate(gCoeff * BP + LP);        // unit delay (state variable)

            // Selects which filter type this function will output.
            double out = 0.0;
            if_constexpr (FilterType == SVFLowpass)
                out = LP;
            else if_constexpr (FilterType == SVFBandpass)
                out = BP;
            else if_constexpr (FilterType == SVFHighpass)
                out = HP;
            else if_constexpr (FilterType == SVFUnitGainBandpass)
                out = 2.0 * RCoeff * BP;
            else if_constexpr (FilterType == SVFBandShelving)
                out = in + 2.0 * RCoeff * KCoeff * BP;
            else if_constexpr (FilterType == SVFNotch)
                out = in - 2.0 * RCoeff * BP;
            else if_constexpr (FilterType == SVFAllpass)
                out = in - 4.0 * RCoeff * BP;
            else if_constexpr (FilterType == SVFPeak)
                out = LP - HP;

            output[i] = out;
        }

        this->z1_A = z1_A;
        this->z2_A = z2_A;
    }

    void VAStateVariableFilter::process(float gain, const float *input, float *output, unsigned count)
    {
        switch (filterType) {
        case SVFLowpass:
            processInternally<SVFLowpass>(gain, input, output, count);
            break;
        case SVFBandpass:
            processInternally<SVFBandpass>(gain, input, output, count);
            break;
        case SVFHighpass:
            processInternally<SVFHighpass>(gain, input, output, count);
            break;
        case SVFUnitGainBandpass:
            processInternally<SVFUnitGainBandpass>(gain, input, output, count);
            break;
        case SVFBandShelving:
            processInternally<SVFBandShelving>(gain, input, output, count);
            break;
        case SVFNotch:
            processInternally<SVFNotch>(gain, input, output, count);
            break;
        case SVFAllpass:
            processInternally<SVFAllpass>(gain, input, output, count);
            break;
        case SVFPeak:
            processInternally<SVFPeak>(gain, input, output, count);
            break;
        default:
            for (unsigned i = 0; i < count; ++i)
                output[i] = gain * input[i];
        }
    }

    void VAStateVariableFilter::clear()
    {
        z1_A = 0;
        z2_A = 0;
    }

    std::complex<double> VAStateVariableFilter::calcTransfer(double freq) const
    {
        double w = 2 * M_PI * freq;
        double wc = 2 * M_PI * cutoffFreq;

        switch (filterType) {
        case SVFLowpass:
            return calcTransferLowpass(w, wc, RCoeff);
        case SVFBandpass:
            return calcTransferBandpass(w, wc, RCoeff);
        case SVFHighpass:
            return calcTransferHighpass(w, wc, RCoeff);
        case SVFUnitGainBandpass:
            return calcTransferUnitGainBandpass(w, wc, RCoeff);
        case SVFBandShelving:
            return calcTransferBandShelving(w, wc, RCoeff, shelfGain);
        case SVFNotch:
            return calcTransferNotch(w, wc, RCoeff);
        case SVFAllpass:
            return calcTransferAllpass(w, wc, RCoeff);
        case SVFPeak:
            return calcTransferPeak(w, wc, RCoeff);
        default:
            return 0.0;
        }
    }

    //==============================================================================

    std::complex<double> VAStateVariableFilter::calcTransferLowpass(double w, double wc, double r)
    {
        std::complex<double> s = w * std::complex<double>(0, 1);
        return (wc * wc) / (s * s + 2.0 * r * wc * s + wc * wc);
    }

    std::complex<double> VAStateVariableFilter::calcTransferBandpass(double w, double wc, double r)
    {
        std::complex<double> s = w * std::complex<double>(0, 1);
        return (wc * s) / (s * s + 2.0 * r * wc * s + wc * wc);
    }

    std::complex<double> VAStateVariableFilter::calcTransferHighpass(double w, double wc, double r)
    {
        std::complex<double> s = w * std::complex<double>(0, 1);
        return (s * s) / (s * s + 2.0 * r * wc * s + wc * wc);
    }

    std::complex<double> VAStateVariableFilter::calcTransferUnitGainBandpass(double w, double wc, double r)
    {
        return 2.0 * r * calcTransferBandpass(w, wc, r);
    }

    std::complex<double> VAStateVariableFilter::calcTransferBandShelving(double w, double wc, double r, double k)
    {
        return 1.0 + k * calcTransferUnitGainBandpass(w, wc, r);
    }

    std::complex<double> VAStateVariableFilter::calcTransferNotch(double w, double wc, double r)
    {
        return calcTransferBandShelving(w, wc, r, -1.0);
    }

    std::complex<double> VAStateVariableFilter::calcTransferAllpass(double w, double wc, double r)
    {
        return calcTransferBandShelving(w, wc, r, -2.0);
    }

    std::complex<double> VAStateVariableFilter::calcTransferPeak(double w, double wc, double r)
    {
        std::complex<double> s = w * std::complex<double>(0, 1);
        return (wc * wc - s * s) / (s * s + 2.0 * r * wc * s + wc * wc);
    }

    void Rezonateur::init(double samplerate)
    {
        allocateWorkBuffers(3 * MaximumOversampling);

        int mode = LowpassMode;
        int ftype = getFilterTypeForMode(mode);

        fMode = mode;
        fOversampling = 1;

        for (unsigned i = 0; i < 3; ++i)
            fFilterGains[i] = 1.0;

        const double cutoffs[] = {300.0, 1800.0, 7600.0};
        const double q = 10.0;

        for (unsigned i = 0; i < 3; ++i) {
            VAStateVariableFilter &filter = fFilters[i];
            filter.setSampleRate(samplerate);
            filter.setFilterType(ftype);
            filter.setCutoffFreq(fFilterCutoffFreqs[i] = cutoffs[i]);
            filter.setQ(fFilterQ[i] = q);
        }
    }

    void Rezonateur::setFilterMode(int mode)
    {
        int ftype = getFilterTypeForMode(mode);

        fMode = mode;

        for (unsigned i = 0; i < 3; ++i) {
            VAStateVariableFilter &filter = fFilters[i];
            filter.setFilterType(ftype);
            filter.clear();
        }
    }

    void Rezonateur::setFilterGain(unsigned nth, float gain)
    {
        assert(nth < 3);
        fFilterGains[nth] = gain;
    }

    void Rezonateur::setFilterCutoff(unsigned nth, float cutoff)
    {
        assert(nth < 3);
        fFilters[nth].setCutoffFreq((fFilterCutoffFreqs[nth] = cutoff) / fOversampling);
    }

    void Rezonateur::setFilterEmph(unsigned nth, float emph)
    {
        assert(nth < 3);
        fFilters[nth].setQ(fFilterQ[nth] = emph);
    }

    int Rezonateur::getFilterMode() const
    {
        return fMode;
    }

    float Rezonateur::getFilterGain(unsigned nth) const
    {
        assert(nth < 3);
        return fFilterGains[nth];
    }

    float Rezonateur::getFilterCutoff(unsigned nth) const
    {
        assert(nth < 3);
        return fFilterCutoffFreqs[nth];
    }

    float Rezonateur::getFilterEmph(unsigned nth) const
    {
        assert(nth < 3);
        return fFilterQ[nth];
    }

    unsigned Rezonateur::getOversampling() const
    {
        return fOversampling;
    }

    void Rezonateur::setOversampling(unsigned oversampling)
    {
        switch (oversampling) {
        default:
            assert(false);
            /* fall through */
        case 1:
            oversampling = 1;
            if (fOversampling == oversampling)
                return;
            break;
        case 2:
            oversampling = 2;
            if (fOversampling == oversampling)
                return;
            fOversampler2x.reset();
            break;
        case 4:
            oversampling = 4;
            if (fOversampling == oversampling)
                return;
            fOversampler4x.reset();
            break;
        case 8:
            oversampling = 8;
            if (fOversampling == oversampling)
                return;
            fOversampler8x.reset();
            break;
        }

        fOversampling = oversampling;

        for (unsigned b = 0; b < 3; ++b) {
            VAStateVariableFilter &filter = fFilters[b];
            filter.setCutoffFreq(fFilterCutoffFreqs[b] / oversampling);
            filter.clear();
        }
    }

    void Rezonateur::process(const float *input, float *output, unsigned count)
    {
        switch (fOversampling) {
        default:
            assert(false);
            /* fall through */
        case 1:
            DSP::NoOversampler noOversampler;
            processOversampled(noOversampler, input, output, count);
            break;
        case 2:
            processOversampled(fOversampler2x, input, output, count);
            break;
        case 4:
            processOversampled(fOversampler4x, input, output, count);
            break;
        case 8:
            processOversampled(fOversampler8x, input, output, count);
            break;
        }
    }

    template <class Oversampler> void Rezonateur::processOversampled(Oversampler &oversampler, const float *input, float *output, unsigned count)
    {
        while (count > 0) {
            unsigned current = (count < sBufferLimit) ? count : sBufferLimit;
            processWithinBufferLimit(oversampler, input, output, current);
            input += current;
            output += current;
            count -= current;
        }
    }

    template <class Oversampler> void Rezonateur::processWithinBufferLimit(Oversampler &oversampler, const float *input, float *output, unsigned count)
    {
        constexpr unsigned ratio = Oversampler::Ratio;

        float filterGains[3];
        getEffectiveFilterGains(filterGains);

        float *accum = getWorkBuffer(0 * MaximumOversampling);
        float *filterOutput = getWorkBuffer(1 * MaximumOversampling);

        ///
        if (ratio > 1) {
            float *filterInput = getWorkBuffer(2 * MaximumOversampling);
            for (unsigned i = 0; i < count; ++i) {
                filterInput[i * ratio] = oversampler.upsample(input[i]);
                for (unsigned o = 1; o < ratio; ++o)
                    filterInput[i * ratio + o] = oversampler.uppad(o);
            }
            input = filterInput;
        }

        ///
        fFilters[0].process(filterGains[0], input, accum, count * ratio);
        for (unsigned b = 1; b < 3; ++b) {
            fFilters[b].process(filterGains[b], input, filterOutput, count * ratio);
            for (unsigned i = 0; i < count * ratio; ++i)
                accum[i] += filterOutput[i];
        }

        ///
        for (unsigned i = 0; i < count; ++i) {
            output[i] = oversampler.downsample(accum[i * ratio]);
            for (unsigned o = 1; o < ratio; ++o)
                oversampler.downstore(accum[i * ratio + o]);
        }
    }

    void Rezonateur::getEffectiveFilterGains(float gains[3]) const
    {
        // must invert the middle filter in bandpass mode
        gains[0] = fFilterGains[0];
        gains[1] = (fMode != BandpassMode) ? fFilterGains[1] : -fFilterGains[1];
        gains[2] = fFilterGains[2];
    }

    int Rezonateur::getFilterTypeForMode(int mode)
    {
        switch (mode) {
        default:
            assert(false);
            /* fall through */
        case LowpassMode:
            return SVFLowpass;
        case BandpassMode:
        case BandpassNotchMode:
            return SVFBandpass;
        case HighpassMode:
            return SVFHighpass;
        }
    }

    double Rezonateur::getResponseGain(double f) const
    {
        std::complex<double> h = 0;
        float filterGains[3];
        getEffectiveFilterGains(filterGains);

        for (unsigned i = 0; i < 3; ++i) {
            double g = filterGains[i];
            h += g * fFilters[i].calcTransfer(f);
        }

        return std::abs(h);
    }

    void Rezonateur::allocateWorkBuffers(unsigned count)
    {
        fWorkBuffers.reset(new float[count * sBufferLimit]);
        fNumWorkBuffers = count;
    }

    float *Rezonateur::getWorkBuffer(unsigned index)
    {
        assert(index < fNumWorkBuffers);
        return &fWorkBuffers[index * sBufferLimit];
    }

    constexpr unsigned Rezonateur::sBufferLimit;
}