#pragma once
#include "Undenormal.hpp"
#include "SoundObject.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

namespace Filters
{

    //////////////////////////////////////////////////////////////////////////////////////////
    // RBJ
    //////////////////////////////////////////////////////////////////////////////////////////

    enum class FilterType
    {
        LowPass = 1,
        HighPass,
        BandPass,
        Notch,
        AllPass,
        Peaking,
        LowShelf,
        HighShelf,
        OnePoleZeroLP,
        OnePoleZeroHP,
    };

    struct Parameters
    {
        FilterType filterType;
        DspFloatType fs;
        DspFloatType f0;
        DspFloatType Q;
        DspFloatType dBGain;
    };

    class RBJBiquadFilter : public FilterProcessor
    {
    private:
        FilterType mfilterType;

        Parameters mparams;

        // coefficients
        DspFloatType ma0, ma1, ma2, mb0, mb1, mb2;

        // buffers
        DspFloatType mx_z1, mx_z2, my_z1, my_z2;
        // prev x,y delayline
        // delayline x_prev;
        // delayline y_prev;

        void calculateCoeffs();

    public:
        RBJBiquadFilter(){};
        ~RBJBiquadFilter(){};
        void setParams(const Parameters &params);
        Parameters getParams();
        DspFloatType process(DspFloatType x);
        void setFilterType(FilterType type)
        {
            mparams.filterType = type;
            setParams(mparams);
        }
        void setCutoff(DspFloatType f)
        {
            mparams.f0 = f;
            setParams(mparams);
        }
        void setQ(DspFloatType q)
        {
            mparams.Q = q;
            setParams(mparams);
        }
        void setGain(DspFloatType g)
        {
            mparams.dBGain = g;
            setParams(mparams);
        }
        enum
        {
            PORT_TYPE,
            PORT_CUTOFF,
            PORT_Q,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v)
        {
            switch (port)
            {
            case PORT_TYPE:
                setFilterType((FilterType)v);
                break;
            case PORT_CUTOFF:
                setCutoff(v);
                break;
            case PORT_Q:
                setQ(v);
                break;
            case PORT_GAIN:
                setGain(v);
                break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
        {
            return process(I);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            Undenormal denormal;
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType x = in[i];
                DspFloatType y = (mb0 / ma0) * x + (mb1 / ma0) * mx_z1 + (mb2 / ma0) * mx_z2 - (ma1 / ma0) * my_z1 - (ma2 / ma0) * my_z2;
                mx_z2 = mx_z1;
                mx_z1 = x;
                my_z2 = my_z1;
                my_z1 = y;
                out[i] = y;
            }            
        }
        void morph(RBJBiquadFilter &other, DspFloatType f = 0.5)
        {
            ma0 = ma0 + f * (other.ma0 - ma0);
            ma1 = ma1 + f * (other.ma1 - ma1);
            ma2 = ma2 + f * (other.ma2 - ma2);
            mb0 = mb0 + f * (other.mb0 - mb0);
            mb1 = mb1 + f * (other.mb1 - mb1);
            mb2 = mb2 + f * (other.mb2 - mb2);
            mparams.f0 = mparams.f0 + f * (other.mparams.f0 - mparams.f0);
            mparams.Q = mparams.Q + f * (other.mparams.Q - mparams.Q);
            mparams.dBGain = mparams.dBGain + f * (other.mparams.dBGain - mparams.dBGain);
        }
    };

    inline void RBJBiquadFilter::setParams(const Parameters &params)
    {
        mparams = params;
        calculateCoeffs();
    }

    inline Parameters RBJBiquadFilter::getParams()
    {
        return mparams;
    }

    inline void RBJBiquadFilter::calculateCoeffs()
    {
        DspFloatType omega0 = 2.0f * M_PI * (mparams.f0 / mparams.fs);
        DspFloatType alpha = std::sin(omega0) / (2.0 * mparams.Q);
        DspFloatType A = std::pow(10, mparams.dBGain / 40.0);
        switch (mparams.filterType)
        {
        case FilterType::LowPass:
        {
            ma0 = 1.0 + alpha;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha;
            mb0 = (1.0 - std::cos(omega0)) / 2.0;
            mb1 = 1.0 - std::cos(omega0);
            mb2 = (1.0 - std::cos(omega0)) / 2.0;
            break;
        }
        case FilterType::HighPass:
        {
            ma0 = 1.0 + alpha;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha;
            mb0 = (1.0 + std::cos(omega0)) / 2.0;
            mb1 = -(1.0 + std::cos(omega0));
            mb2 = (1.0 + std::cos(omega0)) / 2.0;
            break;
        }
        case FilterType::BandPass:
        {
            ma0 = 1.0 + alpha;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha;
            mb0 = alpha;
            mb1 = 0;
            mb2 = -alpha;
            break;
        }
        case FilterType::Notch:
        {
            ma0 = 1.0 + alpha;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha;
            mb0 = 1.0;
            mb1 = -2.0 * std::cos(omega0);
            mb2 = 1.0;
            break;
        }
        case FilterType::AllPass:
        {
            ma0 = 1.0 + alpha;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha;
            mb0 = 1.0 - alpha;
            mb1 = -2.0 * std::cos(omega0);
            mb2 = 1.0 + alpha;
            break;
        }
        case FilterType::Peaking:
        {
            ma0 = 1.0 + alpha / A;
            ma1 = -2.0 * std::cos(omega0);
            ma2 = 1.0 - alpha / A;
            mb0 = 1.0 + alpha * A;
            mb1 = -2.0 * std::cos(omega0);
            mb2 = 1.0 - alpha * A;
            break;
        }
        case FilterType::LowShelf:
        {
            ma0 = (A + 1.0) + (A - 1.0) * std::cos(omega0) + 2.0 * std::sqrt(A) * alpha;
            ma1 = -2.0 * ((A - 1.0) + (A + 1.0) * std::cos(omega0));
            ma2 = (A + 1.0) + (A - 1.0) * std::cos(omega0) - 2.0 * std::sqrt(A) * alpha;
            mb0 = A * ((A + 1.0) - (A - 1.0) * std::cos(omega0) + 2.0 * std::sqrt(A) * alpha);
            mb1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * std::cos(omega0));
            mb2 = A * ((A + 1.0) - (A - 1.0) * std::cos(omega0) - 2.0 * std::sqrt(A) * alpha);
            break;
        }
        case FilterType::HighShelf:
        {
            ma0 = (A + 1.0) - (A - 1.0) * std::cos(omega0) + 2.0 * std::sqrt(A) * alpha;
            ma1 = 2.0 * ((A - 1.0) - (A + 1.0) * std::cos(omega0));
            ma2 = (A + 1.0) - (A - 1.0) * std::cos(omega0) - 2.0 * std::sqrt(A) * alpha;
            mb0 = A * ((A + 1.0) + (A - 1.0) * std::cos(omega0) + 2.0 * std::sqrt(A) * alpha);
            mb1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * std::cos(omega0));
            mb2 = A * ((A + 1.0) + (A - 1.0) * std::cos(omega0) - 2.0 * std::sqrt(A) * alpha);
            break;
        }
        default:
            break;
        }
    }

    inline DspFloatType RBJBiquadFilter::process(DspFloatType x)
    {
        Undenormal denormal;
        DspFloatType y = (mb0 / ma0) * x + (mb1 / ma0) * mx_z1 + (mb2 / ma0) * mx_z2 - (ma1 / ma0) * my_z1 - (ma2 / ma0) * my_z2;

        mx_z2 = mx_z1;
        mx_z1 = x;

        my_z2 = my_z1;
        my_z1 = y;

        return y;
    }

}