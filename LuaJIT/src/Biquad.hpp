#pragma once
#include "Undenormal.hpp"
#include "SoundObject.hpp"

///////////////////////////////////////////////////
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

namespace Filters::Biquad
{
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
        double fs;
        double f0;
        double Q;
        double dBGain;
    };

    class Biquad : public FilterProcessor
    {
    private:
        FilterType mfilterType;
        
        Parameters mparams;

        // coefficients
        double ma0, ma1, ma2, mb0, mb1, mb2;

        // buffers
        double mx_z1, mx_z2, my_z1, my_z2;
        // prev x,y delayline
        // delayline x_prev;
        // delayline y_prev;
        
        void calculateCoeffs();
        
    public:
        Biquad(){};
        ~Biquad(){};
        void setParams(const Parameters& params);
        Parameters getParams();
        double process(double x);

        double Tick(double I, double A=1, double X=0, double Y=0) {
            return process(I);
        }
        void morph(Biquad & other, float f = 0.5)
        {
            ma0 = ma0 + f*(other.ma0 - ma0);
            ma1 = ma1 + f*(other.ma1 - ma1);
            ma2 = ma2 + f*(other.ma2 - ma2);
            mb0 = mb0 + f*(other.mb0 - mb0);
            mb1 = mb1 + f*(other.mb1 - mb1);
            mb2 = mb2 + f*(other.mb2 - mb2);
            m_param.f0 = m_param.f0 + f*(other.m_param.f0 - m_param.f0);
            m_param.Q = m_param.Q + f*(other.m_param.Q - m_param.Q);
            m_param.dbGain = m_param.dbGain + f*(other.m_param.dbGain - m_param.dbGain);
        }
    };



    inline void Biquad::setParams(const Parameters& params)
    {
        mparams = params;
        calculateCoeffs();
    }

    inline Parameters Biquad::getParams()
    {
        return mparams;
    }

    inline void Biquad::calculateCoeffs()
    {
        double omega0 = 2.0f * M_PI * (mparams.f0 / mparams.fs);
        double alpha = std::sin(omega0) / (2.0 * mparams.Q);
        double A = std::pow(10, mparams.dBGain / 40.0);
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

    inline double Biquad::process(double x)
    {    
        Undenormal denormal;
        double y = (mb0 / ma0) * x + (mb1 / ma0) * mx_z1 + (mb2 / ma0) * mx_z2 - (ma1 / ma0) * my_z1 - (ma2 / ma0) * my_z2;

        mx_z2 = mx_z1;
        mx_z1 = x;

        my_z2 = my_z1;
        my_z1 = y;

        return y;
    }

}
