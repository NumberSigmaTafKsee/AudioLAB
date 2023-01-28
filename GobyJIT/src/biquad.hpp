#pragma once


#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "Vc/vectorclass.h"
#include "Vc/vectormath_lib.h"

namespace Filters::biquad
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
        
        void calculateCoeffs();
        
    public:
        Biquad() : FilterProcessor() {};
        ~Biquad(){};
        void setParams(const Parameters& params);
        Parameters getParams();
        double process(double x);  

        double Tick(double I, double A=1, double X=0, double Y=0)  
        {
            return A*process(I);
        }
    };


    void Biquad::setParams(const Parameters& params)
    {
        mparams = params;
        calculateCoeffs();
    }

    Parameters Biquad::getParams()
    {
        return mparams;
    }

    void Biquad::calculateCoeffs()
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
        if(ma0 != 0.0f)
        {
            mb0 /= ma0;
            mb1 /= ma0;
            mb2 /= ma0;
            ma1 /= ma0;
            ma2 /= ma0;
        }
    }

    double Biquad::process(double x)
    {
        Undenormal denormals;
        double y = mb0 * x + mb1 * mx_z1 + mb2 * mx_z2 - ma1 * my_z1 - ma2 * my_z2;

        mx_z2 = mx_z1;
        mx_z1 = x;

        my_z2 = my_z1;
        my_z1 = y;

        return y;
    }


    template<typename SIMD>
    class VecBiquad
    {
    private:
        FilterType mfilterType;
        
        Parameters mparams;

        // coefficients
        double ma0, ma1, ma2, mb0, mb1, mb2;

        // buffers
        SIMD mx_z1, mx_z2, my_z1, my_z2;
        
        void calculateCoeffs();
        
    public:
        VecBiquad(){};
        ~VecBiquad(){};
        void setParams(const Parameters& params);
        Parameters getParams();
        SIMD Tick(SIMD x);    
        void  ProcessBuffer(size_t n, float * input, float * output);
    };

    template<typename T>
    void VecBiquad<T>::setParams(const Parameters& params)
    {
        mparams = params;
        calculateCoeffs();
    }

    template<typename T>
    Parameters VecBiquad<T>::getParams()
    {
        return mparams;
    }

    template<typename T>
    void VecBiquad<T>::calculateCoeffs()
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
        if(ma0 != 0.0f)
        {
            mb0 /= ma0;
            mb1 /= ma0;
            mb2 /= ma0;
            ma1 /= ma0;
            ma2 /= ma0;
        }
    }

    template<typename SIMD>
    SIMD VecBiquad<SIMD>::Tick(SIMD x)
    {
        Undenormal denormals;
        SIMD y = mb0 * x + mb1 * mx_z1 + mb2 * mx_z2 - ma1 * my_z1 - ma2 * my_z2;

        mx_z2 = mx_z1;
        mx_z1 = x;

        my_z2 = my_z1;
        my_z1 = y;

        return y;
    }


    template<typename SIMD>
    void VecBiquad<SIMD>::ProcessBuffer(size_t n, float * in, float * out)
    {
        SIMD r,rout;
        for(size_t i = 0; i < n; i += 8)
        {
            r.load(in+i);
            rout = Tick(r);
            rout.store(out+i);
        }
    }

    using VecBiquad4 = VecBiquad<Vec4f>;
    using VecBiquad8 = VecBiquad<Vec8f>;

}