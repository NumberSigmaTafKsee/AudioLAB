// https://github.com/StudioSixPlusOne/rack-modules/blob/master/src/dsp/HardLimiter.h
/*
 * Copyright (c) 2020 Dave French <contact/dot/dave/dot/french3/at/googlemail/dot/com>
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program (see COPYING); if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <utility>

#include "simd/vector.hpp"
#include "simd/functions.hpp"
#include "simd/sse_mathfun.h"
#include "simd/sse_mathfun_extension.h"
#include "digital.hpp"
#include "LookupTable.h"

using namespace rack;

using float_4 = ::rack::simd::float_4;

namespace sspo
{
    ///
    /// !A fixed parameter limiter
    /// Based on description given in Prikle, Designing Audio Effects 2nd
    ///
    struct Compressor
    {
        Compressor()
        {
            divider.setDivision (divFreq);
        }

        void calcCoeffs()
        {
            attackCoeff = simd::exp (TC / (sampleRate * attackTime));
            releaseCoeff = simd::exp (TC / (sampleRate * releaseTimes));
        }

        void setSampleRate (const DspFloatType sr)
        {
            sampleRate = sr / divFreq;
            calcCoeffs();
        }

        void setTimes (const DspFloatType attack, const DspFloatType release)
        {
            attackTime = attack;
            releaseTimes = release;
            calcCoeffs();
        }
        DspFloatType G = 0.0;
        DspFloatType process (const DspFloatType in)
        {
            if (divider.process())
            {
                //envelope follower
                auto rectIn = std::abs (in);
                currentEnv = rectIn > lastEnv
                                 ? attackCoeff * (lastEnv - rectIn) + rectIn
                                 : releaseCoeff * (lastEnv - rectIn) + rectIn;
                currentEnv = std::max (currentEnv, 0.00000000001f);
                lastEnv = currentEnv;

                //            auto dn = 20.0f * lookup.log10 (currentEnv);
                auto dn = 20.0f * simd::log10 (currentEnv);
                //Hard knee compression
                auto yndB = dn <= threshold ? dn : threshold + ((dn - threshold) / ratio);
                auto gndB = yndB - dn;
                G = simd::pow (10.0f, gndB / 20.0f);
            }

            return in * G;
        }

        DspFloatType attackTime{ 0.0001f };
        DspFloatType releaseTimes{ 0.025f };
        DspFloatType ratio{ 10.5f };
        DspFloatType threshold{ -0.0f }; //dB

    private:
        DspFloatType attackCoeff{ 0.0f };
        DspFloatType releaseCoeff{ 0.0f };
        DspFloatType lastEnv{ 0.0f };
        DspFloatType currentEnv{ 0.0f };
        DspFloatType sampleRate{ 1.0f };
        static constexpr int divFreq = 4;
        dsp::ClockDivider divider;

        static constexpr DspFloatType TC{ -0.9996723408f }; // { std::log (0.368f); } //capacitor discharge to 36.8%
    };

    inline DspFloatType saturate (DspFloatType in, DspFloatType max = 1.0f, DspFloatType kneeWidth = 0.05)
    {
        auto ret = 0.0f;
        if (std::abs (in) < (max - kneeWidth))
        {
            ret = in;
        }
        else
        {
            if (std::abs (in) < max)
            {
                ret = in > 0.0f
                          ? in - ((simd::pow (in - max + (kneeWidth / 2.0), 2.0)) / (2.0 * kneeWidth))
                          : in + ((simd::pow (in + max - (kneeWidth / 2.0), 2.0)) / (2.0 * kneeWidth));
                //ret = rack::math::clamp (ret, -max, max);
            }
            else
                ret = in > 0.0f
                          ? max
                          : -max;
        }

        return ret;
    }

    inline DspFloatType voltageSaturate (DspFloatType in)
    {
        return saturate (in, 11.7f, 0.5f);
    }

    inline float_4 voltageSaturate (float_4 in)
    {
        float_4 ret;
        for (auto i = 0; i < 4; ++i)
            ret[i] = voltageSaturate (in[i]);

        return ret;
    }

    struct Saturator
    {
        DspFloatType max = 1.0f;
        DspFloatType kneeWidth = 0.05f;

        Saturator() = default;

        Saturator (const DspFloatType limit, const DspFloatType knee)
        {
            max = limit;
            kneeWidth = knee;
        }

        DspFloatType process (DspFloatType in) const
        {
            return saturate (in, max, kneeWidth);
        }
    };