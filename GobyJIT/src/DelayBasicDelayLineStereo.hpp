/*
  ==============================================================================

    BasicDelayLine.h
    New version of BasicDelayLine that also accepts an external feedback path
    as input, and which uses fractional delay
    
    Created: 14 Sep 2014 5:15:53pm
    Author:  Keith Hearne

  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "Helper.h"
#include "DelayBasicDelayLine.hpp"

namespace Delay
{

    struct StereoBasicDelayLine : public StereoFXProcessor
    {
        BasicDelayLine delay[2];

        StereoBasicDelayLine(const int sr = 44100, const float d_ms = 0.0f, const float feedback = 0.0f, const float mixLevel = 0.5f)
        : delay[0](sr,d_ms,feedback,mixLevel),
          delay[1](sr,d_ms,feedback,mixLevel)
          {

          }
          
        double Tick(double iL, double iR, double &L, double &R, double A=1, double X=1, double Y=1)
        {
            L = delay[0].Tick(iL,A,X,Y);
            R = delay[1].Tick(iR,A,X,Y);
            return 0.5*(L+R);

        }
        void ProcessBlock(size_t n, float ** in, float ** out) {
            delay[0].ProcessBlock(n,in[0],out[0]);
            delay[1].ProcessBlock(n,in[1],out[1]);
        }
    };
}