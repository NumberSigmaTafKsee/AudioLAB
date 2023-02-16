#pragma once

#include <cmath>
#include "fxdsp.hpp"

namespace FXDSP
{
    template<typename T>
    T midiNoteToFrequency(unsigned note)
    {
        return std::pow(2.0, ((note - 69.0)/12.)) * 440.0;
    }

    template<typename T>
    unsigned frequencyToMidiNote(T f)
    {
        return (unsigned)(69 + (12 * std::log2(f / 440.0)));
    }
}