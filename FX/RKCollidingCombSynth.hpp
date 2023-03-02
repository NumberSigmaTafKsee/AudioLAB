// need to remove all the plugin kaka
// (c) 2020 Takamitsu Endo
//
// This file is part of CollidingCombSynth.
//
// CollidingCombSynth is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CollidingCombSynth is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CollidingCombSynth.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include "RKConstants.h"
#include "RKSmoothFilter.hpp"
//#include "../../common/dsp/somemath.hpp"
//#include "../parameter.hpp"
#include "RKDelay.hpp"
#include "RKEnvelope.hpp"
#include "RKiir.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <random>
#include <cstddef>
#include <cstdint>
#include <string>

enum class NoteState { active, release, rest };

constexpr uint16_t nDelay = 24;
constexpr uint16_t nComb = 8;

#define NOTE_PROCESS_INFO_SMOOTHER(METHOD)                                               \
  lowpassCutoff.METHOD(pv[ID::lowpassCutoff]->getFloat());                               \
  highpassCutoff.METHOD(pv[ID::highpassCutoff]->getFloat());                             \
  noiseGain.METHOD(pv[ID::exciterGain]->getFloat());                                     \
  propagation.METHOD(pv[ID::propagation]->getFloat());

struct NoteProcessInfo {
  std::minstd_rand rngNoise{0};
  std::minstd_rand rngComb{0};
  std::minstd_rand rngString{0};
  std::minstd_rand rngUnison{0};

  ExpSmoother<float> lowpassCutoff;
  ExpSmoother<float> highpassCutoff;
  ExpSmoother<float> noiseGain;
  ExpSmoother<float> propagation;

  void reset(/*GlobalParameter &param*/)
  {
    /*
    using ID = ParameterID::ID;
    auto &pv = param.value;    
    rngNoise.seed(pv[ID::seedNoise]->getInt());
    rngComb.seed(pv[ID::seedComb]->getInt());
    rngString.seed(pv[ID::seedString]->getInt());
    rngUnison.seed(pv[ID::seedUnison]->getInt());
    NOTE_PROCESS_INFO_SMOOTHER(reset);
    */    
  }

  void setParameters( /*GlobalParameter &param */)
  {
    /*
    using ID = ParameterID::ID;
    auto &pv = param.value;
    NOTE_PROCESS_INFO_SMOOTHER(push);
    */    
  }

  void process()
  {
    lowpassCutoff.process();
    highpassCutoff.process();
    noiseGain.process();
  }
};

#define NOTE_CLASS(INSTRSET)                                                             \
  class Note_##INSTRSET {                                                                \
  public:                                                                                \
    NoteState state = NoteState::rest;                                                   \
                                                                                         \
    int32_t id = -1;                                                                     \
    float velocity = 0;                                                                  \
    float noteFreq = 1;                                                                  \
    float pan = 0.5f;                                                                    \
    float gain = 0;                                                                      \
                                                                                         \
    bool isCompressorOn = true;                                                          \
    int32_t releaseCounter = 0;                                                          \
    float releaseLength = 0;                                                             \
                                                                                         \
    ADNoise noise;                                                                       \
    EMAFilter<float> exciterLowpass;                                                     \
    AttackGate<float> gate;                                                              \
    std::array<ShortComb<float>, nComb> comb;                                            \
    KsHat<float, nDelay> cymbal;                                                         \
    ExpADSREnvelopeP<float> cymbalLowpassEnvelope;                                       \
    DCKiller<float> dcKiller;                                                            \
    EasyCompressor<float> compressor;                                                    \
                                                                                         \
    void setup(float sampleRate);                                                        \        
    void release(float sampleRate);                                                      \
    void rest();                                                                         \
    bool isAttacking();                                                                  \
    float getGain();                                                                     \
    std::array<float, 2> process(float sampleRate, NoteProcessInfo &info);               \
  };

#ifdef USE_VECTORCLASS
NOTE_CLASS(AVX512)
NOTE_CLASS(AVX2)
NOTE_CLASS(AVX)
#else
NOTE_CLASS(Plain)
#endif


inline float calcMasterPitch(
  int32_t octave, int32_t semi, int32_t milli, float bend, float equalTemperament)
{
  return equalTemperament * octave + semi + milli / 1000.0f + (bend - 0.5f) * 4.0f;
}

inline float
notePitchToFrequency(float notePitch, float equalTemperament = 12.0f, float a4Hz = 440.0f)
{
  return a4Hz * powf(2.0f, (notePitch - 69.0f) / equalTemperament);
}

inline float calcNotePitch(float notePitch, float equalTemperament = 12.0f)
{
  return powf(2.0f, (notePitch - 69.0f) / equalTemperament);
}

// Fast approximation of EMAFilter::cutoffToP().
// x is normalized frequency. Range is [0, 0.5).
template<typename T> inline T cutoffToPApprox(T x)
{
  return (T(-0.0004930424721520979) + T(2.9650003823571467) * x
          + T(1.8250079928630534) * x * x)
    / (T(0.4649282844668299) + T(1.8754712250208077) * x + T(3.7307819672604023) * x * x)
    + T(0.0010604699447733293);
}

// Fast approximation of OnePoleHighpass::setCutoff().
// x is normalized frequency. Range is [0, 0.5).
template<typename T> T onepoleHighpassPoleApprox(T x)
{
  return (T(5.476984559402437) + T(-13.572160512877103) * x
          + T(9.553503352240815) * x * x)
    / (T(5.479174921068828) + T(20.63587495820437) * x + T(36.02138404173517) * x * x);
}

void NOTE_NAME::setup(float sampleRate)
{
  cymbalLowpassEnvelope.setup(sampleRate);
  cymbal.setup(sampleRate);
  for (auto &cmb : comb) cmb.reset();
}

void NOTE_NAME::noteOn(
  int32_t noteId,
  float notePitch,
  float velocity,
  float pan,
  float sampleRate,
  NoteProcessInfo &info,
  GlobalParameter &param)
{
  using ID = ParameterID::ID;
  auto &pv = param.value;

  state = NoteState::active;
  id = noteId;

  this->velocity = velocity;
  this->pan = pan;
  gain = 1.0f;

  const float eqTemp = pv[ID::equalTemperament]->getFloat() + 1;
  const auto semitone = int32_t(pv[ID::semitone]->getInt()) - 120;
  const auto octave = eqTemp * (int32_t(pv[ID::octave]->getInt()) - 12);
  const auto milli = 0.001f * (int32_t(pv[ID::milli]->getInt()) - 1000);
  const float a4Hz = pv[ID::pitchA4Hz]->getFloat() + 100;
  const auto pitch = calcNotePitch(octave + semitone + milli + notePitch, eqTemp);
  const auto frequency = a4Hz * pitch;

  noise.reset(
    sampleRate, pv[ID::exciterAttack]->getFloat(), pv[ID::exciterDecay]->getFloat(),
    frequency, pv[ID::exciterNoiseMix]->getFloat());
  exciterLowpass.reset();
  exciterLowpass.setCutoff(sampleRate, pv[ID::exciterLowpassCutoff]->getFloat());
  gate.reset(sampleRate, pv[ID::exciterAttack]->getFloat());

  releaseLength = 0.01f * sampleRate;
  releaseCounter = int32_t(releaseLength);

  for (size_t idx = 0; idx < nComb; ++idx) {
    const auto combTime = pv[ID::combTime0 + idx]->getFloat();
    const auto spread = combTime * pv[ID::randomComb]->getFloat();
    std::uniform_real_distribution<float> distCombTime(
      combTime - spread, combTime + spread);
    comb[idx].setTime(sampleRate, distCombTime(info.rngComb));
  }

  for (size_t idx = 0; idx < nDelay; ++idx) {
    const auto freq = pitch * pv[ID::frequency0 + idx]->getFloat();
    const auto spread
      = (freq - float(Scales::frequency.getMin())) * pv[ID::randomFrequency]->getFloat();

    auto distLower = freq - spread;
    auto distUpper = freq + spread;
    if (distLower > distUpper) std::swap(distLower, distUpper);
    std::uniform_real_distribution<float> distFreq(0.0f, 1.0f);
    auto freqValue = distLower + (distUpper - distLower) * distFreq(info.rngString);
    cymbal.string[idx].delay.setTime(sampleRate, 1.0f / freqValue);
  }
  cymbal.trigger(pv[ID::distance]->getFloat(), pv[ID::connection]->getInt());

  cymbalLowpassEnvelope.prepare(
    sampleRate, pv[ID::lowpassA]->getFloat(), pv[ID::lowpassD]->getFloat(),
    pv[ID::lowpassS]->getFloat(), pv[ID::lowpassR]->getFloat());

  dcKiller.reset();

  isCompressorOn = pv[ID::compressor]->getInt();
  compressor.prepare(
    sampleRate, pv[ID::compressorTime]->getFloat(),
    pv[ID::compressorThreshold]->getFloat());
  compressor.reset();
}

void NOTE_NAME::release(float sampleRate)
{
  if (state == NoteState::rest) return;
  state = NoteState::release;
  cymbalLowpassEnvelope.release(sampleRate);
}

void NOTE_NAME::rest() { state = NoteState::rest; }

bool NOTE_NAME::isAttacking() { return cymbalLowpassEnvelope.isAttacking(); }

float NOTE_NAME::getGain() { return gain; }

std::array<float, 2> NOTE_NAME::process(float sampleRate, NoteProcessInfo &info)
{
  float sig = noise.isTerminated
    ? 0
    : info.noiseGain.getValue() * exciterLowpass.process(noise.process(info.rngNoise));

  for (auto &cmb : comb) sig -= cmb.process(sig);
  sig *= gate.process();

  float lpEnv = cymbalLowpassEnvelope.process(sampleRate);
  gain = velocity * lpEnv; // Used to determin most quiet note.

  cymbal.kp = cutoffToPApprox(lpEnv * info.lowpassCutoff.getValue() / sampleRate);
  cymbal.b1 = onepoleHighpassPoleApprox(info.highpassCutoff.getValue() / sampleRate);
  sig = dcKiller.process(cymbal.process(sig, info.propagation.getValue()));

  if (isCompressorOn) sig = compressor.process(sig);

  if (cymbalLowpassEnvelope.isTerminated()) {
    --releaseCounter;
    sig *= releaseCounter / releaseLength;
    if (releaseCounter <= 0) state = NoteState::rest;
  }

  sig *= velocity;
  return {(1.0f - pan) * sig, pan * sig};
}

