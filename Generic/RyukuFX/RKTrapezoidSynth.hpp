// (c) 2019-2020 Takamitsu Endo
//
// This file is part of TrapezoidSynth.
//
// TrapezoidSynth is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TrapezoidSynth is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TrapezoidSynth.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include "RKConstants.h"
#include "RKSmoothFilter.hpp"

#include "RKTrapezoidEnvelope.hpp"
#include "RKTrapezoidIIR.hpp"
#include "RKTrapezoidOsc.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <vector>


struct NoteInfo {
  int32_t id;
  float frequency;
};

template<typename Sample> class TpzMono {
public:
  const static int32_t rngPitchDriftSeed = 987654321;

  const std::array<Sample, 9> octaveTable{0.0625, 0.125, 0.25, 0.5, 1.0,
                                          2.0,    4.0,   8.0,  16.0};

  Sample feedbackBuffer = 0;
  Sample tpzOsc2Buffer = 0;
  PTRTrapezoidOsc tpzOsc1{44100, 0};
  PTRTrapezoidOsc tpzOsc2{44100, 0};
  Random<float> rngPitchDrift{rngPitchDriftSeed};
  DecimationLowpass<Sample> decimationLP;

  SerialZDF1Pole<Sample> filter;

  LFO<Sample> lfo;

  ADSREnvelope<
    Sample,
    TableCurve<Sample, EnvelopeCurveType::attack, 128>,
    TableCurve<Sample, EnvelopeCurveType::decay, 128>,
    TableCurve<Sample, EnvelopeCurveType::decay, 128>>
    gainEnvelope;

  ADSREnvelope<
    Sample,
    TableCurve<Sample, EnvelopeCurveType::attack, 128>,
    TableCurve<Sample, EnvelopeCurveType::decay, 128>,
    TableCurve<Sample, EnvelopeCurveType::decay, 128>>
    filterEnvelope;

  PolyExpEnvelope<double> modEnvelope1;
  PolyExpEnvelope<double> modEnvelope2;

  AMPitchShiter<Sample> shifter1;
  AMPitchShiter<Sample> shifter2;

  Sample noteFreq = 0;
  Sample normalizedKey = 0;

  LinearSmootherLocal<Sample> interpOctave;
  LinearSmootherLocal<Sample> interpOsc1Pitch;
  LinearSmootherLocal<Sample> interpOsc2Pitch;
  LinearSmoother<Sample> interpOsc1Slope;
  LinearSmoother<Sample> interpOsc1PulseWidth;
  LinearSmoother<Sample> interpOsc2Slope;
  LinearSmoother<Sample> interpOsc2PulseWidth;
  LinearSmoother<Sample> interpPhaseMod;
  LinearSmoother<Sample> interpOscMix;
  LinearSmoother<Sample> interpPitchDrift;
  LinearSmoother<Sample> interpFeedback;
  LinearSmoother<Sample> interpFilterCutoff;
  LinearSmoother<Sample> interpFilterFeedback;
  LinearSmoother<Sample> interpFilterSaturation;
  LinearSmoother<Sample> interpFilterEnvToCutoff;
  LinearSmoother<Sample> interpFilterKeyToCutoff;
  LinearSmoother<Sample> interpOscMixToFilterCutoff;
  LinearSmoother<Sample> interpMod1EnvToPhaseMod;
  LinearSmoother<Sample> interpMod2EnvToFeedback;
  LinearSmoother<Sample> interpMod2EnvToLFOFrequency;
  LinearSmoother<Sample> interpModEnv2ToOsc2Slope;
  LinearSmoother<Sample> interpMod2EnvToShifter1;
  LinearSmoother<Sample> interpLFOFrequency;
  LinearSmoother<Sample> interpLFOShape;
  LinearSmoother<Sample> interpLFOToPitch;
  LinearSmoother<Sample> interpLFOToSlope;
  LinearSmoother<Sample> interpLFOToPulseWidth;
  LinearSmoother<Sample> interpLFOToCutoff;
  LinearSmoother<Sample> interpShifter1Pitch;
  LinearSmoother<Sample> interpShifter1Gain;
  LinearSmoother<Sample> interpShifter2Pitch;
  LinearSmoother<Sample> interpShifter2Gain;

  void setup(Sample sampleRate);
  void reset(GlobalParameter &param);
  void startup();
  void setParameters(Sample tempo, GlobalParameter &param);
  void
  noteOn(bool wasResting, Sample frequency, Sample normalizedKey, GlobalParameter &param);
  void noteOff(Sample frequency);
  void release(bool resetPitch);
  Sample process(const size_t bufferSize);

private:
  Sample getOctave(GlobalParameter &param);
  Sample getOsc1Pitch(GlobalParameter &param);
  Sample getOsc2Pitch(GlobalParameter &param);
};



template<typename Sample> void TpzMono<Sample>::setup(Sample sampleRate)
{
  interpOctave.setSampleRate(sampleRate);
  interpOctave.setTime(Sample(0.001));
  interpOsc1Pitch.setSampleRate(sampleRate);
  interpOsc2Pitch.setSampleRate(sampleRate);
  tpzOsc1.sampleRate = 8 * sampleRate;
  tpzOsc2.sampleRate = 8 * sampleRate;
  lfo.setup(sampleRate);
  gainEnvelope.setup(sampleRate);
  filterEnvelope.setup(sampleRate);
  modEnvelope1.setup(sampleRate);
  modEnvelope2.setup(sampleRate);
  filter.setup(8 * sampleRate);
  shifter1.setup(sampleRate);
  shifter2.setup(sampleRate);
}

template<typename Sample> void TpzMono<Sample>::reset(GlobalParameter &param)
{
  decimationLP.reset();
  filter.reset();
  shifter1.reset();
  shifter2.reset();
  gainEnvelope.terminate();
  filterEnvelope.terminate();

  noteFreq = 0;
  normalizedKey = 0;

  ASSIGN_NOTE_PARAMETER(reset);

  interpLFOFrequency.reset(1.0f);

  startup();
}

template<typename Sample> void TpzMono<Sample>::startup()
{
  feedbackBuffer = 0;
  tpzOsc2Buffer = 0;
  tpzOsc1.reset();
  tpzOsc2.reset();
  rngPitchDrift.seed = rngPitchDriftSeed;
  lfo.reset();
}

template<typename Sample>
void TpzMono<Sample>::setParameters(Sample tempo, GlobalParameter &param)
{
  ASSIGN_NOTE_PARAMETER(push);

  float lfoFreq = param.value[ParameterID::lfoFrequency]->getFloat();
  if (param.value[ParameterID::lfoTempoSync]->getInt()) {
    lfoFreq = lfoFreq * tempo / 240.0f;
  }
  interpLFOFrequency.push(lfoFreq);

  filter.setOrder(param.value[ParameterID::filterOrder]->getInt());

  switch (param.value[ParameterID::lfoType]->getInt()) {
    default:
    case 0: // sin
      lfo.type = LFOType::sin;
      break;

    case 1: // saw
      lfo.type = LFOType::saw;
      break;

    case 2: // pulse
      lfo.type = LFOType::pulse;
      break;

    case 3: // noise
      lfo.type = LFOType::noise;
      break;
  }

  gainEnvelope.set(
    param.value[ParameterID::gainA]->getFloat(),
    param.value[ParameterID::gainD]->getFloat(),
    param.value[ParameterID::gainS]->getFloat(),
    param.value[ParameterID::gainR]->getFloat(),
    param.value[ParameterID::gainCurve]->getFloat());
  filterEnvelope.set(
    param.value[ParameterID::filterA]->getFloat(),
    param.value[ParameterID::filterD]->getFloat(),
    param.value[ParameterID::filterS]->getFloat(),
    param.value[ParameterID::filterR]->getFloat(),
    param.value[ParameterID::filterCurve]->getFloat());
}

template<typename Sample>
void TpzMono<Sample>::noteOn(
  bool wasResting, Sample frequency, Sample normalizedKey, GlobalParameter &param)
{
  noteFreq = frequency;

  int32_t pitchSlideType
    = wasResting ? param.value[ParameterID::pitchSlideType]->getInt() : -1;
  switch (pitchSlideType) {
    case 1: { // Sustain
      interpOctave.reset(getOctave(param));
      interpOsc1Pitch.reset(getOsc1Pitch(param));
      interpOsc2Pitch.reset(getOsc2Pitch(param));
    } break;

    case 2: // Reset to 0
      interpOsc1Pitch.reset(0);
      interpOsc2Pitch.reset(0);
      break;

    default:
      break;
  }

  this->normalizedKey = normalizedKey;
  if (
    param.value[ParameterID::gainEnvRetrigger]->getInt() || gainEnvelope.isTerminated()
    || gainEnvelope.isReleasing()) {
    gainEnvelope.reset(
      param.value[ParameterID::gainA]->getFloat(),
      param.value[ParameterID::gainD]->getFloat(),
      param.value[ParameterID::gainS]->getFloat(),
      param.value[ParameterID::gainR]->getFloat(),
      param.value[ParameterID::gainCurve]->getFloat());
  }

  if (
    param.value[ParameterID::filterEnvRetrigger]->getInt()
    || filterEnvelope.isTerminated() || filterEnvelope.isReleasing()) {
    filterEnvelope.reset(
      param.value[ParameterID::filterA]->getFloat(),
      param.value[ParameterID::filterD]->getFloat(),
      param.value[ParameterID::filterS]->getFloat(),
      param.value[ParameterID::filterR]->getFloat(),
      param.value[ParameterID::filterCurve]->getFloat());
  }

  if (param.value[ParameterID::modEnv1Retrigger]->getInt() || wasResting)
    modEnvelope1.reset(
      param.value[ParameterID::modEnv1Attack]->getFloat(),
      param.value[ParameterID::modEnv1Curve]->getFloat());

  if (param.value[ParameterID::modEnv2Retrigger]->getInt() || wasResting)
    modEnvelope2.reset(
      param.value[ParameterID::modEnv2Attack]->getFloat(),
      param.value[ParameterID::modEnv2Curve]->getFloat());
}

template<typename Sample> void TpzMono<Sample>::noteOff(Sample frequency)
{
  noteFreq = frequency;
}

template<typename Sample> void TpzMono<Sample>::release(bool resetPitch)
{
  if (resetPitch) noteFreq = 0;

  gainEnvelope.release();
  filterEnvelope.release();
}

template<typename Sample> Sample TpzMono<Sample>::process(const size_t bufferSize)
{
  if (gainEnvelope.isTerminated()) return 0;

  const auto modEnv2Sig = float(modEnvelope2.process());
  lfo.setFreq(
    interpLFOFrequency.process() + modEnv2Sig * interpMod2EnvToLFOFrequency.process());
  lfo.pw = interpLFOShape.process();
  const auto lfoSig = lfo.process();

  const auto filterEnv = filterEnvelope.process()
    + interpOscMixToFilterCutoff.process() * (1.0f + feedbackBuffer);
  const auto cutoff = interpFilterCutoff.process()
    + interpFilterKeyToCutoff.process() * noteFreq
    + Sample(19800)
      * (interpLFOToCutoff.process() * lfoSig + interpFilterEnvToCutoff.process() * filterEnv);
  filter.setCutoff(clamp(cutoff, Sample(20), Sample(8 * 20000)));
  filter.feedback = interpFilterFeedback.process();
  filter.saturation = interpFilterSaturation.process();

  const auto octave = interpOctave.process();
  tpzOsc1.setFreq(
    octave
    * interpOsc1Pitch.process()
    * (1.0f
      + lfoSig * interpLFOToPitch.process()
      + interpPitchDrift.process() * rngPitchDrift.process()));
  tpzOsc1.setSlope(interpOsc1Slope.process() + lfoSig * interpLFOToSlope.process());
  tpzOsc1.setPulseWidth(
    interpOsc1PulseWidth.process() + lfoSig * interpLFOToPulseWidth.process());

  const auto osc2Freq = octave * interpOsc2Pitch.process();
  tpzOsc2.setFreq(osc2Freq);
  tpzOsc2.setSlope(
    interpOsc2Slope.process() + +modEnv2Sig * interpModEnv2ToOsc2Slope.process());
  tpzOsc2.setPulseWidth(interpOsc2PulseWidth.process());

  const auto oscMix = interpOscMix.process();
  const auto modEnv1Sig = float(modEnvelope1.process());
  tpzOsc1.addPhase(
    feedbackBuffer
      * (interpMod2EnvToFeedback.process() * modEnv2Sig + interpFeedback.process())
    + tpzOsc2Buffer
      * (interpMod1EnvToPhaseMod.process() * modEnv1Sig + interpPhaseMod.process()));
  for (size_t i = 0; i < 8; ++i) {
    tpzOsc2Buffer = tpzOsc2.process();
    const float osc1Sig = tpzOsc1.process();
    feedbackBuffer = decimationLP.process(
      filter.process(osc1Sig + oscMix * (tpzOsc2Buffer - osc1Sig)));
  }

  shifter1.setShift(
    osc2Freq * interpShifter1Pitch.process()
    + modEnv2Sig * interpMod2EnvToShifter1.process());
  shifter2.setShift(osc2Freq * interpShifter2Pitch.process());

  return gainEnvelope.process()
    * (feedbackBuffer + interpShifter1Gain.process() * shifter1.process(feedbackBuffer)
       + interpShifter2Gain.process() * shifter2.process(feedbackBuffer));
}

template<typename Sample> Sample TpzMono<Sample>::getOctave(GlobalParameter &param)
{
  int32_t index = 4
    + int32_t(std::floor(
      param.value[ParameterID::octave]->getFloat()
      + filterEnvelope.getValue()
        * param.value[ParameterID::filterEnvToOctave]->getFloat()));
  if (index < 0) index = 0;
  if (size_t(index) >= octaveTable.size()) index = int32_t(octaveTable.size()) - 1;
  return octaveTable[index];
}

template<typename Sample> Sample TpzMono<Sample>::getOsc1Pitch(GlobalParameter &param)
{
  return noteFreq
    * paramToPitch(
           param.value[ParameterID::osc1Semi]->getFloat(),
           param.value[ParameterID::osc1Cent]->getFloat(),
           param.value[ParameterID::pitchBend]->getFloat());
}

template<typename Sample> Sample TpzMono<Sample>::getOsc2Pitch(GlobalParameter &param)
{
  return noteFreq * param.value[ParameterID::osc2Overtone]->getInt()
    * paramToPitch(
           param.value[ParameterID::osc2Semi]->getFloat(),
           param.value[ParameterID::osc2Cent]->getFloat(),
           param.value[ParameterID::pitchBend]->getFloat());
}
