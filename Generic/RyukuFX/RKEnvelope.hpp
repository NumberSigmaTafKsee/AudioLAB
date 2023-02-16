#pragma once

#include "RKConstants.h"
#include "RKSmoothFilter.hpp"

#include <cfloat>
#include <limits>
#include <random>



template<typename Sample> class ExpADSREnvelopeP {
public:
  void setup(Sample sampleRate) { tailLength = uint32_t(0.01 * sampleRate); }

  void prepare(
    Sample sampleRate,
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime)
  {
    state = State::attack;
    sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    atk = int32_t(sampleRate * attackTime);
    decTime = decayTime;
    relTime = releaseTime;
    pController.setCutoff(sampleRate, Sample(1) / attackTime);
  }

  void reset()
  {
    state = State::terminated;
    pController.reset(0);
  }

  void set(
    Sample sampleRate,
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime)
  {
    switch (state) {
      case State::attack:
        atk = int32_t(sampleRate * attackTime);
        // Fall through.

      case State::decay:
        decTime = decayTime;
        sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
        // Fall through.

      case State::release:
        relTime = releaseTime;

      default:
        break;
    }

    if (state == State::attack)
      pController.setCutoff(sampleRate, Sample(1) / attackTime);
    else if (state == State::decay)
      pController.setCutoff(sampleRate, Sample(1) / decayTime);
    else if (state == State::release)
      pController.setCutoff(sampleRate, Sample(1) / releaseTime);
  }

  void release(Sample sampleRate)
  {
    state = State::release;
    pController.setCutoff(sampleRate, Sample(1) / relTime);
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  Sample process(Sample sampleRate)
  {
    switch (state) {
      case State::attack: {
        value = pController.process(Sample(1));
        --atk;
        if (atk == 0) {
          state = State::decay;
          pController.setCutoff(sampleRate, Sample(1) / decTime);
        }
      } break;

      case State::decay:
        value = pController.process(sustain);
        break;

      case State::release:
        value = pController.process(0);
        if (value < threshold) {
          value = threshold;
          state = State::tail;
          tailCounter = tailLength;
        }
        break;

      case State::tail:
        --tailCounter;
        value = threshold * tailCounter / float(tailLength);
        if (tailCounter == 0) {
          state = State::terminated;
          pController.reset(0);
        } else {
          pController.reset(value);
        }
        break;

      default:
        return 0;
    }
    return value;
  }

private:
  enum class State : int32_t { attack, decay, release, tail, terminated };
  const Sample threshold = Sample(1e-5);

  uint32_t tailLength = 32;
  uint32_t tailCounter = tailLength;

  EMAFilter<Sample> pController;
  State state = State::terminated;
  uint32_t atk = 0;
  Sample decTime = 0;
  Sample relTime = 0;
  Sample sustain = 1;
  Sample value = 0;
};

template<typename Sample> inline Sample cosinterp(Sample t)
{
  return Sample(0.5) * (Sample(1) - std::cos(Sample(pi) * t));
}

template<typename Sample> class ExpDecayCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = Sample(1);
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha = std::pow(threshold, Sample(1) / (seconds * sampleRate));
  }

  bool isTerminated() { return value <= threshold; }

  Sample process()
  {
    if (value <= threshold) return Sample(0);
    value *= alpha;
    return value - threshold;
  }

protected:
  const Sample threshold = Sample(1e-5);
  Sample value = 0;
  Sample alpha = 0;
};

template<typename Sample> class ExpAttackCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = threshold;
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha = std::pow(Sample(1.0) / threshold, Sample(1.0) / (seconds * sampleRate));
  }

  bool isTerminated() { return value >= Sample(1); }

  Sample process()
  {
    value *= alpha;
    if (value >= Sample(1)) return Sample(1 - threshold);
    return value - threshold;
  }

protected:
  const Sample threshold = Sample(1e-5);
  Sample value = 0;
  Sample alpha = 0;
};

// 1 - ExpDecayCurve.process();
template<typename Sample> class NegativeExpAttackCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = Sample(1);
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha = std::pow(threshold, Sample(1) / (seconds * sampleRate));
  }

  bool isTerminated() { return value <= threshold; }

  Sample process()
  {
    if (value <= threshold) return Sample(1 - threshold);
    value *= alpha;
    return Sample(1 - threshold) - value;
  }

protected:
  const Sample threshold = Sample(1e-5);
  Sample value = 0;
  Sample alpha = 0;
};

template<typename Sample> class ExpADSREnvelope {
public:
  void setup(Sample sampleRate)
  {
    this->sampleRate = sampleRate;
    declickLength = int32_t(0.001 * sampleRate);
  }

  Sample adaptTime(Sample seconds, Sample noteFreq)
  {
    const Sample cycle = Sample(1) / noteFreq;
    return seconds < cycle ? cycle : seconds;
  }

  void reset(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq,
    Sample curve)
  {
    if (declickCounter >= declickLength || state == State::terminated) declickCounter = 0;
    state = State::attack;

    sustain = std::max<Sample>(Sample(0.0), std::min<Sample>(sustainLevel, Sample(1.0)));

    offset = value;
    range = Sample(1.0) - value;

    this->curve = curve;

    attackTime = adaptTime(attackTime, noteFreq);
    atk.reset(sampleRate, attackTime);
    atkNeg.reset(sampleRate, attackTime);
    dec.reset(sampleRate, decayTime);
    rel.reset(sampleRate, adaptTime(releaseTime, noteFreq));
  }

  void set(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    switch (state) {
      case State::attack:
        attackTime = adaptTime(attackTime, noteFreq);
        atk.set(sampleRate, attackTime);
        atkNeg.set(sampleRate, attackTime);
        // Fall through.

      case State::decay:
        dec.set(sampleRate, decayTime);
        // Fall through.

      case State::sustain:
        sustain
          = std::max<Sample>(Sample(0.0), std::min<Sample>(sustainLevel, Sample(1.0)));
        // Fall through.

      case State::release:
        rel.set(sampleRate, adaptTime(releaseTime, noteFreq));
        // Fall through.

      default:
        break;
    }
  }

  void release()
  {
    range = value;
    state = State::release;
  }

  void terminate()
  {
    state = State::terminated;
    value = 0;
    offset = 0;
    range = 1;
    sustain = 1;
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  inline Sample declickIn(Sample input)
  {
    if (declickCounter >= declickLength) return input;
    declickCounter += 1;
    return input * cosinterp<Sample>(declickCounter / (Sample)declickLength);
  }

  Sample process()
  {
    switch (state) {
      case State::attack: {
        const auto atkPos = atk.process();
        const auto atkMix = atkPos + curve * (atkNeg.process() - atkPos);
        value = range * declickIn(atkMix) + offset;
        if (atk.isTerminated()) {
          state = State::decay;
          range = Sample(1.0) - sustain;
        }
      } break;

      case State::decay:
        value = range * declickIn(dec.process()) + sustain;
        if (value <= sustain) state = State::sustain;
        break;

      case State::sustain:
        value = declickIn(sustain);
        break;

      case State::release:
        value = range * declickIn(rel.process());
        if (rel.isTerminated()) state = State::terminated;
        break;

      default:
        return 0;
    }
    return value;
  }

protected:
  enum class State : int32_t { attack, decay, sustain, release, terminated };

  int32_t declickLength;
  int32_t declickCounter = 0;

  ExpAttackCurve<Sample> atk{};
  NegativeExpAttackCurve<Sample> atkNeg{};
  ExpDecayCurve<Sample> dec{};
  ExpDecayCurve<Sample> rel{};

  State state = State::terminated;
  Sample value = 0;
  Sample curve = 0;
  Sample sampleRate = 44100;
  Sample offset = 0;
  Sample range = 1;
  Sample sustain = 1;
};
