#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Pitch
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	const DspFloatType referenceFrequency = 261.626; // C4; frequency at which Rack 1v/octave CVs are zero.
	const DspFloatType referenceSemitone = 60.0; // C4; value of C4 in semitones is arbitrary here, so have it match midi note numbers when rounded to integer.
	const DspFloatType twelfthRootTwo = 1.0594630943592953;
	const DspFloatType logTwelfthRootTwo = logf(1.0594630943592953);


	inline DspFloatType frequencyToSemitone(DspFloatType frequency) {
		return logf(frequency / referenceFrequency) / logTwelfthRootTwo + referenceSemitone;
	}

	inline DspFloatType semitoneToFrequency(DspFloatType semitone) {
		return powf(twelfthRootTwo, semitone - referenceSemitone) * referenceFrequency;
	}

	inline DspFloatType frequencyToCV(DspFloatType frequency) {
		return log2f(frequency / referenceFrequency);
	}

	inline DspFloatType cvToFrequency(DspFloatType cv) {
		return powf(2.0, cv) * referenceFrequency;
	}

	inline DspFloatType cvToSemitone(DspFloatType cv) {
		return frequencyToSemitone(cvToFrequency(cv));
	}

	inline DspFloatType semitoneToCV(DspFloatType semitone) {
		return frequencyToCV(semitoneToFrequency(semitone));
	}
}
