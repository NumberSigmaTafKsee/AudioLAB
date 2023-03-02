#pragma once
#include <cmath>
#include <vector>
#include <cstring>
#include <cstdint>

    constexpr DspFloatType kTwoPi = 2.0*M_PI;

    const DspFloatType GUI_Q_MIN = 1.0;
	const DspFloatType GUI_Q_MAX = 10.0;
	const DspFloatType SVF_Q_SLOPE = (25.0 - 0.707) / (GUI_Q_MAX - GUI_Q_MIN);
	const DspFloatType KORG35_Q_SLOPE = (2.0 - 0.707) / (GUI_Q_MAX - GUI_Q_MIN);
	const DspFloatType MOOG_Q_SLOPE = (4.0 - 0.0) / (GUI_Q_MAX - GUI_Q_MIN);
	const DspFloatType DIODE_Q_SLOPE = (17.0 - 0.0) / (GUI_Q_MAX - GUI_Q_MIN);
	const DspFloatType HSYNC_MOD_SLOPE = 3.0;

    enum { LPF1, LPF2, LPF3, LPF4, HPF1, HPF2, HPF3, HPF4, BPF2, BPF4, BSF2, BSF4, 
		APF1, APF2, ANM_LPF1, ANM_LPF2, ANM_LPF3, ANM_LPF4, NUM_FILTER_OUTPUTS };

    inline DspFloatType semitonesBetweenFreqs(DspFloatType startFrequency, DspFloatType endFrequency)
	{
		return log2(endFrequency / startFrequency)*12.0;
	}

	inline void boundValue(DspFloatType& value, DspFloatType minValue, DspFloatType maxValue)
	{
		const DspFloatType t = value < minValue ? minValue : value;
		value = t > maxValue ? maxValue : t;
	}
	inline void mapDoubleValue(DspFloatType& value, DspFloatType min, DspFloatType minMap, DspFloatType slope)
	{
		// --- bound to limits
		value = minMap + slope * (value - min);
	}

	inline void mapDoubleValue(DspFloatType& value, DspFloatType min, DspFloatType max, DspFloatType minMap, DspFloatType maxMap)
	{
		// --- bound to limits
		boundValue(value, min, max);
		DspFloatType mapped = ((value - min) / (max - min)) * (maxMap - minMap) + minMap;
		value = mapped;
	}	
    struct FilterOutput
	{
        DspFloatType filter[NUM_FILTER_OUTPUTS];
		FilterOutput() { clearData(); }		
		void clearData()
		{
			memset(&filter[0], 0, sizeof(DspFloatType) * NUM_FILTER_OUTPUTS);
		}
	};

   class IFilterBase
	{
		virtual bool reset(DspFloatType _sampleRate) = 0;
		virtual bool update() = 0;
		virtual FilterOutput* process(DspFloatType xn) = 0;
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q) = 0;
	};

    
	struct VA1Coeffs
	{
		// --- filter coefficients
		DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		DspFloatType beta = 1.0;			///< beta value, not used
	};

	
    //@{
	/**
	\ingroup Constants-Enums
	Constants for the synth filter objects and cores
	*/
	const DspFloatType freqModLow = 20.0;
	const DspFloatType freqModHigh = 18000.0; // -- this is reduced from 20480.0 only for self oscillation at upper frequency range
	const DspFloatType freqModSemitoneRange = semitonesBetweenFreqs(freqModLow, freqModHigh);
	const uint32_t FILTER_AUDIO_INPUTS = 2;
	const uint32_t FILTER_AUDIO_OUTPUTS = 2;
	enum class FilterModel { kFirstOrder, kSVF, kKorg35, kMoog, kDiode };
	enum { FLT1, FLT2, FLT3, FLT4 };
	const int MOOG_SUBFILTERS = 4;
	const int DIODE_SUBFILTERS = 4;
	const int KORG_SUBFILTERS = 3;
	enum class VAFilterAlgorithm {
		kBypassFilter, kLPF1, kHPF1, kAPF1, kSVF_LP, kSVF_HP, kSVF_BP, kSVF_BS, kKorg35_LP, kKorg35_HP, kMoog_LP1, kMoog_LP2, kMoog_LP3, kMoog_LP4, kDiode_LP4
	};
	enum class BQFilterAlgorithm {
		kBypassFilter, k1PLPF, k1PHPF, kLPF2, kHPF2
	};
	//@}


	/**
	\struct FilterParameters
	\ingroup SynthStructures
	\brief
	Custom parameter structure for the any of the synth filters. This structure is designed to take care
	of both VA and biquad filter parameter requirements. Notable members:
	- filterIndex the selected index from a GUI control that the user toggles
	- modKnobValue[4] the four mod knob values on the range [0, 1]
	\author Will Pirkle http://www.willpirkle.com
	\remark This object is included in Designing Software Synthesizer Plugins in C++ 2nd Ed. by Will Pirkle
	\version Revision : 1.0
	\date Date : 2021 / 05 / 02
	*/
	struct FilterParameters
	{
		// --- use with strongly typed enums
		int32_t filterIndex = 0;	///< filter index in GUI control

		DspFloatType fc = 1000.0;					///< parameter fc
		DspFloatType Q = 1.0;						///< parameter Q
		DspFloatType filterOutputGain_dB = 0.0;	///< parameter output gain in dB
		DspFloatType filterDrive = 1.0;			///< parameter drive (distortion)
		DspFloatType bassGainComp = 0.0;			///< 0.0 = no bass compensation, 1.0 = restore all bass 
		bool analogFGN = true;				///< use analog FGN filters; adds to CPU load

		// --- key tracking
		bool enableKeyTrack = false;		///< key track flag
		DspFloatType keyTrackSemis = 0.0;			///< key tracking ratio in semitones

		// --- Mod Knobs and core support
		DspFloatType modKnobValue[4] = { 0.5, 0.0, 0.0, 0.0 }; ///< mod knobs
		uint32_t moduleIndex = 0;	///< module identifier
	};


	class VA1Filter : public IFilterBase
	{
	public:
		// --- constructor/destructor
		VA1Filter();
		virtual ~VA1Filter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);

		// --- set coeffs directly, bypassing coeff calculation
		void setAlpha(DspFloatType _alpha) { coeffs.alpha = _alpha; }
		void setBeta(DspFloatType _beta) { coeffs.beta = _beta; }
		void setCoeffs(VA1Coeffs& _coeffs) {
			coeffs = _coeffs;
		}

		void copyCoeffs(VA1Filter& destination) { 
			destination.setCoeffs(coeffs);
		}
		
		// --- added for MOOG & K35, need access to this output value, scaled by beta
		DspFloatType getFBOutput() { return coeffs.beta * sn; }

	protected:
		FilterOutput output;
		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;

		// --- state storage
		DspFloatType sn = 0.0;						///< state variables

		// --- filter coefficients
		VA1Coeffs coeffs;
	};

	struct VASVFCoeffs
	{
		// --- filter coefficients
		DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		DspFloatType rho = 1.0;			///< beta value, not used
		DspFloatType sigma = 1.0;			///< beta value, not used
		DspFloatType alpha0 = 1.0;			///< beta value, not used
	};

	// --- makes all filter outputs: LPF1, LPF1A, HPF1, APF1
	class VASVFilter : public IFilterBase
	{
	public:
		// --- constructor/destructor
		VASVFilter();
		virtual ~VASVFilter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);

		// --- set coeffs directly, bypassing coeff calculation
		void setCoeffs(VASVFCoeffs& _coeffs) {
			coeffs = _coeffs;
		}

		void copyCoeffs(VASVFilter& destination) {
			destination.setCoeffs(coeffs);
		}

	protected:
		FilterOutput output;
		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;
		DspFloatType Q = 0.0;

		// --- state storage
		DspFloatType integrator_z[2];						///< state variables

		// --- filter coefficients
		VASVFCoeffs coeffs;
	};

	struct VAKorg35Coeffs
	{
		// --- filter coefficients
		DspFloatType K = 1.0;			///< beta value, not used
		DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		DspFloatType alpha0 = 1.0;			///< beta value, not used
		DspFloatType g = 1.0;			///< beta value, not used
	};

	// --- makes both LPF and HPF (DspFloatType filter)
	class VAKorg35Filter : public IFilterBase
	{
	public:
		// --- constructor/destructor
		VAKorg35Filter();
		virtual ~VAKorg35Filter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);

		// --- set coeffs directly, bypassing coeff calculation
		void setCoeffs(VAKorg35Coeffs& _coeffs) {
			coeffs = _coeffs;

			// --- three sync-tuned filters
			for (uint32_t i = 0; i < KORG_SUBFILTERS; i++)
			{
				lpfVAFilters[i].setAlpha(coeffs.alpha);
				hpfVAFilters[i].setAlpha(coeffs.alpha);
			}

			// --- set filter beta values
			DspFloatType deno = 1.0 + coeffs.g;

			lpfVAFilters[FLT2].setBeta((coeffs.K * (1.0 - coeffs.alpha)) / deno);
			lpfVAFilters[FLT3].setBeta(-1.0 / deno);

			hpfVAFilters[FLT2].setBeta(-coeffs.alpha / deno);
			hpfVAFilters[FLT3].setBeta(1.0 / deno);
		//	hpfVAFilters[FLT3].setBeta(lpfVAFilters[FLT3].getBeta);
		}

		void copyCoeffs(VAKorg35Filter& destination) {
			destination.setCoeffs(coeffs);
		}

	protected:
		FilterOutput output;
		VA1Filter lpfVAFilters[KORG_SUBFILTERS];
		VA1Filter hpfVAFilters[KORG_SUBFILTERS];
		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;

		// --- filter coefficients
		VAKorg35Coeffs coeffs;

		//DspFloatType K = 0.0;
		//DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		//DspFloatType alpha0 = 0.0;		///< input scalar, correct delay-free loop
	};

	struct VAMoogCoeffs
	{
		// --- filter coefficients
		DspFloatType K = 1.0;			///< beta value, not used
		DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		DspFloatType alpha0 = 1.0;			///< beta value, not used
		DspFloatType sigma = 1.0;			///< beta value, not used
		DspFloatType bassComp = 1.0;			///< beta value, not used
		DspFloatType g = 1.0;			///< beta value, not used

		// --- these are to minimize repeat calculations for left/right pairs
		DspFloatType subFilterBeta[MOOG_SUBFILTERS] = { 0.0, 0.0, 0.0, 0.0 };
	};

	// --- makes both LPF and HPF (DspFloatType filter)
	class VAMoogFilter : public IFilterBase
	{
	public:
		// --- constructor/destructor
		VAMoogFilter();
		virtual ~VAMoogFilter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);

		// --- set coeffs directly, bypassing coeff calculation
		void setCoeffs(const VAMoogCoeffs& _coeffs) {
			coeffs = _coeffs;

			// --- four sync-tuned filters
			for (uint32_t i = 0; i < MOOG_SUBFILTERS; i++)
			{
				// --- set alpha directly
				subFilter[i].setAlpha(coeffs.alpha);
				subFilterFGN[i].setAlpha(coeffs.alpha);

				// --- set beta directly
				subFilter[i].setBeta(coeffs.subFilterBeta[i]);
				subFilterFGN[i].setBeta(coeffs.subFilterBeta[i]);
			}
		}

		void copyCoeffs(VAMoogFilter& destination) {
			destination.setCoeffs(coeffs);
		}

	protected:
		FilterOutput output;
		VA1Filter subFilter[MOOG_SUBFILTERS];
		VA1Filter subFilterFGN[MOOG_SUBFILTERS];
		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;

		// --- filter coefficients
		VAMoogCoeffs coeffs;
	};

	struct DiodeVA1Coeffs
	{
		// --- filter coefficients
		DspFloatType alpha0 = 0.0;		///< input scalar, correct delay-free loop
		DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
		DspFloatType beta = 1.0;			///< beta value, not used
		DspFloatType gamma = 1.0;			///< beta value, not used
		DspFloatType delta = 1.0;			///< beta value, not used
		DspFloatType epsilon = 1.0;		///< beta value, not used
	};

	class VADiodeSubFilter : public IFilterBase
	{
	public:
		// --- constructor/destructor
		VADiodeSubFilter();
		virtual ~VADiodeSubFilter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);

		void setCoeffs(const DiodeVA1Coeffs& _coeffs) { coeffs = _coeffs; }
		void copyCoeffs(VADiodeSubFilter& destination) {
			destination.setCoeffs(coeffs);
		}

		void setFBInput(DspFloatType _feedbackIn) { feedbackIn = _feedbackIn; }
		DspFloatType getFBOutput() { return coeffs.beta * (sn + feedbackIn*coeffs.delta); }

	protected:
		DiodeVA1Coeffs coeffs;
		FilterOutput output;
		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;

		// --- state storage
		DspFloatType sn = 0.0;						///< state variables
		DspFloatType feedbackIn = 0.0;
	};
	
	struct VADiodeCoeffs
	{
		// --- filter coefficients

		DspFloatType alpha0 = 0.0;		///< input scalar, correct delay-free loop
		DspFloatType gamma = 0.0;		///< input scalar, correct delay-free loop
		DspFloatType beta[4] = { 0.0, 0.0, 0.0, 0.0 };
		DiodeVA1Coeffs diodeCoeffs[4];
		DspFloatType bassComp = 0.0;		// --- increase for MORE bass
		DspFloatType alpha1 = 1.0;		// --- FGN amp correction
		DspFloatType K = 1.0;		// --- 
	};

	// --- makes both LPF and HPF (DspFloatType filter)
	class VADiodeFilter
	{
	public:
		// --- constructor/destructor
		VADiodeFilter();
		virtual ~VADiodeFilter() {}

		// --- these match SynthModule names
		virtual bool reset(DspFloatType _sampleRate);
		virtual bool update();
		virtual FilterOutput* process(DspFloatType xn); 
		virtual void setFilterParams(DspFloatType _fc, DspFloatType _Q);
	
		void setCoeffs(const VADiodeCoeffs& _coeffs) {
			coeffs = _coeffs; 

			// --- update subfilter coeffs
			for (uint32_t i = 0; i < DIODE_SUBFILTERS; i++)
			{
				subFilter[i].setCoeffs(coeffs.diodeCoeffs[i]);
				subFilterFGN[i].setCoeffs(coeffs.diodeCoeffs[i]);
			}
		}

		void copyCoeffs(VADiodeFilter& destination) {
			destination.setCoeffs(coeffs);
		}

	protected:
		FilterOutput output;
		VADiodeSubFilter subFilter[DIODE_SUBFILTERS];
		VADiodeSubFilter subFilterFGN[DIODE_SUBFILTERS];

		DspFloatType sampleRate = 44100.0;				///< current sample rate
		DspFloatType halfSamplePeriod = 1.0;
		DspFloatType fc = 0.0;

		// --- filter coefficients
		VADiodeCoeffs coeffs;
	};



