
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: diodeladder.dsp
// Name: DiodeLadder
// Author: Christopher Arndt
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------

//#include "DiodeLadder.hpp"

#include <utility>
#include <cmath>
#include <algorithm>


namespace Analog::Filters::DiodeLadder 
{
	class DiodeLadder::BasicDsp {
	public:
		virtual ~BasicDsp() {}
	};

	//------------------------------------------------------------------------------
	// Begin the Faust code section

	namespace {

	template <class T> inline T min(T a, T b) { return (a < b) ? a : b; }
	template <class T> inline T max(T a, T b) { return (a > b) ? a : b; }

	class Meta {
	public:
		// dummy
		void declare(...) {}
	};

	class UI {
	public:
		// dummy
		void openHorizontalBox(...) {}
		void openVerticalBox(...) {}
		void closeBox(...) {}
		void declare(...) {}
		void addButton(...) {}
		void addCheckButton(...) {}
		void addVerticalSlider(...) {}
		void addHorizontalSlider(...) {}
		void addVerticalBargraph(...) {}
		void addHorizontalBargraph(...) {}
	};

	typedef DiodeLadder::BasicDsp dsp;

	} // namespace

	#define FAUSTPP_VIRTUAL // do not declare any methods virtual
	#define FAUSTPP_PRIVATE public // do not hide any members
	#define FAUSTPP_PROTECTED public // do not hide any members

	// define the DSP in the anonymous namespace
	#define FAUSTPP_BEGIN_NAMESPACE namespace {
	#define FAUSTPP_END_NAMESPACE }


	#if defined(__GNUC__)
	#   pragma GCC diagnostic push
	#   pragma GCC diagnostic ignored "-Wunused-parameter"
	#endif

	#ifndef FAUSTPP_PRIVATE
	#   define FAUSTPP_PRIVATE private
	#endif
	#ifndef FAUSTPP_PROTECTED
	#   define FAUSTPP_PROTECTED protected
	#endif
	#ifndef FAUSTPP_VIRTUAL
	#   define FAUSTPP_VIRTUAL virtual
	#endif

	#ifndef FAUSTPP_BEGIN_NAMESPACE
	#   define FAUSTPP_BEGIN_NAMESPACE
	#endif
	#ifndef FAUSTPP_END_NAMESPACE
	#   define FAUSTPP_END_NAMESPACE
	#endif

	FAUSTPP_BEGIN_NAMESPACE

	#ifndef FAUSTDspFloatType
	#define FAUSTDspFloatType DspFloatType
	#endif 

	static DspFloatType mydsp_faustpower2_f(DspFloatType value) {
		return (value * value);
	}
	static DspFloatType mydsp_faustpower10_f(DspFloatType value) {
		return (((((((((value * value) * value) * value) * value) * value) * value) * value) * value) * value);
	}
	static DspFloatType mydsp_faustpower3_f(DspFloatType value) {
		return ((value * value) * value);
	}
	static DspFloatType mydsp_faustpower4_f(DspFloatType value) {
		return (((value * value) * value) * value);
	}

	#ifndef FAUSTCLASS 
	#define FAUSTCLASS mydsp
	#endif

	#ifdef __APPLE__ 
	#define exp10f __exp10f
	#define exp10 __exp10
	#endif

	class mydsp : public dsp {
		
	FAUSTPP_PRIVATE:
		
		int fSampleRate;
		DspFloatType fConst0;
		FAUSTDspFloatType fHslider0;
		DspFloatType fRec5[2];
		FAUSTDspFloatType fHslider1;
		DspFloatType fRec1[2];
		DspFloatType fRec2[2];
		DspFloatType fRec3[2];
		DspFloatType fRec4[2];
		
	public:
		
		void metadata(Meta* m) { 
			m->declare("../../faust/diodeladder.dsp/diodeLadder:author", "Eric Tarr");
			m->declare("../../faust/diodeladder.dsp/diodeLadder:license", "MIT-style STK-4.3 license");
			m->declare("author", "Christopher Arndt");
			m->declare("description", "FAUST Diode Ladder 24 dB LPF");
			m->declare("filename", "diodeladder.dsp");
			m->declare("license", "MIT-style STK-4.3 license");
			m->declare("maths.lib/author", "GRAME");
			m->declare("maths.lib/copyright", "GRAME");
			m->declare("maths.lib/license", "LGPL with exception");
			m->declare("maths.lib/name", "Faust Math Library");
			m->declare("maths.lib/version", "2.3");
			m->declare("misceffects.lib/name", "Misc Effects Library");
			m->declare("misceffects.lib/version", "2.0");
			m->declare("name", "DiodeLadder");
			m->declare("platform.lib/name", "Generic Platform Library");
			m->declare("platform.lib/version", "0.1");
			m->declare("signals.lib/name", "Faust Signal Routing Library");
			m->declare("signals.lib/version", "0.0");
		}

		FAUSTPP_VIRTUAL int getNumInputs() {
			return 1;
		}
		FAUSTPP_VIRTUAL int getNumOutputs() {
			return 1;
		}
		FAUSTPP_VIRTUAL int getInputRate(int channel) {
			int rate;
			switch ((channel)) {
				case 0: {
					rate = 1;
					break;
				}
				default: {
					rate = -1;
					break;
				}
			}
			return rate;
		}
		FAUSTPP_VIRTUAL int getOutputRate(int channel) {
			int rate;
			switch ((channel)) {
				case 0: {
					rate = 1;
					break;
				}
				default: {
					rate = -1;
					break;
				}
			}
			return rate;
		}
		
		static void classInit(int sample_rate) {
		}
		
		FAUSTPP_VIRTUAL void instanceConstants(int sample_rate) {
			fSampleRate = sample_rate;
			fConst0 = (3.14159274f / std::min<DspFloatType>(192000.0f, std::max<DspFloatType>(1.0f, DspFloatType(fSampleRate))));
		}
		
		FAUSTPP_VIRTUAL void instanceResetUserInterface() {
			fHslider0 = FAUSTDspFloatType(20000.0f);
			fHslider1 = FAUSTDspFloatType(1.0f);
		}
		
		FAUSTPP_VIRTUAL void instanceClear() {
			for (int l0 = 0; (l0 < 2); l0 = (l0 + 1)) {
				fRec5[l0] = 0.0f;
			}
			for (int l1 = 0; (l1 < 2); l1 = (l1 + 1)) {
				fRec1[l1] = 0.0f;
			}
			for (int l2 = 0; (l2 < 2); l2 = (l2 + 1)) {
				fRec2[l2] = 0.0f;
			}
			for (int l3 = 0; (l3 < 2); l3 = (l3 + 1)) {
				fRec3[l3] = 0.0f;
			}
			for (int l4 = 0; (l4 < 2); l4 = (l4 + 1)) {
				fRec4[l4] = 0.0f;
			}
		}
		
		FAUSTPP_VIRTUAL void init(int sample_rate) {
			classInit(sample_rate);
			instanceInit(sample_rate);
		}
		FAUSTPP_VIRTUAL void instanceInit(int sample_rate) {
			instanceConstants(sample_rate);
			instanceResetUserInterface();
			instanceClear();
		}
		
		FAUSTPP_VIRTUAL mydsp* clone() {
			return new mydsp();
		}
		
		FAUSTPP_VIRTUAL int getSampleRate() {
			return fSampleRate;
		}
		
		FAUSTPP_VIRTUAL void buildUserInterface(UI* ui_interface) {
			ui_interface->openVerticalBox("DiodeLadder");
			ui_interface->declare(&fHslider0, "0", "");
			ui_interface->declare(&fHslider0, "abbrev", "cutoff");
			ui_interface->declare(&fHslider0, "scale", "log");
			ui_interface->declare(&fHslider0, "style", "knob");
			ui_interface->declare(&fHslider0, "symbol", "cutoff");
			ui_interface->declare(&fHslider0, "unit", "hz");
			ui_interface->addHorizontalSlider("Cutoff frequency", &fHslider0, 20000.0f, 20.0f, 20000.0f, 0.100000001f);
			ui_interface->declare(&fHslider1, "1", "");
			ui_interface->declare(&fHslider1, "abbrev", "q");
			ui_interface->declare(&fHslider1, "style", "knob");
			ui_interface->declare(&fHslider1, "symbol", "q");
			ui_interface->addHorizontalSlider("Q", &fHslider1, 1.0f, 0.707199991f, 25.0f, 0.00999999978f);
			ui_interface->closeBox();
		}
		
		FAUSTPP_VIRTUAL void compute(int count, FAUSTDspFloatType** inputs, FAUSTDspFloatType** outputs) {
			FAUSTDspFloatType* input0 = inputs[0];
			FAUSTDspFloatType* output0 = outputs[0];
			DspFloatType fSlow0 = (0.00100000005f * DspFloatType(fHslider0));
			DspFloatType fSlow1 = (DspFloatType(fHslider1) + -0.707000017f);
			DspFloatType fSlow2 = (0.00514551532f * fSlow1);
			for (int i = 0; (i < count); i = (i + 1)) {
				fRec5[0] = (fSlow0 + (0.999000013f * fRec5[1]));
				DspFloatType fTemp0 = std::tan((fConst0 * fRec5[0]));
				DspFloatType fTemp1 = std::max<DspFloatType>(-1.0f, std::min<DspFloatType>(1.0f, (100.0f * DspFloatType(input0[i]))));
				DspFloatType fTemp2 = (17.0f - (9.69999981f * mydsp_faustpower10_f((0.0f - (0.333333343f * (1.0f - (std::log10(fRec5[0]) + -0.30103001f)))))));
				DspFloatType fTemp3 = (fTemp0 + 1.0f);
				DspFloatType fTemp4 = ((0.5f * ((fRec1[1] * fTemp0) / fTemp3)) + fRec2[1]);
				DspFloatType fTemp5 = ((fTemp0 * (1.0f - (0.25f * (fTemp0 / fTemp3)))) + 1.0f);
				DspFloatType fTemp6 = ((fTemp0 * fTemp4) / fTemp5);
				DspFloatType fTemp7 = (0.5f * fTemp6);
				DspFloatType fTemp8 = (fTemp7 + fRec3[1]);
				DspFloatType fTemp9 = ((fTemp0 * (1.0f - (0.25f * (fTemp0 / fTemp5)))) + 1.0f);
				DspFloatType fTemp10 = ((fTemp0 * fTemp8) / fTemp9);
				DspFloatType fTemp11 = (fTemp10 + fRec4[1]);
				DspFloatType fTemp12 = (fTemp5 * fTemp9);
				DspFloatType fTemp13 = ((fTemp0 * (1.0f - (0.5f * (fTemp0 / fTemp9)))) + 1.0f);
				DspFloatType fTemp14 = mydsp_faustpower2_f(fTemp0);
				DspFloatType fTemp15 = (fTemp3 * fTemp5);
				DspFloatType fTemp16 = ((fTemp0 * ((((((1.5f * (fTemp1 * (1.0f - (0.333333343f * mydsp_faustpower2_f(fTemp1))))) + (fSlow1 * ((fTemp2 * (((0.0f - ((0.0205820613f * fTemp6) + (0.0411641225f * fRec1[1]))) - (0.0205820613f * fTemp10)) - (0.00514551532f * ((mydsp_faustpower3_f(fTemp0) * fTemp11) / (fTemp12 * fTemp13))))) / fTemp3))) * ((0.5f * (fTemp14 / (fTemp9 * fTemp13))) + 1.0f)) / ((fSlow2 * ((mydsp_faustpower4_f(fTemp0) * fTemp2) / ((fTemp15 * fTemp9) * fTemp13))) + 1.0f)) + ((fTemp8 + (0.5f * ((fTemp0 * fTemp11) / fTemp13))) / fTemp9)) - fRec4[1])) / fTemp3);
				DspFloatType fTemp17 = ((fTemp0 * ((0.5f * (((fRec4[1] + fTemp16) * ((0.25f * (fTemp14 / fTemp12)) + 1.0f)) + ((fTemp4 + (0.5f * fTemp10)) / fTemp5))) - fRec3[1])) / fTemp3);
				DspFloatType fTemp18 = ((fTemp0 * ((0.5f * (((fRec3[1] + fTemp17) * ((0.25f * (fTemp14 / fTemp15)) + 1.0f)) + ((fRec1[1] + fTemp7) / fTemp3))) - fRec2[1])) / fTemp3);
				DspFloatType fTemp19 = ((fTemp0 * ((0.5f * (fRec2[1] + fTemp18)) - fRec1[1])) / fTemp3);
				DspFloatType fRec0 = (fRec1[1] + fTemp19);
				fRec1[0] = (fRec1[1] + (2.0f * fTemp19));
				fRec2[0] = (fRec2[1] + (2.0f * fTemp18));
				fRec3[0] = (fRec3[1] + (2.0f * fTemp17));
				fRec4[0] = (fRec4[1] + (2.0f * fTemp16));
				output0[i] = FAUSTDspFloatType(fRec0);
				fRec5[1] = fRec5[0];
				fRec1[1] = fRec1[0];
				fRec2[1] = fRec2[0];
				fRec3[1] = fRec3[0];
				fRec4[1] = fRec4[0];
			}
		}

	};
	FAUSTPP_END_NAMESPACE


	#if defined(__GNUC__)
	#   pragma GCC diagnostic pop
	#endif



	//------------------------------------------------------------------------------
	// End the Faust code section




	DiodeLadder::DiodeLadder()
	{

		mydsp *dsp = new mydsp;
		fDsp.reset(dsp);
		dsp->instanceResetUserInterface();

	}

	DiodeLadder::~DiodeLadder()
	{
	}

	void DiodeLadder::init(DspFloatType sample_rate)
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.classInit(sample_rate);
		dsp.instanceConstants(sample_rate);
		clear();

	}

	void DiodeLadder::clear() noexcept
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.instanceClear();

	}

	void DiodeLadder::process(
		const DspFloatType *in0,
		DspFloatType *out0,
		unsigned count) noexcept
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		DspFloatType *inputs[] = {
			const_cast<DspFloatType *>(in0),
		};
		DspFloatType *outputs[] = {
			out0,
		};
		dsp.compute(count, inputs, outputs);

	}

	const char *DiodeLadder::parameter_label(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return "Cutoff frequency";
		
		case 1:
			return "Q";
		
		default:
			return 0;
		}
	}

	const char *DiodeLadder::parameter_short_label(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return "cutoff";
		
		case 1:
			return "q";
		
		default:
			return 0;
		}
	}

	const char *DiodeLadder::parameter_symbol(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return "cutoff";
		
		case 1:
			return "q";
		
		default:
			return 0;
		}
	}

	const char *DiodeLadder::parameter_unit(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return "hz";
		
		case 1:
			return "";
		
		default:
			return 0;
		}
	}

	const ParameterRange *DiodeLadder::parameter_range(unsigned index) noexcept
	{
		switch (index) {
		
		case 0: {
			static const ParameterRange range = { 20000, 20, 20000 };
			return &range;
		}
		
		case 1: {
			static const ParameterRange range = { 1, 0.70719999, 25 };
			return &range;
		}
		
		default:
			return 0;
		}
	}

	bool DiodeLadder::parameter_is_trigger(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool DiodeLadder::parameter_is_boolean(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool DiodeLadder::parameter_is_integer(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool DiodeLadder::parameter_is_logarithmic(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return true;
		
		default:
			return false;
		}
	}

	DspFloatType DiodeLadder::get_parameter(unsigned index) const noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		switch (index) {
		
		case 0:
			return dsp.fHslider0;
		
		case 1:
			return dsp.fHslider1;
		
		default:
			(void)dsp;
			return 0;
		}
	}

	void DiodeLadder::set_parameter(unsigned index, DspFloatType value) noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		switch (index) {
		
		case 0:
			dsp.fHslider0 = value;
			break;
		
		case 1:
			dsp.fHslider1 = value;
			break;
		
		default:
			(void)dsp;
			(void)value;
			break;
		}
	}


	DspFloatType DiodeLadder::get_cutoff() const noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		return dsp.fHslider0;
	}

	DspFloatType DiodeLadder::get_q() const noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		return dsp.fHslider1;
	}


	void DiodeLadder::set_cutoff(DspFloatType value) noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.fHslider0 = value;
	}

	void DiodeLadder::set_q(DspFloatType value) noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.fHslider1 = value;
	}
}