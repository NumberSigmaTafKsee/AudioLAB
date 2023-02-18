
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: oberheim.dsp
// Name: Oberheim
// Author: Christopher Arndt
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------

#include <utility>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <math.h>

namespace Analog::Filters::Oberheim
{
	class Oberheim::BasicDsp {
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

	typedef Oberheim::BasicDsp dsp;

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


	static DspFloatType mydsp_faustpower2_f(DspFloatType value) {
		return (value * value);
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
		DspFloatType fHslider0;
		DspFloatType fRec6[2];
		DspFloatType fHslider1;
		DspFloatType fRec4[2];
		DspFloatType fRec5[2];
		
	public:
		
		void metadata(Meta* m) { 
			m->declare("../../faust/oberheim.dsp/oberheim:author", "Eric Tarr");
			m->declare("../../faust/oberheim.dsp/oberheim:license", "MIT-style STK-4.3 license");
			m->declare("author", "Christopher Arndt");
			m->declare("description", "Oberheim multi-mode, state-variable filter");
			m->declare("filename", "oberheim.dsp");
			m->declare("license", "MIT-style STK-4.3 license");
			m->declare("maths.lib/author", "GRAME");
			m->declare("maths.lib/copyright", "GRAME");
			m->declare("maths.lib/license", "LGPL with exception");
			m->declare("maths.lib/name", "Faust Math Library");
			m->declare("maths.lib/version", "2.3");
			m->declare("misceffects.lib/name", "Misc Effects Library");
			m->declare("misceffects.lib/version", "2.0");
			m->declare("name", "Oberheim");
			m->declare("platform.lib/name", "Generic Platform Library");
			m->declare("platform.lib/version", "0.1");
			m->declare("signals.lib/name", "Faust Signal Routing Library");
			m->declare("signals.lib/version", "0.0");
		}

		FAUSTPP_VIRTUAL int getNumInputs() {
			return 1;
		}
		FAUSTPP_VIRTUAL int getNumOutputs() {
			return 4;
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
				case 1: {
					rate = 1;
					break;
				}
				case 2: {
					rate = 1;
					break;
				}
				case 3: {
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
			fHslider0 = DspFloatType(20000.0f);
			fHslider1 = DspFloatType(1.0f);
		}
		
		FAUSTPP_VIRTUAL void instanceClear() {
			for (int l0 = 0; (l0 < 2); l0 = (l0 + 1)) {
				fRec6[l0] = 0.0f;
			}
			for (int l1 = 0; (l1 < 2); l1 = (l1 + 1)) {
				fRec4[l1] = 0.0f;
			}
			for (int l2 = 0; (l2 < 2); l2 = (l2 + 1)) {
				fRec5[l2] = 0.0f;
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
			ui_interface->openVerticalBox("Oberheim");
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
			ui_interface->addHorizontalSlider("Q", &fHslider1, 1.0f, 0.5f, 10.0f, 0.00999999978f);
			ui_interface->closeBox();
		}
		
		FAUSTPP_VIRTUAL void compute(int count, DspFloatType** inputs, DspFloatType** outputs) {
			DspFloatType* input0 = inputs[0];
			DspFloatType* output0 = outputs[0];
			DspFloatType* output1 = outputs[1];
			DspFloatType* output2 = outputs[2];
			DspFloatType* output3 = outputs[3];
			DspFloatType fSlow0 = (0.00100000005f * DspFloatType(fHslider0));
			DspFloatType fSlow1 = (1.0f / DspFloatType(fHslider1));
			for (int i = 0; (i < count); i = (i + 1)) {
				fRec6[0] = (fSlow0 + (0.999000013f * fRec6[1]));
				DspFloatType fTemp0 = std::tan((fConst0 * fRec6[0]));
				DspFloatType fTemp1 = (fSlow1 + fTemp0);
				DspFloatType fTemp2 = (DspFloatType(input0[i]) - (fRec4[1] + (fRec5[1] * fTemp1)));
				DspFloatType fTemp3 = ((fTemp0 * fTemp1) + 1.0f);
				DspFloatType fTemp4 = ((fTemp0 * fTemp2) / fTemp3);
				DspFloatType fTemp5 = std::max<DspFloatType>(-1.0f, std::min<DspFloatType>(1.0f, (fRec5[1] + fTemp4)));
				DspFloatType fTemp6 = (1.0f - (0.333333343f * mydsp_faustpower2_f(fTemp5)));
				DspFloatType fTemp7 = ((fTemp0 * fTemp5) * fTemp6);
				DspFloatType fRec0 = (fRec4[1] + fTemp7);
				DspFloatType fTemp8 = (fTemp2 / fTemp3);
				DspFloatType fRec1 = fTemp8;
				DspFloatType fTemp9 = (fTemp5 * fTemp6);
				DspFloatType fRec2 = fTemp9;
				DspFloatType fRec3 = (fTemp7 + (fRec4[1] + fTemp8));
				fRec4[0] = (fRec4[1] + (2.0f * fTemp7));
				fRec5[0] = (fTemp4 + fTemp9);
				output0[i] = DspFloatType(fRec3);
				output1[i] = DspFloatType(fRec2);
				output2[i] = DspFloatType(fRec1);
				output3[i] = DspFloatType(fRec0);
				fRec6[1] = fRec6[0];
				fRec4[1] = fRec4[0];
				fRec5[1] = fRec5[0];
			}
		}

	};


	#if defined(__GNUC__)
	#   pragma GCC diagnostic pop
	#endif



	//------------------------------------------------------------------------------
	// End the Faust code section




	Oberheim::Oberheim()
	{

		mydsp *dsp = new mydsp;
		fDsp.reset(dsp);
		dsp->instanceResetUserInterface();

	}

	Oberheim::~Oberheim()
	{
	}

	void Oberheim::init(DspFloatType sample_rate)
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.classInit(sample_rate);
		dsp.instanceConstants(sample_rate);
		clear();

	}

	void Oberheim::clear() noexcept
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.instanceClear();

	}

	void Oberheim::process(
		const DspFloatType *in0,
		DspFloatType *out0,DspFloatType *out1,DspFloatType *out2,DspFloatType *out3,
		unsigned count) noexcept
	{

		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		DspFloatType *inputs[] = {
			const_cast<DspFloatType *>(in0),
		};
		DspFloatType *outputs[] = {
			out0,out1,out2,out3,
		};
		dsp.compute(count, inputs, outputs);

	}

	const char *Oberheim::parameter_label(unsigned index) noexcept
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

	const char *Oberheim::parameter_short_label(unsigned index) noexcept
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

	const char *Oberheim::parameter_symbol(unsigned index) noexcept
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

	const char *Oberheim::parameter_unit(unsigned index) noexcept
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

	const ParameterRange *Oberheim::parameter_range(unsigned index) noexcept
	{
		switch (index) {
		
		case 0: {
			static const ParameterRange range = { 20000, 20, 20000 };
			return &range;
		}
		
		case 1: {
			static const ParameterRange range = { 1, 0.5, 10 };
			return &range;
		}
		
		default:
			return 0;
		}
	}

	bool Oberheim::parameter_is_trigger(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool Oberheim::parameter_is_boolean(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool Oberheim::parameter_is_integer(unsigned index) noexcept
	{
		switch (index) {
		
		default:
			return false;
		}
	}

	bool Oberheim::parameter_is_logarithmic(unsigned index) noexcept
	{
		switch (index) {
		
		case 0:
			return true;
		
		default:
			return false;
		}
	}

	DspFloatType Oberheim::get_parameter(unsigned index) const noexcept
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

	void Oberheim::set_parameter(unsigned index, DspFloatType value) noexcept
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


	DspFloatType Oberheim::get_cutoff() const noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		return dsp.fHslider0;
	}

	DspFloatType Oberheim::get_q() const noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		return dsp.fHslider1;
	}


	void Oberheim::set_cutoff(DspFloatType value) noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.fHslider0 = value;
	}

	void Oberheim::set_q(DspFloatType value) noexcept
	{
		mydsp &dsp = static_cast<mydsp &>(*fDsp);
		dsp.fHslider1 = value;
	}

}


