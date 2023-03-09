#pragma once

#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <math.h>
#include <utility>
#include <cmath>


class Fverb {
public:
    Fverb();
    ~Fverb();

    void init(DspFloatType sample_rate);
    void clear() noexcept;

    void process(const DspFloatType *in0, const DspFloatType *in1, DspFloatType *out0, DspFloatType *out1,
                 unsigned count) noexcept;

    enum { NumInputs = 2 };
    enum { NumOutputs = 2 };
    enum { NumActives = 11 };
    enum { NumPassives = 0 };
    enum { NumParameters = 11 };

    enum Parameter {
        p_predelay,
        p_input,
        p_input_lowpass,
        p_input_highpass,
        p_input_diffusion_1,
        p_input_diffusion_2,
        p_tail_density,
        p_decay,
        p_damping,
        p_mod_frequency,
        p_mod_depth,

    };

    struct ParameterRange {
        DspFloatType init;
        DspFloatType min;
        DspFloatType max;
    };

    static const char *parameter_label(unsigned index) noexcept;
    static const char *parameter_short_label(unsigned index) noexcept;
    static const char *parameter_symbol(unsigned index) noexcept;
    static const char *parameter_unit(unsigned index) noexcept;
    static const ParameterRange *parameter_range(unsigned index) noexcept;
    static bool parameter_is_trigger(unsigned index) noexcept;
    static bool parameter_is_boolean(unsigned index) noexcept;
    static bool parameter_is_integer(unsigned index) noexcept;
    static bool parameter_is_logarithmic(unsigned index) noexcept;

    DspFloatType get_parameter(unsigned index) const noexcept;
    void set_parameter(unsigned index, DspFloatType value) noexcept;


    DspFloatType get_predelay() const noexcept;

    DspFloatType get_input() const noexcept;

    DspFloatType get_input_lowpass() const noexcept;

    DspFloatType get_input_highpass() const noexcept;

    DspFloatType get_input_diffusion_1() const noexcept;

    DspFloatType get_input_diffusion_2() const noexcept;

    DspFloatType get_tail_density() const noexcept;

    DspFloatType get_decay() const noexcept;

    DspFloatType get_damping() const noexcept;

    DspFloatType get_mod_frequency() const noexcept;

    DspFloatType get_mod_depth() const noexcept;


    void set_predelay(DspFloatType value) noexcept;

    void set_input(DspFloatType value) noexcept;

    void set_input_lowpass(DspFloatType value) noexcept;

    void set_input_highpass(DspFloatType value) noexcept;

    void set_input_diffusion_1(DspFloatType value) noexcept;

    void set_input_diffusion_2(DspFloatType value) noexcept;

    void set_tail_density(DspFloatType value) noexcept;

    void set_decay(DspFloatType value) noexcept;

    void set_damping(DspFloatType value) noexcept;

    void set_mod_frequency(DspFloatType value) noexcept;

    void set_mod_depth(DspFloatType value) noexcept;


public:
    class BasicDsp;

private:
    std::unique_ptr<BasicDsp> fDsp;
};



class Fverb::BasicDsp {
public:
    virtual ~BasicDsp()
    {
    }
};

//------------------------------------------------------------------------------
// Begin the Faust code section

namespace {

template <class T> inline T min(T a, T b)
{
    return (a < b) ? a : b;
}
template <class T> inline T max(T a, T b)
{
    return (a > b) ? a : b;
}

class Meta {
public:
    // dummy
    void declare(...)
    {
    }
};

class UI {
public:
    // dummy
    void openHorizontalBox(...)
    {
    }
    void openVerticalBox(...)
    {
    }
    void closeBox(...)
    {
    }
    void declare(...)
    {
    }
    void addButton(...)
    {
    }
    void addCheckButton(...)
    {
    }
    void addVerticalSlider(...)
    {
    }
    void addHorizontalSlider(...)
    {
    }
    void addVerticalBargraph(...)
    {
    }
    void addHorizontalBargraph(...)
    {
    }
};

typedef Fverb::BasicDsp dsp;

}  // namespace

#define FAUSTPP_VIRTUAL  // do not declare any methods virtual
#define FAUSTPP_PRIVATE public  // do not hide any members
#define FAUSTPP_PROTECTED public  // do not hide any members

// define the DSP in the anonymous namespace
#define FAUSTPP_BEGIN_NAMESPACE namespace {
#define FAUSTPP_END_NAMESPACE }


#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifndef FAUSTPP_PRIVATE
#define FAUSTPP_PRIVATE private
#endif
#ifndef FAUSTPP_PROTECTED
#define FAUSTPP_PROTECTED protected
#endif
#ifndef FAUSTPP_VIRTUAL
#define FAUSTPP_VIRTUAL virtual
#endif

#ifndef FAUSTPP_BEGIN_NAMESPACE
#define FAUSTPP_BEGIN_NAMESPACE
#endif
#ifndef FAUSTPP_END_NAMESPACE
#define FAUSTPP_END_NAMESPACE
#endif

FAUSTPP_BEGIN_NAMESPACE

#ifndef FAUSTFLOAT
#define FAUSTFLOAT DspFloatType
#endif

FAUSTPP_END_NAMESPACE

FAUSTPP_BEGIN_NAMESPACE

class mydspSIG0 {

    FAUSTPP_PRIVATE :

        int iVec6[2];
    int iRec34[2];

public:
    int getNumInputsmydspSIG0()
    {
        return 0;
    }
    int getNumOutputsmydspSIG0()
    {
        return 1;
    }

    void instanceInitmydspSIG0(int sample_rate)
    {
        for (int l21 = 0; (l21 < 2); l21 = (l21 + 1)) {
            iVec6[l21] = 0;
        }
        for (int l22 = 0; (l22 < 2); l22 = (l22 + 1)) {
            iRec34[l22] = 0;
        }
    }

    void fillmydspSIG0(int count, DspFloatType *table)
    {
        for (int i1 = 0; (i1 < count); i1 = (i1 + 1)) {
            iVec6[0] = 1;
            iRec34[0] = ((iVec6[1] + iRec34[1]) % 65536);
            table[i1] = std::sin((9.58738019e-05f * DspFloatType(iRec34[0])));
            iVec6[1] = iVec6[0];
            iRec34[1] = iRec34[0];
        }
    }
};

inline static mydspSIG0 *newmydspSIG0()
{
    return (mydspSIG0 *)new mydspSIG0();
}
inline static void deletemydspSIG0(mydspSIG0 *dsp)
{
    delete dsp;
}

static DspFloatType ftbl0mydspSIG0[65536];

#ifndef FAUSTCLASS
#define FAUSTCLASS mydsp
#endif

#ifdef __APPLE__
#define exp10f __exp10f
#define exp10 __exp10
#endif

class mydsp : public dsp {

    FAUSTPP_PRIVATE :

        int fSampleRate;
    DspFloatType fConst0;
    DspFloatType fConst1;
    FAUSTFLOAT fHslider0;
    DspFloatType fConst2;
    DspFloatType fConst3;
    DspFloatType fRec8[2];
    FAUSTFLOAT fHslider1;
    DspFloatType fRec22[2];
    int IOTA;
    DspFloatType fVec0[131072];
    DspFloatType fConst4;
    FAUSTFLOAT fHslider2;
    DspFloatType fRec23[2];
    DspFloatType fConst5;
    FAUSTFLOAT fHslider3;
    DspFloatType fRec24[2];
    DspFloatType fRec21[2];
    FAUSTFLOAT fHslider4;
    DspFloatType fRec25[2];
    DspFloatType fRec20[2];
    FAUSTFLOAT fHslider5;
    DspFloatType fRec26[2];
    DspFloatType fVec1[1024];
    int iConst6;
    DspFloatType fRec18[2];
    DspFloatType fVec2[1024];
    int iConst7;
    DspFloatType fRec16[2];
    FAUSTFLOAT fHslider6;
    DspFloatType fRec27[2];
    DspFloatType fVec3[4096];
    int iConst8;
    DspFloatType fRec14[2];
    DspFloatType fVec4[2048];
    int iConst9;
    DspFloatType fRec12[2];
    int iConst10;
    FAUSTFLOAT fHslider7;
    DspFloatType fRec28[2];
    DspFloatType fVec5[131072];
    FAUSTFLOAT fHslider8;
    DspFloatType fRec33[2];
    FAUSTFLOAT fHslider9;
    DspFloatType fRec36[2];
    DspFloatType fRec35[2];
    DspFloatType fConst11;
    DspFloatType fConst12;
    DspFloatType fRec29[2];
    DspFloatType fRec30[2];
    int iRec31[2];
    int iRec32[2];
    DspFloatType fRec10[2];
    DspFloatType fVec7[32768];
    int iConst13;
    FAUSTFLOAT fHslider10;
    DspFloatType fRec37[2];
    DspFloatType fRec9[2];
    DspFloatType fVec8[32768];
    int iConst14;
    DspFloatType fRec6[2];
    DspFloatType fRec0[32768];
    DspFloatType fRec1[16384];
    DspFloatType fRec2[32768];
    DspFloatType fVec9[131072];
    DspFloatType fRec52[2];
    DspFloatType fRec51[2];
    DspFloatType fVec10[1024];
    int iConst15;
    DspFloatType fRec49[2];
    DspFloatType fVec11[1024];
    int iConst16;
    DspFloatType fRec47[2];
    DspFloatType fVec12[4096];
    int iConst17;
    DspFloatType fRec45[2];
    DspFloatType fVec13[2048];
    int iConst18;
    DspFloatType fRec43[2];
    int iConst19;
    DspFloatType fVec14[131072];
    DspFloatType fRec53[2];
    DspFloatType fRec54[2];
    int iRec55[2];
    int iRec56[2];
    DspFloatType fRec41[2];
    DspFloatType fVec15[32768];
    int iConst20;
    DspFloatType fRec40[2];
    DspFloatType fVec16[16384];
    int iConst21;
    DspFloatType fRec38[2];
    DspFloatType fRec3[32768];
    DspFloatType fRec4[8192];
    DspFloatType fRec5[32768];
    int iConst22;
    int iConst23;
    int iConst24;
    int iConst25;
    int iConst26;
    int iConst27;
    int iConst28;
    int iConst29;
    int iConst30;
    int iConst31;
    int iConst32;
    int iConst33;
    int iConst34;
    int iConst35;

public:
    void metadata(Meta *m)
    {
        m->declare("author", "Jean Pierre Cimalando");
        m->declare("basics.lib/name", "Faust Basic Element Library");
        m->declare("basics.lib/version", "0.2");
        m->declare("compile_options",
                   "-a /home/jpc/.cache/faustpp/1278494-md.cpp -lang cpp -es 1 "
                   "-single -ftz 0");
        m->declare("delays.lib/name", "Faust Delay Library");
        m->declare("delays.lib/version", "0.1");
        m->declare("filename", "fverb.dsp");
        m->declare("filters.lib/allpass_comb:author", "Julius O. Smith III");
        m->declare("filters.lib/allpass_comb:copyright",
                   "Copyright (C) 2003-2019 by Julius O. Smith III "
                   "<jos@ccrma.stanford.edu>");
        m->declare("filters.lib/allpass_comb:license",
                   "MIT-style STK-4.3 license");
        m->declare("filters.lib/fir:author", "Julius O. Smith III");
        m->declare("filters.lib/fir:copyright",
                   "Copyright (C) 2003-2019 by Julius O. Smith III "
                   "<jos@ccrma.stanford.edu>");
        m->declare("filters.lib/fir:license", "MIT-style STK-4.3 license");
        m->declare("filters.lib/iir:author", "Julius O. Smith III");
        m->declare("filters.lib/iir:copyright",
                   "Copyright (C) 2003-2019 by Julius O. Smith III "
                   "<jos@ccrma.stanford.edu>");
        m->declare("filters.lib/iir:license", "MIT-style STK-4.3 license");
        m->declare("filters.lib/lowpass0_highpass1",
                   "Copyright (C) 2003-2019 by Julius O. Smith III "
                   "<jos@ccrma.stanford.edu>");
        m->declare("filters.lib/name", "Faust Filters Library");
        m->declare("filters.lib/version", "0.3");
        m->declare("license", "BSD-2-Clause");
        m->declare("maths.lib/author", "GRAME");
        m->declare("maths.lib/copyright", "GRAME");
        m->declare("maths.lib/license", "LGPL with exception");
        m->declare("maths.lib/name", "Faust Math Library");
        m->declare("maths.lib/version", "2.5");
        m->declare("name", "fverb");
        m->declare("oscillators.lib/name", "Faust Oscillator Library");
        m->declare("oscillators.lib/version", "0.1");
        m->declare("platform.lib/name", "Generic Platform Library");
        m->declare("platform.lib/version", "0.2");
        m->declare("signals.lib/name", "Faust Signal Routing Library");
        m->declare("signals.lib/version", "0.1");
        m->declare("version", "0.5");
    }

    FAUSTPP_VIRTUAL int getNumInputs()
    {
        return 2;
    }
    FAUSTPP_VIRTUAL int getNumOutputs()
    {
        return 2;
    }

    static void classInit(int sample_rate)
    {
        mydspSIG0 *sig0 = newmydspSIG0();
        sig0->instanceInitmydspSIG0(sample_rate);
        sig0->fillmydspSIG0(65536, ftbl0mydspSIG0);
        deletemydspSIG0(sig0);
    }

    FAUSTPP_VIRTUAL void instanceConstants(int sample_rate)
    {
        fSampleRate = sample_rate;
        fConst0 = std::min<DspFloatType>(192000.0f, std::max<DspFloatType>(1.0f, DspFloatType(fSampleRate)));
        fConst1 = (0.441000015f / fConst0);
        fConst2 = (44.0999985f / fConst0);
        fConst3 = (1.0f - fConst2);
        fConst4 = (0.0441000015f / fConst0);
        fConst5 = (1.0f / fConst0);
        iConst6 =
            std::min<int>(65536, std::max<int>(0, (int((0.00462820474f * fConst0)) + -1)));
        iConst7 =
            std::min<int>(65536, std::max<int>(0, (int((0.00370316859f * fConst0)) + -1)));
        iConst8 =
            std::min<int>(65536, std::max<int>(0, (int((0.013116831f * fConst0)) + -1)));
        iConst9 =
            std::min<int>(65536, std::max<int>(0, (int((0.00902825873f * fConst0)) + -1)));
        iConst10 =
            (std::min<int>(65536, std::max<int>(0, int((0.106280029f * fConst0)))) + 1);
        fConst11 = (1.0f / DspFloatType(int((0.00999999978f * fConst0))));
        fConst12 = (0.0f - fConst11);
        iConst13 = std::min<int>(65536, std::max<int>(0, int((0.141695514f * fConst0))));
        iConst14 =
            std::min<int>(65536, std::max<int>(0, (int((0.0892443135f * fConst0)) + -1)));
        iConst15 =
            std::min<int>(65536, std::max<int>(0, (int((0.00491448538f * fConst0)) + -1)));
        iConst16 =
            std::min<int>(65536, std::max<int>(0, (int((0.00348745007f * fConst0)) + -1)));
        iConst17 =
            std::min<int>(65536, std::max<int>(0, (int((0.0123527432f * fConst0)) + -1)));
        iConst18 =
            std::min<int>(65536, std::max<int>(0, (int((0.00958670769f * fConst0)) + -1)));
        iConst19 =
            (std::min<int>(65536, std::max<int>(0, int((0.124995798f * fConst0)))) + 1);
        iConst20 = std::min<int>(65536, std::max<int>(0, int((0.149625346f * fConst0))));
        iConst21 =
            std::min<int>(65536, std::max<int>(0, (int((0.0604818389f * fConst0)) + -1)));
        iConst22 =
            std::min<int>(65536, std::max<int>(0, int((0.00893787201f * fConst0))));
        iConst23 = std::min<int>(65536, std::max<int>(0, int((0.099929437f * fConst0))));
        iConst24 = std::min<int>(65536, std::max<int>(0, int((0.067067638f * fConst0))));
        iConst25 =
            std::min<int>(65536, std::max<int>(0, int((0.0642787516f * fConst0))));
        iConst26 =
            std::min<int>(65536, std::max<int>(0, int((0.0668660328f * fConst0))));
        iConst27 =
            std::min<int>(65536, std::max<int>(0, int((0.0062833908f * fConst0))));
        iConst28 =
            std::min<int>(65536, std::max<int>(0, int((0.0358186886f * fConst0))));
        iConst29 =
            std::min<int>(65536, std::max<int>(0, int((0.0118611604f * fConst0))));
        iConst30 = std::min<int>(65536, std::max<int>(0, int((0.121870905f * fConst0))));
        iConst31 =
            std::min<int>(65536, std::max<int>(0, int((0.0898155272f * fConst0))));
        iConst32 = std::min<int>(65536, std::max<int>(0, int((0.041262053f * fConst0))));
        iConst33 = std::min<int>(65536, std::max<int>(0, int((0.070931755f * fConst0))));
        iConst34 =
            std::min<int>(65536, std::max<int>(0, int((0.0112563418f * fConst0))));
        iConst35 =
            std::min<int>(65536, std::max<int>(0, int((0.00406572362f * fConst0))));
    }

    FAUSTPP_VIRTUAL void instanceResetUserInterface()
    {
        fHslider0 = FAUSTFLOAT(50.0f);
        fHslider1 = FAUSTFLOAT(100.0f);
        fHslider2 = FAUSTFLOAT(0.0f);
        fHslider3 = FAUSTFLOAT(10000.0f);
        fHslider4 = FAUSTFLOAT(100.0f);
        fHslider5 = FAUSTFLOAT(75.0f);
        fHslider6 = FAUSTFLOAT(62.5f);
        fHslider7 = FAUSTFLOAT(70.0f);
        fHslider8 = FAUSTFLOAT(0.5f);
        fHslider9 = FAUSTFLOAT(1.0f);
        fHslider10 = FAUSTFLOAT(5500.0f);
    }

    FAUSTPP_VIRTUAL void instanceClear()
    {
        for (int l0 = 0; (l0 < 2); l0 = (l0 + 1)) {
            fRec8[l0] = 0.0f;
        }
        for (int l1 = 0; (l1 < 2); l1 = (l1 + 1)) {
            fRec22[l1] = 0.0f;
        }
        IOTA = 0;
        for (int l2 = 0; (l2 < 131072); l2 = (l2 + 1)) {
            fVec0[l2] = 0.0f;
        }
        for (int l3 = 0; (l3 < 2); l3 = (l3 + 1)) {
            fRec23[l3] = 0.0f;
        }
        for (int l4 = 0; (l4 < 2); l4 = (l4 + 1)) {
            fRec24[l4] = 0.0f;
        }
        for (int l5 = 0; (l5 < 2); l5 = (l5 + 1)) {
            fRec21[l5] = 0.0f;
        }
        for (int l6 = 0; (l6 < 2); l6 = (l6 + 1)) {
            fRec25[l6] = 0.0f;
        }
        for (int l7 = 0; (l7 < 2); l7 = (l7 + 1)) {
            fRec20[l7] = 0.0f;
        }
        for (int l8 = 0; (l8 < 2); l8 = (l8 + 1)) {
            fRec26[l8] = 0.0f;
        }
        for (int l9 = 0; (l9 < 1024); l9 = (l9 + 1)) {
            fVec1[l9] = 0.0f;
        }
        for (int l10 = 0; (l10 < 2); l10 = (l10 + 1)) {
            fRec18[l10] = 0.0f;
        }
        for (int l11 = 0; (l11 < 1024); l11 = (l11 + 1)) {
            fVec2[l11] = 0.0f;
        }
        for (int l12 = 0; (l12 < 2); l12 = (l12 + 1)) {
            fRec16[l12] = 0.0f;
        }
        for (int l13 = 0; (l13 < 2); l13 = (l13 + 1)) {
            fRec27[l13] = 0.0f;
        }
        for (int l14 = 0; (l14 < 4096); l14 = (l14 + 1)) {
            fVec3[l14] = 0.0f;
        }
        for (int l15 = 0; (l15 < 2); l15 = (l15 + 1)) {
            fRec14[l15] = 0.0f;
        }
        for (int l16 = 0; (l16 < 2048); l16 = (l16 + 1)) {
            fVec4[l16] = 0.0f;
        }
        for (int l17 = 0; (l17 < 2); l17 = (l17 + 1)) {
            fRec12[l17] = 0.0f;
        }
        for (int l18 = 0; (l18 < 2); l18 = (l18 + 1)) {
            fRec28[l18] = 0.0f;
        }
        for (int l19 = 0; (l19 < 131072); l19 = (l19 + 1)) {
            fVec5[l19] = 0.0f;
        }
        for (int l20 = 0; (l20 < 2); l20 = (l20 + 1)) {
            fRec33[l20] = 0.0f;
        }
        for (int l23 = 0; (l23 < 2); l23 = (l23 + 1)) {
            fRec36[l23] = 0.0f;
        }
        for (int l24 = 0; (l24 < 2); l24 = (l24 + 1)) {
            fRec35[l24] = 0.0f;
        }
        for (int l25 = 0; (l25 < 2); l25 = (l25 + 1)) {
            fRec29[l25] = 0.0f;
        }
        for (int l26 = 0; (l26 < 2); l26 = (l26 + 1)) {
            fRec30[l26] = 0.0f;
        }
        for (int l27 = 0; (l27 < 2); l27 = (l27 + 1)) {
            iRec31[l27] = 0;
        }
        for (int l28 = 0; (l28 < 2); l28 = (l28 + 1)) {
            iRec32[l28] = 0;
        }
        for (int l29 = 0; (l29 < 2); l29 = (l29 + 1)) {
            fRec10[l29] = 0.0f;
        }
        for (int l30 = 0; (l30 < 32768); l30 = (l30 + 1)) {
            fVec7[l30] = 0.0f;
        }
        for (int l31 = 0; (l31 < 2); l31 = (l31 + 1)) {
            fRec37[l31] = 0.0f;
        }
        for (int l32 = 0; (l32 < 2); l32 = (l32 + 1)) {
            fRec9[l32] = 0.0f;
        }
        for (int l33 = 0; (l33 < 32768); l33 = (l33 + 1)) {
            fVec8[l33] = 0.0f;
        }
        for (int l34 = 0; (l34 < 2); l34 = (l34 + 1)) {
            fRec6[l34] = 0.0f;
        }
        for (int l35 = 0; (l35 < 32768); l35 = (l35 + 1)) {
            fRec0[l35] = 0.0f;
        }
        for (int l36 = 0; (l36 < 16384); l36 = (l36 + 1)) {
            fRec1[l36] = 0.0f;
        }
        for (int l37 = 0; (l37 < 32768); l37 = (l37 + 1)) {
            fRec2[l37] = 0.0f;
        }
        for (int l38 = 0; (l38 < 131072); l38 = (l38 + 1)) {
            fVec9[l38] = 0.0f;
        }
        for (int l39 = 0; (l39 < 2); l39 = (l39 + 1)) {
            fRec52[l39] = 0.0f;
        }
        for (int l40 = 0; (l40 < 2); l40 = (l40 + 1)) {
            fRec51[l40] = 0.0f;
        }
        for (int l41 = 0; (l41 < 1024); l41 = (l41 + 1)) {
            fVec10[l41] = 0.0f;
        }
        for (int l42 = 0; (l42 < 2); l42 = (l42 + 1)) {
            fRec49[l42] = 0.0f;
        }
        for (int l43 = 0; (l43 < 1024); l43 = (l43 + 1)) {
            fVec11[l43] = 0.0f;
        }
        for (int l44 = 0; (l44 < 2); l44 = (l44 + 1)) {
            fRec47[l44] = 0.0f;
        }
        for (int l45 = 0; (l45 < 4096); l45 = (l45 + 1)) {
            fVec12[l45] = 0.0f;
        }
        for (int l46 = 0; (l46 < 2); l46 = (l46 + 1)) {
            fRec45[l46] = 0.0f;
        }
        for (int l47 = 0; (l47 < 2048); l47 = (l47 + 1)) {
            fVec13[l47] = 0.0f;
        }
        for (int l48 = 0; (l48 < 2); l48 = (l48 + 1)) {
            fRec43[l48] = 0.0f;
        }
        for (int l49 = 0; (l49 < 131072); l49 = (l49 + 1)) {
            fVec14[l49] = 0.0f;
        }
        for (int l50 = 0; (l50 < 2); l50 = (l50 + 1)) {
            fRec53[l50] = 0.0f;
        }
        for (int l51 = 0; (l51 < 2); l51 = (l51 + 1)) {
            fRec54[l51] = 0.0f;
        }
        for (int l52 = 0; (l52 < 2); l52 = (l52 + 1)) {
            iRec55[l52] = 0;
        }
        for (int l53 = 0; (l53 < 2); l53 = (l53 + 1)) {
            iRec56[l53] = 0;
        }
        for (int l54 = 0; (l54 < 2); l54 = (l54 + 1)) {
            fRec41[l54] = 0.0f;
        }
        for (int l55 = 0; (l55 < 32768); l55 = (l55 + 1)) {
            fVec15[l55] = 0.0f;
        }
        for (int l56 = 0; (l56 < 2); l56 = (l56 + 1)) {
            fRec40[l56] = 0.0f;
        }
        for (int l57 = 0; (l57 < 16384); l57 = (l57 + 1)) {
            fVec16[l57] = 0.0f;
        }
        for (int l58 = 0; (l58 < 2); l58 = (l58 + 1)) {
            fRec38[l58] = 0.0f;
        }
        for (int l59 = 0; (l59 < 32768); l59 = (l59 + 1)) {
            fRec3[l59] = 0.0f;
        }
        for (int l60 = 0; (l60 < 8192); l60 = (l60 + 1)) {
            fRec4[l60] = 0.0f;
        }
        for (int l61 = 0; (l61 < 32768); l61 = (l61 + 1)) {
            fRec5[l61] = 0.0f;
        }
    }

    FAUSTPP_VIRTUAL void init(int sample_rate)
    {
        classInit(sample_rate);
        instanceInit(sample_rate);
    }
    FAUSTPP_VIRTUAL void instanceInit(int sample_rate)
    {
        instanceConstants(sample_rate);
        instanceResetUserInterface();
        instanceClear();
    }

    FAUSTPP_VIRTUAL mydsp *clone()
    {
        return new mydsp();
    }

    FAUSTPP_VIRTUAL int getSampleRate()
    {
        return fSampleRate;
    }

    FAUSTPP_VIRTUAL void buildUserInterface(UI *ui_interface)
    {
        ui_interface->openVerticalBox("fverb");
        ui_interface->declare(&fHslider2, "01", "");
        ui_interface->declare(&fHslider2, "symbol", "predelay");
        ui_interface->declare(&fHslider2, "unit", "ms");
        ui_interface->addHorizontalSlider("Predelay", &fHslider2,
                                          FAUSTFLOAT(0.0f), FAUSTFLOAT(0.0f),
                                          FAUSTFLOAT(300.0f), FAUSTFLOAT(1.0f));
        ui_interface->declare(&fHslider1, "02", "");
        ui_interface->declare(&fHslider1, "symbol", "input");
        ui_interface->declare(&fHslider1, "unit", "%");
        ui_interface->addHorizontalSlider("Input amount", &fHslider1, FAUSTFLOAT(100.0f),
                                          FAUSTFLOAT(0.0f), FAUSTFLOAT(100.0f),
                                          FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider3, "03", "");
        ui_interface->declare(&fHslider3, "scale", "log");
        ui_interface->declare(&fHslider3, "symbol", "input_lowpass");
        ui_interface->declare(&fHslider3, "unit", "Hz");
        ui_interface->addHorizontalSlider("Input low-pass cutoff", &fHslider3,
                                          FAUSTFLOAT(10000.0f), FAUSTFLOAT(1.0f),
                                          FAUSTFLOAT(20000.0f), FAUSTFLOAT(1.0f));
        ui_interface->declare(&fHslider4, "04", "");
        ui_interface->declare(&fHslider4, "scale", "log");
        ui_interface->declare(&fHslider4, "symbol", "input_highpass");
        ui_interface->declare(&fHslider4, "unit", "Hz");
        ui_interface->addHorizontalSlider("Input high-pass cutoff", &fHslider4,
                                          FAUSTFLOAT(100.0f), FAUSTFLOAT(1.0f),
                                          FAUSTFLOAT(1000.0f), FAUSTFLOAT(1.0f));
        ui_interface->declare(&fHslider5, "05", "");
        ui_interface->declare(&fHslider5, "symbol", "input_diffusion_1");
        ui_interface->declare(&fHslider5, "unit", "%");
        ui_interface->addHorizontalSlider("Input diffusion 1", &fHslider5,
                                          FAUSTFLOAT(75.0f), FAUSTFLOAT(0.0f),
                                          FAUSTFLOAT(100.0f), FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider6, "06", "");
        ui_interface->declare(&fHslider6, "symbol", "input_diffusion_2");
        ui_interface->declare(&fHslider6, "unit", "%");
        ui_interface->addHorizontalSlider("Input diffusion 2", &fHslider6,
                                          FAUSTFLOAT(62.5f), FAUSTFLOAT(0.0f),
                                          FAUSTFLOAT(100.0f), FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider7, "07", "");
        ui_interface->declare(&fHslider7, "symbol", "tail_density");
        ui_interface->declare(&fHslider7, "unit", "%");
        ui_interface->addHorizontalSlider("Tail density", &fHslider7, FAUSTFLOAT(70.0f),
                                          FAUSTFLOAT(0.0f), FAUSTFLOAT(100.0f),
                                          FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider0, "08", "");
        ui_interface->declare(&fHslider0, "symbol", "decay");
        ui_interface->declare(&fHslider0, "unit", "%");
        ui_interface->addHorizontalSlider("Decay", &fHslider0, FAUSTFLOAT(50.0f),
                                          FAUSTFLOAT(0.0f), FAUSTFLOAT(100.0f),
                                          FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider10, "09", "");
        ui_interface->declare(&fHslider10, "scale", "log");
        ui_interface->declare(&fHslider10, "symbol", "damping");
        ui_interface->declare(&fHslider10, "unit", "Hz");
        ui_interface->addHorizontalSlider("Damping", &fHslider10,
                                          FAUSTFLOAT(5500.0f), FAUSTFLOAT(10.0f),
                                          FAUSTFLOAT(20000.0f), FAUSTFLOAT(1.0f));
        ui_interface->declare(&fHslider9, "10", "");
        ui_interface->declare(&fHslider9, "symbol", "mod_frequency");
        ui_interface->declare(&fHslider9, "unit", "Hz");
        ui_interface->addHorizontalSlider("Modulator frequency", &fHslider9,
                                          FAUSTFLOAT(1.0f), FAUSTFLOAT(0.00999999978f),
                                          FAUSTFLOAT(4.0f), FAUSTFLOAT(0.00999999978f));
        ui_interface->declare(&fHslider8, "11", "");
        ui_interface->declare(&fHslider8, "symbol", "mod_depth");
        ui_interface->declare(&fHslider8, "unit", "ms");
        ui_interface->addHorizontalSlider("Modulator depth", &fHslider8,
                                          FAUSTFLOAT(0.5f), FAUSTFLOAT(0.0f),
                                          FAUSTFLOAT(10.0f), FAUSTFLOAT(0.100000001f));
        ui_interface->closeBox();
    }

    FAUSTPP_VIRTUAL void compute(int count, FAUSTFLOAT **inputs, FAUSTFLOAT **outputs)
    {
        FAUSTFLOAT *input0 = inputs[0];
        FAUSTFLOAT *input1 = inputs[1];
        FAUSTFLOAT *output0 = outputs[0];
        FAUSTFLOAT *output1 = outputs[1];
        DspFloatType fSlow0 = (fConst1 * DspFloatType(fHslider0));
        DspFloatType fSlow1 = (fConst1 * DspFloatType(fHslider1));
        DspFloatType fSlow2 = (fConst4 * DspFloatType(fHslider2));
        DspFloatType fSlow3 =
            (fConst2 * std::exp((fConst5 * (0.0f - (6.28318548f * DspFloatType(fHslider3))))));
        DspFloatType fSlow4 =
            (fConst2 * std::exp((fConst5 * (0.0f - (6.28318548f * DspFloatType(fHslider4))))));
        DspFloatType fSlow5 = (fConst1 * DspFloatType(fHslider5));
        DspFloatType fSlow6 = (fConst1 * DspFloatType(fHslider6));
        DspFloatType fSlow7 = (fConst1 * DspFloatType(fHslider7));
        DspFloatType fSlow8 = (fConst4 * DspFloatType(fHslider8));
        DspFloatType fSlow9 = (fConst2 * DspFloatType(fHslider9));
        DspFloatType fSlow10 =
            (fConst2 * std::exp((fConst5 * (0.0f - (6.28318548f * DspFloatType(fHslider10))))));
        for (int i0 = 0; (i0 < count); i0 = (i0 + 1)) {
            fRec8[0] = (fSlow0 + (fConst3 * fRec8[1]));
            DspFloatType fTemp0 =
                std::min<DspFloatType>(0.5f, std::max<DspFloatType>(0.25f, (fRec8[0] + 0.150000006f)));
            fRec22[0] = (fSlow1 + (fConst3 * fRec22[1]));
            fVec0[(IOTA & 131071)] = (DspFloatType(input1[i0]) * fRec22[0]);
            fRec23[0] = (fSlow2 + (fConst3 * fRec23[1]));
            int iTemp1 =
                std::min<int>(65536, std::max<int>(0, int((fConst0 * fRec23[0]))));
            fRec24[0] = (fSlow3 + (fConst3 * fRec24[1]));
            fRec21[0] = (fVec0[((IOTA - iTemp1) & 131071)] + (fRec24[0] * fRec21[1]));
            DspFloatType fTemp2 = (1.0f - fRec24[0]);
            fRec25[0] = (fSlow4 + (fConst3 * fRec25[1]));
            fRec20[0] = ((fRec21[0] * fTemp2) + (fRec25[0] * fRec20[1]));
            DspFloatType fTemp3 = (fRec25[0] + 1.0f);
            DspFloatType fTemp4 = (0.0f - (0.5f * fTemp3));
            fRec26[0] = (fSlow5 + (fConst3 * fRec26[1]));
            DspFloatType fTemp5 = (((0.5f * (fRec20[0] * fTemp3)) + (fRec20[1] * fTemp4)) -
                            (fRec26[0] * fRec18[1]));
            fVec1[(IOTA & 1023)] = fTemp5;
            fRec18[0] = fVec1[((IOTA - iConst6) & 1023)];
            DspFloatType fRec19 = (fRec26[0] * fTemp5);
            DspFloatType fTemp6 = ((fRec19 + fRec18[1]) - (fRec26[0] * fRec16[1]));
            fVec2[(IOTA & 1023)] = fTemp6;
            fRec16[0] = fVec2[((IOTA - iConst7) & 1023)];
            DspFloatType fRec17 = (fRec26[0] * fTemp6);
            fRec27[0] = (fSlow6 + (fConst3 * fRec27[1]));
            DspFloatType fTemp7 = ((fRec17 + fRec16[1]) - (fRec27[0] * fRec14[1]));
            fVec3[(IOTA & 4095)] = fTemp7;
            fRec14[0] = fVec3[((IOTA - iConst8) & 4095)];
            DspFloatType fRec15 = (fRec27[0] * fTemp7);
            DspFloatType fTemp8 = ((fRec15 + fRec14[1]) - (fRec27[0] * fRec12[1]));
            fVec4[(IOTA & 2047)] = fTemp8;
            fRec12[0] = fVec4[((IOTA - iConst9) & 2047)];
            DspFloatType fRec13 = (fRec27[0] * fTemp8);
            fRec28[0] = (fSlow7 + (fConst3 * fRec28[1]));
            DspFloatType fTemp9 = (fRec12[1] + ((fRec8[0] * fRec3[((IOTA - iConst10) & 32767)]) +
                                         (fRec13 + (fRec28[0] * fRec10[1]))));
            fVec5[(IOTA & 131071)] = fTemp9;
            fRec33[0] = (fSlow8 + (fConst3 * fRec33[1]));
            fRec36[0] = (fSlow9 + (fConst3 * fRec36[1]));
            DspFloatType fTemp10 = (fRec35[1] + (fConst5 * fRec36[0]));
            fRec35[0] = (fTemp10 - DspFloatType(int(fTemp10)));
            int iTemp11 =
                (int((fConst0 *
                      ((fRec33[0] *
                        ftbl0mydspSIG0[int((65536.0f * (fRec35[0] + (0.25f - DspFloatType(int((fRec35[0] + 0.25f)))))))]) +
                       0.0305097271f))) +
                 -1);
            DspFloatType fThen1 =
                (((fRec30[1] == 1.0f) & (iTemp11 != iRec32[1])) ? fConst12 : 0.0f);
            DspFloatType fThen3 =
                (((fRec30[1] == 0.0f) & (iTemp11 != iRec31[1])) ? fConst11 : fThen1);
            DspFloatType fElse3 = (((fRec30[1] > 0.0f) & (fRec30[1] < 1.0f)) ? fRec29[1] : 0.0f);
            DspFloatType fTemp12 = ((fRec29[1] != 0.0f) ? fElse3 : fThen3);
            fRec29[0] = fTemp12;
            fRec30[0] =
                std::max<DspFloatType>(0.0f, std::min<DspFloatType>(1.0f, (fRec30[1] + fTemp12)));
            iRec31[0] =
                (((fRec30[1] >= 1.0f) & (iRec32[1] != iTemp11)) ? iTemp11 : iRec31[1]);
            iRec32[0] =
                (((fRec30[1] <= 0.0f) & (iRec31[1] != iTemp11)) ? iTemp11 : iRec32[1]);
            DspFloatType fTemp13 =
                fVec5[((IOTA - std::min<int>(65536, std::max<int>(0, iRec31[0]))) & 131071)];
            fRec10[0] =
                (fTemp13 +
                 (fRec30[0] *
                  (fVec5[((IOTA - std::min<int>(65536, std::max<int>(0, iRec32[0]))) & 131071)] -
                   fTemp13)));
            DspFloatType fRec11 = (0.0f - (fRec28[0] * fTemp9));
            DspFloatType fTemp14 = (fRec11 + fRec10[1]);
            fVec7[(IOTA & 32767)] = fTemp14;
            fRec37[0] = (fSlow10 + (fConst3 * fRec37[1]));
            fRec9[0] = (fVec7[((IOTA - iConst13) & 32767)] + (fRec37[0] * fRec9[1]));
            DspFloatType fTemp15 = (1.0f - fRec37[0]);
            DspFloatType fTemp16 = ((fTemp0 * fRec6[1]) + ((fRec8[0] * fRec9[0]) * fTemp15));
            fVec8[(IOTA & 32767)] = fTemp16;
            fRec6[0] = fVec8[((IOTA - iConst14) & 32767)];
            DspFloatType fRec7 = (0.0f - (fTemp0 * fTemp16));
            fRec0[(IOTA & 32767)] = (fRec7 + fRec6[1]);
            fRec1[(IOTA & 16383)] = (fRec9[0] * fTemp15);
            fRec2[(IOTA & 32767)] = fTemp14;
            fVec9[(IOTA & 131071)] = (DspFloatType(input0[i0]) * fRec22[0]);
            fRec52[0] = (fVec9[((IOTA - iTemp1) & 131071)] + (fRec24[0] * fRec52[1]));
            fRec51[0] = ((fTemp2 * fRec52[0]) + (fRec25[0] * fRec51[1]));
            DspFloatType fTemp17 = (((0.5f * (fRec51[0] * fTemp3)) + (fTemp4 * fRec51[1])) -
                             (fRec26[0] * fRec49[1]));
            fVec10[(IOTA & 1023)] = fTemp17;
            fRec49[0] = fVec10[((IOTA - iConst15) & 1023)];
            DspFloatType fRec50 = (fRec26[0] * fTemp17);
            DspFloatType fTemp18 = ((fRec50 + fRec49[1]) - (fRec26[0] * fRec47[1]));
            fVec11[(IOTA & 1023)] = fTemp18;
            fRec47[0] = fVec11[((IOTA - iConst16) & 1023)];
            DspFloatType fRec48 = (fRec26[0] * fTemp18);
            DspFloatType fTemp19 = ((fRec48 + fRec47[1]) - (fRec27[0] * fRec45[1]));
            fVec12[(IOTA & 4095)] = fTemp19;
            fRec45[0] = fVec12[((IOTA - iConst17) & 4095)];
            DspFloatType fRec46 = (fRec27[0] * fTemp19);
            DspFloatType fTemp20 = ((fRec46 + fRec45[1]) - (fRec27[0] * fRec43[1]));
            fVec13[(IOTA & 2047)] = fTemp20;
            fRec43[0] = fVec13[((IOTA - iConst18) & 2047)];
            DspFloatType fRec44 = (fRec27[0] * fTemp20);
            DspFloatType fTemp21 =
                (fRec43[1] + ((fRec8[0] * fRec0[((IOTA - iConst19) & 32767)]) +
                              (fRec44 + (fRec28[0] * fRec41[1]))));
            fVec14[(IOTA & 131071)] = fTemp21;
            int iTemp22 =
                (int((fConst0 * ((fRec33[0] * ftbl0mydspSIG0[int((65536.0f * fRec35[0]))]) +
                                 0.025603978f))) +
                 -1);
            DspFloatType fThen7 =
                (((fRec54[1] == 1.0f) & (iTemp22 != iRec56[1])) ? fConst12 : 0.0f);
            DspFloatType fThen9 =
                (((fRec54[1] == 0.0f) & (iTemp22 != iRec55[1])) ? fConst11 : fThen7);
            DspFloatType fElse9 = (((fRec54[1] > 0.0f) & (fRec54[1] < 1.0f)) ? fRec53[1] : 0.0f);
            DspFloatType fTemp23 = ((fRec53[1] != 0.0f) ? fElse9 : fThen9);
            fRec53[0] = fTemp23;
            fRec54[0] =
                std::max<DspFloatType>(0.0f, std::min<DspFloatType>(1.0f, (fRec54[1] + fTemp23)));
            iRec55[0] =
                (((fRec54[1] >= 1.0f) & (iRec56[1] != iTemp22)) ? iTemp22 : iRec55[1]);
            iRec56[0] =
                (((fRec54[1] <= 0.0f) & (iRec55[1] != iTemp22)) ? iTemp22 : iRec56[1]);
            DspFloatType fTemp24 =
                fVec14[((IOTA - std::min<int>(65536, std::max<int>(0, iRec55[0]))) & 131071)];
            fRec41[0] =
                (fTemp24 +
                 (fRec54[0] *
                  (fVec14[((IOTA - std::min<int>(65536, std::max<int>(0, iRec56[0]))) & 131071)] -
                   fTemp24)));
            DspFloatType fRec42 = (0.0f - (fRec28[0] * fTemp21));
            DspFloatType fTemp25 = (fRec42 + fRec41[1]);
            fVec15[(IOTA & 32767)] = fTemp25;
            fRec40[0] = (fVec15[((IOTA - iConst20) & 32767)] + (fRec37[0] * fRec40[1]));
            DspFloatType fTemp26 = ((fTemp0 * fRec38[1]) + ((fRec8[0] * fTemp15) * fRec40[0]));
            fVec16[(IOTA & 16383)] = fTemp26;
            fRec38[0] = fVec16[((IOTA - iConst21) & 16383)];
            DspFloatType fRec39 = (0.0f - (fTemp0 * fTemp26));
            fRec3[(IOTA & 32767)] = (fRec39 + fRec38[1]);
            fRec4[(IOTA & 8191)] = (fTemp15 * fRec40[0]);
            fRec5[(IOTA & 32767)] = fTemp25;
            output0[i0] = FAUSTFLOAT(
                (0.600000024f *
                 (((fRec2[((IOTA - iConst22) & 32767)] + fRec2[((IOTA - iConst23) & 32767)]) +
                   fRec0[((IOTA - iConst24) & 32767)]) -
                  (((fRec1[((IOTA - iConst25) & 16383)] + fRec5[((IOTA - iConst26) & 32767)]) +
                    fRec4[((IOTA - iConst27) & 8191)]) +
                   fRec3[((IOTA - iConst28) & 32767)]))));
            output1[i0] = FAUSTFLOAT(
                (0.600000024f *
                 (((fRec5[((IOTA - iConst29) & 32767)] + fRec5[((IOTA - iConst30) & 32767)]) +
                   fRec3[((IOTA - iConst31) & 32767)]) -
                  (((fRec4[((IOTA - iConst32) & 8191)] + fRec2[((IOTA - iConst33) & 32767)]) +
                    fRec1[((IOTA - iConst34) & 16383)]) +
                   fRec0[((IOTA - iConst35) & 32767)]))));
            fRec8[1] = fRec8[0];
            fRec22[1] = fRec22[0];
            IOTA = (IOTA + 1);
            fRec23[1] = fRec23[0];
            fRec24[1] = fRec24[0];
            fRec21[1] = fRec21[0];
            fRec25[1] = fRec25[0];
            fRec20[1] = fRec20[0];
            fRec26[1] = fRec26[0];
            fRec18[1] = fRec18[0];
            fRec16[1] = fRec16[0];
            fRec27[1] = fRec27[0];
            fRec14[1] = fRec14[0];
            fRec12[1] = fRec12[0];
            fRec28[1] = fRec28[0];
            fRec33[1] = fRec33[0];
            fRec36[1] = fRec36[0];
            fRec35[1] = fRec35[0];
            fRec29[1] = fRec29[0];
            fRec30[1] = fRec30[0];
            iRec31[1] = iRec31[0];
            iRec32[1] = iRec32[0];
            fRec10[1] = fRec10[0];
            fRec37[1] = fRec37[0];
            fRec9[1] = fRec9[0];
            fRec6[1] = fRec6[0];
            fRec52[1] = fRec52[0];
            fRec51[1] = fRec51[0];
            fRec49[1] = fRec49[0];
            fRec47[1] = fRec47[0];
            fRec45[1] = fRec45[0];
            fRec43[1] = fRec43[0];
            fRec53[1] = fRec53[0];
            fRec54[1] = fRec54[0];
            iRec55[1] = iRec55[0];
            iRec56[1] = iRec56[0];
            fRec41[1] = fRec41[0];
            fRec40[1] = fRec40[0];
            fRec38[1] = fRec38[0];
        }
    }
};
FAUSTPP_END_NAMESPACE

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------
// End the Faust code section


inline Fverb::Fverb()
{

    mydsp *dsp = new mydsp;
    fDsp.reset(dsp);
    dsp->instanceResetUserInterface();
}

inline Fverb::~Fverb()
{
}

inline void Fverb::init(DspFloatType sample_rate)
{

    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.classInit(sample_rate);
    dsp.instanceConstants(sample_rate);
    clear();
}

inline void Fverb::clear() noexcept
{

    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.instanceClear();
}

inline void Fverb::process(const DspFloatType *in0, const DspFloatType *in1, DspFloatType *out0,
                    DspFloatType *out1, unsigned count) noexcept
{

    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    DspFloatType *inputs[] = {
        const_cast<DspFloatType *>(in0),
        const_cast<DspFloatType *>(in1),
    };
    DspFloatType *outputs[] = {
        out0,
        out1,
    };
    dsp.compute(count, inputs, outputs);
}

inline const char *Fverb::parameter_label(unsigned index) noexcept
{
    switch (index) {

    case 0:
        return "Predelay";

    case 1:
        return "Input amount";

    case 2:
        return "Input low-pass cutoff";

    case 3:
        return "Input high-pass cutoff";

    case 4:
        return "Input diffusion 1";

    case 5:
        return "Input diffusion 2";

    case 6:
        return "Tail density";

    case 7:
        return "Decay";

    case 8:
        return "Damping";

    case 9:
        return "Modulator frequency";

    case 10:
        return "Modulator depth";

    default:
        return 0;
    }
}

inline const char *Fverb::parameter_short_label(unsigned index) noexcept
{
    switch (index) {

    case 0:
        return "";

    case 1:
        return "";

    case 2:
        return "";

    case 3:
        return "";

    case 4:
        return "";

    case 5:
        return "";

    case 6:
        return "";

    case 7:
        return "";

    case 8:
        return "";

    case 9:
        return "";

    case 10:
        return "";

    default:
        return 0;
    }
}

inline const char *Fverb::parameter_symbol(unsigned index) noexcept
{
    switch (index) {

    case 0:
        return "predelay";

    case 1:
        return "input";

    case 2:
        return "input_lowpass";

    case 3:
        return "input_highpass";

    case 4:
        return "input_diffusion_1";

    case 5:
        return "input_diffusion_2";

    case 6:
        return "tail_density";

    case 7:
        return "decay";

    case 8:
        return "damping";

    case 9:
        return "mod_frequency";

    case 10:
        return "mod_depth";

    default:
        return 0;
    }
}

inline const char *Fverb::parameter_unit(unsigned index) noexcept
{
    switch (index) {

    case 0:
        return "ms";

    case 1:
        return "%";

    case 2:
        return "Hz";

    case 3:
        return "Hz";

    case 4:
        return "%";

    case 5:
        return "%";

    case 6:
        return "%";

    case 7:
        return "%";

    case 8:
        return "Hz";

    case 9:
        return "Hz";

    case 10:
        return "ms";

    default:
        return 0;
    }
}

inline const Fverb::ParameterRange *Fverb::parameter_range(unsigned index) noexcept
{
    switch (index) {

    case 0: {
        static const ParameterRange range = {0, 0, 300};
        return &range;
    }

    case 1: {
        static const ParameterRange range = {100, 0, 100};
        return &range;
    }

    case 2: {
        static const ParameterRange range = {10000, 1, 20000};
        return &range;
    }

    case 3: {
        static const ParameterRange range = {100, 1, 1000};
        return &range;
    }

    case 4: {
        static const ParameterRange range = {75, 0, 100};
        return &range;
    }

    case 5: {
        static const ParameterRange range = {62.5, 0, 100};
        return &range;
    }

    case 6: {
        static const ParameterRange range = {70, 0, 100};
        return &range;
    }

    case 7: {
        static const ParameterRange range = {50, 0, 100};
        return &range;
    }

    case 8: {
        static const ParameterRange range = {5500, 10, 20000};
        return &range;
    }

    case 9: {
        static const ParameterRange range = {1, 0.0099999998, 4};
        return &range;
    }

    case 10: {
        static const ParameterRange range = {0.5, 0, 10};
        return &range;
    }

    default:
        return 0;
    }
}

inline bool Fverb::parameter_is_trigger(unsigned index) noexcept
{
    switch (index) {

    default:
        return false;
    }
}

inline bool Fverb::parameter_is_boolean(unsigned index) noexcept
{
    switch (index) {

    default:
        return false;
    }
}

inline bool Fverb::parameter_is_integer(unsigned index) noexcept
{
    switch (index) {

    default:
        return false;
    }
}

inline bool Fverb::parameter_is_logarithmic(unsigned index) noexcept
{
    switch (index) {

    case 2:
        return true;

    case 3:
        return true;

    case 8:
        return true;

    default:
        return false;
    }
}

inline DspFloatType Fverb::get_parameter(unsigned index) const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    switch (index) {

    case 0:
        return dsp.fHslider2;

    case 1:
        return dsp.fHslider1;

    case 2:
        return dsp.fHslider3;

    case 3:
        return dsp.fHslider4;

    case 4:
        return dsp.fHslider5;

    case 5:
        return dsp.fHslider6;

    case 6:
        return dsp.fHslider7;

    case 7:
        return dsp.fHslider0;

    case 8:
        return dsp.fHslider10;

    case 9:
        return dsp.fHslider9;

    case 10:
        return dsp.fHslider8;

    default:
        (void)dsp;
        return 0;
    }
}

inline void Fverb::set_parameter(unsigned index, DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    switch (index) {

    case 0:
        dsp.fHslider2 = value;
        break;

    case 1:
        dsp.fHslider1 = value;
        break;

    case 2:
        dsp.fHslider3 = value;
        break;

    case 3:
        dsp.fHslider4 = value;
        break;

    case 4:
        dsp.fHslider5 = value;
        break;

    case 5:
        dsp.fHslider6 = value;
        break;

    case 6:
        dsp.fHslider7 = value;
        break;

    case 7:
        dsp.fHslider0 = value;
        break;

    case 8:
        dsp.fHslider10 = value;
        break;

    case 9:
        dsp.fHslider9 = value;
        break;

    case 10:
        dsp.fHslider8 = value;
        break;

    default:
        (void)dsp;
        (void)value;
        break;
    }
}


inline DspFloatType Fverb::get_predelay() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider2;
}

inline DspFloatType Fverb::get_input() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider1;
}

inline DspFloatType Fverb::get_input_lowpass() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider3;
}

inline DspFloatType Fverb::get_input_highpass() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider4;
}

inline DspFloatType Fverb::get_input_diffusion_1() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider5;
}

inline DspFloatType Fverb::get_input_diffusion_2() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider6;
}

inline DspFloatType Fverb::get_tail_density() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider7;
}

inline DspFloatType Fverb::get_decay() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider0;
}

inline DspFloatType Fverb::get_damping() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider10;
}

inline DspFloatType Fverb::get_mod_frequency() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider9;
}

inline DspFloatType Fverb::get_mod_depth() const noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    return dsp.fHslider8;
}


inline void Fverb::set_predelay(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider2 = value;
}

inline void Fverb::set_input(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider1 = value;
}

inline void Fverb::set_input_lowpass(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider3 = value;
}

inline void Fverb::set_input_highpass(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider4 = value;
}

inline void Fverb::set_input_diffusion_1(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider5 = value;
}

inline void Fverb::set_input_diffusion_2(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider6 = value;
}

inline void Fverb::set_tail_density(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider7 = value;
}

inline void Fverb::set_decay(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider0 = value;
}

inline void Fverb::set_damping(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider10 = value;
}

inline void Fverb::set_mod_frequency(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider9 = value;
}

inline void Fverb::set_mod_depth(DspFloatType value) noexcept
{
    mydsp &dsp = static_cast<mydsp &>(*fDsp);
    dsp.fHslider8 = value;
}