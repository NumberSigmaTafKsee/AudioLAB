#pragma once

#include <memory>
#include <utility>
#include <cmath>
#include <algorithm>
#include <math.h>
#include <vector>



namespace SoundAlchemy::MoogLadder
{

    template<typename T> class mydsp;

    template<typename T>
    struct ParameterRange {
        T init;
        T min;
        T max;
    };

    template<class T>
    class TMoogLadder {
    public:
        TMoogLadder();
        ~TMoogLadder();

        void init(T sample_rate);
        void clear() noexcept;

        void process(
            const T *in0,
            T *out0,
            unsigned count) noexcept;

        enum { NumInputs = 1 };
        enum { NumOutputs = 1 };
        enum { NumActives = 2 };
        enum { NumPassives = 0 };
        enum { NumParameters = 2 };

        enum Parameter {
            p_cutoff,
            p_q,
            
        };

        
        static const char *parameter_label(unsigned index) noexcept;
        static const char *parameter_short_label(unsigned index) noexcept;
        static const char *parameter_symbol(unsigned index) noexcept;
        static const char *parameter_unit(unsigned index) noexcept;
        static const ParameterRange<T> *parameter_range(unsigned index) noexcept;
        static bool parameter_is_trigger(unsigned index) noexcept;
        static bool parameter_is_boolean(unsigned index) noexcept;
        static bool parameter_is_integer(unsigned index) noexcept;
        static bool parameter_is_logarithmic(unsigned index) noexcept;

        T get_parameter(unsigned index) const noexcept;
        void set_parameter(unsigned index, T value) noexcept;

        T get_cutoff() const noexcept;        
        T get_q() const noexcept;
                
        void set_cutoff(T value) noexcept;        
        void set_q(T value) noexcept;

        void Process(const size_t n, const std::vector<T> & in, std::vector<T> & out) {
            process(in.data(),out.data(),n);
        }

        T Tick(T input);

        mydsp<T> * fDsp;
    };

    
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


    #define FAUSTPP_VIRTUAL // do not declare any methods virtual
    #define FAUSTPP_PRIVATE public // do not hide any members
    #define FAUSTPP_PROTECTED public // do not hide any members

    // define the DSP in the anonymous namespace
    #define FAUSTPP_BEGIN_NAMESPACE namespace {
    #define FAUSTPP_END_NAMESPACE }


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

    #ifndef FAUSTFLOAT
    #define FAUSTFLOAT float
    #endif 


    template<typename T>
    static T mydsp_faustpower2_f(T value) {
        return (value * value);
    }
    template<typename T>
    static T mydsp_faustpower3_f(T value) {
        return ((value * value) * value);
    }
    template<typename T>
    static T mydsp_faustpower4_f(T value) {
        return (((value * value) * value) * value);
    }

    #ifndef FAUSTCLASS 
    #define FAUSTCLASS mydsp<T>
    #endif

    #ifdef __APPLE__ 
    #define exp10f __exp10f
    #define exp10 __exp10
    #endif

    template<typename T>
    class mydsp {
        
    FAUSTPP_PRIVATE:
        
        int fSampleRate;
        T fConst0;
        T fHslider0;
        T fRec5[2];
        T fHslider1;
        T fRec1[2];
        T fRec2[2];
        T fRec3[2];
        T fRec4[2];
        
    public:
        
        void metadata(Meta* m) { 
            m->declare("../../faust/moogladder.dsp/moogLadder:author", "Eric Tarr");
            m->declare("../../faust/moogladder.dsp/moogLadder:license", "MIT-style STK-4.3 license");
            m->declare("author", "Christopher Arndt");
            m->declare("description", "FAUST Moog Ladder 24 dB LPF");
            m->declare("filename", "moogladder.dsp");
            m->declare("license", "MIT-style STK-4.3 license");
            m->declare("maths.lib/author", "GRAME");
            m->declare("maths.lib/copyright", "GRAME");
            m->declare("maths.lib/license", "LGPL with exception");
            m->declare("maths.lib/name", "Faust Math Library");
            m->declare("maths.lib/version", "2.3");
            m->declare("name", "TMoogLadder");
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
            fConst0 = (3.14159274f / std::min<T>(192000.0f, std::max<T>(1.0f, T(fSampleRate))));
        }
        
        FAUSTPP_VIRTUAL void instanceResetUserInterface() {
            fHslider0 = T(20000.0f);
            fHslider1 = T(1.0f);
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
        
        FAUSTPP_VIRTUAL mydsp<T>* clone() {
            return new mydsp<T>();
        }
        
        FAUSTPP_VIRTUAL int getSampleRate() {
            return fSampleRate;
        }
        
        FAUSTPP_VIRTUAL void buildUserInterface(UI* ui_interface) {
            ui_interface->openVerticalBox("TMoogLadder");
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
        
        FAUSTPP_VIRTUAL void compute(int count, T** inputs, T** outputs) {
            T* input0 = inputs[0];
            T* output0 = outputs[0];
            T fSlow0 = (0.00100000005f * T(fHslider0));
            T fSlow1 = (0.0411641225f * (T(fHslider1) + -0.707000017f));
            for (int i = 0; (i < count); i = (i + 1)) {
                fRec5[0] = (fSlow0 + (0.999000013f * fRec5[1]));
                T fTemp0 = std::tan((fConst0 * fRec5[0]));
                T fTemp1 = (3.9000001f - (0.899999976f * std::pow((0.0f - (0.333333343f * (1.0f - (std::log10(fRec5[0]) + -0.30103001f)))), 0.200000003f)));
                T fTemp2 = (fTemp0 + 1.0f);
                T fTemp3 = ((fTemp0 * (((T(input0[i]) - (fSlow1 * (fTemp1 * (((fRec1[1] + (fRec2[1] * fTemp0)) + (mydsp_faustpower2_f(fTemp0) * fRec3[1])) + (mydsp_faustpower3_f(fTemp0) * fRec4[1]))))) / ((fSlow1 * (fTemp1 * mydsp_faustpower4_f(fTemp0))) + 1.0f)) - fRec4[1])) / fTemp2);
                T fTemp4 = ((fTemp0 * ((fRec4[1] + fTemp3) - fRec3[1])) / fTemp2);
                T fTemp5 = ((fTemp0 * ((fRec3[1] + fTemp4) - fRec2[1])) / fTemp2);
                T fTemp6 = ((fTemp0 * ((fRec2[1] + fTemp5) - fRec1[1])) / fTemp2);
                T fRec0 = (fRec1[1] + fTemp6);
                fRec1[0] = (fRec1[1] + (2.0f * fTemp6));
                fRec2[0] = (fRec2[1] + (2.0f * fTemp5));
                fRec3[0] = (fRec3[1] + (2.0f * fTemp4));
                fRec4[0] = (fRec4[1] + (2.0f * fTemp3));
                output0[i] = T(fRec0);
                fRec5[1] = fRec5[0];
                fRec1[1] = fRec1[0];
                fRec2[1] = fRec2[0];
                fRec3[1] = fRec3[0];
                fRec4[1] = fRec4[0];
            }
        }

        FAUSTPP_VIRTUAL void tick(T& inputs, T& outputs) {
            T fSlow0 = (0.00100000005f * T(fHslider0));
            T fSlow1 = (0.0411641225f * (T(fHslider1) + -0.707000017f));
            fRec5[0] = (fSlow0 + (0.999000013f * fRec5[1]));
            T fTemp0 = std::tan((fConst0 * fRec5[0]));
            T fTemp1 = (3.9000001f - (0.899999976f * std::pow((0.0f - (0.333333343f * (1.0f - (std::log10(fRec5[0]) + -0.30103001f)))), 0.200000003f)));
            T fTemp2 = (fTemp0 + 1.0f);
            T fTemp3 = ((fTemp0 * (((inputs - (fSlow1 * (fTemp1 * (((fRec1[1] + (fRec2[1] * fTemp0)) + (mydsp_faustpower2_f(fTemp0) * fRec3[1])) + (mydsp_faustpower3_f(fTemp0) * fRec4[1]))))) / ((fSlow1 * (fTemp1 * mydsp_faustpower4_f(fTemp0))) + 1.0f)) - fRec4[1])) / fTemp2);
            T fTemp4 = ((fTemp0 * ((fRec4[1] + fTemp3) - fRec3[1])) / fTemp2);
            T fTemp5 = ((fTemp0 * ((fRec3[1] + fTemp4) - fRec2[1])) / fTemp2);
            T fTemp6 = ((fTemp0 * ((fRec2[1] + fTemp5) - fRec1[1])) / fTemp2);
            T fRec0 = (fRec1[1] + fTemp6);
            fRec1[0] = (fRec1[1] + (2.0f * fTemp6));
            fRec2[0] = (fRec2[1] + (2.0f * fTemp5));
            fRec3[0] = (fRec3[1] + (2.0f * fTemp4));
            fRec4[0] = (fRec4[1] + (2.0f * fTemp3));
            outputs = T(fRec0);
            fRec5[1] = fRec5[0];
            fRec1[1] = fRec1[0];
            fRec2[1] = fRec2[0];
            fRec3[1] = fRec3[0];
            fRec4[1] = fRec4[0];
        }

    };
    
    template<typename T>
    TMoogLadder<T>::TMoogLadder()
    {
        fDsp = new mydsp<T>;        
        fDsp->instanceResetUserInterface();
    }

    template<typename T>
    TMoogLadder<T>::~TMoogLadder()
    {
        delete fDsp;
    }

    template<typename T>
    void TMoogLadder<T>::init(T sample_rate)
    {

        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        dsp.classInit(sample_rate);
        dsp.instanceConstants(sample_rate);
        clear();

    }

    template<typename T>
    void TMoogLadder<T>::clear() noexcept
    {

        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        dsp.instanceClear();

    }

    template<typename T>
    void TMoogLadder<T>::process(
        const T *in0,
        T *out0,
        unsigned count) noexcept
    {

        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        T *inputs[] = {
            const_cast<T *>(in0),
        };
        T *outputs[] = {
            out0,
        };
        dsp.compute(count, inputs, outputs);

    }

    template<typename T>
    const char *TMoogLadder<T>::parameter_label(unsigned index) noexcept
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

    template<typename T>
    const char *TMoogLadder<T>::parameter_short_label(unsigned index) noexcept
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

    template<typename T>
    const char *TMoogLadder<T>::parameter_symbol(unsigned index) noexcept
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

    template<typename T>
    const char *TMoogLadder<T>::parameter_unit(unsigned index) noexcept
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

    template<typename T>
    const ParameterRange<T> *TMoogLadder<T>::parameter_range(unsigned index) noexcept
    {
        switch (index) {
        
        case 0: {
            static const ParameterRange<T> range = { 20000, 20, 20000 };
            return &range;
        }
        
        case 1: {
            static const ParameterRange<T> range = { 1, 0.70719999, 25 };
            return &range;
        }
        
        default:
            return 0;
        }
    }

    template<typename T>
    bool TMoogLadder<T>::parameter_is_trigger(unsigned index) noexcept
    {
        switch (index) {
        
        default:
            return false;
        }
    }

    template<typename T>
    bool TMoogLadder<T>::parameter_is_boolean(unsigned index) noexcept
    {
        switch (index) {
        
        default:
            return false;
        }
    }

    template<typename T>
    bool TMoogLadder<T>::parameter_is_integer(unsigned index) noexcept
    {
        switch (index) {
        
        default:
            return false;
        }
    }

    template<typename T>
    bool TMoogLadder<T>::parameter_is_logarithmic(unsigned index) noexcept
    {
        switch (index) {
        
        case 0:
            return true;
        
        default:
            return false;
        }
    }

    template<typename T>
    T TMoogLadder<T>::get_parameter(unsigned index) const noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
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

    template<typename T>
    void TMoogLadder<T>::set_parameter(unsigned index, T value) noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
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


    template<typename T>
    T TMoogLadder<T>::get_cutoff() const noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        return dsp.fHslider0;
    }

    template<typename T>
    T TMoogLadder<T>::get_q() const noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        return dsp.fHslider1;
    }

    template<typename T>
    void TMoogLadder<T>::set_cutoff(T value) noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        dsp.fHslider0 = value;
    }

    template<typename T>
    void TMoogLadder<T>::set_q(T value) noexcept
    {
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        dsp.fHslider1 = value;
    }

    template<typename T>
    T TMoogLadder<T>::Tick(T input) { 
        T output=0.0; 
        mydsp<T> &dsp = static_cast<mydsp<T> &>(*fDsp);
        dsp.tick(input,output); 
        return output; 
    }
};





