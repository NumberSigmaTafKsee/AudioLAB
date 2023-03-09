
/////////////////////////////////////////////////////////////////////////////////
// AudioTK
/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <ATK/Core/TypedBaseFilter.h>

#include <ATK/Adaptive/BlockLMSFilter.h>
#include <ATK/Adaptive/LMSFilter.h>
#include <ATK/Adaptive/RLSFilter.h>

#include <ATK/Delay/FeedbackDelayNetworkFilter.h>
#include <ATK/Delay/FixedDelayLineFilter.h>
#include <ATK/Delay/HadamardMixture.h>
#include <ATK/Delay/HouseholderMixture.h>
#include <ATK/Delay/MultipleUniversalFixedDelayLineFilter.h>
#include <ATK/Delay/UniversalFixedDelayLineFilter.h>
#include <ATK/Delay/UniversalVariableDelayLineFilter.h>
#include <ATK/Delay/VariableDelayLineFilter.h>

#include <ATK/Distortion/DiodeClipperFilter.h>
#include <ATK/Distortion/HalfTanhShaperFilter.h>
#include <ATK/Distortion/SD1OverdriveFilter.h>
#include <ATK/Distortion/SimpleOverdriveFilter.h>
#include <ATK/Distortion/TS9OverdriveFilter.h>
#include <ATK/Distortion/TanhShaperFilter.h>

//#include <ATK/Distortion/WaveshaperFilter.h>

#include <ATK/Dynamic/AttackReleaseFilter.h>
#include <ATK/Dynamic/AttackReleaseHysteresisFilter.h>
#include <ATK/Dynamic/GainFilter.h>
#include <ATK/Dynamic/GainColoredCompressorFilter.h>
#include <ATK/Dynamic/GainColoredExpanderFilter.h>
#include <ATK/Dynamic/GainCompressorFilter.h>
#include <ATK/Dynamic/GainExpanderFilter.h>
#include <ATK/Dynamic/GainLimiterFilter.h>
#include <ATK/Dynamic/GainMaxColoredExpanderFilter.h>
#include <ATK/Dynamic/GainMaxCompressorFilter.h>
#include <ATK/Dynamic/GainMaxExpanderFilter.h>
#include <ATK/Dynamic/GainSwellFilter.h>
#include <ATK/Dynamic/PowerFilter.h>
#include <ATK/Dynamic/RelativePowerFilter.h>

#include <ATK/EQ/BesselFilter.h>
#include <ATK/EQ/ButterworthFilter.h>
#include <ATK/EQ/ChamberlinFilter.h>
#include <ATK/EQ/Chebyshev1Filter.h>
#include <ATK/EQ/Chebyshev2Filter.h>
#include <ATK/EQ/CustomFIRFilter.h>
#include <ATK/EQ/CustomIIRFilter.h>
#include <ATK/EQ/FIRFilter.h>
#include <ATK/EQ/IIRFilter.h>
#include <ATK/EQ/LinkwitzRileyFilter.h>
#include <ATK/EQ/PedalToneStackFilter.h>
#include <ATK/EQ/RIAAFilter.h>
#include <ATK/EQ/RemezBasedFilter.h>
#include <ATK/EQ/RobertBristowJohnsonFilter.h>
#include <ATK/EQ/FourthOrderFilter.h>
#include <ATK/EQ/SecondOrderFilter.h>
#include <ATK/EQ/SecondOrderSVFFilter.h>
#include <ATK/EQ/SimpleIIRFilter.h>
#include <ATK/EQ/TimeVaryingIIRFilter.h>
#include <ATK/EQ/TimeVaryingSecondOrderFilter.h>
#include <ATK/EQ/TimeVaryingSecondOrderSVFFilter.h>
#include <ATK/EQ/ToneStackFilter.h>

#include <ATK/Preamplifier/TriodeFilter.h>
#include <ATK/Preamplifier/Triode2Filter.h>
#include <ATK/Preamplifier/FollowerTransistorClassAFilter.h>
#include <ATK/Preamplifier/TransistorClassAFilter.h>
#include <ATK/Preamplifier/Triode2Filter.h>
#include <ATK/Preamplifier/TriodeFilter.h>
#include <ATK/Preamplifier/KorenTriodeFunction.h>
#include <ATK/Preamplifier/EnhancedKorenTriodeFunction.h>
#include <ATK/Preamplifier/LeachTriodeFunction.h>
#include <ATK/Preamplifier/MunroPiazzaTriodeFunction.h>
#include <ATK/Preamplifier/DempwolfTriodeFunction.h>

#include <ATK/Special/ConvolutionFilter.h>

#include <ATK/Reverberation/AllPassReverbFilter.h>
#include <ATK/Reverberation/LowPassReverbFilter.h>

#include <ATK/Tools/ApplyGainFilter.h>
#include <ATK/Tools/BufferFilter.h>
#include <ATK/Tools/CachedCosinusGeneratorFilter.h>
#include <ATK/Tools/CachedSinusGeneratorFilter.h>
#include <ATK/Tools/DecimationFilter.h>
#include <ATK/Tools/DerivativeFilter.h>
#include <ATK/Tools/DryWetFilter.h>
#include <ATK/Tools/MaxFilter.h>
#include <ATK/Tools/MSFilter.h>
#include <ATK/Tools/MuteSoloBufferFilter.h>
#include <ATK/Tools/OffsetVolumeFilter.h>
#include <ATK/Tools/OneMinusFilter.h>
#include <ATK/Tools/OversamplingFilter.h>
#include <ATK/Tools/PanFilter.h>
#include <ATK/Tools/SinusGeneratorFilter.h>
#include <ATK/Tools/SumFilter.h>
#include <ATK/Tools/TanFilter.h>
#include <ATK/Tools/VolumeFilter.h>
#include <ATK/Tools/WhiteNoiseGeneratorFilter.h>

namespace AudioTK
{

    template <typename DataType>
    class ATKMonoInputFilter : public ATK::TypedBaseFilter<DataType> {
    public:
        using ATK::TypedBaseFilter<DataType>::outputs;

        explicit ATKMonoInputFilter() : ATK::TypedBaseFilter<DataType>(0, 1){}

        void set_input(DataType * inputs, int size) {
            mInputs = inputs;
            mSize = size;
        }

    protected:
        DataType * mInputs = nullptr;
        int mSize = 0;
        
        virtual void process_impl(int64_t size) const {
            for (int64_t i = 0; i < size; ++i) {
                outputs[0][i] = mInputs[i];
            }            
        }
    };

    template <typename DataType>
    class ATKInputFilter : public ATK::TypedBaseFilter<DataType> {
    public:
        using ATK::TypedBaseFilter<DataType>::outputs;

        explicit ATKInputFilter(int channels) : ATK::TypedBaseFilter<DataType>(0, channels), mChannels(channels) {}

        void set_inputs(DataType ** inputs, int size) {
            mInputs = inputs;
            mSize = size;
        }

    protected:
        DataType** mInputs = nullptr;
        int mSize = 0;
        int mChannels = 1;

        virtual void process_impl(int64_t size) const {
            for (int c = 0; c < mChannels; ++c) {
                for (int64_t i = 0; i < size; ++i) {
                    outputs[c][i] = mInputs[c][i];
                }
            }
        }
    };

    template <typename DataType>
    class ATKMonoOutputFilter : public ATK::TypedBaseFilter<DataType> {
    public:
        using ATK::TypedBaseFilter<DataType>::converted_inputs;

        explicit ATKMonoOutputFilter() : ATK::TypedBaseFilter<DataType>(1,0) {}

        void set_output(DataType* outputs, int size) {
            mOutputs = outputs;
            mSize = size;
        }

    protected:
        DataType * mOutputs = nullptr;
        int mSize = 0;
        

        virtual void process_impl(int64_t size) const {
            for (int64_t i = 0; i < size; ++i) {
                mOutputs[i] = converted_inputs[0][i];
            }            
        }
    };

    template <typename DataType>
    class ATKOutputFilter : public ATK::TypedBaseFilter<DataType> {
    public:
        using ATK::TypedBaseFilter<DataType>::converted_inputs;

        explicit ATKOutputFilter(int channels) : ATK::TypedBaseFilter<DataType>(channels, 0), mChannels(channels) {}

        void set_outputs(DataType** outputs, int size) {
            mOutputs = outputs;
            mSize = size;
        }

    protected:
        DataType** mOutputs = nullptr;
        int mSize = 0;
        int mChannels = 1;

        virtual void process_impl(int64_t size) const {
            for (int c = 0; c < mChannels; ++c) {
                for (int64_t i = 0; i < size; ++i) {
                    mOutputs[c][i] = converted_inputs[c][i];
                }
            }
        }
    };

    
    using InType = ATKInputFilter<DspFloatType>;
    using OutType = ATKOutputFilter<DspFloatType>;
    using MonoInType = ATKMonoInputFilter<DspFloatType>;
    using MonoOutType = ATKMonoOutputFilter<DspFloatType>;
    using BaseFilterType = ATK::TypedBaseFilter<DspFloatType>;
    using SumType = ATK::SumFilter<DspFloatType>;
    using VolumeType = ATK::VolumeFilter<DspFloatType>;
    using OversamplingType_2 = ATK::OversamplingFilter<DspFloatType, ATK::Oversampling6points5order_2<DspFloatType>>;
    using OversamplingType_4 = ATK::OversamplingFilter<DspFloatType, ATK::Oversampling6points5order_4<DspFloatType>>;
    using OversamplingType_8 = ATK::OversamplingFilter<DspFloatType, ATK::Oversampling6points5order_8<DspFloatType>>;
    using DownsamplingType = ATK::DecimationFilter<DspFloatType>;
    using LowPassType = ATK::IIRFilter<ATK::ButterworthLowPassCoefficients<DspFloatType>>;
    using ToneStackType = ATK::IIRFilter<ATK::ToneStackCoefficients<DspFloatType>>;
    //using CabType = ATK::ConvolutionFilter<DspFloatType>;
    
    struct ATKMonoMultiFilterProcessor : public MonoFXProcessor
    {
    private:
        std::unique_ptr<MonoInType> mInputs;
        std::unique_ptr<MonoOutType> mOutputs;
        std::vector<ATK::BaseFilter*> filters;
        std::unique_ptr<OversamplingType_2> mOversample;
        std::unique_ptr<DownsamplingType> mDecimator;
        DspFloatType sampleRate;
    public:

        ATKMonoMultiFilterProcessor(DspFloatType sampleRate=44100.0) : MonoFXProcessor()
        {
            this->sampleRate = sampleRate;
            mInputs.reset(new MonoInType());
            mInputs->set_output_sampling_rate(sampleRate);

            mOutputs.reset(new MonoOutType());
            mOutputs->set_input_sampling_rate(sampleRate);        

            mOversample.reset(new OversamplingType_2(1));
            mOversample->set_input_sampling_rate(sampleRate);
            mOversample->set_output_sampling_rate(2*sampleRate);

            mDecimator.reset(new DownsamplingType(1));
            mDecimator->set_input_sampling_rate(2*sampleRate);
            mDecimator->set_output_sampling_rate(sampleRate);
            
            filters.push_back(mOversample.get());
        }
        ~ATKMonoMultiFilterProcessor() {

        }

        void addFilter(ATK::BaseFilter * p) {
            
        }
        void connect() 
        {
            
            filters.push_back(mDecimator.get());
            

            auto p = filters.begin();
            auto last = p;
            (*p)->set_input_port(0, mInputs.get(), 0);            
            p++;

            while(p != filters.end())
            {
                (*p)->set_input_port(0,(*last),0);                
                last = p++;
            }
            mOutputs->set_input_port(0,(*last),0);            
        }
        void ProcessBlock(size_t nFrames, DspFloatType * inputs, DspFloatType * outputs)
        {        
            mInputs->set_input(inputs, nFrames);
            mOutputs->set_output(outputs, nFrames);
            mOutputs->process(nFrames);
        }
        void ProcessInplace(size_t nFrames, DspFloatType * buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };
    
    struct ATKMultiFilterProcessor : public StereoFXProcessor
    {
    private:
        std::unique_ptr<InType> mInputs;
        std::unique_ptr<OutType> mOutputs;
        std::vector<ATK::BaseFilter*> filters;
        std::unique_ptr<OversamplingType_2> mOversample;
        std::unique_ptr<DownsamplingType> mDecimator;
        DspFloatType sampleRate;
    
    public:

        ATKMultiFilterProcessor(DspFloatType sampleRate=44100.0) : StereoFXProcessor()
        {
            this->sampleRate = sampleRate;
            mInputs.reset(new InType(2));
            mInputs->set_output_sampling_rate(sampleRate);

            mOutputs.reset(new OutType(2));
            mOutputs->set_input_sampling_rate(sampleRate);        

            mOversample.reset(new OversamplingType_2(2));
            mOversample->set_input_sampling_rate(sampleRate);
            mOversample->set_output_sampling_rate(2*sampleRate);

            mDecimator.reset(new DownsamplingType(2));
            mDecimator->set_input_sampling_rate(2*sampleRate);
            mDecimator->set_output_sampling_rate(sampleRate);

            
            filters.push_back(mOversample.get());
        }
        ~ATKMultiFilterProcessor() {

        }

        void addFilter(ATK::BaseFilter * p) {
            
        }
        void connect() 
        {
            
            filters.push_back(mDecimator.get());
            

            auto p = filters.begin();
            auto last = p;
            (*p)->set_input_port(0, mInputs.get(), 0);
            (*p)->set_input_port(1, mInputs.get(), 1);
            p++;

            while(p != filters.end())
            {
                (*p)->set_input_port(0,(*last),0);
                (*p)->set_input_port(1,(*last),1);
                last = p++;
            }
            mOutputs->set_input_port(0,(*last),0);
            mOutputs->set_input_port(1,(*last),1);
        }
        void ProcessBlock(size_t nFrames, DspFloatType ** inputs, DspFloatType ** outputs)
        {        
            mInputs->set_inputs(inputs, nFrames);
            mOutputs->set_outputs(outputs, nFrames);
            mOutputs->process(nFrames);
        }
        void ProcessInplace(size_t nFrames, DspFloatType ** buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };

    
    // process a filter
    struct ATKMonoFilterProcessor : public MonoFXProcessor
    {
    private:
        std::unique_ptr<MonoInType> mInputs;
        std::unique_ptr<MonoOutType> mOutputs;
        ATK::BaseFilter * filter;
        std::unique_ptr<OversamplingType_2> mOversample;
        std::unique_ptr<DownsamplingType> mDecimator;
        DspFloatType sampleRate;
    public:
        ATKMonoFilterProcessor(ATK::BaseFilter* f,DspFloatType sampleRate=44100) : MonoFXProcessor()
        {
            this->sampleRate = sampleRate;
            mInputs.reset(new MonoInType());
            //mInputs->set_input_sampling_rate(sampleRate);
            mInputs->set_output_sampling_rate(sampleRate);

            mOutputs.reset(new MonoOutType());
            mOutputs->set_input_sampling_rate(sampleRate);        
            //mOutputs->set_output_sampling_rate(sampleRate);        
            
            mOversample.reset(new OversamplingType_2(1));
            mOversample->set_input_sampling_rate(sampleRate);
            mOversample->set_output_sampling_rate(2*sampleRate);

            mDecimator.reset(new DownsamplingType);
            mDecimator->set_input_sampling_rate(2*sampleRate);
            mDecimator->set_output_sampling_rate(sampleRate);
            
            filter = f;
            filter->set_input_sampling_rate(sampleRate*2);        
            //filter->set_output_sampling_rate(sampleRate);        
            connect();
        }
        ~ATKMonoFilterProcessor() {

        }

        void setFilter(ATK::BaseFilter * p)
        {
            filter = p;
        }

        void connect() 
        {               
            mOversample->set_input_port(0, mInputs.get(), 0);                        
            filter->set_input_port(0, mOversample.get(), 0);                    
            mDecimator->set_input_port(0, filter, 0);                                      
            mOutputs->set_input_port(0,mDecimator.get(),0);            
        }
        void ProcessBlock(size_t nFrames, DspFloatType * inputs, DspFloatType * outputs)
        {                    
            mInputs->set_input(inputs, nFrames);
            mOutputs->set_output(outputs, nFrames);
            mOutputs->process(nFrames);
        }
        void ProcessInplace(size_t nFrames, DspFloatType * buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };

    struct ATKMonoFilter : public MonoFXProcessor
    {
        ATKMonoFilterProcessor * p;
        ATKMonoFilter() : MonoFXProcessor()
        {
            p = NULL;
        }
        ~ATKMonoFilter() {
            if(p) delete p;
        }
        void setFilter(ATK::BaseFilter * f) {
            if(p) delete p;
            p = new ATKMonoFilterProcessor(f);
        }
        void ProcessBlock(size_t nFrames, DspFloatType * inputs, DspFloatType * outputs)
        {        
            if(p != NULL)
                p->ProcessBlock(nFrames,inputs,outputs);
        }
        void ProcessInplace(size_t nFrames, DspFloatType * buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };   


    // process a filter
    struct ATKFilterProcessor : public StereoFXProcessor
    {
    private:
        std::unique_ptr<InType> mInputs;
        std::unique_ptr<OutType> mOutputs;
        ATK::BaseFilter * filter;
        std::unique_ptr<OversamplingType_2> mOversample;
        std::unique_ptr<DownsamplingType> mDecimator;
        DspFloatType sampleRate;
    public:
        ATKFilterProcessor(ATK::BaseFilter* f,DspFloatType sampleRate=44100) : StereoFXProcessor()
        {
            this->sampleRate = sampleRate;
            mInputs.reset(new InType(2));
            //mInputs->set_input_sampling_rate(sampleRate);
            mInputs->set_output_sampling_rate(sampleRate);

            mOutputs.reset(new OutType(2));
            mOutputs->set_input_sampling_rate(sampleRate);        
            //mOutputs->set_output_sampling_rate(sampleRate);        
            
            mOversample.reset(new OversamplingType_2(2));
            mOversample->set_input_sampling_rate(sampleRate);
            mOversample->set_output_sampling_rate(2*sampleRate);

            mDecimator.reset(new DownsamplingType);
            mDecimator->set_input_sampling_rate(2*sampleRate);
            mDecimator->set_output_sampling_rate(sampleRate);
            
            filter = f;
            filter->set_input_sampling_rate(sampleRate*2);        
            //filter->set_output_sampling_rate(sampleRate);        
            connect();
        }
        ~ATKFilterProcessor() {

        }

        void setFilter(ATK::BaseFilter * p)
        {
            filter = p;
        }

        void connect() 
        {   
            
            mOversample->set_input_port(0, mInputs.get(), 0);
            mOversample->set_input_port(1, mInputs.get(), 1);
            
            filter->set_input_port(0, mOversample.get(), 0);
            filter->set_input_port(1, mOversample.get(), 1);
            
            mDecimator->set_input_port(0, filter, 0);
            mDecimator->set_input_port(1, filter, 1);
                                      
            mOutputs->set_input_port(0,mDecimator.get(),0);
            mOutputs->set_input_port(1,mDecimator.get(),1);
        }
        void ProcessBlock(size_t nFrames, DspFloatType ** inputs, DspFloatType ** outputs)
        {                    
            mInputs->set_inputs(inputs, nFrames);
            mOutputs->set_outputs(outputs, nFrames);
            mOutputs->process(nFrames);
        }
        void ProcessInplace(size_t nFrames, DspFloatType ** buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };
    
    struct ATKFilter : public StereoFXProcessor
    {
        ATKFilterProcessor * p;
        ATKFilter() : StereoFXProcessor()
        {
            p = NULL;
        }
        ~ATKFilter() {
            if(p) delete p;
        }
        void setFilter(ATK::BaseFilter * f) {
            if(p) delete p;
            p = new ATKFilterProcessor(f);
        }
        void ProcessBlock(size_t nFrames, DspFloatType ** inputs, DspFloatType ** outputs)
        {        
            if(p != NULL)
                p->ProcessBlock(nFrames,inputs,outputs);
        }
        void ProcessInplace(size_t nFrames, DspFloatType ** buffer) {
			ProcessBlock(nFrames,buffer,buffer);
		}
    };   
}
