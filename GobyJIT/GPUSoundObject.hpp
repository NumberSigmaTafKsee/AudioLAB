// FX Object
#pragma once


#include <memory>
#include <list>
#include <vector>
#include <map>
#include <random>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cassert>

//#define DSPFLOATDOUBLE 1
typedef float DspFloatType;

#include "Undenormal.hpp"
#include "StdNoise.hpp"
#include "MusicFunctions.hpp"
#include "FX/ClipFunctions.hpp"



extern Default noise;
extern DSPType sampleRate;
extern DSPType invSampleRate;

namespace AudioDSP
{
    struct Random
    {
        Random() { srand(time(NULL)); }

        DSPType      frand() { return ((DSPType)rand()/(DSPType)RAND_MAX); }
        uint64_t    randint(int min, int max) { return round((max-min)*frand())+min; }
        bool        flip(DSPType prob) { return frand() < prob; }
        uint64_t    random(int mod) { return rand() % mod; }    
    };

    enum ObjectType
    {    
        
        GPU_PARAMETER2_PROCESSOR,     // F(X,Y)
        GPU_STEREOSPLITTER_PROCESSOR,
        
        GPU_MONO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
        GPU_MONO_FUNCTION_PROCESSOR,       // Input    
        GPU_MONO_OSCILLATOR_PROCESSOR,     // Generator no input
        GPU_MONO_FILTER_PROCESSOR,         // DSPType Input
        GPU_MONO_AMPLIFIER_PROCESSOR,    
        GPU_MONO_FX_PROCESSOR,    
        GPU_MONO_CASCADE_PROCESSOR,
        GPU_MONO_MIXER_PROCESSOR,
        GPU_MONO_MORPHER_PROCESSOR,
        GPU_MONO_OPERATOR_PROCESSOR,
        GPU_MONO_OVERSAMPLE_PROCESSOR,    // Resampler up/fx/down
        GPU_MONO_UPSAMPLE_PROCESSOR,      // resample up
        GPU_MONO_DOWNSAMPLE_PROCESSOR,    // resample down

        GPU_STEREO_FX_PROCESSOR,
        GPU_STEREO_OVERSAMPLER_PROCESSOR,    // Resampler up/fx/down
        GPU_STEREO_UPSAMPLER_PROCESSOR,      // resample up
        GPU_STEREO_DOWNSAMPLER_PROCESSOR,    // resample down
        GPU_STEREO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
        GPU_STEREO_FUNCTION_PROCESSOR,       // Input
        GPU_STEREO_PARAMETER2_PROCESSOR,     // F(X,Y)
        GPU_STEREO_OSCILLATOR_PROCESSOR,     // Generator no input
        GPU_STEREO_FILTER_PROCESSOR,         // DSPType Input
        GPU_STEREO_AMPLIFIER_PROCESSOR,
        GPU_STEREO_CASCADE_PROCESSOR,
        GPU_STEREO_OPERATOR_PROCESSOR,

        GPU_MONO_SIGNAL_SOURCE_PROCESSOR,
        GPU_STEREO_SIGNAL_SOURCE_PROCESSOR,
        GPU_MONO_SIGNAL_SINK_PROCESSOR,
        GPU_STEREO_SIGNAL_SINK_PROCESSOR,
        GPU_INTERLEAVE_PROCESSOR,
        GPU_DEINERLEAVE_PROCESSOR,
        GPU_FILTER_BANK_PROCESSOR,
        GPU_SPECTRUM_PROCESSOR,    
    };


    template<typename DSPType>
    struct GPUSoundProcessor
    {    
        
        DSPType preGain = 1;
        DSPType postGain = 1;

        virtual ObjectType getType() const = 0;

        // i do not want any kind of complicated data structure
        // just a simple function to set the port value    
        virtual void setPort(int port, DSPType value) {
            printf("No port %d\n",port);
        }
        virtual void setPort2(int port, DSPType a, DSPType b) {
            printf("No port %d\n",port);
        }
        virtual void setPortV(int port, const std::vector<DSPType> & v) {
            printf("No port %d\n",port);
        }
        virtual DSPType getPort(int port) {
            printf("No port %d\n",port);
            return 0;
        }
        virtual DSPType getPort2(int port, DSPType v) {
            printf("No port %d\n",port);
            return 0;
        }
        virtual void getPortV(int port, std::vector<DSPType> & v) {
            printf("No port %d\n",port);
        }
        virtual void printPortMap() {
            printf("No ports\n");
        }
        virtual void randomize() {

        }
        bool loadPreset(const char * filename) {
            return false;
        }
        bool savePreset(const char * filename) {
            return false;
        }    
    };

    template<typename DSPType>
    struct GPUPort
    {
        int port;
        DSPType value;
        GPUSoundProcessor<DSPType> * p;
    };

    template<typename DSPType>
    struct GPUPorts
    {
        std::list<std::shared_ptr<GPUPort<DSPType>> ports;

        using PortMap = std::map<std::string,GPUPort<DSPType>*>;
        PortMap portmap;

        GPUPorts() {

        }
        void addPort(const std::string & name, GPUPort<DSPType> * p) {
            ports.push_back( std::shared_ptr<Port>(p, [](GPUPort<DSPType> * p){ delete p; }));
            portmap[name] = p;
        }
        void Run() {
            auto i = ports.begin();
            while(i != ports.end())
            {
                auto port = *i;
                port->p->setPort(port->port,port->value);
                i++;
            }
        }
    };

    template<typename DSPType>
    struct GPUMonoProcessor : public GPUSoundProcessor<DSPType>
    {
        
        GPUMonoProcessor() : GPUSoundProcessortemplate<typename DSPType>()
        {

        }

            
        void InplaceProcess(size_t n, DSPType * buffer) {
            ProcessBlock(n,buffer,buffer);
        }

        virtual DSPType Tick(DSPType I=1, DSPType A=1, DSPType X=0, DSPType Y=0)
        {
            assert(1==0);
        }
        virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) 
        {
            assert(1==0);
        }

    };


    template<typename DSPType>
    struct GPUMonoCascadeProcessor : public GPUMonoProcessor<DSPType>
    {
        std::list<MonoProcessor*> procs;

        GPUMonoCascadeProcessor() : GPUMonoProcessor<DSPType>()
        {

        }

        ObjectType getType() const {
            return GPU_MONO_CASCADE_PROCESSOR;
        }

        void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) {
            auto i = procs.begin();
            memcpy(outputs,inputs,n*sizeof(DSPType));
            for(; i != procs.end(); i++)    
            {
                switch((*i)->getType())
                {
                case MONO_GENERATOR_PROCESSOR:
                case MONO_FUNCTION_PROCESSOR:
                case MONO_FILTER_PROCESSOR:
                case MONO_AMPLIFIER_PROCESSOR:
                    for(size_t j = 0; j < n; j++)
                        outputs[j] = (*i)->Tick(outputs[j]);
                    break;
                case MONO_OSCILLATOR_PROCESSOR:
                    for(size_t j = 0; j < n; j++)
                        outputs[j] *= (*i)->Tick(outputs[j]);
                    break;
                case MONO_FX_PROCESSOR:
                    (*i)->ProcessBlock(n,outputs,outputs);
                    break;
                }
            }
        }
    };

    template<typename DSPType>
    struct StereoProcessor : public SoundProcessor<DSPType>
    {
        
        DSPType pan;

        GPUStereoProcessor() : GPUSoundProcessor<DSPType>() {
            pan = 0.5;
        }

        
        virtual void ProcessBlock(size_t n, DSPType ** inputs, DSPType ** outputs) 
        {
            assert(1==0);
        }

        virtual DSPType Tick(DSPType IL, DSPType IR, DSPType &L, DSPType &R, DSPType A=1, DSPType X=0, DSPType Y=0)
        {
            assert(1==0);
        }   
        void InplaceProcess(size_t n, DSPType ** buffer) {
            ProcessBlock(n,buffer,buffer);
        }

        void Run(size_t n, DSPType ** inputs, DSPType ** outputs) {
            ProcessBlock(n,inputs,outputs);    
        }
        
    };

    template<typename DSPType>
    struct GPUStereoCascadeProcessor : public GPUStereoProcessor<DSPType>
    {
        std::list<GPUStereoProcessor<DSPType>*> procs;

        GPUStereoCascadeProcessor() : GPUStereoProcessor<DSPType>()
        {

        }
        virtual ObjectType getType() const { return GPU_STEREO_CASCADE_PROCESSOR; }

        void ProcessBlock(size_t n, DSPType ** inputs, DSPType ** outputs) {
            auto i = procs.begin();
            memcpy(outputs[0],inputs[0],n*sizeof(DSPType));
            memcpy(outputs[1],inputs[1],n*sizeof(DSPType));
            for(; i != procs.end(); i++)    
            {
                switch((*i)->getType())
                {            
                case STEREO_GENERATOR_PROCESSOR:
                case STEREO_FUNCTION_PROCESSOR:
                case STEREO_FILTER_PROCESSOR:
                case STEREO_AMPLIFIER_PROCESSOR:
                    for(size_t j = 0; j < n; j++)
                    {
                        DSPType L,R;
                        (*i)->Tick(outputs[0][j],outputs[1][j],L,R);
                        outputs[0][j] = L;
                        outputs[1][j] = R;
                    }
                    break;
                case STEREO_OSCILLATOR_PROCESSOR:
                    for(size_t j = 0; j < n; j++)
                    {
                        DSPType L,R;
                        (*i)->Tick(outputs[0][j],outputs[1][j],L,R);
                        outputs[0][j] *= L;
                        outputs[1][j] *= R;
                    }
                    break;
                case STEREO_FX_PROCESSOR:
                    (*i)->ProcessBlock(n,outputs,outputs);
                    break;
                }
            }
        }
    };

    // FX Processors should always be a stream and should not use Tick
    // It is a defect that they have Tick and will be eventually removed

    template<typename DSPType>
    struct GPUMonoFXProcessor : public GPUMonoProcessor<DSPType>
    { 
        GPUMonoFXProcessor() : GPUMonoProcessor<DSPType>() {
    
        }

        ObjectType getType() const {
            return GPU_MONO_FX_PROCESSOR;
        }
        
        virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) = 0;
    };

    template<typename DSPType>
    struct GPUStereoFXProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUStereoFXProcessor() : GPUStereoProcessor<DSPType>() {
    
        }
        
        ObjectType getType() const {
            return GPU_STEREO_FX_PROCESSOR;
        }        

        virtual void ProcessBlock(size_t n, DSPType ** inputs, DSPType ** outputs) = 0;
    };

    // Envelopes, LFOs, Noise
    // Difference between generator and function is small
    // A function takes input (I=input)
    // generator does not take input (I=index or unused)
    // F(I,A,X,Y)
    // G(Index,A,X,Y)
    // I must always be used by a function and it must be a function of I
    // a Parameter is a function of multiple variables without modulation Parametric(X,Y...)
    // I may or may not be used by a generator
    // A is almost always amplitude or gain and may not be used
    // X and Y are variable parameter depend on the unit and may not be used at all

    template<typename DSPType>
    struct GPUGeneratorProcessor : public GPUMonoProcessor<DSPType>
    {
        GPUGeneratorProcessor() : GPUMonoProcessor<DSPType>() 
        {
            
        }
        ObjectType getType() const {
            return GPU_MONO_GENERATOR_PROCESSOR;
        }
        virtual DSPType Tick(DSPType I=0, DSPType A=0, DSPType X=0, DSPType Y=0) = 0;

        void Generate(size_t n, DSPType * output) {
            for(size_t i = 0; i < n; i++)
                output[i] = Tick();
        }
        virtual void ProcessBlock(size_t n, DSPType * in, DSPType * out)
        {
            for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
        }
    };

    template<typename DSPType>
    struct GPUMixerProcessor : public GPUSoundProcessor<DSPType>
    {
        GPUMixerProcessor() : GPUSoundProcessor<DSPType>()
        {

        }
        ObjectType getType() const {
            return GPU_MONO_MIXER_PROCESSOR;
        }

        virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType *) {

        }
        virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType **) {
            
        }

    };

    template<typename DSPType>
    struct GPUFunctionProcessor : public GPUMonoProcessor<DSPType>
    {
        GPUFunctionProcessor() : GPUMonoProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return GPU_MONO_FUNCTION_PROCESSOR;
        }

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=0, DSPType Y=0) = 0;

        void ProcessBlock(size_t n, DSPType * input, DSPType * output) {
            for(size_t i = 0; i < n; i++)
                output[i] = Tick(input[i]);
        }
    };

    template<typename DSPType>
    struct GPUParameter2Processor : public GPUSoundProcessor<DSPType>
    {
        GPUParameter2Processor() : GPUSoundProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return GPU_PARAMETER2_PROCESSOR;
        }

        virtual DSPType Tick(DSPType a, DSPType b) = 0;

        void ProcessBlock(size_t n, DSPType * x, DSPType * y, DSPType * output) {
            for(size_t i = 0; i < n; i++)
                output[i] = Tick(x[i],y[i]);
        }
    };

    template<typename DSPType>
    struct GPUStereoSplitterProcessor : public GPUSoundProcessor<DSPType>
    {
        GPUStereoSplitterProcessor() : GPUSoundProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return GPU_STEREOSPLITTER_PROCESSOR;
        }

        virtual DSPType Tick(DSPType in, DSPType &a, DSPType &b) = 0;

        void ProcessBlock(size_t n, DSPType * in, DSPType * a, DSPType * b) {
            for(size_t i = 0; i < n; i++)
            {   
                DSPType x = a[i];
                DSPType y = b[i];
                Tick(in[i],x,y);
            }
        }
    };

    template<typename DSPType>
    struct GPUOscillatorProcessor : public GPUMonoProcessor<DSPType>
    {
        /*
        std::vector<OscillatorProcessor*> slaves;    
        int    m_waveform = 0;
        DSPType m_morph = 0;
        DSPType m_freq  = 440.0f;
        DSPType m_phase = 0;
        DSPType m_index = 1;
        DSPType m_gain = 1;
        DSPType m_fm = 0;
        DSPType m_pm = 0;
        DSPType m_fenv = 1;
        DSPType m_penv = 1;  
        DSPType m_drift = 0;  
        DSPType m_mod = 1;
        DSPType m_cmod = 1;
        DSPType m_env = 1;
        DSPType m_lfo = 1;
        DSPType m_pwm = 0.5;
        */

        GPUOscillatorProcessor() : GPUMonoProcessor<DSPType>() 
        {
            
        }
        ObjectType getType() const {
            return GPU_MONO_OSCILLATOR_PROCESSOR;
        }

        /*
        virtual void  setWaveform(int waveform) { m_waveform = waveform; }
        virtual void  setFrequency(DSPType f) { m_freq = f; }
        virtual void  setFrequencyCV(DSPType f) { m_freq = cv2freq(f); }
        virtual void  setPWM(DSPType p) { m_pwm = p; }
        virtual void  setPhase(DSPType p) { m_phase =p; }
        virtual void  setGain(DSPType a) { m_gain = a; }
        virtual void  setIndex(DSPType i) { m_index = i; }
        virtual void  setFM(DSPType f) { m_fm = f; }
        virtual void  setPM(DSPType p) { m_pm = p; }
        virtual void  setFreqEnv(DSPType e) { m_fenv = e; }
        virtual void  setPhaseEnv(DSPType e) { m_penv = e; }
        virtual void  setLFO(DSPType r) { m_lfo = r; }    
        virtual void  setModulator(DSPType r) { m_mod = r; }
        virtual void  setCModulator(DSPType r) { m_cmod = r; }
        virtual void  setEnvelope(DSPType e) { m_env = e; }
        virtual void  setDrift(DSPType d) { m_drift = d; }
        virtual void  sync() { m_phase = 0; }

        virtual void  randomize() {

        }
        */
        virtual DSPType Tick(DSPType I=0, DSPType A=1, DSPType X=0, DSPType Y=0) = 0;    

        virtual void ProcessBlock(size_t n, DSPType * in, DSPType * out)
        {
            for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
        }
    };

    template<typename DSPType>
    struct GPUFilterProcessor : public MonoProcessor<DSPType>
    {
        GPUFilterProcessor() : GPUMonoProcessor<DSPType>() {
    
        }
        ObjectType getType() const {
            return GPU_MONO_FILTER_PROCESSOR;
        }
        /*
        DSPType m_cutoff = 1.0f;
        DSPType m_resonance = 0.5f;
        DSPType m_q = 0.5;
        DSPType m_bw = 1;
        DSPType m_gain = 1;
        DSPType m_slope = 1;
        DSPType m_A = 1;
        DSPType m_X = 0;
        DSPType m_Y = 0;
        DSPType m_cMin = -1.0f;
        DSPType m_cMax = 1.0f;
        DSPType m_dc   = 0.0f;
        DSPType m_distDB = 0.0f;
        DSPType m_pre = 1.0f;
        DSPType m_post = 1.0f;

        virtual void updateTick(DSPType A = 1, DSPType X = 0, DSPType Y = 0) {
            m_A = A;
            m_X = X;
            m_Y = Y;
        }
        virtual void setA(DSPType A) {
            m_A = A;
        }
        virtual void setX(DSPType X) {
            m_X = X;        
        }
        virtual void setY(DSPType Y) {
            m_Y = Y;
        }
        virtual void setDC(DSPType dc) {
            m_dc = dc;
        }
        virtual void setMin(DSPType min) {
            m_cMin = min;
        }
        virtual void setMax(DSPType max) {
            m_cMax = max;
        }
        virtual void setDistDB(DSPType db) {
            m_distDB = db;
        }
        virtual void setPre(DSPType p) {
            m_pre = p;
        }
        virtual void setPost(DSPType p) {
            m_post = p;
        }
        virtual void setCutoff(DSPType c) {
            m_cutoff = c;
        }
        virtual void setResonance(DSPType r) {
            m_resonance = r;
        }
        virtual void setCutoffCV(DSPType f) {
            m_cutoff = cv2freq(f);
        }
        virtual void setResonanceCV(DSPType r) {
            m_resonance  = r/10.0;
        }
        virtual void setQ(DSPType q) {
            m_q = q;
        }
        virtual void setBW(DSPType bw) {
            m_bw = bw;
        }
        virtual void setSlope(DSPType s) {
            m_slope = s;
        }
        virtual void setGain(DSPType g) {
            m_gain = g;
        }
        */
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X =0, DSPType Y=0)=0;

        // cascade
        DSPType Run(DSPType I, DSPType A=1, DSPType X=0, DSPType Y=0) {
            DSPType x = Tick(I,A,X,Y);
            return x;
        }
    
        // cascade
        void ProcessBlock(size_t numSamples, DSPType * inputs, DSPType * outputs)
        {
            for(size_t i = 0; i < numSamples; i++)
            {
                outputs[i] = Run(inputs[i]);
            }        
        }
    };


    template<typename DSPType>
    struct GPUAmplifierProcessor : public GPUMonoProcessor<DSPType>
    {
        
        GPUAmplifierProcessor() : GPUMonoProcessor<DSPType>()
        {
        
        }
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1)=0;

        ObjectType getType() const {
            return GPU_MONO_AMPLIFIER_PROCESSOR;
        }
        
        // cascade
        void ProcessBlock(size_t numSamples, DSPType * inputs, DSPType * outputs)
        {
            for(size_t i = 0; i < numSamples; i++)
            {
                outputs[i] = Tick(inputs[i]);
            }        
        }
    };

    // StereoOscillator
    // StereoGenerator
    // StereoFunction

    template<typename DSPType>
    struct GPUStereoOscillatorProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUOscillatorProcessor<DSPType> * osc;
        
        GPUStereoOscillatorProcessor(OscillatorProcessor * o) : GPUStereoProcessor<DSPType>() {
            osc = o;
            pan = 0.5;
        }
        ObjectType getType() const {
            return GPU_STEREO_OSCILLATOR_PROCESSOR;
        }
        void setPan(DSPType p ) {
            pan = p;
        }
        void ProcessBlock(size_t n, DSPType ** out)
        {
            DSPType tick;
            for(size_t i = 0; i < n; i++)
            {
                tick = osc->Tick();
                out[0][i] = tick * sin((1-pan)*M_PI/2);
                out[1][i] = tick * cos(pan*M_PI/2);
            }
        }
        DSPType Tick(DSPType iL, DSPType iR, DSPType & L, DSPType & R, DSPType A=1, DSPType X=1, DSPType Y=1)
        {
            DSPType r = osc->Tick(0.5*(iL+iR),A,X,Y);
            L = r * sin((1-pan)*M_PI/2);;
            R = r * cos(pan*M_PI/2);;
            return 0.5*(L+R);
        }
    };

    template<typename DSPType>
    struct GPUStereoGeneratorProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUGeneratorProcessor<DSPType> * osc;
        

        GPUStereoGeneratorProcessor(GeneratorProcessor * o) : GPUStereoProcessor<DSPType>() {
            osc = o;
            pan = 0.5;
        }
        ObjectType getType() const {
            return GPU_STEREO_GENERATOR_PROCESSOR;
        }
        void setPan(DSPType p ) {
            pan = p;
        }
        DSPType Tick(DSPType iL, DSPType iR, DSPType & L, DSPType & R, DSPType A=1, DSPType X=1, DSPType Y=1)
        {
            DSPType r = osc->Tick(0.5*(iL+iR),A,X,Y);
            L = r * sin((1-pan)*M_PI/2);
            R = r * cos(pan*M_PI/2);
            return 0.5*(iL+iR);
        }
        void ProcessBlock(size_t n, DSPType ** out)
        {
            DSPType tick;
            for(size_t i = 0; i < n; i++)
            {
                tick = osc->Tick();
                out[0][i] = tick * sin((1-pan)*M_PI/2);
                out[1][i] = tick * cos(pan*M_PI/2);
            }
        }
    };


    template<typename DSPType>
    struct GPUStereoFunctionProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUFunctionProcessor<DSPType> * filter[2];
        GPUStereoFunctionProcessor(GPUFunctionProcessor<DSPType> * L, GPUFunctionProcessor<DSPType> * R) 
        : GPUStereoProcessor<DSPType>()
        {
            filter[0] = L;
            filter[1] = R;
        }
        ObjectType getType() const {
            return GPU_STEREO_FUNCTION_PROCESSOR;
        }
        void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) {
            for(size_t i = 0; i < 2; i++) filter[i]->ProcessBlock(n,in[i],out[i]);
        }
        DSPType Tick(DSPType iL, DSPType iR, DSPType & L, DSPType & R, DSPType A=1, DSPType X=1, DSPType Y=1)
        {
            
            L = filter[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
            R = filter[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
            return 0.5*(L+R);
        }
    };


    template<typename DSPType>
    struct GPUStereoFilterProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUFilterProcessor<DSPType> * filter[2];
        GPUStereoFilterProcessor(GPUFilterProcessor<DSPType> * L, GPUFilterProcessor<DSPType> * R) 
        : GPUStereoProcessor<DSPType>()
        {
            filter[0] = L;
            filter[1] = R;
        }
        ObjectType getType() const {
            return GPU_STEREO_FILTER_PROCESSOR;
        }
        void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) {
            for(size_t i = 0; i < 2; i++) filter[i]->ProcessBlock(n,in[i],out[i]);
        }
        DSPType Tick(DSPType iL, DSPType iR, DSPType & L, DSPType & R, DSPType A=1, DSPType X=1, DSPType Y=1)
        {        
            L = filter[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
            R = filter[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
            return 0.5*(L+R);
        }
    };


    template<typename DSPType>
    struct GPUStereoAmplifierProcessor : public GPUStereoProcessor<DSPType>
    {
        GPUAmplifierProcessor<DSPType> * amp[2];
        GPUStereoAmplifierProcessor(GPUAmplifierProcessor<DSPType> * L, GPUAmplifierProcessor<DSPType> * R) 
        : GPUStereoProcessor<DSPType>()
        {
            amp[0] = L;
            amp[1] = R;
        }
        ObjectType getType() const {
            return GPU_STEREO_AMPLIFIER_PROCESSOR;
        }
        void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) {
            for(size_t i = 0; i < 2; i++) amp[i]->ProcessBlock(n,in[i],out[i]);
        }
        DSPType Tick(DSPType iL, DSPType iR, DSPType & L, DSPType & R, DSPType A=0, DSPType X=0, DSPType Y=0)
        {        
            L = amp[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
            R = amp[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
            return 0.5*(L+R);
        }
    };



    template<typename DSPType, int N>
    struct GPUFilterBankProcessor : public GPUFilterProcessor<DSPType>
    {
        GPUFilterProcessor<DSPType> * taps[N];

        GPUFilterBankProcessor() : GPUFilterProcessor<DSPType>()
        {
            memset(taps,0,N*sizeof(GPUFilterProcessor<DSPType>*));
        }
        ObjectType getType() const {
            return GPU_FILTER_BANK_PROCESSOR;
        }
        void setTap(size_t t, FilterProcessor * m) {
            taps[t] = m;
        }
        DSPType Tick(DSPType I, DSPType A=1, DSPType X=0, DSPType Y=0)
        {
            DSPType r = taps[0]->Tick(I,A,X,Y);
            for(size_t i = 0; i < N; i++) {
                if(taps[i]) r = taps[i]->Tick(I,A,X,Y);
            }
            return r;
        }
        void ProcessBlock(size_t numSamples, DSPType * in, DSPType * out)
        {
            if(taps[0]) taps[0]->ProcessBlock(numSamples,in,out);
            for(size_t i = 1; i < numSamples; i++)
                if(taps[i]) taps[i]->ProcessBlock(numSamples,out,out);            
        }
    };

    // these are all variously different
    // FFT : input, output, bins, magnitude, phase
    // Convolution: ProcessBlock, ProcessComplexBlock
    template<typename DSPType>
    struct GPUSpectrumProcessor : public GPUSoundProcessor<DSPType>
    {
    virtual ObjectType getType() const { return GPU_SPECTRUM_PROCESSOR; }   
    }; 


    template<typename DSPType, class BASE>
    struct GPUParameter2ProcessorPlugin : public BASE, public GPUParameter2Processor<DSPType>
    {

        GPUParameter2ProcessorPlugin() : BASE(), GPUParameter2Processor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType X, DSPType Y) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUAmplifierProcessorPlugin : public BASE, public GPUAmplifierProcessor<DSPType>
    {

        GPUAmplifierProcessorPlugin() : BASE(), GPUAmplifierProcessor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUFilterProcessorPlugin : public BASE, public GPUFilterProcessor<DSPType>
    {

        GPUFilterProcessorPlugin() : BASE(),GPUFilterProcessor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUOscillatorProcessorPlugin : public BASE, public GPUOscillatorProcessor<DSPType>
    {

        GPUOscillatorProcessorPlugin() : BASE(), GPUOscillatorProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUFunctionProcessorPlugin : public BASE, public GPUFunctionProcessor<DSPType>
    {

        GPUFunctionProcessorPlugin() : BASE(),GPUFunctionProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUGeneratorProcessorPlugin : public BASE, public GPUGeneratorProcessor<DSPType>
    {

        GPUGeneratorProcessorPlugin() : BASE(), GPUGeneratorProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct GPUMonoProcessorPlugin : public BASE, public GPUMonoProcessor<DSPType>
    {

        GPUMonoProcessorPlugin() : BASE(), GPUMonoProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, TYPE * in, TYPE * out) = 0;
    };

    template<typename DSPType,class BASE>
    struct GPUStereoProcessorPlugin : public BASE, public GPUStereoProcessor<DSPType>
    {

        GPUStereoProcessorPlugin() : BASE(), GPUStereoProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

    };

    template<typename DSPType,class BASE>
    struct GPUMonoFXProcessorPlugin : public BASE, public GPUMonoFXProcessor<DSPType>
    {

        GPUMonoFXProcessorPlugin() : BASE(),GPUFXProcessor<DSPType>()
        {

        };
        
        virtual void ProcessBlock(size_t n, DSPType * in, DSPType * out) = 0;
    };

    template<typename DSPType,class BASE>
    struct GPUStereoFXProcessorPlugin : public BASE, public GPUStereoFXProcessor<DSPType>
    {

        GPUStereoFXProcessorPlugin() : BASE(), GPUStereoFXProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

    };



    template<typename DSPType,class BASE>
    struct GPUClassProcessor : public BASE
    {
        GPUClassProcessor() : BASE()
        {

        }    
    };



    // signal comes from outside (an audio stream or file)// it should be normalized [-1,1] and is assumed continuous and periodic
    template<typename DSPType>
    struct GPUSignalSourceProcessor : public GPUSoundProcessor<DSPType>
    {

    };


    // a signal sink records data (to a buffer or file)
    template<typename DSPType>
    struct GPUSignalSinkProcessor : public GPUSoundProcessor<DSPType>
    {

    };

    template<typename DSPType,class BASE>
    struct GPUStereoSignalSourceProcessor : public BASE, public GPUSignalSourceProcessor<DSPType>
    {

        GPUStereoSignalSourceProcessor() : BASE(), GPUSignalSourceProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

    };

    template<typename DSPType,class BASE>
    struct GPUStereoSignalSinkProcessor : public BASE, public GPUSignalSinkProcessor<DSPType>
    {

        GPUStereoSignalSinkProcessor() : BASE(), GPUSignalSinkProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in) = 0;
    };


    template<typename DSPType>
    struct GPUMonoOversampleProcessor : public GPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return GPU_MONO_OVERSAMPLE_PROCESSOR;
        }

        
    };

    template<typename DSPType>
    struct GPUMonoUpsampleProcessor : public GPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return GPU_MONO_UPSAMPLE_PROCESSOR;
        }
    };

    template<typename DSPType>
    struct GPUMonoDownsampleProcessor : public GPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return GPU_MONO_DOWNSAMPLE_PROCESSOR;
        }
    };
}
