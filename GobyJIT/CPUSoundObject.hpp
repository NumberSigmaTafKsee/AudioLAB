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
        
        CPU_PARAMETER2_PROCESSOR,     // F(X,Y)
        CPU_STEREOSPLITTER_PROCESSOR,
        
        CPU_MONO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
        CPU_MONO_FUNCTION_PROCESSOR,       // Input    
        CPU_MONO_OSCILLATOR_PROCESSOR,     // Generator no input
        CPU_MONO_FILTER_PROCESSOR,         // DSPType Input
        CPU_MONO_AMPLIFIER_PROCESSOR,    
        CPU_MONO_FX_PROCESSOR,    
        CPU_MONO_CASCADE_PROCESSOR,
        CPU_MONO_MIXER_PROCESSOR,
        CPU_MONO_MORPHER_PROCESSOR,
        CPU_MONO_OPERATOR_PROCESSOR,
        CPU_MONO_OVERSAMPLE_PROCESSOR,    // Resampler up/fx/down
        CPU_MONO_UPSAMPLE_PROCESSOR,      // resample up
        CPU_MONO_DOWNSAMPLE_PROCESSOR,    // resample down

        CPU_STEREO_FX_PROCESSOR,
        CPU_STEREO_OVERSAMPLER_PROCESSOR,    // Resampler up/fx/down
        CPU_STEREO_UPSAMPLER_PROCESSOR,      // resample up
        CPU_STEREO_DOWNSAMPLER_PROCESSOR,    // resample down
        CPU_STEREO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
        CPU_STEREO_FUNCTION_PROCESSOR,       // Input
        CPU_STEREO_PARAMETER2_PROCESSOR,     // F(X,Y)
        CPU_STEREO_OSCILLATOR_PROCESSOR,     // Generator no input
        CPU_STEREO_FILTER_PROCESSOR,         // DSPType Input
        CPU_STEREO_AMPLIFIER_PROCESSOR,
        CPU_STEREO_CASCADE_PROCESSOR,
        CPU_STEREO_OPERATOR_PROCESSOR,

        CPU_MONO_SIGNAL_SOURCE_PROCESSOR,
        CPU_STEREO_SIGNAL_SOURCE_PROCESSOR,
        CPU_MONO_SIGNAL_SINK_PROCESSOR,
        CPU_STEREO_SIGNAL_SINK_PROCESSOR,
        CPU_INTERLEAVE_PROCESSOR,
        CPU_DEINERLEAVE_PROCESSOR,
        CPU_FILTER_BANK_PROCESSOR,
        CPU_SPECTRUM_PROCESSOR,    
    };


    template<typename DSPType>
    struct CPUSoundProcessor
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
    struct CPUPort
    {
        int port;
        DSPType value;
        CPUSoundProcessor<DSPType> * p;
    };

    template<typename DSPType>
    struct CPUPorts
    {
        std::list<std::shared_ptr<CPUPort<DSPType>> ports;

        using PortMap = std::map<std::string,CPUPort<DSPType>*>;
        PortMap portmap;

        CPUPorts() {

        }
        void addPort(const std::string & name, CPUPort<DSPType> * p) {
            ports.push_back( std::shared_ptr<Port>(p, [](CPUPort<DSPType> * p){ delete p; }));
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
    struct CPUMonoProcessor : public CPUSoundProcessor<DSPType>
    {
        
        CPUMonoProcessor() : CPUSoundProcessortemplate<typename DSPType>()
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
    struct CPUMonoCascadeProcessor : public CPUMonoProcessor<DSPType>
    {
        std::list<MonoProcessor*> procs;

        CPUMonoCascadeProcessor() : CPUMonoProcessor<DSPType>()
        {

        }

        ObjectType getType() const {
            return CPU_MONO_CASCADE_PROCESSOR;
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

        CPUStereoProcessor() : CPUSoundProcessor<DSPType>() {
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
    struct CPUStereoCascadeProcessor : public CPUStereoProcessor<DSPType>
    {
        std::list<CPUStereoProcessor<DSPType>*> procs;

        CPUStereoCascadeProcessor() : CPUStereoProcessor<DSPType>()
        {

        }
        virtual ObjectType getType() const { return CPU_STEREO_CASCADE_PROCESSOR; }

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
    struct CPUMonoFXProcessor : public CPUMonoProcessor<DSPType>
    { 
        CPUMonoFXProcessor() : CPUMonoProcessor<DSPType>() {
    
        }

        ObjectType getType() const {
            return CPU_MONO_FX_PROCESSOR;
        }
        
        virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) = 0;
    };

    template<typename DSPType>
    struct CPUStereoFXProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUStereoFXProcessor() : CPUStereoProcessor<DSPType>() {
    
        }
        
        ObjectType getType() const {
            return CPU_STEREO_FX_PROCESSOR;
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
    struct CPUGeneratorProcessor : public CPUMonoProcessor<DSPType>
    {
        CPUGeneratorProcessor() : CPUMonoProcessor<DSPType>() 
        {
            
        }
        ObjectType getType() const {
            return CPU_MONO_GENERATOR_PROCESSOR;
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
    struct CPUMixerProcessor : public CPUSoundProcessor<DSPType>
    {
        CPUMixerProcessor() : CPUSoundProcessor<DSPType>()
        {

        }
        ObjectType getType() const {
            return CPU_MONO_MIXER_PROCESSOR;
        }

        virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType *) {

        }
        virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType **) {
            
        }

    };

    template<typename DSPType>
    struct CPUFunctionProcessor : public CPUMonoProcessor<DSPType>
    {
        CPUFunctionProcessor() : CPUMonoProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return CPU_MONO_FUNCTION_PROCESSOR;
        }

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=0, DSPType Y=0) = 0;

        void ProcessBlock(size_t n, DSPType * input, DSPType * output) {
            for(size_t i = 0; i < n; i++)
                output[i] = Tick(input[i]);
        }
    };

    template<typename DSPType>
    struct CPUParameter2Processor : public CPUSoundProcessor<DSPType>
    {
        CPUParameter2Processor() : CPUSoundProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return CPU_PARAMETER2_PROCESSOR;
        }

        virtual DSPType Tick(DSPType a, DSPType b) = 0;

        void ProcessBlock(size_t n, DSPType * x, DSPType * y, DSPType * output) {
            for(size_t i = 0; i < n; i++)
                output[i] = Tick(x[i],y[i]);
        }
    };

    template<typename DSPType>
    struct CPUStereoSplitterProcessor : public CPUSoundProcessor<DSPType>
    {
        CPUStereoSplitterProcessor() : CPUSoundProcessor<DSPType>() 
        {
    
        }

        ObjectType getType() const {
            return CPU_STEREOSPLITTER_PROCESSOR;
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
    struct CPUOscillatorProcessor : public CPUMonoProcessor<DSPType>
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

        CPUOscillatorProcessor() : CPUMonoProcessor<DSPType>() 
        {
            
        }
        ObjectType getType() const {
            return CPU_MONO_OSCILLATOR_PROCESSOR;
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
    struct CPUFilterProcessor : public MonoProcessor<DSPType>
    {
        CPUFilterProcessor() : CPUMonoProcessor<DSPType>() {
    
        }
        ObjectType getType() const {
            return CPU_MONO_FILTER_PROCESSOR;
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
    struct CPUAmplifierProcessor : public CPUMonoProcessor<DSPType>
    {
        
        CPUAmplifierProcessor() : CPUMonoProcessor<DSPType>()
        {
        
        }
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1)=0;

        ObjectType getType() const {
            return CPU_MONO_AMPLIFIER_PROCESSOR;
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
    struct CPUStereoOscillatorProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUOscillatorProcessor<DSPType> * osc;
        
        CPUStereoOscillatorProcessor(OscillatorProcessor * o) : CPUStereoProcessor<DSPType>() {
            osc = o;
            pan = 0.5;
        }
        ObjectType getType() const {
            return CPU_STEREO_OSCILLATOR_PROCESSOR;
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
    struct CPUStereoGeneratorProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUGeneratorProcessor<DSPType> * osc;
        

        CPUStereoGeneratorProcessor(GeneratorProcessor * o) : CPUStereoProcessor<DSPType>() {
            osc = o;
            pan = 0.5;
        }
        ObjectType getType() const {
            return CPU_STEREO_GENERATOR_PROCESSOR;
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
    struct CPUStereoFunctionProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUFunctionProcessor<DSPType> * filter[2];
        CPUStereoFunctionProcessor(CPUFunctionProcessor<DSPType> * L, CPUFunctionProcessor<DSPType> * R) 
        : CPUStereoProcessor<DSPType>()
        {
            filter[0] = L;
            filter[1] = R;
        }
        ObjectType getType() const {
            return CPU_STEREO_FUNCTION_PROCESSOR;
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
    struct CPUStereoFilterProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUFilterProcessor<DSPType> * filter[2];
        CPUStereoFilterProcessor(CPUFilterProcessor<DSPType> * L, CPUFilterProcessor<DSPType> * R) 
        : CPUStereoProcessor<DSPType>()
        {
            filter[0] = L;
            filter[1] = R;
        }
        ObjectType getType() const {
            return CPU_STEREO_FILTER_PROCESSOR;
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
    struct CPUStereoAmplifierProcessor : public CPUStereoProcessor<DSPType>
    {
        CPUAmplifierProcessor<DSPType> * amp[2];
        CPUStereoAmplifierProcessor(CPUAmplifierProcessor<DSPType> * L, CPUAmplifierProcessor<DSPType> * R) 
        : CPUStereoProcessor<DSPType>()
        {
            amp[0] = L;
            amp[1] = R;
        }
        ObjectType getType() const {
            return CPU_STEREO_AMPLIFIER_PROCESSOR;
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
    struct CPUFilterBankProcessor : public CPUFilterProcessor<DSPType>
    {
        CPUFilterProcessor<DSPType> * taps[N];

        CPUFilterBankProcessor() : CPUFilterProcessor<DSPType>()
        {
            memset(taps,0,N*sizeof(CPUFilterProcessor<DSPType>*));
        }
        ObjectType getType() const {
            return CPU_FILTER_BANK_PROCESSOR;
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
    struct CPUSpectrumProcessor : public CPUSoundProcessor<DSPType>
    {
    virtual ObjectType getType() const { return CPU_SPECTRUM_PROCESSOR; }   
    }; 


    template<typename DSPType, class BASE>
    struct CPUParameter2ProcessorPlugin : public BASE, public CPUParameter2Processor<DSPType>
    {

        CPUParameter2ProcessorPlugin() : BASE(), CPUParameter2Processor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType X, DSPType Y) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUAmplifierProcessorPlugin : public BASE, public CPUAmplifierProcessor<DSPType>
    {

        CPUAmplifierProcessorPlugin() : BASE(), CPUAmplifierProcessor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUFilterProcessorPlugin : public BASE, public CPUFilterProcessor<DSPType>
    {

        CPUFilterProcessorPlugin() : BASE(),CPUFilterProcessor<DSPType>()
        {

        };
        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUOscillatorProcessorPlugin : public BASE, public CPUOscillatorProcessor<DSPType>
    {

        CPUOscillatorProcessorPlugin() : BASE(), CPUOscillatorProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUFunctionProcessorPlugin : public BASE, public CPUFunctionProcessor<DSPType>
    {

        CPUFunctionProcessorPlugin() : BASE(),CPUFunctionProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUGeneratorProcessorPlugin : public BASE, public CPUGeneratorProcessor<DSPType>
    {

        CPUGeneratorProcessorPlugin() : BASE(), CPUGeneratorProcessor<DSPType>()
        {

        };

        virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
    };

    template<typename DSPType, class BASE>
    struct CPUMonoProcessorPlugin : public BASE, public CPUMonoProcessor<DSPType>
    {

        CPUMonoProcessorPlugin() : BASE(), CPUMonoProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, TYPE * in, TYPE * out) = 0;
    };

    template<typename DSPType,class BASE>
    struct CPUStereoProcessorPlugin : public BASE, public CPUStereoProcessor<DSPType>
    {

        CPUStereoProcessorPlugin() : BASE(), CPUStereoProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

    };

    template<typename DSPType,class BASE>
    struct CPUMonoFXProcessorPlugin : public BASE, public CPUMonoFXProcessor<DSPType>
    {

        CPUMonoFXProcessorPlugin() : BASE(),CPUFXProcessor<DSPType>()
        {

        };
        
        virtual void ProcessBlock(size_t n, DSPType * in, DSPType * out) = 0;
    };

    template<typename DSPType,class BASE>
    struct CPUStereoFXProcessorPlugin : public BASE, public CPUStereoFXProcessor<DSPType>
    {

        CPUStereoFXProcessorPlugin() : BASE(), CPUStereoFXProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

    };



    template<typename DSPType,class BASE>
    struct CPUClassProcessor : public BASE
    {
        CPUClassProcessor() : BASE()
        {

        }    
    };



    // signal comes from outside (an audio stream or file)// it should be normalized [-1,1] and is assumed continuous and periodic
    template<typename DSPType>
    struct CPUSignalSourceProcessor : public CPUSoundProcessor<DSPType>
    {

    };


    // a signal sink records data (to a buffer or file)
    template<typename DSPType>
    struct CPUSignalSinkProcessor : public CPUSoundProcessor<DSPType>
    {

    };

    template<typename DSPType,class BASE>
    struct CPUStereoSignalSourceProcessor : public BASE, public CPUSignalSourceProcessor<DSPType>
    {

        CPUStereoSignalSourceProcessor() : BASE(), CPUSignalSourceProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

    };

    template<typename DSPType,class BASE>
    struct CPUStereoSignalSinkProcessor : public BASE, public CPUSignalSinkProcessor<DSPType>
    {

        CPUStereoSignalSinkProcessor() : BASE(), CPUSignalSinkProcessor<DSPType>()
        {

        };

        virtual void ProcessBlock(size_t n, DSPType ** in) = 0;
    };


    template<typename DSPType>
    struct CPUMonoOversampleProcessor : public CPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return CPU_MONO_OVERSAMPLE_PROCESSOR;
        }

        
    };

    template<typename DSPType>
    struct CPUMonoUpsampleProcessor : public CPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return CPU_MONO_UPSAMPLE_PROCESSOR;
        }
    };

    template<typename DSPType>
    struct CPUMonoDownsampleProcessor : public CPUSoundProcessor<DSPType>
    {
        ObjectType getType() const {
            return CPU_MONO_DOWNSAMPLE_PROCESSOR;
        }
    };
};
