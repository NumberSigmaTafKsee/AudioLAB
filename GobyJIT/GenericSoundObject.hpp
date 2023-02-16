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




struct Random
{
    Random(){ srand(time(NULL)); }
    DSPType     frand() { return ((DSPType)rand()/(DSPType)RAND_MAX); }
    DSPType     rand() { return ((DSPType)rand()/(DSPType)RAND_MAX); }
    uint64_t    randint(int min, int max) { return round((max-min)*frand())+min; }
    bool        flip(DSPType prob) { return frand() < prob; }
    uint64_t    random(int mod) { return rand() % mod; }    
};

enum ObjectType
{    
    
    GS_PARAMETER2_PROCESSOR,     // F(X,Y)
    GS_STEREOSPLITTER_PROCESSOR,
    
    GS_MONO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
    GS_MONO_FUNCTION_PROCESSOR,       // Input    
    GS_MONO_OSCILLATOR_PROCESSOR,     // Generator no input
    GS_MONO_FILTER_PROCESSOR,         // DSPType Input
    GS_MONO_AMPLIFIER_PROCESSOR,    
    GS_MONO_FX_PROCESSOR,    
    GS_MONO_CASCADE_PROCESSOR,
    GS_MONO_MIXER_PROCESSOR,
    GS_MONO_MORPHER_PROCESSOR,
    GS_MONO_OPERATOR_PROCESSOR,
    GS_MONO_OVERSAMPLE_PROCESSOR,    // Resampler up/fx/down
    GS_MONO_UPSAMPLE_PROCESSOR,      // resample up
    GS_MONO_DOWNSAMPLE_PROCESSOR,    // resample down

    GS_STEREO_FX_PROCESSOR,
    GS_STEREO_OVERSAMPLER_PROCESSOR,    // Resampler up/fx/down
    GS_STEREO_UPSAMPLER_PROCESSOR,      // resample up
    GS_STEREO_DOWNSAMPLER_PROCESSOR,    // resample down
    GS_STEREO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
    GS_STEREO_FUNCTION_PROCESSOR,       // Input
    GS_STEREO_PARAMETER2_PROCESSOR,     // F(X,Y)
    GS_STEREO_OSCILLATOR_PROCESSOR,     // Generator no input
    GS_STEREO_FILTER_PROCESSOR,         // DSPType Input
    GS_STEREO_AMPLIFIER_PROCESSOR,
    GS_STEREO_CASCADE_PROCESSOR,
    GS_STEREO_OPERATOR_PROCESSOR,

    GS_MONO_SIGNAL_SOURCE_PROCESSOR,
    GS_STEREO_SIGNAL_SOURCE_PROCESSOR,
    GS_MONO_SIGNAL_SINK_PROCESSOR,
    GS_STEREO_SIGNAL_SINK_PROCESSOR,
    GS_INTERLEAVE_PROCESSOR,
    GS_DEINERLEAVE_PROCESSOR,
    GS_FILTER_BANK_PROCESSOR,
    GS_SPECTRUM_PROCESSOR,    
};


template<typename DSPType>
struct GSSoundProcessor
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
struct GSPort
{
    int port;
    DSPType value;
    GSSoundProcessor<DSPType> * p;
};

template<typename DSPType>
struct GSPorts
{
    std::list<std::shared_ptr<GSPort<DSPType>> ports;

    using PortMap = std::map<std::string,GSPort<DSPType>*>;
    PortMap portmap;

    GSPorts() {

    }
    void addPort(const std::string & name, GSPort<DSPType> * p) {
        ports.push_back( std::shared_ptr<Port>(p, [](GSPort<DSPType> * p){ delete p; }));
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
struct GSMonoProcessor : public GSSoundProcessor<DSPType>
{
    
    GSMonoProcessor() : GSSoundProcessortemplate<typename DSPType>()
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
struct GSMonoCascadeProcessor : public GSMonoProcessor<DSPType>
{
    std::list<MonoProcessor*> procs;

    GSMonoCascadeProcessor() : GSMonoProcessor<DSPType>()
    {

    }

    ObjectType getType() const {
        return GS_MONO_CASCADE_PROCESSOR;
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

    GSStereoProcessor() : GSSoundProcessor<DSPType>() {
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
struct GSStereoCascadeProcessor : public GSStereoProcessor<DSPType>
{
    std::list<GSStereoProcessor<DSPType>*> procs;

    GSStereoCascadeProcessor() : GSStereoProcessor<DSPType>()
    {

    }
    virtual ObjectType getType() const { return GS_STEREO_CASCADE_PROCESSOR; }

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
struct GSMonoFXProcessor : public GSMonoProcessor<DSPType>
{ 
    GSMonoFXProcessor() : GSMonoProcessor<DSPType>() {
 
    }

    ObjectType getType() const {
        return GS_MONO_FX_PROCESSOR;
    }
    
    virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) = 0;
};

template<typename DSPType>
struct GSStereoFXProcessor : public GSStereoProcessor<DSPType>
{
    GSStereoFXProcessor() : GSStereoProcessor<DSPType>() {
 
    }
    
    ObjectType getType() const {
        return GS_STEREO_FX_PROCESSOR;
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
struct GSGeneratorProcessor : public GSMonoProcessor<DSPType>
{
    GSGeneratorProcessor() : GSMonoProcessor<DSPType>() 
    {
        
    }
    ObjectType getType() const {
        return GS_MONO_GENERATOR_PROCESSOR;
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
struct GSMixerProcessor : public GSSoundProcessor<DSPType>
{
    GSMixerProcessor() : GSSoundProcessor<DSPType>()
    {

    }
    ObjectType getType() const {
        return GS_MONO_MIXER_PROCESSOR;
    }

    virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType *) {

    }
    virtual void ProcessBlock(size_t n, size_t block, DSPType **, DSPType **) {
        
    }

};

template<typename DSPType>
struct GSFunctionProcessor : public GSMonoProcessor<DSPType>
{
    GSFunctionProcessor() : GSMonoProcessor<DSPType>() 
    {
 
    }

    ObjectType getType() const {
        return GS_MONO_FUNCTION_PROCESSOR;
    }

    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=0, DSPType Y=0) = 0;

    void ProcessBlock(size_t n, DSPType * input, DSPType * output) {
        for(size_t i = 0; i < n; i++)
            output[i] = Tick(input[i]);
    }
};

template<typename DSPType>
struct GSParameter2Processor : public GSSoundProcessor<DSPType>
{
    GSParameter2Processor() : GSSoundProcessor<DSPType>() 
    {
 
    }

    ObjectType getType() const {
        return GS_PARAMETER2_PROCESSOR;
    }

    virtual DSPType Tick(DSPType a, DSPType b) = 0;

    void ProcessBlock(size_t n, DSPType * x, DSPType * y, DSPType * output) {
        for(size_t i = 0; i < n; i++)
            output[i] = Tick(x[i],y[i]);
    }
};

template<typename DSPType>
struct GSStereoSplitterProcessor : public GSSoundProcessor<DSPType>
{
    GSStereoSplitterProcessor() : GSSoundProcessor<DSPType>() 
    {
 
    }

    ObjectType getType() const {
        return GS_STEREOSPLITTER_PROCESSOR;
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
struct GSOscillatorProcessor : public GSMonoProcessor<DSPType>
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

    GSOscillatorProcessor() : GSMonoProcessor<DSPType>() 
    {
        
    }
    ObjectType getType() const {
        return GS_MONO_OSCILLATOR_PROCESSOR;
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
struct GSFilterProcessor : public MonoProcessor<DSPType>
{
    GSFilterProcessor() : GSMonoProcessor<DSPType>() {
 
    }
    ObjectType getType() const {
        return GS_MONO_FILTER_PROCESSOR;
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
struct GSAmplifierProcessor : public GSMonoProcessor<DSPType>
{
    
    GSAmplifierProcessor() : GSMonoProcessor<DSPType>()
    {
    
    }
    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1)=0;

    ObjectType getType() const {
        return GS_MONO_AMPLIFIER_PROCESSOR;
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
struct GSStereoOscillatorProcessor : public GSStereoProcessor<DSPType>
{
    GSOscillatorProcessor<DSPType> * osc;
    
    GSStereoOscillatorProcessor(OscillatorProcessor * o) : GSStereoProcessor<DSPType>() {
        osc = o;
        pan = 0.5;
    }
    ObjectType getType() const {
        return GS_STEREO_OSCILLATOR_PROCESSOR;
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
struct GSStereoGeneratorProcessor : public GSStereoProcessor<DSPType>
{
    GSGeneratorProcessor<DSPType> * osc;
    

    GSStereoGeneratorProcessor(GeneratorProcessor * o) : GSStereoProcessor<DSPType>() {
        osc = o;
        pan = 0.5;
    }
    ObjectType getType() const {
        return GS_STEREO_GENERATOR_PROCESSOR;
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
struct GSStereoFunctionProcessor : public GSStereoProcessor<DSPType>
{
    GSFunctionProcessor<DSPType> * filter[2];
    GSStereoFunctionProcessor(GSFunctionProcessor<DSPType> * L, GSFunctionProcessor<DSPType> * R) 
    : GSStereoProcessor<DSPType>()
    {
        filter[0] = L;
        filter[1] = R;
    }
    ObjectType getType() const {
        return GS_STEREO_FUNCTION_PROCESSOR;
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
struct GSStereoFilterProcessor : public GSStereoProcessor<DSPType>
{
    GSFilterProcessor<DSPType> * filter[2];
    GSStereoFilterProcessor(GSFilterProcessor<DSPType> * L, GSFilterProcessor<DSPType> * R) 
    : GSStereoProcessor<DSPType>()
    {
        filter[0] = L;
        filter[1] = R;
    }
    ObjectType getType() const {
        return GS_STEREO_FILTER_PROCESSOR;
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
struct GSStereoAmplifierProcessor : public GSStereoProcessor<DSPType>
{
    GSAmplifierProcessor<DSPType> * amp[2];
    GSStereoAmplifierProcessor(GSAmplifierProcessor<DSPType> * L, GSAmplifierProcessor<DSPType> * R) 
    : GSStereoProcessor<DSPType>()
    {
        amp[0] = L;
        amp[1] = R;
    }
    ObjectType getType() const {
        return GS_STEREO_AMPLIFIER_PROCESSOR;
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
struct GSFilterBankProcessor : public GSFilterProcessor<DSPType>
{
    GSFilterProcessor<DSPType> * taps[N];

    GSFilterBankProcessor() : GSFilterProcessor<DSPType>()
    {
        memset(taps,0,N*sizeof(GSFilterProcessor<DSPType>*));
    }
    ObjectType getType() const {
        return GS_FILTER_BANK_PROCESSOR;
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
struct GSSpectrumProcessor : public GSSoundProcessor<DSPType>
{
   virtual ObjectType getType() const { return GS_SPECTRUM_PROCESSOR; }   
}; 


template<typename DSPType, class BASE>
struct GSParameter2ProcessorPlugin : public BASE, public GSParameter2Processor<DSPType>
{

    GSParameter2ProcessorPlugin() : BASE(), GSParameter2Processor<DSPType>()
    {

    };
    virtual DSPType Tick(DSPType X, DSPType Y) = 0;
};

template<typename DSPType, class BASE>
struct GSAmplifierProcessorPlugin : public BASE, public GSAmplifierProcessor<DSPType>
{

    GSAmplifierProcessorPlugin() : BASE(), GSAmplifierProcessor<DSPType>()
    {

    };
    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
};

template<typename DSPType, class BASE>
struct GSFilterProcessorPlugin : public BASE, public GSFilterProcessor<DSPType>
{

    GSFilterProcessorPlugin() : BASE(),GSFilterProcessor<DSPType>()
    {

    };
    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
};

template<typename DSPType, class BASE>
struct GSOscillatorProcessorPlugin : public BASE, public GSOscillatorProcessor<DSPType>
{

    GSOscillatorProcessorPlugin() : BASE(), GSOscillatorProcessor<DSPType>()
    {

    };

    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
};

template<typename DSPType, class BASE>
struct GSFunctionProcessorPlugin : public BASE, public GSFunctionProcessor<DSPType>
{

    GSFunctionProcessorPlugin() : BASE(),GSFunctionProcessor<DSPType>()
    {

    };

    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
};

template<typename DSPType, class BASE>
struct GSGeneratorProcessorPlugin : public BASE, public GSGeneratorProcessor<DSPType>
{

    GSGeneratorProcessorPlugin() : BASE(), GSGeneratorProcessor<DSPType>()
    {

    };

    virtual DSPType Tick(DSPType I, DSPType A=1, DSPType X=1, DSPType Y=1) = 0;
};

template<typename DSPType, class BASE>
struct GSMonoProcessorPlugin : public BASE, public GSMonoProcessor<DSPType>
{

    GSMonoProcessorPlugin() : BASE(), GSMonoProcessor<DSPType>()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE * in, TYPE * out) = 0;
};

template<typename DSPType,class BASE>
struct GSStereoProcessorPlugin : public BASE, public GSStereoProcessor<DSPType>
{

    GSStereoProcessorPlugin() : BASE(), GSStereoProcessor<DSPType>()
    {

    };

    virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

};

template<typename DSPType,class BASE>
struct GSMonoFXProcessorPlugin : public BASE, public GSMonoFXProcessor<DSPType>
{

    GSMonoFXProcessorPlugin() : BASE(),GSFXProcessor<DSPType>()
    {

    };
    
    virtual void ProcessBlock(size_t n, DSPType * in, DSPType * out) = 0;
};

template<typename DSPType,class BASE>
struct GSStereoFXProcessorPlugin : public BASE, public GSStereoFXProcessor<DSPType>
{

    GSStereoFXProcessorPlugin() : BASE(), GSStereoFXProcessor<DSPType>()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

};



template<typename DSPType,class BASE>
struct GSClassProcessor : public BASE
{
    GSClassProcessor() : BASE()
    {

    }    
};



// signal comes from outside (an audio stream or file)// it should be normalized [-1,1] and is assumed continuous and periodic
template<typename DSPType>
struct GSSignalSourceProcessor : public GSSoundProcessor<DSPType>
{

};


// a signal sink records data (to a buffer or file)
template<typename DSPType>
struct GSSignalSinkProcessor : public GSSoundProcessor<DSPType>
{

};

template<typename DSPType,class BASE>
struct GSStereoSignalSourceProcessor : public BASE, public GSSignalSourceProcessor<DSPType>
{

    GSStereoSignalSourceProcessor() : BASE(), GSSignalSourceProcessor<DSPType>()
    {

    };

    virtual void ProcessBlock(size_t n, DSPType ** in, DSPType ** out) = 0;

};

template<typename DSPType,class BASE>
struct GSStereoSignalSinkProcessor : public BASE, public GSSignalSinkProcessor<DSPType>
{

    GSStereoSignalSinkProcessor() : BASE(), GSSignalSinkProcessor<DSPType>()
    {

    };

    virtual void ProcessBlock(size_t n, DSPType ** in) = 0;
};


template<typename DSPType>
struct GSMonoOversampleProcessor : public GSSoundProcessor<DSPType>
{
    ObjectType getType() const {
        return GS_MONO_OVERSAMPLE_PROCESSOR;
    }

    
};

template<typename DSPType>
struct GSMonoUpsampleProcessor : public GSSoundProcessor<DSPType>
{
    ObjectType getType() const {
        return GS_MONO_UPSAMPLE_PROCESSOR;
    }
};

template<typename DSPType>
struct GSMonoDownsampleProcessor : public GSSoundProcessor<DSPType>
{
    ObjectType getType() const {
        return GS_MONO_DOWNSAMPLE_PROCESSOR;
    }
};


