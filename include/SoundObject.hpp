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
#ifdef DSPFLOATDOUBLE
typedef double DspFloatType;
#else
typedef float DspFloatType;
#endif

#include "Undenormal.hpp"
#include "StdNoise.hpp"
#include "MusicFunctions.hpp"
#include "FX/ClipFunctions.hpp"


/*
extern Default noise;
extern DspFloatType sampleRate;
extern DspFloatType invSampleRate;
*/
struct Random
{
    Random() {  }

    static void seed() { srand(time(NULL)); }
    DspFloatType      frand() { return ((DspFloatType)::rand()/(DspFloatType)RAND_MAX); }
    DspFloatType      rand() { return ((DspFloatType)::rand()/(DspFloatType)RAND_MAX); }
    uint64_t    randint(int min, int max) { return round((max-min)*frand())+min; }
    bool        flip(DspFloatType prob) { return frand() < prob; }
    uint64_t    random(int mod) { return ::rand() % mod; }    
};

enum ObjectType
{    
    
    PARAMETER2_PROCESSOR,     // F(X,Y)
    STEREOSPLITTER_PROCESSOR,
    
    MONO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
    MONO_FUNCTION_PROCESSOR,       // Input    
    MONO_OSCILLATOR_PROCESSOR,     // Generator no input
    MONO_FILTER_PROCESSOR,         // DspFloatType Input
    MONO_AMPLIFIER_PROCESSOR,    
    MONO_FX_PROCESSOR,    
    MONO_CASCADE_PROCESSOR,
    MONO_MIXER_PROCESSOR,
    MONO_MORPHER_PROCESSOR,
    MONO_OPERATOR_PROCESSOR,
    MONO_OVERSAMPLE_PROCESSOR,    // Resampler up/fx/down
    MONO_UPSAMPLE_PROCESSOR,      // resample up
    MONO_DOWNSAMPLE_PROCESSOR,    // resample down

    STEREO_FX_PROCESSOR,
    STEREO_OVERSAMPLER_PROCESSOR,    // Resampler up/fx/down
    STEREO_UPSAMPLER_PROCESSOR,      // resample up
    STEREO_DOWNSAMPLER_PROCESSOR,    // resample down
    STEREO_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
    STEREO_FUNCTION_PROCESSOR,       // Input
    STEREO_PARAMETER2_PROCESSOR,     // F(X,Y)
    STEREO_OSCILLATOR_PROCESSOR,     // Generator no input
    STEREO_FILTER_PROCESSOR,         // DspFloatType Input
    STEREO_AMPLIFIER_PROCESSOR,
    STEREO_CASCADE_PROCESSOR,
    STEREO_OPERATOR_PROCESSOR,

    MONO_SIGNAL_SOURCE_PROCESSOR,
    STEREO_SIGNAL_SOURCE_PROCESSOR,
    MONO_SIGNAL_SINK_PROCESSOR,
    STEREO_SIGNAL_SINK_PROCESSOR,
    INTERLEAVE_PROCESSOR,
    DEINERLEAVE_PROCESSOR,
    FILTER_BANK_PROCESSOR,
    SPECTRUM_PROCESSOR,    
};


// LuaJIT Functions
struct SoundProcessor
{    
    
    DspFloatType preGain = 1;
    DspFloatType postGain = 1;

    virtual ObjectType getType() const = 0;

    // i do not want any kind of complicated data structure
    // just a simple function to set the port value    
    virtual void setPort(int port, DspFloatType value) {
        printf("No port %d\n",port);
    }
    virtual void setPort2(int port, DspFloatType a, DspFloatType b) {
        printf("No port %d\n",port);
    }
    virtual void setPortV(int port, const std::vector<DspFloatType> & v) {
        printf("No port %d\n",port);
    }
    virtual DspFloatType getPort(int port) {
        printf("No port %d\n",port);
        return 0;
    }
    virtual DspFloatType getPort2(int port, DspFloatType v) {
        printf("No port %d\n",port);
        return 0;
    }
    virtual void getPortV(int port, std::vector<DspFloatType> & v) {
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

struct Port
{
    int port;
    DspFloatType value;
    SoundProcessor * p;
};

struct Ports
{
    std::list<std::shared_ptr<Port>> ports;

    using PortMap = std::map<std::string,Port*>;
    PortMap portmap;

    Ports() {

    }
    void addPort(const std::string & name, Port * p) {
        ports.push_back( std::shared_ptr<Port>(p, [](Port * p){ delete p; }));
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

struct MonoProcessor : public SoundProcessor
{
    
    MonoProcessor() : SoundProcessor()
    {

    }

    
    //virtual void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs, DspFloatType *A, DspFloatType *X, DspFloatType *Y) = 0;
        
    void InplaceProcess(size_t n, DspFloatType * buffer) {
        ProcessBlock(n,buffer,buffer);
    }

    virtual DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {
        assert(1==0);
    }
    virtual void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) 
    {
        assert(1==0);
    }

};


struct MonoCascadeProcessor : public MonoProcessor
{
    std::list<MonoProcessor*> procs;

    MonoCascadeProcessor() : MonoProcessor()
    {

    }

    ObjectType getType() const {
        return MONO_CASCADE_PROCESSOR;
    }

    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        auto i = procs.begin();
        memcpy(outputs,inputs,n*sizeof(DspFloatType));
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

struct StereoProcessor : public SoundProcessor
{
    
    DspFloatType pan;

    StereoProcessor() : SoundProcessor() {        
        pan = 0.5;
    }

    
    virtual void ProcessBlock(size_t n, DspFloatType ** inputs, DspFloatType ** outputs) 
    {
        assert(1==0);
    }

    virtual DspFloatType Tick(DspFloatType IL, DspFloatType IR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {
        assert(1==0);
    }   
    void InplaceProcess(size_t n, DspFloatType ** buffer) {
        ProcessBlock(n,buffer,buffer);
    }

    void Run(size_t n, DspFloatType ** inputs, DspFloatType ** outputs) {
        ProcessBlock(n,inputs,outputs);    
    }
    
};

struct StereoCascadeProcessor : public StereoProcessor
{
    std::list<StereoProcessor*> procs;

    StereoCascadeProcessor() : StereoProcessor()
    {

    }
    virtual ObjectType getType() const { return STEREO_CASCADE_PROCESSOR; }

    void ProcessBlock(size_t n, DspFloatType ** inputs, DspFloatType ** outputs) {
        auto i = procs.begin();
        memcpy(outputs[0],inputs[0],n*sizeof(DspFloatType));
        memcpy(outputs[1],inputs[1],n*sizeof(DspFloatType));
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
                    DspFloatType L,R;
                    (*i)->Tick(outputs[0][j],outputs[1][j],L,R);
                    outputs[0][j] = L;
                    outputs[1][j] = R;
                }
                break;
            case STEREO_OSCILLATOR_PROCESSOR:
                for(size_t j = 0; j < n; j++)
                {
                    DspFloatType L,R;
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
struct MonoFXProcessor : public MonoProcessor
{ 
    MonoFXProcessor() : MonoProcessor() {
 
    }

    ObjectType getType() const {
        return MONO_FX_PROCESSOR;
    }
    
    virtual void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) = 0;
};

struct StereoFXProcessor : public StereoProcessor
{
    StereoFXProcessor() : StereoProcessor() {
 
    }
    
    ObjectType getType() const {
        return STEREO_FX_PROCESSOR;
    }        

    virtual void ProcessBlock(size_t n, DspFloatType ** inputs, DspFloatType ** outputs) = 0;
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

struct GeneratorProcessor : public MonoProcessor
{
    GeneratorProcessor() : MonoProcessor() 
    {
        
    }
    ObjectType getType() const {
        return MONO_GENERATOR_PROCESSOR;
    }
    virtual DspFloatType Tick(DspFloatType I=0, DspFloatType A=0, DspFloatType X=0, DspFloatType Y=0) = 0;

    void Generate(size_t n, DspFloatType * output) {
        for(size_t i = 0; i < n; i++)
            output[i] = Tick();
    }
    virtual void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
};

struct MixerProcessor : public SoundProcessor
{
    MixerProcessor() : SoundProcessor()
    {

    }
    ObjectType getType() const {
        return MONO_MIXER_PROCESSOR;
    }

    virtual void ProcessBlock(size_t n, size_t block, DspFloatType **, DspFloatType *) {

    }
    virtual void ProcessBlock(size_t n, size_t block, DspFloatType **, DspFloatType **) {
        
    }

};

struct FunctionProcessor : public MonoProcessor
{
    FunctionProcessor() : MonoProcessor() 
    {
 
    }

    ObjectType getType() const {
        return MONO_FUNCTION_PROCESSOR;
    }

    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) = 0;

    void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            output[i] = Tick(input[i]);
    }
};

struct Parameter2Processor : public SoundProcessor
{
    Parameter2Processor() : SoundProcessor() 
    {
 
    }

    ObjectType getType() const {
        return PARAMETER2_PROCESSOR;
    }

    virtual DspFloatType Tick(DspFloatType a, DspFloatType b) = 0;

    void ProcessBlock(size_t n, DspFloatType * x, DspFloatType * y, DspFloatType * output) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            output[i] = Tick(x[i],y[i]);
    }
};

struct StereoSplitterProcessor : public SoundProcessor
{
    StereoSplitterProcessor() : SoundProcessor() 
    {
 
    }

    ObjectType getType() const {
        return STEREOSPLITTER_PROCESSOR;
    }

    virtual DspFloatType Tick(DspFloatType in, DspFloatType &a, DspFloatType &b) = 0;

    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * a, DspFloatType * b) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {   
            DspFloatType x = a[i];
            DspFloatType y = b[i];
            Tick(in[i],x,y);
        }
    }
};

struct OscillatorProcessor : public MonoProcessor
{    
    OscillatorProcessor() : MonoProcessor() 
    {
        
    }
    ObjectType getType() const {
        return MONO_OSCILLATOR_PROCESSOR;
    }

    virtual DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) = 0;    

    virtual void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
};

struct FilterProcessor : public MonoProcessor
{
    FilterProcessor() : MonoProcessor() {
 
    }
    ObjectType getType() const {
        return MONO_FILTER_PROCESSOR;
    }
    
    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X =0, DspFloatType Y=0)=0;

    // cascade
    DspFloatType Run(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
        DspFloatType x = Tick(I,A,X,Y);
        return x;
    }
 
    // cascade
    void ProcessBlock(size_t numSamples, DspFloatType * inputs, DspFloatType * outputs)
    {
        #pragma omp simd
        for(size_t i = 0; i < numSamples; i++)
        {
            outputs[i] = Run(inputs[i]);
        }        
    }
};

struct AmplifierProcessor : public MonoProcessor
{
    
    AmplifierProcessor() : MonoProcessor()
    {
    
    }
    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)=0;

    ObjectType getType() const {
        return MONO_AMPLIFIER_PROCESSOR;
    }
    
    // cascade
    void ProcessBlock(size_t numSamples, DspFloatType * inputs, DspFloatType * outputs)
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

struct StereoOscillatorProcessor : public StereoProcessor
{
    OscillatorProcessor * osc;
    
    StereoOscillatorProcessor(OscillatorProcessor * o) : StereoProcessor() {
        osc = o;
        pan = 0.5;
    }
    ObjectType getType() const {
        return STEREO_OSCILLATOR_PROCESSOR;
    }
    void setPan(DspFloatType p ) {
        pan = p;
    }
    void ProcessBlock(size_t n, DspFloatType ** out)
    {
        DspFloatType tick;
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            tick = osc->Tick();
            out[0][i] = tick * sin((1-pan)*M_PI/2);
            out[1][i] = tick * cos(pan*M_PI/2);
        }
    }
    DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & L, DspFloatType & R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
    {
        DspFloatType r = osc->Tick(0.5*(iL+iR),A,X,Y);
        L = r * sin((1-pan)*M_PI/2);;
        R = r * cos(pan*M_PI/2);;
        return 0.5*(L+R);
    }
};

struct StereoGeneratorProcessor : public StereoProcessor
{
    GeneratorProcessor * osc;
    

    StereoGeneratorProcessor(GeneratorProcessor * o) : StereoProcessor() {
        osc = o;
        pan = 0.5;
    }
    ObjectType getType() const {
        return STEREO_GENERATOR_PROCESSOR;
    }
    void setPan(DspFloatType p ) {
        pan = p;
    }
    DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & L, DspFloatType & R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
    {
        DspFloatType r = osc->Tick(0.5*(iL+iR),A,X,Y);
        L = r * sin((1-pan)*M_PI/2);
        R = r * cos(pan*M_PI/2);
        return 0.5*(iL+iR);
    }
    void ProcessBlock(size_t n, DspFloatType ** out)
    {
        DspFloatType tick;
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
        {
            tick = osc->Tick();
            out[0][i] = tick * sin((1-pan)*M_PI/2);
            out[1][i] = tick * cos(pan*M_PI/2);
        }
    }
};

struct StereoFunctionProcessor : public StereoProcessor
{
    FunctionProcessor * filter[2];
    StereoFunctionProcessor(FunctionProcessor * L, FunctionProcessor * R) 
    : StereoProcessor()
    {
        filter[0] = L;
        filter[1] = R;
    }
    ObjectType getType() const {
        return STEREO_FUNCTION_PROCESSOR;
    }
    void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
        for(size_t i = 0; i < 2; i++) filter[i]->ProcessBlock(n,in[i],out[i]);
    }
    DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & L, DspFloatType & R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
    {
        
        L = filter[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
        R = filter[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
        return 0.5*(L+R);
    }
};

struct StereoFilterProcessor : public StereoProcessor
{
    FilterProcessor * filter[2];
    StereoFilterProcessor(FilterProcessor * L, FilterProcessor * R) 
    : StereoProcessor()
    {
        filter[0] = L;
        filter[1] = R;
    }
    ObjectType getType() const {
        return STEREO_FILTER_PROCESSOR;
    }
    void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
        for(size_t i = 0; i < 2; i++) filter[i]->ProcessBlock(n,in[i],out[i]);
    }
    DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & L, DspFloatType & R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
    {        
        L = filter[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
        R = filter[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
        return 0.5*(L+R);
    }
};

struct StereoAmplifierProcessor : public StereoProcessor
{
    AmplifierProcessor * amp[2];
    StereoAmplifierProcessor(AmplifierProcessor * L, AmplifierProcessor * R) 
    : StereoProcessor()
    {
        amp[0] = L;
        amp[1] = R;
    }
    ObjectType getType() const {
        return STEREO_AMPLIFIER_PROCESSOR;
    }
    void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
        for(size_t i = 0; i < 2; i++) amp[i]->ProcessBlock(n,in[i],out[i]);
    }
    DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & L, DspFloatType & R, DspFloatType A=0, DspFloatType X=0, DspFloatType Y=0)
    {        
        L = amp[0]->Tick(iL,A,X,Y) * sin((1-pan)*M_PI/2);;
        R = amp[1]->Tick(iR,A,X,Y) * cos(pan*M_PI/2);;
        return 0.5*(L+R);
    }
};



template<int N>
struct FilterBankProcessor : public FilterProcessor
{
    FilterProcessor * taps[N];

    FilterBankProcessor() : FilterProcessor()
    {
        memset(taps,0,N*sizeof(FilterProcessor*));
    }
    ObjectType getType() const {
        return FILTER_BANK_PROCESSOR;
    }
    void setTap(size_t t, FilterProcessor * m) {
        taps[t] = m;
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {
        DspFloatType r = taps[0]->Tick(I,A,X,Y);
        for(size_t i = 0; i < N; i++) {
            if(taps[i]) r = taps[i]->Tick(I,A,X,Y);
        }
        return r;
    }
    void ProcessBlock(size_t numSamples, DspFloatType * in, DspFloatType * out)
    {
        if(taps[0]) taps[0]->ProcessBlock(numSamples,in,out);
        for(size_t i = 1; i < numSamples; i++)
            if(taps[i]) taps[i]->ProcessBlock(numSamples,out,out);            
    }
};

// these are all variously different
// FFT : input, output, bins, magnitude, phase
// Convolution: ProcessBlock, ProcessComplexBlock
struct SpectrumProcessor : public SoundProcessor
{
   virtual ObjectType getType() const { return SPECTRUM_PROCESSOR; }   
}; 


template<class BASE>
struct Parameter2ProcessorPlugin : public BASE, public Parameter2Processor
{

    Parameter2ProcessorPlugin() : BASE(),Parameter2Processor()
    {

    };
    virtual DspFloatType Tick(DspFloatType X, DspFloatType Y) = 0;
};

template<class BASE>
struct AmplifierProcessorPlugin : public BASE, public AmplifierProcessor
{

    AmplifierProcessorPlugin() : BASE(),AmplifierProcessor()
    {

    };
    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) = 0;
};

template<class BASE>
struct FilterProcessorPlugin : public BASE, public FilterProcessor
{

    FilterProcessorPlugin() : BASE(),FilterProcessor()
    {

    };
    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) = 0;
};

template<class BASE>
struct OscillatorProcessorPlugin : public BASE, public OscillatorProcessor
{

    OscillatorProcessorPlugin() : BASE(),OscillatorProcessor()
    {

    };

    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) = 0;
};

template<class BASE>
struct FunctionProcessorPlugin : public BASE, public FunctionProcessor
{

    FunctionProcessorPlugin() : BASE(),FunctionProcessor()
    {

    };

    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) = 0;
};

template<class BASE>
struct GeneratorProcessorPlugin : public BASE, public GeneratorProcessor
{

    GeneratorProcessorPlugin() : BASE(),GeneratorProcessor()
    {

    };

    virtual DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) = 0;
};

template<class BASE, typename TYPE=DspFloatType>
struct MonoProcessorPlugin : public BASE, public MonoProcessor
{

    MonoProcessorPlugin() : BASE(),MonoProcessor()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE * in, TYPE * out) = 0;
};

template<class BASE, typename TYPE=DspFloatType>
struct StereoProcessorPlugin : public BASE, public StereoProcessor
{

    StereoProcessorPlugin() : BASE(),StereoProcessor()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

};

template<class BASE, typename TYPE=DspFloatType>
struct MonoFXProcessorPlugin : public BASE, public MonoFXProcessor
{

    MonoFXProcessorPlugin() : BASE(),MonoFXProcessor()
    {

    };
    
    virtual void ProcessBlock(size_t n, TYPE * in, TYPE * out) = 0;
};

template<class BASE, typename TYPE=DspFloatType>
struct StereoFXProcessorPlugin : public BASE, public StereoFXProcessor
{

    StereoFXProcessorPlugin() : BASE(),StereoFXProcessor()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

};



template<typename BASE>
struct ClassProcessor : public BASE
{
    ClassProcessor() : BASE()
    {

    }    
};


/////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////
using SoundProcessorList = std::list<SoundProcessor*>;
using SoundProcessorVector = std::vector<SoundProcessor*>;
using SoundProcessorNameMap= std::map<std::string,SoundProcessor*>;

using MonoProcessorList = std::list<MonoProcessor*>;
using MonoProcessorVector = std::vector<MonoProcessor*>;
using MonoProcessorNameMap= std::map<std::string,MonoProcessor*>;

using StereoProcessorList = std::list<StereoProcessor*>;
using StereoProcessorVector = std::vector<StereoProcessor*>;
using StereoProcessorNameMap= std::map<std::string,StereoProcessor*>;

// signal comes from outside (an audio stream or file)// it should be normalized [-1,1] and is assumed continuous and periodic
struct SignalSourceProcessor : public SoundProcessor
{

};
// a signal sink records data (to a buffer or file)
struct SignalSinkProcessor : public SoundProcessor
{

};

template<class BASE, typename TYPE=DspFloatType>
struct StereoSignalSourceProcessor : public BASE, public SignalSourceProcessor
{

    StereoSignalSourceProcessor() : BASE(),SignalSourceProcessor()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE ** in, TYPE ** out) = 0;

};

template<class BASE, typename TYPE=DspFloatType>
struct StereoSignalSinkProcessor : public BASE, public SignalSinkProcessor
{

    StereoSignalSinkProcessor() : BASE(),SignalSinkProcessor()
    {

    };

    virtual void ProcessBlock(size_t n, TYPE ** in) = 0;
};


struct MonoOversampleProcessor : public SoundProcessor
{
    ObjectType getType() const {
        return MONO_OVERSAMPLE_PROCESSOR;
    }

    
};
struct MonoUpsampleProcessor : public SoundProcessor
{
    ObjectType getType() const {
        return MONO_UPSAMPLE_PROCESSOR;
    }
};
struct MonoDownsampleProcessor : public SoundProcessor
{
    ObjectType getType() const {
        return MONO_DOWNSAMPLE_PROCESSOR;
    }
};

struct OperatorProcessor : public SoundProcessor
{
    OperatorProcessor() : SoundProcessor()
    {

    }    
};
struct MonoOperatorProcessor : public SoundProcessor
{   
    MonoProcessorList list;

    MonoOperatorProcessor() : SoundProcessor()
    {

    }    
    ObjectType getType() const {
        return MONO_OPERATOR_PROCESSOR;
    }
    void GetSamples(MonoProcessor * p, size_t n, DspFloatType * input, DspFloatType * output) {
        memcpy(output,input,n*sizeof(DspFloatType));
        switch(p->getType())
        {
            case MONO_OSCILLATOR_PROCESSOR:
            case MONO_GENERATOR_PROCESSOR:
            {
                for(size_t i = 0; i < n; i++)
                    output[i] = p->Tick();
                break;
            }
            case MONO_FILTER_PROCESSOR:
            case MONO_FUNCTION_PROCESSOR:
            case MONO_AMPLIFIER_PROCESSOR:
                {
                for(size_t i = 0; i < n; i++)
                    output[i] = p->Tick(input[i]);
                break;
            }
            case MONO_FX_PROCESSOR:
                p->ProcessBlock(n,input,output);
                break;
        }
    }
};

struct MonoCascadeOperatorProcessor : public MonoOperatorProcessor
{
    MonoCascadeOperatorProcessor() : MonoOperatorProcessor()
    {

    }
    
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,outputs,outputs);            
        }
    }
};

struct MonoMultiplierOperatorProcessor : public MonoOperatorProcessor
{
    MonoMultiplierOperatorProcessor() : MonoOperatorProcessor()
    {

    }    
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,temp.data(),temp.data());
            for(size_t j = 0; j < n; j++)
                outputs[j] *= temp[j];
        }
    }
};

struct MonoAdderOperatorProcessor : public MonoOperatorProcessor
{
    MonoAdderOperatorProcessor() : MonoOperatorProcessor()
    {

    }    
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,temp.data(),temp.data());
            for(size_t j = 0; j < n; j++)
                outputs[j] += temp[j];
        }
    }
};

struct MonoSubtractorOperatorProcessor : public MonoOperatorProcessor
{
    MonoSubtractorOperatorProcessor() : MonoOperatorProcessor()
    {

    }   
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,temp.data(),temp.data());
            for(size_t j = 0; j < n; j++)
                outputs[j] -= temp[j];
        }
    }
};

struct MonoModulatorOperatorProcessor : public MonoOperatorProcessor
{
    MonoModulatorOperatorProcessor() : MonoOperatorProcessor()
    {

    }
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,temp.data(),temp.data());
            for(size_t j = 0; j < n; j++)
                outputs[j] = fmod(outputs[j],temp[j]);
        }
    }
};

struct MonoBlendOperatorProcessor : public MonoOperatorProcessor
{
    DspFloatType blend = 0.5;
    MonoBlendOperatorProcessor() : MonoOperatorProcessor()
    {

    }
    void ProcessBlock(size_t n, DspFloatType * inputs, DspFloatType * outputs) {
        std::vector<DspFloatType> temp(n);
        auto i = list.begin();
        GetSamples(*i,n,inputs,outputs);
        memcpy(outputs,temp.data(),n*sizeof(DspFloatType));
        for(i++; i != list.end(); i++)
        {
            GetSamples(*i,n,temp.data(),temp.data());
            #pragma omp simd
            for(size_t j = 0; j < n; j++)
                outputs[j] = outputs[j] + blend*(temp[j]-outputs[j]);
        }
    }
};


