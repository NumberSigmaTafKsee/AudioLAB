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

#include "vector2d.hpp"
#include "vector4d.hpp"
#include "vector4f.hpp"
#include "vector8f.hpp"

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
        
        CPU_GENERATOR_PROCESSOR,      // Generator LFO,Envelope, Noise no input
        CPU_FUNCTION_PROCESSOR,       // Input    
        CPU_OSCILLATOR_PROCESSOR,     // Generator no input
        CPU_FILTER_PROCESSOR,         // DSPType Input
        CPU_AMPLIFIER_PROCESSOR,    
        CPU_FX_PROCESSOR,    
        
        CPU_CASCADE_PROCESSOR,
        CPU_MIXER_PROCESSOR,
        CPU_MORPHER_PROCESSOR,        
        CPU_OVERSAMPLE_PROCESSOR,    // Resampler up/fx/down
        CPU_UPSAMPLE_PROCESSOR,      // resample up
        CPU_DOWNSAMPLE_PROCESSOR,    // resample down

        
        CPU_SIGNAL_SOURCE_PROCESSOR,        
        CPU_SIGNAL_SINK_PROCESSOR,        
        CPU_INTERLEAVE_PROCESSOR,
        CPU_DEINTERLEAVE_PROCESSOR,
        CPU_FILTER_BANK_PROCESSOR,
        CPU_SPECTRUM_PROCESSOR,   

        CPU_WAVEGUIDE_PROCESSOR,
        CPU_WAVE_DIGITAL_FILTER_PROCESSOR,
        CPU_RTSPICE_PROCESSOR, 

        CPU_FFT_FILTER_PROCESSOR,
        CPU_WAVELET_FILTER_PROCESSOR,
        CPU_CONVOLUTION_FILTER_PROCESSOR,
        CPU_AUTOCORRELATION_FILTER_PROCESSOR,
        CPU_CROSSCORRELATION_FILTER_PROCESSOR,
        CPU_FIR_FILTER_PROCESSOR,
        CPU_IIR_FILTER_PROCESSOR,
        CPU_RESAMPLER_PROCESSOR,
        CPU_FFT_DIGITAL_FILTER_PROCESSOR,
    };


    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUSoundProcessor
    {    
        
        using SIMD = DSPType::SIMD;
        SIMD voices = 0;
        Float preGain  = 1.0;
        Float postGain = 1.0;

        virtual ObjectType getType() const = 0;

        // i do not want any kind of complicated data structure
        // just a simple function to set the port value    
        virtual void setPort(int port, Float value) {
            printf("No port %d\n",port);
        }
        virtual void setPort2(int port, Float a, Float b) {
            printf("No port %d\n",port);
        }
        virtual void setPortV(int port, const std::vector<Float> & v) {
            printf("No port %d\n",port);
        }
        virtual Float getPort(int port) {
            printf("No port %d\n",port);
            return 0;
        }
        virtual Float getPort2(int port, Float v) {
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

        virtual SIMD Tick(SIMD I=1, Float A=1, Float X=0, Float Y=0)
        {
            assert(1==0);
        }
        virtual void ProcessBlock(size_t n, Float * inputs, Float * outputs) 
        {
            for(size_t i = 0; i < n; i +=  bump)
            {
                SIMD r;
                r.load_a(inputs+i);
                r = Tick(r);
                r.store_a(outputs+i);
            }
        }
        void InplaceProcess(size_t n, Float * buffer) {
            ProcessBlock(n,buffer,buffer);
        }
        virtual void setVoices(SIMD v) {
            voices = v;
        }

    };

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUPort
    {
        using SIMD = DSPType::SIMD;
        int port;
        Float value;
        CPUSoundProcessor<DSPType,Float,bump> * p;
    };

    template<typename DSPType=Vector4f, typename Float=float,int bump=4>
    struct CPUPorts
    {
        using SIMD = DSPType::SIMD;
        std::list<std::shared_ptr<CPUPort<DSPType,Float,bump>> ports;
        using PortMap = std::map<std::string,CPUPort<DSPType,Float,bump>*>;
        PortMap portmap;

        CPUPorts() {

        }
        void addPort(const std::string & name, CPUPort<DSPType,Float,bump> * p) {
            ports.push_back( std::shared_ptr<Port>(p, [](CPUPort<DSPType,Float,bump> * p){ delete p; }));
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

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUMonoProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        using SIMD = DSPType::SIMD;

        CPUMonoProcessor() : CPUSoundProcessortemplate<DSPType,Float,bump>()
        {

        }        
        

    };


    
    // FX Processors should always be a stream and should not use Tick
    // It is a defect that they have Tick and will be eventually removed

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUFXProcessor : public CPUSoundProcessor<DSPType>
    { 
        using SIMD = DSPType::SIMD;

        CPUFXProcessor() : CPUSoundProcessor<DSPType,Float,bump>() {
    
        
        ObjectType getType() const {
            return CPU_FX_PROCESSOR;
        }            
    };
    
    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUGeneratorProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        CPUGeneratorProcessor() : CPUSoundProcessor<DSPType,Float,bump>() 
        {
            
        }
        ObjectType getType() const {
            return CPU_GENERATOR_PROCESSOR;
        }
        virtual SIMD Tick(SIMD I=0, Float A=0, Float X=0, Float Y=0) = 0;

        void Generate(size_t n, Float * output) {
            for(size_t i = 0; i < n; i+=bump)
            {
                SIMD o;
                o = Tick();
                o.store_a(output+i);
            }                
        }        
    };

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUMixerProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        using SIMD = DSPType::SIMD;
        CPUMixerProcessor() : CPUSoundProcessor<DSPType,Float,bump>()
        {

        }
        ObjectType getType() const {
            return CPU_MIXER_PROCESSOR;
        }

        virtual void ProcessBlock(size_t n, size_t block, Float **, Float *) {

        }
        virtual void ProcessBlock(size_t n, size_t block, Float **, Float **) {
            
        }

    };

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUFunctionProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        using SIMD = DSPType::SIMD;
        CPUFunctionProcessor() : CPUSoundProcessor<DSPType>() 
        {
    
        }
        ObjectType getType() const {
            return CPU_FUNCTION_PROCESSOR;
        }        
    };

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUParameter2Processor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        CPUParameter2Processor() : CPUSoundProcessor<DSPType,Float,bump>() 
        {
    
        }

        ObjectType getType() const {
            return CPU_PARAMETER2_PROCESSOR;
        }

        virtual SIMD Tick(SIMD a, SIMD b) = 0;

        void ProcessBlock(size_t n, Float * x, Float * y, Float * output) {
            for(size_t i = 0; i < n; i+=bump)
            {
                SIMD o1,o2,o;
                o1.load_a(x+i);
                o2.load_a(y+i);
                o = Tick(o1,o2);
                o.store_a(output+i);
            }
        }
    };
    
    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUOscillatorProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {        
        CPUOscillatorProcessor() : CPUSoundProcessor<DSPType,Float,bump>() 
        {
            
        }
        ObjectType getType() const {
            return CPU_OSCILLATOR_PROCESSOR;
        }
    };

    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUFilterProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {
        CPUFilterProcessor() : CPUSoundProcessor<DSPType,Float,bump>() {
    
        }
        ObjectType getType() const {
            return CPU_FILTER_PROCESSOR;
        }                
    };


    template<typename DSPType=Vector4f, typename Float=float, int bump=4>
    struct CPUAmplifierProcessor : public CPUSoundProcessor<DSPType,Float,bump>
    {        
        CPUAmplifierProcessor() : CPUSoundProcessor<DSPType,Float,bump>()
        {
        
        }
        ObjectType getType() const {
            return CPU_AMPLIFIER_PROCESSOR;
        }                    
    };


};
