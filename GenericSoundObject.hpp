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

    virtual void InplaceProcess(size_t n, DSPType * buffer) {
        ProcessBlock(n,buffer,buffer);
    }

    virtual DSPType Tick(DSPType I=1, DSPType A=1, DSPType X=0, DSPType Y=0)
    {
        assert(1==0);
    }
    virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) 
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            inputs[i] = Tick(outputs[i]);
    }
};

// VC Vector and AVec
template<typename DSPType, typename SIMDType, int bump>
struct VCSoundProcessor : public GSSoundProcessor<DSPType>
{
    VCSoundProcessor() : GSSoundProcessor<DSPType>
    {
    
    }

    virtual void InplaceProcess(size_t n, DSPType * buffer) {
        ProcessBlock(n,buffer,buffer);
    }

    virtual SIMDType Tick(SIMDType I=1, DSPType A=1, DSPType X=0, DSPType Y=0)
    {
        assert(1==0);
    }
    virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) 
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i+=bump)
        {
            SIMDType ri;
            ri.load_a(inputs+i);
            ri = Tick(ri);
            ri.store_a(outputs+i);            
        }
    }
}

// stereo = matrix
template<typename DSPType>
struct VectorSoundProcessor : public GSSoundProcessor<DSPType>
{

};