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
#include "ClipFunctions.hpp"

struct Random
{
    Random(){ srand(time(NULL)); }
    DspFloatType     frand() { return ((DspFloatType)rand()/(DspFloatType)RAND_MAX); }
    DspFloatType     rand() { return ((DspFloatType)rand()/(DspFloatType)RAND_MAX); }
    uint64_t    randint(int min, int max) { return round((max-min)*frand())+min; }
    bool        flip(DspFloatType prob) { return frand() < prob; }
    uint64_t    random(int mod) { return std::rand() % mod; }    
};



template<typename DSPType>
struct GSSoundProcessor
{    
    
    DSPType preGain = 1;
    DSPType postGain = 1;
    
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
        return I;
    }
    virtual void ProcessBlock(size_t n, DSPType * inputs, DSPType * outputs) 
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            inputs[i] = Tick(outputs[i]);
    }    
    virtual void ProcessInplace(size_t n, DSPType * in) {
        ProcessBlock(n,in,in);
    }
};

