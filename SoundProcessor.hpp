
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


typedef float DspFloatType;

#include "Undenormal.hpp"
#include "StdNoise.hpp"
#include "MusicFunctions.hpp"
#include "ClipFunctions.hpp"

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
    
    virtual DspFloatType Tick(DspFloatType I, DspFloatType A, DspFloatType X, DspFloatType Y) {
		assert(1==0);
	}
	virtual DspFloatType Tick(DspFloatType IL, DspFloatType IR, DspFloatType& oL, DspFloatType& oR, DspFloatType A, DspFloatType X, DspFloatType Y) {
		assert(1==0);
	}
	virtual void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		assert(1==0);
	}
	virtual void ProcessSIMD(size_t n, DspFloatType ** in, DspFloatType ** out) {
		assert(1==0);
	}
	virtual void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
		ProcessSIMD(n,in,out);
	}
	virtual void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
		ProcessSIMD(n,in,out);
	}
	virtual void ProcessInplace(size_t n, DspFloatType * in) {
		ProcessSIMD(n,in,in);
	}
	virtual void ProcessInplace(size_t n, DspFloatType ** in) {
		ProcessSIMD(n,in,in);
	}
};


