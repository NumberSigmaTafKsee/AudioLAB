#pragma once 
#include <vector>
#include <iostream>
#include "WaveTable.hpp"

struct FourierWave : public OscillatorProcessor
{
    WaveTableOsc wave;
    DspFloatType        *F;
    DspFloatType        *P;
    DspFloatType        *A;
    size_t        n;
    DspFloatType         sr;

    FourierWave(DspFloatType sample_rate, int wave_size=1024) : OscillatorProcessor() {
        n = sample_rate/2;
        F = new DspFloatType[n];
        P = new DspFloatType[n];
        A = new DspFloatType[n];
        sr= sample_rate;

        std::vector<DspFloatType> sine(wave_size);
        MakeSineTable(sine,wave_size,1.0f,sr);
        wave.addWaveTable(wave_size,sine,sr/2.0f);
    }
    ~FourierWave() {
        if(F) delete [] F;
        if(A) delete [] A;
    }

    void SetSaw(DspFloatType f) {
        size_t max_harmonics = (sr/2)/f;
        memset(F,0,n*sizeof(DspFloatType));
        memset(A,0,n*sizeof(DspFloatType));   
        #pragma omp simd     
        for(size_t i = 1; i < max_harmonics; i++)
        {
            F[i-1] = i*f;
            A[i-1] = 1.0f/(DspFloatType)i;
        }
    }
    void SetSquare(DspFloatType f) {
        size_t max_harmonics = (sr/2)/f;
        memset(F,0,n*sizeof(DspFloatType));
        memset(A,0,n*sizeof(DspFloatType));        
        #pragma omp simd
        for(size_t i = 1; i < max_harmonics; i++)
        {                
            F[i-1] = i*f;
            if(i % 2 == 0) A[i-1] = 0;
            else A[i-1] = 1.0f/(DspFloatType)i;
        }
    }
    void SetTriangle(DspFloatType f) {
        size_t max_harmonics = (sr/2)/f;
        memset(F,0,n*sizeof(DspFloatType));
        memset(A,0,n*sizeof(DspFloatType));        
        #pragma omp simd
        for(size_t i = 1; i < max_harmonics; i++)
        {                
            F[i-1] = i*f;
            if(i % 2 == 0) A[i-1] = 0;
            else A[i-1] = std::pow(-1.0f,(i-1)/2.0f)/(DspFloatType)(i*i);
        }
    }
    void SetSine(DspFloatType f)
    {
        memset(F,0,n*sizeof(DspFloatType));
        memset(A,0,n*sizeof(DspFloatType));        
        F[0] = f;
        A[0] = 1.0f;
    }
    DspFloatType Tick(DspFloatType Index=1, DspFloatType G=1, DspFloatType FM=1, DspFloatType PM=1) {
        DspFloatType r = 0.0f;
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            if(F[i] == 0.0f) break;       
            wave.SetPhase(P[i]);                         
            P[i] += F[i]/sr;
            if(P[i] >= 1.0f) P[i] -= 1.0f;                
            r += Index*wave.Tick();
        }
        return G*r;
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
    }
    void ProcessInplace(size_t n, DspFloatType * in) {
        ProcessBlock(n,in,in);
    }
};

