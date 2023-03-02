//
//  BMRoundRobinFilter.c
//  BMAudioFilters
//
//  Created by Hans on 23/3/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "BMRoundRobinFilter.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void BMRoundRobinFilter_processBufferStereo(BMRoundRobinFilter* f, const float* inL, const float* inR, float* outR, float* outL, size_t numSamples){
    BMMultiLevelBiquad_processBufferStereo(&f->bqf, inL, inR, outL, outR, numSamples);
}

void BMRoundRobinFilter_processBufferMono(BMRoundRobinFilter* f, const float* input, float* output, size_t numSamples){
    BMMultiLevelBiquad_processBufferMono(&f->bqf, input, output, numSamples);
}

void BMRoundRobinFilter_init(BMRoundRobinFilter* f, float sampleRate, size_t numBands, float minFc, float maxFc, float gainRange_db, bool stereo){
    f->numBands = numBands;
    f->minFc = minFc;
    f->maxFc = maxFc;
    f->gainRange_db = gainRange_db;
    
    BMMultiLevelBiquad_init(&f->bqf, numBands, sampleRate, stereo, true);
    
    // initialize the round robin filter with a random setting
    BMRoundRobinFilter_newNote(f);
}


inline void BMRoundRobinFilter_newNote(BMRoundRobinFilter* f){
    float fc = f->minFc;
    
    // find the space between individual peaks in the filters
    // spacing them evenly in linear frequency
    float peakSpacing = (f->maxFc - f->minFc)/(float)(f->numBands - 1);
    
    for (size_t i=0; i<f->numBands; i++) {
        float gain_db = (((float)rand()/(float)RAND_MAX)-0.5)*2.0*f->gainRange_db;
        
        // note that Q = peakSpacing
        BMMultiLevelBiquad_setBell(&(f->bqf), fc, peakSpacing, gain_db, i);
        
        // move the fc over for the next peak
        fc += peakSpacing;
    }
}

void BMRoundRobinFilter_destroy(BMRoundRobinFilter* f){
    BMMultiLevelBiquad_destroy(&f->bqf);
}
    

#ifdef __cplusplus
}
#endif
