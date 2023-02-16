//
//  BMVelocityFilter.c
//  VelocityFilter
//
//  Created by Hans on 16/3/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "BMVelocityFilter.h"
#include "Constants.h"


#ifdef __cplusplus
extern "C" {
#endif

void BMVelocityFilter_init(BMVelocityFilter* f, float sampleRate, bool stereo){
    BMMultiLevelBiquad_init(&f->bqf, 1, sampleRate, stereo, true);
    
    f->centreVelocity = 100.0;
    
    // initialise the filter to allpass
    float fcForMIDINote = 1000.0;
    BMVelocityFilter_newNote(f, 100, 0, &fcForMIDINote);
}



void BMVelocityFilter_processBufferStereo(BMVelocityFilter* f, const float* inL, const float* inR, float* outR, float* outL, size_t numSamples){
    BMMultiLevelBiquad_processBufferStereo(&f->bqf, inL, inR, outL, outR, numSamples);
}



void BMVelocityFilter_processBufferMono(BMVelocityFilter* f, const float* input, float* output, size_t numSamples){
    BMMultiLevelBiquad_processBufferMono(&f->bqf, input, output, numSamples);
}



void BMVelocityFilter_newNote(BMVelocityFilter* f, float velocity, size_t MIDINoteNumber, float* fcForMIDINote){
    assert(velocity >= 0.0 && velocity <= 127.0);
    
    
    /*
     *  three cases: velovity > centre
     *               velocity < centre
     *               velocity == centre
     */
    // velocity > centre => boost gain
    float gain;
    if (velocity > f->centreVelocity)
        gain = f->maxGainDb*(velocity-f->centreVelocity)/(127.0 - f->centreVelocity);
    
    // velocity < centre => cut gain
    else if (velocity < f->centreVelocity)
        gain = f->minGainDb*(f->centreVelocity - velocity)/(f->centreVelocity);
    
    // velocity == centre => unity gain
    else gain = 1.0;
    
    
    // convert to decibels
    float gain_db = BM_GAIN_TO_DB(gain);
    
    
    // set the high shelf filter gain and cutoff frequency,
    // looking up the fc from the table in fcForMIDINote
    BMMultiLevelBiquad_setHighShelf(&f->bqf,
                                    fcForMIDINote[MIDINoteNumber],
                                    gain_db,
                                    1);
}




void BMVelocityFilter_setVelocityGainRange(BMVelocityFilter* f, float minGainDb, float maxGainDb, float centreVelocity){
    assert(maxGainDb>=0.0 && minGainDb<=0.0);
    assert(centreVelocity < 127.0 && centreVelocity > 0.0);
    
    f->minGainDb = minGainDb;
    f->maxGainDb = maxGainDb;
}

void BMVelocityFilter_destroy(BMVelocityFilter* f){
    BMMultiLevelBiquad_destroy(&f->bqf);
}

#ifdef __cplusplus
}
#endif