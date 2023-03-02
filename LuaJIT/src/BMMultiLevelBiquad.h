//
//  BMMultiLevelBiquad.h
//  VelocityFilter
//
//  Created by Hans on 14/3/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#ifndef BMMultiLevelBiquad_h
#define BMMultiLevelBiquad_h

#include <stdio.h>
#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BMMultiLevelBiquad {
    // dynamic memory
    vDSP_biquadm_Setup multiChannelFilterSetup;
    vDSP_biquad_Setup singleChannelFilterSetup;
    float* monoDelays;
    double* coefficients;
    float* coefficients_float;
    
    // static memory
    float gain;
    size_t numLevels;
    size_t numChannels;
    double sampleRate;
    bool needsUpdate, useRealTimeUpdate, useBiquadm;
    
} BMMultiLevelBiquad;



// process a stereo buffer of samples
void BMMultiLevelBiquad_processBufferStereo(BMMultiLevelBiquad* bqf, const float* inL, const float* inR, float* outL, float* outR, size_t numSamples);

// process a mono buffer of samples
void BMMultiLevelBiquad_processBufferMono(BMMultiLevelBiquad* bqf, const float* input, float* output, size_t numSamples);

/*
 * init must be called once before using the filter.  To change the number of
 * levels in the fitler, call destroy first, then call this function with
 * the new number of levels
 *
 * monoRealTimeUpdate: If you are updating coefficients of a MONO filter in
 *                     realtime, set this to true. Processing of audio is
 *                     slightly slower, but updates can happen in realtime.
 *                     This setting has no effect on stereo filters.
 *                     This setting has no effect if the OS does not support
 *                     realtime updates of vDSP_biquadm filter coefficients.
 *
 */
void BMMultiLevelBiquad_init(BMMultiLevelBiquad* bqf, size_t numLevels, float sampleRate, bool stereo, bool monoRealTimeUpdate);

// free up memory objects
void BMMultiLevelBiquad_destroy(BMMultiLevelBiquad* bqf);


// set a bell-shape filter at on the specified level in both channels
// and update filter settings
void BMMultiLevelBiquad_setBell(BMMultiLevelBiquad* bqf, float fc, float bandwidth, float gain_db, size_t level);


// set a high shelf filter at on the specified level in both
// channels and update filter settings
void BMMultiLevelBiquad_setHighShelf(BMMultiLevelBiquad* bqf, float fc, float gain_db, size_t level);


void BMMultiLevelBiquad_setLowPass12db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level);


//void BMMultiLevelBiquad_setHighPass12db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level);


void BMMultiLevelBiquad_setHighPass6db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level);


// Calling this sets the filter coefficients at 'level' to bypass.
// Note that the filter still processes through the bypassed section
// but the output is the same as the input.
void BMMultiLevelBiquad_setBypass(BMMultiLevelBiquad* bqf, size_t level);


// set gain in db
void BMMultiLevelBiquad_setGain(BMMultiLevelBiquad* bqf, float gain_db);


/*
 * frequency: an an array specifying frequencies at which we want to evaluate
 * the transfer function magnitude of the filter.
 *
 * magnitude: an array for storing the result
 *
 * fs: the sampling rate
 *
 * length: the number of elements in frequency and magnitude
 *
 */
void BMMultiLevelBiquad_tfMagVector(BMMultiLevelBiquad* bqf, const float *frequency, float *magnitude, float fs, size_t length);


// returns the total group delay of all levels of the filter at the
// specified frequency.
double BMMultiLevelBiquad_groupDelay(BMMultiLevelBiquad* bqf, double freq, double sampleRate);

#ifdef __cplusplus
}
#endif

#endif /* BMMultiLevelBiquad_h */
