//
//  BMMultiLevelBiquad.c
//  VelocityFilter
//
//  Created by Hans on 14/3/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "BMMultiLevelBiquad.h"
#include "Constants.h"
#include "BMComplexMath.h"
#include "BMGetOSVersion.h"
#include <stdlib.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * function declarations for use within this file
 */

// evaluate the combined transfer function of all levels of the filter
// at the frequency specified by z
DSPDoubleComplex BMMultiLevelBiquad_tfEval(BMMultiLevelBiquad* bqf, DSPDoubleComplex z);


// Destroy the current filter setup and create a new one.
// This function must be called only from the audio processing thread.
void BMMultiLevelBiquad_recreate(BMMultiLevelBiquad* bqf);


// this is a thread-safe way to update filter coefficients
void BMMultiLevelBiquad_enqueueUpdate(BMMultiLevelBiquad* bqf);


// this function updates the filter immediately and is safe to call
// only from the audio thread. It changes the filter coefficients in realtime
// if possible; otherwise it calls _recreate.
void BMMultiLevelBiquad_updateNow(BMMultiLevelBiquad* bqf);


// returns true if the operating system supports vDSP_biquadm_SetCoefficentsSingle()
bool BMMultiLevelBiquad_OSSupportsRealtimeUpdate();

/* end internal function declarations */



void BMMultiLevelBiquad_processBufferStereo(BMMultiLevelBiquad* bqf, const float* inL, const float* inR, float* outL, float* outR, size_t numSamples){
    
    // this function is only for two channel filtering
    assert(bqf->numChannels == 2);
    
    // update filter coefficients if necessary
    if (bqf->needsUpdate) BMMultiLevelBiquad_updateNow(bqf);
    
    // link the two input buffers into a single multidimensional array
    const float* twoChannelInput [2];
    twoChannelInput[0] = inL;
    twoChannelInput[1] = inR;
    
    // link the two input buffers into a single multidimensional array
    float* twoChannelOutput [2];
    twoChannelOutput[0] = outL;
    twoChannelOutput[1] = outR;
    
    // apply a multilevel biquad filter to both channels
    vDSP_biquadm(bqf->multiChannelFilterSetup, (const float* _Nonnull * _Nonnull)twoChannelInput, 1, twoChannelOutput, 1, numSamples);
    
    
    // apply gain to both channels
    vDSP_vsmul(outL, 1, &bqf->gain, outL, 1, numSamples);
    vDSP_vsmul(outR, 1, &bqf->gain, outR, 1, numSamples);
}




void BMMultiLevelBiquad_processBufferMono(BMMultiLevelBiquad* bqf, const float* input, float* output, size_t numSamples){
    
    // this function is only for single channel filtering
    assert(bqf->numChannels == 1);
    
    // update filter coefficients if necessary
    if (bqf->needsUpdate) BMMultiLevelBiquad_updateNow(bqf);
    
    // if using the multiChannel filter for single channel processing
    if(bqf->useBiquadm){
        // biquadm requires arrays of pointers as input and output
        const float* inputP [1];
        float* outputP [1];
        inputP[0]  = input;
        outputP[0] = output;
        
        // apply a multiChannel biquad filter
        vDSP_biquadm(bqf->multiChannelFilterSetup, (const float* _Nonnull * _Nonnull)inputP, 1, outputP, 1, numSamples);
        
        
        // if using the single channel filter
    } else {
        vDSP_biquad(bqf->singleChannelFilterSetup, bqf->monoDelays, input, 1, output, 1, numSamples);
    }
    
    // apply gain
    vDSP_vsmul(output, 1, &bqf->gain, output, 1, numSamples);
}


// Find out if the OS supports vDSP_biquadm updates in realtime
bool BMMultiLevelBiquad_OSSupportsRealtimeUpdate(){
    
    bool OSSupportsRealtimeUpdate = false;
    
    // iOS >= 9.0 (build# 13) supports realtime updates
    if (BM_isiOS() && BM_getOSMajorBuildNumber() >= 13)OSSupportsRealtimeUpdate = true;
    
    // Mac OS X >= 10.0 (build#14) supports them too
    if (BM_isMacOS() && BM_getOSMajorBuildNumber() >=14)OSSupportsRealtimeUpdate = true;
    
    return OSSupportsRealtimeUpdate;
}

void BMMultiLevelBiquad_init(BMMultiLevelBiquad* bqf,
                             size_t numLevels,
                             float sampleRate,
                             bool isStereo,
                             bool monoRealTimeUpdate){
    
    bqf->needsUpdate = false;
    bqf->sampleRate = sampleRate;
    bqf->numLevels = numLevels;
    bqf->numChannels = isStereo ? 2 : 1;
    
    // Should we use the multichannel biquad?
    // Even for mono signals, we have to use biquadm if we need realtime
    // update of filter coefficients.
    bqf->useBiquadm = false;
    if (isStereo || monoRealTimeUpdate) bqf->useBiquadm = true;
    
    // We will update in realtime if the OS supports it and we are using
    // vDSP_biquadm
    bqf->useRealTimeUpdate = false;
    if(BMMultiLevelBiquad_OSSupportsRealtimeUpdate() && bqf->useBiquadm)
        bqf->useRealTimeUpdate = true;
    
    
    // Allocate memory for 5 coefficients per filter,
    // 2 filters per level (left and right channels)
    free(bqf->coefficients);
    bqf->coefficients = malloc(numLevels*5*bqf->numChannels*sizeof(double));
    
    // repeat the allocation for floating point coefficients. We need
    // both double and float to support realtime updates
    free(bqf->coefficients_float);
    bqf->coefficients_float = malloc(numLevels*5*bqf->numChannels*sizeof(float));
    
    // start with all levels on bypass
    for (size_t i=0; i<numLevels; i++) {
        BMMultiLevelBiquad_setBypass(bqf, i);
    }
    
    // set 0db of gain
    BMMultiLevelBiquad_setGain(bqf,0.0);
    
    // setup filter struct
    BMMultiLevelBiquad_recreate(bqf);
}

void BMMultiLevelBiquad_setGain(BMMultiLevelBiquad* bqf, float gain_db){
    bqf->gain = BM_DB_TO_GAIN(gain_db);
}

void BMMultiLevelBiquad_queueUpdate(BMMultiLevelBiquad* bqf){
    bqf->needsUpdate = true;
}

inline void BMMultiLevelBiquad_updateNow(BMMultiLevelBiquad* bqf){
    
    // using realtime updates
    if(bqf->useRealTimeUpdate){
        // convert the coefficients to floating point
        for(size_t i=0; i<bqf->numLevels*bqf->numChannels*5; i++)
            bqf->coefficients_float[i] = bqf->coefficients[i];
        
        // update the coefficients
        vDSP_biquadm_SetCoefficientsSingle(bqf->multiChannelFilterSetup, bqf->coefficients_float, 0, 0, bqf->numLevels, bqf->numChannels);
        
        // not using realtime updates
    } else {
        BMMultiLevelBiquad_recreate(bqf);
    }
    
    bqf->needsUpdate = false;
}



inline void BMMultiLevelBiquad_recreate(BMMultiLevelBiquad* bqf){
    // using multichannel vDSP_biquadm
    if (bqf->multiChannelFilterSetup)
        vDSP_biquadm_DestroySetup(bqf->multiChannelFilterSetup);
    
    // using single channel vDSP_biquad
    if(bqf->singleChannelFilterSetup)
        vDSP_biquad_DestroySetup(bqf->singleChannelFilterSetup);
    
    if(bqf->useBiquadm)
        bqf->multiChannelFilterSetup =
        vDSP_biquadm_CreateSetup(bqf->coefficients, bqf->numLevels, bqf->numChannels);
    else
        bqf->singleChannelFilterSetup =
        vDSP_biquad_CreateSetup(bqf->coefficients, bqf->numLevels);
}

void BMMultiLevelBiquad_destroy(BMMultiLevelBiquad* bqf){
    if(bqf->coefficients) free(bqf->coefficients);
    if(bqf->coefficients_float) free(bqf->coefficients_float);
    vDSP_biquadm_DestroySetup(bqf->multiChannelFilterSetup);
}


void BMMultiLevelBiquad_setBypass(BMMultiLevelBiquad* bqf, size_t level){
    assert(level < bqf->numLevels);
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b0 + 2;
        double* a1 = b0 + 3;
        double* a2 = b0 + 4;
        
        *b0 = 1.0;
        *b1 = *b2 = *a1 = *a2 = 0.0;
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}


// based on formula in 2.3.10 of Digital Filters for Everyone by Rusty Allred
void BMMultiLevelBiquad_setHighShelf(BMMultiLevelBiquad* bqf, float fc, float gain_db, size_t level){
    assert(level < bqf->numLevels);
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b0 + 2;
        double* a1 = b0 + 3;
        double* a2 = b0 + 4;
        
        float gainV = BM_DB_TO_GAIN(gain_db);
        
        // if gain is close to 1.0, bypass the filter
        if (fabsf(gain_db) < 0.01){
            *b0 = 1.0;
            *b1 = *b2 = *a1 = *a2 = 0.0;
        }
        
        // if the gain is nontrivial
        else {
            double gamma = tanf(M_PI * fc / bqf->sampleRate);
            double gamma_2 = gamma*gamma;
            double sqrt_gain = sqrtf(gainV);
            
            // conditionally set G
            double G;
            if (gainV > 2.0f) G = gainV * M_SQRT2 * 0.5f;
            else {
                if (gainV >= 0.5f) G = sqrt_gain;
                else G = gainV * M_SQRT2;
            }
            double G_2 = G*G;
            
            // compute reuseable variables
            double g_d = powf((G_2 - 1.0f)/(gainV*gainV - G_2), 0.25f);
            double g_d_2 = g_d*g_d;
            double g_n = g_d * sqrt_gain;
            double g_n_2 = g_n * g_n;
            double sqrt_2_g_d_gamma = M_SQRT2 * g_d * gamma;
            double sqrt_2_g_n_gamma = M_SQRT2 * g_n * gamma;
            double gamma_2_plus_g_d_2 = gamma_2 + g_d_2;
            double gamma_2_plus_g_n_2 = gamma_2 + g_n_2;
            
            double one_over_denominator = 1.0f / (gamma_2_plus_g_d_2 + sqrt_2_g_d_gamma);
            
            *b0 = (gamma_2_plus_g_n_2 + sqrt_2_g_n_gamma) * one_over_denominator;
            *b1 = 2.0f * (gamma_2 - g_n_2) * one_over_denominator;
            *b2 = (gamma_2_plus_g_n_2 - sqrt_2_g_n_gamma) * one_over_denominator;
            
            *a1 = 2.0f * (gamma_2 - g_d_2) * one_over_denominator;
            *a2 = (gamma_2_plus_g_d_2 - sqrt_2_g_d_gamma)*one_over_denominator;
        }
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}



// based on formulae in 2.3.8 in Digital Filters are for Everyone,
// 2nd ed. by Rusty Allred
void BMMultiLevelBiquad_setBell(BMMultiLevelBiquad* bqf, float fc, float bandwidth,float gain_db, size_t level){
    assert(level < bqf->numLevels);
    
    float gainV = BM_DB_TO_GAIN(gain_db);
    
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b0 + 2;
        double* a1 = b0 + 3;
        double* a2 = b0 + 4;
        
        // if gain is close to 1.0, bypass the filter
        if (fabsf(gain_db) < 0.01){
            *b0 = 1.0;
            *b1 = *b2 = *a1 = *a2 = 0.0;
        }
        
        // if the gain is nontrivial
        else {
            double alpha =  tan( (M_PI * bandwidth)   / bqf->sampleRate);
            double beta  = -cos( (M_2_PI * fc) / bqf->sampleRate);
            double oneOverD;
            
            if (gainV < 1.0) {
                oneOverD = 1.0 / (alpha + gainV);
                // feed-forward coefficients
                *b0 = (gainV + alpha*gainV) * oneOverD;
                *b1 = 2.0 * beta * gainV * oneOverD;
                *b2 = (gainV - alpha*gainV) * oneOverD;
                
                // recursive coefficients
                *a1 = 2.0 * beta * gainV * oneOverD;
                *a2 = (gainV - alpha) * oneOverD;
            } else { // gain >= 1
                oneOverD = 1.0 / (alpha + 1.0);
                // feed-forward coefficients
                *b0 = (1.0 + alpha*gainV) * oneOverD;
                *b1 = 2.0 * beta * oneOverD;
                *b2 = (1.0 - alpha*gainV) * oneOverD;
                
                // recursive coefficients
                *a1 = 2.0 * beta * oneOverD;
                *a2 = (1.0 - alpha) * oneOverD;
            }
        }
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}

void BMMultiLevelBiquad_setLowPass12db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level){
    assert(level < bqf->numLevels);
    
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b1 + 1;
        double* a1 = b2 + 1;
        double* a2 = a1 + 1;
        
        // if fc is greater than 99% of the nyquyst frequency, bypass the filter
        if (fc > 0.99*0.5*sampleRate){
            *b0 = 1.0;
            *b1 = *b2 = *a1 = *a2 = 0.0;
        }
        
        // else, the cutoff frequency is non-trivial
        else {
            double gamma = tan(M_PI * fc / sampleRate);
            double gamma_2 = gamma * gamma;
            double gamma_x_sqrt_2 = gamma * M_SQRT2;
            double one_over_denominator = 1.0 / (gamma_2 + gamma_x_sqrt_2 + 1.0);
            
            *b0 = gamma_2 * one_over_denominator;
            *b1 = 2.0 * *b0;
            *b2 = *b0;
            
            *a1 = 2.0 * (gamma_2 - 1.0) * one_over_denominator;
            *a2 = (gamma_2 - gamma_x_sqrt_2 + 1.0) * one_over_denominator;
        }
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}

void BMMultiLevelBiquad_setHighPass12db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level){
    assert(level < bqf->numLevels);
    
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b1 + 1;
        double* a1 = b2 + 1;
        double* a2 = a1 + 1;
        
        // if the cutoff frequency is less than 1 hz, bypass the filter
        if (fc < 1.0){
            *b0 = 1.0;
            *b1 = *b2 = *a1 = *a2 = 0.0;
        }
        
        // else, for useful settings of fc, setup the filter
        else {
            double gamma = tan(M_PI * fc / bqf->sampleRate);
            double gamma_2 = gamma * gamma;
            double gamma_x_sqrt_2 = gamma * M_SQRT2;
            double one_over_denominator = 1.0 / (gamma_2 + gamma_x_sqrt_2 + 1.0);
            
            *b0 = 1.0 * one_over_denominator;
            *b1 = -2.0 * one_over_denominator;
            *b2 = *b0;
            
            *a1 = 2.0 * (gamma_2 - 1.0) * one_over_denominator;
            *a2 = (gamma_2 - gamma_x_sqrt_2 + 1.0) * one_over_denominator;
        }
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}


void BMMultiLevelBiquad_setHighPass6db(BMMultiLevelBiquad* bqf, double fc, double sampleRate, size_t level){
    assert(level < bqf->numLevels);
    
    
    // for left and right channels, set coefficients
    for(size_t i=0; i < bqf->numChannels; i++){
        
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + i*5;
        double* b1 = b0 + 1;
        double* b2 = b1 + 1;
        double* a1 = b2 + 1;
        double* a2 = a1 + 1;
        
        // if the cutoff frequency is less than 10 hz, bypass the filter
        if (fc < 10.0){
            *b0 = 1.0;
            *b1 = *b2 = *a1 = *a2 = 0.0;
        }
        
        // else, for useful settings of fc, setup the filter
        else {
            double gamma = tan(M_PI * fc / bqf->sampleRate);
            double one_over_denominator = 1.0 / (gamma + 1.0);
            
            *b0 = 1.0 * one_over_denominator;
            *b1 = -1.0 * one_over_denominator;
            *b2 = 0.0;
            
            *a1 = (gamma - 1.0) * one_over_denominator;
            *a2 = 0.0;
        }
    }
    
    BMMultiLevelBiquad_queueUpdate(bqf);
}


// evaluate the transfer function of the filter at all levels for the
// frequency specified by the complex number z
inline DSPDoubleComplex BMMultiLevelBiquad_tfEval(BMMultiLevelBiquad* bqf, DSPDoubleComplex z){
    
    DSPDoubleComplex z2 = DSPDoubleComplex_cmul(z, z);
    
    DSPDoubleComplex out = DSPDoubleComplex_init(bqf->gain, 0.0);
    
    
    for (size_t level = 0; level < bqf->numLevels; level++) {
        
        // both channels are the same so we just check the left one
        size_t channel = 0;
        
        double* b0 = bqf->coefficients + level*bqf->numChannels*5 + channel*5;
        double* b1 = b0+1;
        double* b2 = b0+2;
        double* a1 = b0+3;
        double* a2 = b0+4;
        
        
        DSPDoubleComplex numerator =
        DSPDoubleComplex_add3(DSPDoubleComplex_smul(*b0, z2),
                              DSPDoubleComplex_smul(*b1, z),
                              DSPDoubleComplex_init(*b2, 0.0));
        
        DSPDoubleComplex denominator =
        DSPDoubleComplex_add3(z2,
                              DSPDoubleComplex_smul(*a1, z),
                              DSPDoubleComplex_init(*a2, 0.0));
        
        out = DSPDoubleComplex_cmul(out,
                                    DSPDoubleComplex_divide(numerator,
                                                            denominator));
    }
    
    return out;
}

/*
 * frequency: an an array specifying frequencies at which we want to evaluate
 * the transfer function magnitude of the filter
 *
 * magnitude: an array for storing the result
 *
 * fs: the sampling rate
 *
 * length: the number of elements in frequency and magnitude
 *
 */
void BMMultiLevelBiquad_tfMagVector(BMMultiLevelBiquad* bqf, const float *frequency, float *magnitude, float fs, size_t length){
    for (size_t i = 0; i<length; i++) {
        // convert from frequency into z (complex angular velocity)
        DSPDoubleComplex z = DSPDoubleComplex_z(frequency[i], fs);
        // evaluate the transfer function at z and take absolute value
        magnitude[i] = DSPDoubleComplex_abs(BMMultiLevelBiquad_tfEval(bqf,z));
    }
}


// returns the total group delay of all levels of the filter at the
// specified frequency.
//
// uses a cookbook formula for group delay of biquad filters, based on
// the fft derivative method.
double BMMultiLevelBiquad_groupDelay(BMMultiLevelBiquad* bqf, double freq, double sampleRate){
    double delay = 0.0;
    
    for (size_t level=0; level<bqf->numLevels; level++) {
        
        double b0 = bqf->coefficients[5*level];
        double b1 = bqf->coefficients[5*level + 1];
        double b2 = bqf->coefficients[5*level + 2];
        double a1 = bqf->coefficients[5*level + 3];
        double a2 = bqf->coefficients[5*level + 4];
        
        // normalize the feed forward coefficients so that b0=1
        // see: see: http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
        b1 /= b0;
        b2 /= b0;
        b0 = 1.0;
        
        // radian normalised frequency
        double w = freq / (0.5*sampleRate);
        
        // calculate the group delay of the normalized filter using a cookbook formula
        // http://music-dsp.music.columbia.narkive.com/9F6BIvHy/group-delay
        // or
        // http://music.columbia.edu/pipermail/music-dsp/1998-April/053307.html
        //
        //    T(w) =
        //
        //      b1^2 + 2*b2^2 + b1*(1 + 3*b2)*cos(w) + 2*b2*cos(2*w)
        //    --------------------------------------------------------
        //     1 + b1^2 + b2^2 + 2*b1*(1 + b2)*cos(w) + 2*b2*cos(2*w)
        //
        //
        //        a1^2 + 2*a2^2 + a1*(1 + 3*a2)*cos(w) + 2*a2*cos(2*w)
        //    - --------------------------------------------------------
        //        1 + a1^2 + a2^2 + 2*a1*(1 + a2)*cos(w) + 2*a2*cos(2*w)
        //
        //
        //    w is normalized radian frequency and T(w) is measured in sample units.
        
        
        //      b1^2 + 2*b2^2 + b1*(1 + 3*b2)*cos(w) + 2*b2*cos(2*w)
        double num1 = b1*b1 + 2.0*b2*b2 + b1*(1.0 + 3.0*b2)*cos(w) + 2.0*b2*cos(2.0*w);
        //     1 + b1^2 + b2^2 + 2*b1*(1 + b2)*cos(w) + 2*b2*cos(2*w)
        double den1 = 1.0 + b1*b1 + b2*b2 + 2.0*b1*(1.0 + b2)*cos(w) + 2.0*b2*cos(2.0*w);
        double frac1 = num1/den1;
        
        
        //        a1^2 + 2*a2^2 + a1*(1 + 3*a2)*cos(w) + 2*a2*cos(2*w)
        double num2 = a1*a1 + 2.0*a2*a2 + a1*(1.0 + 3.0*a2)*cos(w) + 2.0*a2*cos(2.0*w);
        //        1 + a1^2 + a2^2 + 2*a1*(1 + a2)*cos(w) + 2*a2*cos(2*w)
        double den2 = 1.0 + a1*a1 + a2*a2 + 2.0*a1*(1.0 + a2)*cos(w) + 2.0*a2*cos(2.0*w);
        double frac2 = num2/den2;
        
        // add the delay of the current level to the total delay
        delay += frac1 - frac2;
    }
    
    return delay;
}

#ifdef __cplusplus
}
#endif