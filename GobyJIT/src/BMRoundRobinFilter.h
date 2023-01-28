//
//  BMRoundRobinFilter.h
//  BMAudioFilters
//
//  This struct simulates round robin sampling from a single sample by
//  applying a bank of bell filters with randomised gain.  Each time
//  the functino newNote() is called, the filters are randomised anew,
//  creating subtle tonal differences between notes.
//
//  Created by Hans on 23/3/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#ifndef BMRoundRobinFilter_h
#define BMRoundRobinFilter_h

#include "BMMultiLevelBiquad.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BMRoundRobinFilter{
    BMMultiLevelBiquad bqf;
    float minFc, maxFc, gainRange_db;
    size_t numBands;
} BMRoundRobinFilter;

/*
 * Must be called before using the filter.
 *
 * sampleRate - self explanitory
 * numBands - the number of bell-shaped filters to use for tone shaping
 * minFc - centre frequency of the lowest bell
 * maxFc - centre frequency of the highest bell
 * gainRange_db - the maximum excursion from zero gain for a single bell.
 *                for example, if this is 3db, then each bell will have gain
 *                set randomly between -3 and +3 db.
 * stereo - set true if you need stereo processing.
 *
 */
void BMRoundRobinFilter_init(BMRoundRobinFilter* f, float sampleRate, size_t numBands, float minFc, float maxFc, float gainRange_db, bool stereo);

/*
 * to use this function, you must call _init() with stereo=true
 */
void BMRoundRobinFilter_processBufferStereo(BMRoundRobinFilter* f, const float* inL, const float* inR, float* outR, float* outL, size_t numSamples);

/*
 * to use this function, you must call _init() with stereo=false
 */
void BMRoundRobinFilter_processBufferMono(BMRoundRobinFilter* f, const float* input, float* output, size_t numSamples);

/*
 * Call this function before the start of each new note to re-randomise the
 * filters.
 */
void BMRoundRobinFilter_newNote(BMRoundRobinFilter* f);

/*
 * free the memory used by the filter
 */
void BMRoundRobinFilter_destroy(BMRoundRobinFilter* f);

#ifdef __cplusplus
}
#endif

#endif /* BMRoundRobinFilter_h */
