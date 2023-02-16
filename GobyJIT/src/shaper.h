#pragma once


/*
 * Scale takes input of [0.0-1.0] and scales to [min-max], descale does the
 * reverse. Inverse scaling (ie min > max) is supported.
 * Originally written as helper functions for interfacing with VST.
 */
float linearScale( float in, float min, float max );
float linearDescale( float in, float min, float max );

/*
 * Similar to the above functions, but exponential (base 2 for now)
 * Not to be confused with expoShape below which is a DSP function.
 */
float expoScale( float in, float min, float max );
float expoDescale( float in, float min, float max );

/*
 * Note we don't need a floorDescale (makes no sense since the output has to
 * be 0.0-1.0), we can just use linearDescale instead.
 */
float floorScale( float in, float min, float max );

/*
 * Audio shaping functions.
 */
float expoShape( float in, float amount );
float softClipShape( float in, float amount );
float sineShape( float in, float amount );
float chebyshevShape( float in, float amount );
float chebyshevRec( float in, int depth );
