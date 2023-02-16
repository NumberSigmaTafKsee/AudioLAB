#pragma once 

#include <cmath>
namespace DistortionFunctions {
// First two functions from here:
// http://www.mitpressjournals.org/doi/pdf/10.1162/comj.2009.33.2.85
float DistortionFunctions::distortionOne(float input)
{
    return 3 * input / 2 * ((1 - pow(input, 2) / 3));
}


float DistortionFunctions::distortionTwo(float input)
{
    // if input > 0 sign is minus
    const float sign = input > 0 ? 1.0f : -1.0f;
    
    return sign * (fabs(2 * input) - pow(input, 2));
}



// Functions from MusicDSP: http://musicdsp.org



/*
 Soft saturation: http://musicdsp.org/showone.php?id=42
 
 Type : waveshaper
 References : Posted by Bram de Jong
 
 Notes :
 This only works for positive values of x. param should be in the range 0..1
 */
float DistortionFunctions::softSaturation(float input, float param)
{
    float output;
    
    if (input < param)
    {
        output = input;
    }

    if (input > param)
    {
        output = param + (input - param) / (1 + pow(((input - param) / (1 - param)), 2));
    }


    if (input > 1)
    {
        output = (param + 1) / 2;
    }
    
    return output;
}


// Overdrive, default drive = 0
float DistortionFunctions::overdrive(float input, float drive = 0.0f)
{
     return (((sqrt(input) - input) * drive) + input) * (1 - ((drive / 4) * 0.5));
}


// Clipping, not very cool: http://musicdsp.org/showone.php?id=104
float DistortionFunctions::clipper(float input, float drive)
{
    float sign = (input < 0) ? -1.0f : 1.0f;
    
    return sign * pow(atan (pow (fabs(input), drive)), (1 / drive));
}    
}

