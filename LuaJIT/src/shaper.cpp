#pragma once
#include <cmath>



float linearScale( float in, float min, float max )
{
	float ret;
	if ( min == 0.0f && max == 0.0f )
	{
		ret = 0.0f;
	}
	else if ( min > max )
	{
		ret = min - ( in * ( min - max ) );
	}
	else
	{
		ret = min + ( in * ( max - min ) );
	}
	return ret;
}

float linearDescale( float in, float min, float max )
{
	float ret;
	if ( min == 0.0f && max == 0.0f )
	{
		ret = 0.0f;
	}
	else if ( min > max )
	{
		ret = ( min - in ) / ( min - max );
	}
	else
	{
		ret = ( in - min ) / ( max - min );
	}
	return ret;
}

float expoScale( float in, float min, float max )
{
	// negative log makes no sense...
	if ( min < 0.0f || max < 0.0f )
	{
		return 0.0f;
	}

	// not handling min > max (inverse) case yet

	// figure out how many "octaves" (doublings) it takes to get from min to
	// max
	// we only have log & log10 so we have to do change of base
	// note this uses + instead of * so we can handle min == 0
	float octaves = log( max - min + 1 ) / log( 2.0f );
	return ( min - 1 ) + pow( 2.0f, in * octaves );
}

float expoDescale( float in, float min, float max )
{
	// see above
	if ( min < 0.0f || max < 0.0f )
	{
		return 0.0f;
	}

	// again, not handling min > max (inverse) case yet
	
	// note this was derived by simply inverting the previous function
	float log2 = log( 2.0f );
	return ( log( in - min + 1 ) / log2 ) / ( log( max - min + 1 ) / log2 );
}

float floorScale( float in, float min, float max )
{
	if ( min > max )
	{
		return ceil( linearScale( in, min, max ) );
	}
	else
	{
		return floor( linearScale( in, min, max ) );
	}
}

float expoShape( float in, float amount )
{
	if ( in == 0.0f )
		return in;

	float flip = in < 0.0f ? -1.0f : 1.0f;

	return pow( in * flip, amount ) * flip;
}

float softClipShape( float in, float amount )
{
	return in / ( 1 + fabs( in ) );
}

float sineShape( float in, float amount )
{
	return in * cos( in * amount );
}

float chebyshevShape( float in, float amount )
{
	return chebyshevRec( in, (int)amount );
}

float chebyshevRec( float in, int depth )
{
	if ( depth == 0 )
	{
		return 1.0f;
	}

	// lastval represents C(k-1)
	float lastVal = 1.0f;
	float out = in;
	float temp;

	// depth=1 is the base case
	for( int i = 1; i < depth; i++ )
	{
		temp = out;
		out = ( 2.0f * in * out ) - lastVal;
		lastVal = temp;
	}
	return out;
}