#pragma once
#include <cmath>

namespace Distortion
{
	struct Shaper
	{
	
		static double linearScale( double in, double min, double max )
		{
			double ret;
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

		static double linearDescale( double in, double min, double max )
		{
			double ret;
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

		static double expoScale( double in, double min, double max )
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
			double octaves = log( max - min + 1 ) / log( 2.0f );
			return ( min - 1 ) + pow( 2.0f, in * octaves );
		}

		static double expoDescale( double in, double min, double max )
		{
			// see above
			if ( min < 0.0f || max < 0.0f )
			{
				return 0.0f;
			}

			// again, not handling min > max (inverse) case yet
			
			// note this was derived by simply inverting the previous function
			double log2 = log( 2.0f );
			return ( log( in - min + 1 ) / log2 ) / ( log( max - min + 1 ) / log2 );
		}

		static double floorScale( double in, double min, double max )
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

		static double expoShape( double in, double amount )
		{
			if ( in == 0.0f )
				return in;

			double flip = in < 0.0f ? -1.0f : 1.0f;

			return pow( in * flip, amount ) * flip;
		}

		static double softClipShape( double in, double amount )
		{
			return in / ( 1 + fabs( in ) );
		}

		static double sineShape( double in, double amount )
		{
			return in * cos( in * amount );
		}


		static double chebyshevRec( double in, int depth )
		{
			if ( depth == 0 )
			{
				return 1.0f;
			}

			// lastval represents C(k-1)
			double lastVal = 1.0f;
			double out = in;
			double temp;

			// depth=1 is the base case
			for( int i = 1; i < depth; i++ )
			{
				temp = out;
				out = ( 2.0f * in * out ) - lastVal;
				lastVal = temp;
			}
			return out;
		}


		static double chebyshevShape( double in, double amount )
		{
			return chebyshevRec( in, (int)amount );
		}
	};
}