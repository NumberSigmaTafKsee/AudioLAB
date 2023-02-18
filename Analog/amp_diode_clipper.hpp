#pragma once
#include <cmath>
#include "GenericSoundProcessor.hpp"

namespace FX::Distortion::Diode
{
	// https://github.com/a-carson/DiodeClipper
	Template<typename T>
	class DiodeClipper : public GSSoundProcessor<T>
	{
	public:

		DiodeClipper() : GSSoundProcessor<T>()
		{
			// initialise
			x = 0.0f;
			y = 0.0f;

			// default parameters
			R = 1000.0f;
			C = 33e-9f;
			fs = 48000.0f;
			Is = 2.52e-9f;
			Vt = 25.83e-3f;
			Ni = 1.752f;
		};

		void setSampleRate(T sampleRate)
		{
			fs = sampleRate;
		}

		void setCircuitParameters(T resistance, T capacitance)
		{
			R = resistance;
			C = capacitance;
		}

		void setDiodeParameters(T saturationCurrent, T thermalVoltage, T idealityFactor)
		{
			Is = saturationCurrent;
			Vt = thermalVoltage;
			Ni = idealityFactor;
		}

		void initialise()
		{
			cap = Ni * Vt * acosh((2.0f*fs*R*C + 1.0f) * Ni * Vt / (2.0f * Is * R));
		}

		// Non linear function
		T func(T Y, T p)
		{
			return Y * R * C + (Y + (2.0f * Is * R) * sinh(Y /(Ni * Vt)) - p) / (2.0f * fs);
		}

		// Derivative of Non-linear function
		T dfunc(T Y)
		{
			return R * C + (1 + (2.0f * Is * R / (Ni * Vt)) * cosh(Y / (Ni * Vt))) / (2.0f * fs);
		}

		// Fast approximation of function
		T fastFunc(T Y, T p)
		{
			T sinhy = sinh(Y / (Ni * Vt));
			return Y * R *C + (Y + (2.0f * Is * R) * sinhy - p) / (2.0f*fs);

		}

		// Fast aproximation of derivative
		T fastDfunc(T Y)
		{
			T coshy = cosh(Y / (Ni * Vt));
			return R * C + (1 + (2 * Is * R / (Ni * Vt))*coshy) / (2.0f * fs);
		}

		// Main Process
		T process(T in)
		{
			iter = 0;
			const T p = x + in;		
			y = Ni * Vt * asinh(p / (2 * Is * R));
			T res = func(y, p);
			T J = dfunc(y);
			T step = res / J;
			T cond = fabsf(step);

			while ((cond > tol) && (iter < maxIters))
			{
				// Cap step size if necessary
				if (step > cap)
				{
					step = cap;
				}
				if (step < -1.0f * cap)
				{
					step = -1.0f * cap;
				}

				// Newton step
				y -= step;
				
				// Check argument
				T arg = y / (Ni * Vt);

				// Compute residual and jacobian
				if (fabsf(arg) < 5)
				{
					res = fastFunc(y, p);
					J = fastDfunc(y);
				}
				else
				{
					res = func(y, p);
					J = dfunc(y);
				}

				// Calculate step
				step = res / J;

				iter++;
				cond = fabsf(step);
			}

			// fail safe
			if (y != y)
			{
				y = 0;
			}

			// update state variable
			x =  4.0f * fs * y * R * C - x;
			return y;
		}

		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
		{
			return A*process(I);
		}
		void ProcessSIMD(size_t n, T * in, T * out) {
			#pragma omp simd
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, T * in, T * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, T * buffer) {
			ProcessSIMD(n,buffer,buffer);
		}

	private:

		// Sample Rate
		T fs;

		// state variable
		T x;

		// output variable
		T y;

		// Circuit parameters
		T C;									// capactance
		T R;									// resistance
		T Is;								// saturation current
		T Vt;								// thermal voltage
		T Ni;								// ideality factor

		// Newton raphson parameters
		T cap;
		const T tol = 1e-7;						// tolerance
		const unsigned int maxIters = 25;  // maximum number of iterations
		unsigned int iter = 0;
	};
}


