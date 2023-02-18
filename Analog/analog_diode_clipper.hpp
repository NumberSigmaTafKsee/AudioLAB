#include <cmath>
#include "GenericSoundObject.hpp"

namespace Analog::Distortion::Diode
{
	// https://github.com/a-carson/DiodeClipper
	template<class DSP>
	class DiodeClipper : public GSSoundProcessor<DSP>
	{
	public:

		DiodeClipper() : GSSoundProcessor<DSP>()
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

		void setSampleRate(DSP sampleRate)
		{
			fs = sampleRate;
		}

		void setCircuitParameters(DSP resistance, DSP capacitance)
		{
			R = resistance;
			C = capacitance;
		}

		void setDiodeParameters(DSP saturationCurrent, DSP thermalVoltage, DSP idealityFactor)
		{
			Is = saturationCurrent;
			Vt = thermalVoltage;
			Ni = idealityFactor;
		}
		enum {
			PORT_R,
			PORT_C,
			PORT_IS,
			PORT_VT,
			PORT_NI,
		};
		void setPort(int port, DspFloatType v) {
			switch(port)
			{
				case PORT_R: R = v; break;
				case PORT_C: C = v; break;
				case PORT_IS: Is = v; break;
				case PORT_VT: Vt = v; break;
				case PORT_NI: Ni = v; break;
				default: printf("No port %d\n", port);
			}
		}
		void initialise()
		{
			cap = Ni * Vt * acosh((2.0f*fs*R*C + 1.0f) * Ni * Vt / (2.0f * Is * R));
		}

		// Non linear function
		DSP func(DSP Y, DSP p)
		{
			return Y * R * C + (Y + (2.0f * Is * R) * sinh(Y /(Ni * Vt)) - p) / (2.0f * fs);
		}

		// Derivative of Non-linear function
		DSP dfunc(DSP Y)
		{
			return R * C + (1 + (2.0f * Is * R / (Ni * Vt)) * cosh(Y / (Ni * Vt))) / (2.0f * fs);
		}

		// Fast approximation of function
		DSP fastFunc(DSP Y, DSP p)
		{
			DSP sinhy = sinh(Y / (Ni * Vt));
			return Y * R *C + (Y + (2.0f * Is * R) * sinhy - p) / (2.0f*fs);

		}

		// Fast aproximation of derivative
		DSP fastDfunc(DSP Y)
		{
			DSP coshy = cosh(Y / (Ni * Vt));
			return R * C + (1 + (2 * Is * R / (Ni * Vt))*coshy) / (2.0f * fs);
		}

		// Main Process
		DSP process(DSP in)
		{
			Undenormal denormal;            			
			const DSP p = x + in;		
			y = Ni * Vt * std::asinh(p / (2 * Is * R));
			DSP res = func(y, p);
			DSP J = dfunc(y);
			DSP step = res / J;
			DSP cond = std::fabs(step);
			iter = 0;

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
				DSP arg = y / (Ni * Vt);

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
				cond = std::fabs(step);
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

		void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            Undenormal denormal;            
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {            
				iter = 0;
				const DSP in = in[i];
				const DSP p = x + in;		
				y = Ni * Vt * std::asinh(p / (2 * Is * R));
				DSP res = func(y, p);
				DSP J = dfunc(y);
				DSP step = res / J;
				DSP cond = std::fabs(step);

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
					DSP arg = y / (Ni * Vt);

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
					cond = std::fabs(step);
				}

				// fail safe
				if (y != y)
				{
					y = 0;
				}

				// update state variable
				x =  4.0f * fs * y * R * C - x;
				out[i] = y;
            }
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DSP * out) {
            ProcessSIMD(n,out,out);
        }
	private:

		// Sample Rate
		DSP fs;

		// state variable
		DSP x;

		// output variable
		DSP y;

		// Circuit parameters
		DSP C;									// capactance
		DSP R;									// resistance
		DSP Is;								// saturation current
		DSP Vt;								// thermal voltage
		DSP Ni;								// ideality factor

		// Newton raphson parameters
		DSP cap;
		const DSP tol = 1e-7;						// tolerance
		const unsigned int maxIters = 25;  // maximum number of iterations
		unsigned int iter = 0;
	};
}


