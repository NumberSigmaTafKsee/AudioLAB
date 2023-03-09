#pragma once

#include "GenericSoundObject.hpp"

namespace Analog::MoogFilters::MoogCat
{
	template<typename DSP>
    struct MoogCat : public GSSoundProcessor<DSP>
    {
        
        DSP fc = 800.0f; //cutoff frequency
        DSP res = 0.0f; //resonance
        DSP Fs = 44100.0f; //sampling frequency
        DSP g = tan(M_PI * fc / Fs); //frequency warping factor
        DSP r_cat = 1.064f; //filter type value    

        //vector of states
        static constexpr size_t numStates = 4;
        std::vector<std::array<DSP, numStates>> state;

        MoogCat() : GSSoundProcessor<DSP>() {

        }
        void setCutoff(DSP f) {
            fc = f;        
        }
        void setResonance(DSP r) {
            res = r;
        }

        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,		
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setCutoff(v);
                break;
            case PORT_RESONANCE:
                setResonance(v);
                break;			
            }
        }
        DSP processSample(DSP input, size_t channel) noexcept {
			Undenormal denormal;
            auto& s = state[channel];
            const DSP g_2 = g*g;
            const DSP g_3 = g_2*g;
            const DSP g_4 = g_3*g;
            const DSP r_cat_2 = r_cat*r_cat;
            const DSP val1 = res*g_4*r_cat_2;
            const DSP val2 = g_3*r_cat;
            const DSP val3 = g*r_cat;
            const DSP val4 = g_2*r_cat_2;
            const DSP val5 = 2*val3 + 1;
            const DSP val6 = g_2 + val5;
            const DSP val7 = g*val6;
            const auto den = (4*val1 + g_4 + 4*val2 + 4*val4 + 2*g_2 + 4*val3 + 1);
            DSP out = (g_3*s[0] - g_2*(val5)*s[1] + val7*s[2] - (2*val2 + 4*val4 + g_2 + 4*val3 + 1)*s[3])/den;
            const DSP a = -(val1*4 + g_4 + 4*val2 + 4*val4 - 1)*s[0] + 2*g*(res*val4*4 + val6)*s[1] - 8*g_2*res*r_cat_2*s[2] + 8*g*res*r_cat_2*(2*val3+1)*s[3] + val7*2*input;
            const DSP b = - 2*val7*s[0] + (-val1*4 - g_4 + 4*val4 + 4*val3 + 1)*s[1] + 8*g_3*res*r_cat_2*s[2] - 8*g_2*res*r_cat_2*(2*val3+1)*s[3] - g_2*(val6)*2*input;
            const DSP c = 2*g_2*s[0] - g*(val5)*2*s[1] - (val1*4 + g_4 + 4*g_3*r_cat + 4*val4 - 1)*s[2] + 2*g*(res*val4*4 + val6)*s[3] + 2*g_3*input;
            const DSP d = -2*g_3*s[0] + g_2*(val5)*2*s[1] - val7*2*s[2] + (-val1*4 - g_4 + 4*val4 + 4*val3 + 1)*s[3] - 2*g_4*input;
            s[0] = a/den;
            s[1] = b/den;
            s[2] = c/den;
            s[3] = d/den;
            return out;
        }

        DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1)
        {
            DSP F = fc;
            DSP R = res;			
            setCutoff(F*fabs(X));
            setResonance(R*fabs(Y));
            DSP out = processSample(I,0);
            return A*out;
        }
        void ProcessSIMD(size_t n, DSP * inp, DSP * outp, size_t channel=0)
        {
			Undenormal denormal;
			#pragma omp simd aligned(inp,outp)
			for(size_t i = 0; i < n; i++)
			{
				auto& s = state[channel];
				const DSP input = inp[i];
				const DSP g_2 = g*g;
				const DSP g_3 = g_2*g;
				const DSP g_4 = g_3*g;
				const DSP r_cat_2 = r_cat*r_cat;
				const DSP val1 = res*g_4*r_cat_2;
				const DSP val2 = g_3*r_cat;
				const DSP val3 = g*r_cat;
				const DSP val4 = g_2*r_cat_2;
				const DSP val5 = 2*val3 + 1;
				const DSP val6 = g_2 + val5;
				const DSP val7 = g*val6;
				const auto den = (4*val1 + g_4 + 4*val2 + 4*val4 + 2*g_2 + 4*val3 + 1);
				DSP out = (g_3*s[0] - g_2*(val5)*s[1] + val7*s[2] - (2*val2 + 4*val4 + g_2 + 4*val3 + 1)*s[3])/den;
				const DSP a = -(val1*4 + g_4 + 4*val2 + 4*val4 - 1)*s[0] + 2*g*(res*val4*4 + val6)*s[1] - 8*g_2*res*r_cat_2*s[2] + 8*g*res*r_cat_2*(2*val3+1)*s[3] + val7*2*input;
				const DSP b = - 2*val7*s[0] + (-val1*4 - g_4 + 4*val4 + 4*val3 + 1)*s[1] + 8*g_3*res*r_cat_2*s[2] - 8*g_2*res*r_cat_2*(2*val3+1)*s[3] - g_2*(val6)*2*input;
				const DSP c = 2*g_2*s[0] - g*(val5)*2*s[1] - (val1*4 + g_4 + 4*g_3*r_cat + 4*val4 - 1)*s[2] + 2*g*(res*val4*4 + val6)*s[3] + 2*g_3*input;
				const DSP d = -2*g_3*s[0] + g_2*(val5)*2*s[1] - val7*2*s[2] + (-val1*4 - g_4 + 4*val4 + 4*val3 + 1)*s[3] - 2*g_4*input;
				s[0] = a/den;
				s[1] = b/den;
				s[2] = c/den;
				s[3] = d/den;
				outp[i] = out;
			}
		}
		void ProcessBlock(size_t n, DSP * inp, DSP * outp) {
			ProcessSIMD(n,inp,outp);
		}
		void ProcessInplace(size_t n, DSP * inp) {
			ProcessSIMD(n,inp,inp);
		}		
    };
}
