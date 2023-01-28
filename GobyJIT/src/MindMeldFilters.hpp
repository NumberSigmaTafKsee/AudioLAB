// https://github.com/MarcBoule/MindMeldModular
#pragma once

class FirstOrderCoefficients {
	protected: 

	float b[2];// coefficients b0, b1
	float a;// coefficient a1
	
	
	public: 
	
	void setParameters(bool isHighPass, float nfc) {// normalized freq
		// nfc: normalized cutoff frequency (cutoff frequency / sample rate), must be > 0
		// freq pre-warping with inclusion of M_PI factor; 
		//   avoid tan() if fc is low (< 1102.5 Hz @ 44.1 kHz, since error at this freq is 2 Hz)
		float nfcw = nfc < 0.025f ? float(M_PI) * nfc : std::tan(float(M_PI) * std::min(0.499f, nfc));
		
		// denominator coefficient (same for both LPF and HPF)
		a = (nfcw - 1.0f) / (nfcw + 1.0f);
		
		// numerator coefficients
		float hbcst = 1.0f / (1.0f + nfcw);
		float lbcst = 1.0f - hbcst;// equivalent to: hbcst * nfcw;
		b[0] = isHighPass ? hbcst : lbcst;
		b[1] = isHighPass ? -hbcst : lbcst;
	}		
};


class FirstOrderFilter : public FirstOrderCoefficients {
	float x;
	float y;
	
	public: 
		
	void reset() {
		x = 0.0f;
		y = 0.0f;
	}

	float process(float in) {
		y = b[0] * in + b[1] * x - a * y;
		x = in;
		return y;
	}
};


class FirstOrderStereoFilter : public FirstOrderCoefficients {
	float x[2];
	float y[2];
	
	public: 
		
	void reset() {
		x[0] = 0.0f;
		x[1] = 0.0f;
		y[0] = 0.0f;
		y[1] = 0.0f;
	}
	
	void process(float* out, float* in) {
		y[0] = b[0] * in[0] + b[1] * x[0] - a * y[0];
		y[1] = b[0] * in[1] + b[1] * x[1] - a * y[1];
		x[0] = in[0];
		x[1] = in[1];
		out[0] = y[0];
		out[1] = y[1];
	}
};

class ButterworthSecondOrder {
	float b[3];// coefficients b0, b1 and b2
	float a[3 - 1];// coefficients a1 and a2
	float x[3 - 1];
	float y[3 - 1];
	float midCoef = float(M_SQRT2);
	
	public:
	
	void setMidCoef(float _midCoef) {
		midCoef = _midCoef;
	}
	
	void reset() {
		for (int i = 0; i < 2; i++) {
			x[i] = 0.0f;
			y[i] = 0.0f;
		}
	}

	void setParameters(bool isHighPass, float nfc) {// normalized freq
		// nfc: normalized cutoff frequency (cutoff frequency / sample rate), must be > 0
		// freq pre-warping with inclusion of M_PI factor; 
		//   avoid tan() if fc is low (< 1102.5 Hz @ 44.1 kHz, since error at this freq is 2 Hz)
		float nfcw = nfc < 0.025f ? float(M_PI) * nfc : std::tan(float(M_PI) * std::min(0.499f, nfc));

		// denominator coefficients (same for both LPF and HPF)
		float acst = nfcw * nfcw + nfcw * midCoef + 1.0f;
		a[0] = 2.0f * (nfcw * nfcw - 1.0f) / acst;
		a[1] = (nfcw * nfcw - nfcw * midCoef + 1.0f) / acst;
		
		// numerator coefficients
		float hbcst = 1.0f / acst;
		float lbcst = hbcst * nfcw * nfcw;			
		b[0] = isHighPass ? hbcst : lbcst;
		b[1] = (isHighPass ? -hbcst : lbcst) * 2.0f;
		b[2] = b[0];
	}
	
	float process(float in) {
		float out = b[0] * in + b[1] * x[0] + b[2] * x[1] - a[0] * y[0] - a[1] * y[1];
		x[1] = x[0];
		x[0] = in;
		y[1] = y[0];
		y[0] = out;
		return out;
	}
};


class ButterworthThirdOrder {
	FirstOrderFilter f1;
	ButterworthSecondOrder f2;
	
	public:
	
	ButterworthThirdOrder() {
		f2.setMidCoef(1.0f);
	}
	
	void reset() {
		f1.reset();
		f2.reset();
	}
	
	void setParameters(bool isHighPass, float nfc) {// normalized freq
		f1.setParameters(isHighPass, nfc);
		f2.setParameters(isHighPass, nfc);
	}
	
	float process(float in) {
		return f2.process(f1.process(in));
	}
};


class ButterworthFourthOrder {
	ButterworthSecondOrder f1;
	ButterworthSecondOrder f2;
	
	public:
	
	ButterworthFourthOrder() {
		f1.setMidCoef(0.765367f);
		f2.setMidCoef(1.847759f);
	}
	
	void reset() {
		f1.reset();
		f2.reset();
	}
	
	void setParameters(bool isHighPass, float nfc) {// normalized freq
		f1.setParameters(isHighPass, nfc);
		f2.setParameters(isHighPass, nfc);
	}
	
	float process(float in) {
		return f2.process(f1.process(in));
	}
};


class LinkwitzRileyCoefficients {
	protected:
	
	bool secondOrderFilters = false;// local memory of what is in iirs		
	
	simd::float_4 b[3];// coefficients b0, b1 and b2, where each float of float_4 is LeftLow, LeftHigh, RightLow, RightHigh
	simd::float_4 a[3 - 1];// coefficients a1 and a2, where each float of float_4 is LeftLow, LeftHigh, RightLow, RightHigh


	public: 
	
	void setFilterCutoffs(float nfc, bool _secondOrder) {
		secondOrderFilters = _secondOrder;
		
		// nfc: normalized cutoff frequency (cutoff frequency / sample rate), must be > 0
		// freq pre-warping with inclusion of M_PI factor; 
		//   avoid tan() if fc is low (< 1102.5 Hz @ 44.1 kHz, since error at this freq is 2 Hz)
		float nfcw = nfc < 0.025f ? float(M_PI) * nfc : std::tan(float(M_PI) * std::min(0.499f, nfc));
		
		if (secondOrderFilters) {	
			// denominator coefficients (same for both LPF and HPF)
			float acst = nfcw * nfcw + nfcw * float(M_SQRT2) + 1.0f;
			a[0] = simd::float_4(2.0f * (nfcw * nfcw - 1.0f) / acst);
			a[1] = simd::float_4((nfcw * nfcw - nfcw * float(M_SQRT2) + 1.0f) / acst);
			
			// numerator coefficients
			float hbcst = 1.0f / acst;
			float lbcst = hbcst * nfcw * nfcw;			
			b[0] = simd::float_4(lbcst, hbcst, lbcst, hbcst);
			b[1] = simd::float_4(lbcst, -hbcst, lbcst, -hbcst) * 2.0f;
			b[2] = b[0];
		}
		else {
			// denominator coefficients (same for both LPF and HPF)
			float acst = (nfcw - 1.0f) / (nfcw + 1.0f);
			a[0] = simd::float_4(acst);
			a[1] = simd::float_4(0.0f);
			
			// numerator coefficients
			float hbcst = 1.0f / (1.0f + nfcw);
			float lbcst = 1.0f - hbcst;// equivalent to: hbcst * nfcw;
			b[0] = simd::float_4(lbcst, hbcst, lbcst, hbcst);
			b[1] = simd::float_4(lbcst, -hbcst, lbcst, -hbcst);
			b[2] = simd::float_4(0.0f);
		}
	}
};


class LinkwitzRileyStereoCrossover : public LinkwitzRileyCoefficients {	
	simd::float_4 xS1[3 - 1];
	simd::float_4 yS1[3 - 1];
	simd::float_4 xS2[3 - 1];
	simd::float_4 yS2[3 - 1];
	
	
	public: 
		
	void reset() {
		for (int i = 0; i < 2; i++) {
			xS1[i] = 0.0f;
			yS1[i] = 0.0f;
			xS2[i] = 0.0f;
			yS2[i] = 0.0f;
		}
	}

	simd::float_4 process(float left, float right) {
		// return [0] = left low, left high, right low, [3] = right high
		simd::float_4 in = simd::float_4(left, left, right, right);
		if (!secondOrderFilters) {
			in[0] *= -1.0f;// phase correction needed for first order filters (used to make 2nd order L-R crossover)
			in[2] *= -1.0f;
		}

		// stage 1
		simd::float_4 outS1 = b[0] * in + b[1] * xS1[0] + b[2] * xS1[1] - a[0] * yS1[0] - a[1] * yS1[1];
		xS1[1] = xS1[0];
		xS1[0] = in;
		yS1[1] = yS1[0];
		yS1[0] = outS1;

		// stage 2 (outS1 used as in)
		simd::float_4 outS2 = b[0] * outS1 + b[1] * xS2[0] + b[2] * xS2[1] - a[0] * yS2[0] - a[1] * yS2[1];
		xS2[1] = xS2[0];
		xS2[0] = outS1;
		yS2[1] = yS2[0];
		yS2[0] = outS2;

		return outS2;
	}
};


class LinkwitzRileyStereo8xCrossover : public LinkwitzRileyCoefficients {	
	simd::float_4 xS1[8][3 - 1];
	simd::float_4 yS1[8][3 - 1];
	simd::float_4 xS2[8][3 - 1];
	simd::float_4 yS2[8][3 - 1];
	
	
	public: 
		
	void reset() {
		for (int c = 0; c < 8; c++) {
			for (int i = 0; i < 2; i++) {
				xS1[c][i] = 0.0f;
				yS1[c][i] = 0.0f;
				xS2[c][i] = 0.0f;
				yS2[c][i] = 0.0f;
			}
		}
	}

	simd::float_4 process(float left, float right, int c) {
		// return [0] = left low, left high, right low, [3] = right high
		simd::float_4 in = simd::float_4(left, left, right, right);
		if (!secondOrderFilters) {
			in[0] *= -1.0f;// phase correction needed for first order filters (used to make 2nd order L-R crossover)
			in[2] *= -1.0f;
		}

		
		// stage 1
		simd::float_4 outS1 = b[0] * in + b[1] * xS1[c][0] + b[2] * xS1[c][1] - a[0] * yS1[c][0] - a[1] * yS1[c][1];
		xS1[c][1] = xS1[c][0];
		xS1[c][0] = in;
		yS1[c][1] = yS1[c][0];
		yS1[c][0] = outS1;

		// stage 2 (outS1 used as in)
		simd::float_4 outS2 = b[0] * outS1 + b[1] * xS2[c][0] + b[2] * xS2[c][1] - a[0] * yS2[c][0] - a[1] * yS2[c][1];
		xS2[c][1] = xS2[c][0];
		xS2[c][0] = outS1;
		yS2[c][1] = yS2[c][0];
		yS2[c][0] = outS2;

		return outS2;
	}
};

class QuattroBiQuadCoeff {
	protected: 
	
	// coefficients
	simd::float_4 b0;
	simd::float_4 b1;
	simd::float_4 b2;
	simd::float_4 a1;
	simd::float_4 a2;
	
	
	public: 
	
	
	enum Type {
		LOWSHELF,
		HIGHSHELF,
		PEAK,
	};
	
	
	void setParameters(int i, Type type, float nfc, float V, float Q) {
		// i: eq index (0 to 3),
		// type: type of filter/eq
		// nfc: normalized cutoff frequency (fc/sampleRate)
		// V: linearGain for peak or shelving
		// Q: quality factor
		
		// nfc: normalized cutoff frequency (cutoff frequency / sample rate), must be > 0
		// freq pre-warping with inclusion of M_PI factor; 
		//   avoid tan() if fc is low (< 1102.5 Hz @ 44.1 kHz, since error at this freq is 2 Hz)
		float K = nfc < 0.025f ? float(M_PI) * nfc : std::tan(float(M_PI) * std::min(0.499f, nfc));

		switch (type) {
			case LOWSHELF: {
				float sqrtV = std::sqrt(V);
				Q = std::sqrt(Q) / float(M_SQRT2);
				if (V >= 1.f) {// when V = 1, b0 = 1, a1 = b1, a2 = b2
					float norm = 1.f / (1.f + K / Q + K * K);
					b0[i] = (1.f + sqrtV * K / Q + V * K * K) * norm;
					b1[i] = 2.f * (V * K * K - 1.f) * norm;
					b2[i] = (1.f - sqrtV * K / Q + V * K * K) * norm;
					a1[i] = 2.f * (K * K - 1.f) * norm;
					a2[i] = (1.f - K / Q + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f + K / (Q * sqrtV) + K * K / V);
					b0[i] = (1.f + K / Q + K * K) * norm;
					b1[i] = 2.f * (K * K - 1) * norm;
					b2[i] = (1.f - K / Q + K * K) * norm;
					a1[i] = 2.f * (K * K / V - 1.f) * norm;
					a2[i] = (1.f - K / (Q * sqrtV) + K * K / V) * norm;
				}
			} break;

			case HIGHSHELF: {
				float sqrtV = std::sqrt(V);
				Q = std::sqrt(Q) / float(M_SQRT2);
				if (V >= 1.f) {// when V = 1, b0 = 1, a1 = b1, a2 = b2
					float norm = 1.f / (1.f + K / Q + K * K);
					b0[i] = (V + sqrtV * K / Q + K * K) * norm;
					b1[i] = 2.f * (K * K - V) * norm;
					b2[i] = (V - sqrtV * K / Q + K * K) * norm;
					a1[i] = 2.f * (K * K - 1.f) * norm;
					a2[i] = (1.f - K / Q + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f / V + K / (Q * sqrtV) + K * K);
					b0[i] = (1.f + K / Q + K * K) * norm;
					b1[i] = 2.f * (K * K - 1.f) * norm;
					b2[i] = (1.f - K / Q + K * K) * norm;
					a1[i] = 2.f * (K * K - 1.f / V) * norm;
					a2[i] = (1.f / V - K / (Q * sqrtV) + K * K) * norm;
				}
			} break;

			case PEAK: {
				if (V >= 1.f) {
					float norm = 1.f / (1.f + K / Q + K * K);
					b0[i] = (1.f + K / Q * V + K * K) * norm;
					b1[i] = 2.f * (K * K - 1.f) * norm;
					b2[i] = (1.f - K / Q * V + K * K) * norm;
					a1[i] = b1[i];
					a2[i] = (1.f - K / Q + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f + K / Q / V + K * K);
					b0[i] = (1.f + K / Q + K * K) * norm;
					b1[i] = 2.f * (K * K - 1.f) * norm;
					b2[i] = (1.f - K / Q + K * K) * norm;
					a1[i] = b1[i];
					a2[i] = (1.f - K / Q / V + K * K) * norm;
				}
			} break;

			default: break;
		}
	}


	// add all 4 values in return vector to get total gain (dB) since each float is gain (dB) of one biquad
	simd::float_4 getFrequencyResponse(float f) {
		// Compute sum(b_k z^-k) / sum(a_k z^-k) where z = e^(i s)
		
		float s = 2 * float(M_PI) * f;// s: normalized angular frequency equal to $2 \pi f / f_{sr}$ ($\pi$ is the Nyquist frequency)
		
		simd::float_4 bSum[2] = {b0, 0.0f};
		simd::float_4 aSum[2] = {1.0f, 0.0f};
			
		float p = -1 * s;
		simd::float_4 z[2] = {std::cos(p), std::sin(p)};
		bSum[0] += b1 * z[0];
		bSum[1] += b1 * z[1];
		aSum[0] += a1 * z[0];
		aSum[1] += a1 * z[1];
	
		p = -2 * s;
		z[0] = std::cos(p); z[1] = std::sin(p);
		bSum[0] += b2 * z[0];
		bSum[1] += b2 * z[1];
		aSum[0] += a2 * z[0];
		aSum[1] += a2 * z[1];
		
		simd::float_4 num[2] = {bSum[0] * aSum[0] +	bSum[1] * aSum[1],  bSum[1] * aSum[0] - bSum[0] * aSum[1]};
		simd::float_4 denom = aSum[0] * aSum[0] + aSum[1] * aSum[1];
		simd::float_4 norm = simd::hypot(num[0] / denom,  num[1] / denom);
		return 20.0f * simd::log10(norm);// return in dB
	}
};



// Four stereo biquad filters in pipeline series, where each biquad's parameters can be set separately
class QuattroBiQuad : public QuattroBiQuadCoeff {
	
	// input/output shift registers
	simd::float_4 x0L, x0R;// input is x0[0]
	simd::float_4 x1L, x1R;
	simd::float_4 x2L, x2R;
	simd::float_4 y0L, y0R;// output is y0[3]
	simd::float_4 y1L, y1R;
	simd::float_4 y2L, y2R;
	
	// other
	bool optResetDone;
	int8_t gainsDifferentThanOne; // 4 ls bits are bool bits, when all zero, can bypass y0 math
	
	
	public:


	void reset() {
		x0L = 0.0f;
		x0R = 0.0f;
		x1L = 0.0f;
		x1R = 0.0f;
		x2L = 0.0f;
		x2R = 0.0f;
		y0L = 0.0f;
		y0R = 0.0f;
		y1L = 0.0f;
		y1R = 0.0f;
		y2L = 0.0f;
		y2R = 0.0f;
		gainsDifferentThanOne = 0xF;	
		optResetDone = false;
	}
	
	
	void setParameters(int i, Type type, float f, float V, float Q) {
		// type: type of filter/eq
		// i: eq index (0 to 3),
		// f: normalized frequency (fc/sampleRate)
		// V: linearGain for peak or shelving
		// Q: quality factor
		if (V == 1.0f) {
			gainsDifferentThanOne &= ~(0x1 << i);
		}
		else {
			gainsDifferentThanOne |= (0x1 << i);
		}
		QuattroBiQuadCoeff::setParameters(i, type, f, V, Q);
	}
	

	void process(float* out, float* in) {
		if (gainsDifferentThanOne == 0) {
			if (!optResetDone) {
				x2L = 0.0f;
				x2R = 0.0f;
				x1L = 0.0f;
				x1R = 0.0f;
				y2L = 0.0f;
				y2R = 0.0f;
				y1L = 0.0f;
				y1R = 0.0f;
				optResetDone = true;
			}
			
			x0L[0] = in[0];
			x0L[1] = y0L[0];
			x0L[2] = y0L[1];
			x0L[3] = y0L[2];
			y0L = x0L;
			
			x0R[0] = in[1];
			x0R[1] = y0R[0];
			x0R[2] = y0R[1];
			x0R[3] = y0R[2];
			y0R = x0R;
		}
		else {
			optResetDone = false;
			x2L = x1L;
			x1L = x0L;
			x0L[0] = in[0];
			x0L[1] = y0L[0];
			x0L[2] = y0L[1];
			x0L[3] = y0L[2];
			y2L = y1L;
			y1L = y0L;
			y0L = b0 * x0L + b1 * x1L + b2 * x2L - a1 * y1L - a2 * y2L;// https://en.wikipedia.org/wiki/Infinite_impulse_response

			x2R = x1R;
			x1R = x0R;
			x0R[0] = in[1];
			x0R[1] = y0R[0];
			x0R[2] = y0R[1];
			x0R[3] = y0R[2];
			y2R = y1R;
			y1R = y0R;
			y0R = b0 * x0R + b1 * x1R + b2 * x2R - a1 * y1R - a2 * y2R;// https://en.wikipedia.org/wiki/Infinite_impulse_response
		}
		
		out[0] = y0L[3];
		out[1] = y0R[3];
	}	
};

