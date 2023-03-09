#pragma once

namespace Analog::Moog
{
	template<typename DSP>
	class VAOnePole
	{
	public:

		VAOnePole(DspFloatType sr) : sampleRate(sr)
		{
			Reset();
		}

		void Reset()
		{
			alpha = 1.0;
			beta = 0.0;
			gamma = 1.0;
			delta = 0.0;
			epsilon = 0.0;
			a0 = 1.0;
			feedback = 0.0;
			z1 = 0.0;
		}

		DspFloatType Tick(DspFloatType s)
		{
			Undenormal denormal;
			s = s * gamma + feedback + epsilon * GetFeedbackOutput();
			DspFloatType vn = (a0 * s - z1) * alpha;
			DspFloatType out = vn + z1;
			z1 = vn + out;
			return out;
		}

		void SetFeedback(DspFloatType fb) { feedback = fb; }
		DspFloatType GetFeedbackOutput(){ return beta * (z1 + feedback * delta); }
		void SetAlpha(DspFloatType a) { alpha = a; };
		void SetBeta(DspFloatType b) { beta = b; };

	private:

		DspFloatType sampleRate;
		DspFloatType alpha;
		DspFloatType beta;
		DspFloatType gamma;
		DspFloatType delta;
		DspFloatType epsilon;
		DspFloatType a0;
		DspFloatType feedback;
		DspFloatType z1;
	};

	template<typename DSP>
	class OberheimVariationMoog : public LadderFilterBase
	{

	public:

		OberheimVariationMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			LPF1 = new VAOnePole(sampleRate);
			LPF2 = new VAOnePole(sampleRate);
			LPF3 = new VAOnePole(sampleRate);
			LPF4 = new VAOnePole(sampleRate);

			saturation = 1.0;
			Q = 3.0;

			SetCutoff(1000.f);
			SetResonance(0.1f);
		}

		virtual ~OberheimVariationMoog()
		{
			delete LPF1;
			delete LPF2;
			delete LPF3;
			delete LPF4;
		}

		void Process(size_t n, DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				DspFloatType input = samples[s];

				DspFloatType sigma =
					LPF1->GetFeedbackOutput() +
					LPF2->GetFeedbackOutput() +
					LPF3->GetFeedbackOutput() +
					LPF4->GetFeedbackOutput();

				input *= 1.0 + K;

				// calculate input to first filter
				DspFloatType u = (input - K * sigma) * alpha0;

				u = tanh(saturation * u);

				DspFloatType stage1 = LPF1->Tick(u);
				DspFloatType stage2 = LPF2->Tick(stage1);
				DspFloatType stage3 = LPF3->Tick(stage2);
				DspFloatType stage4 = LPF4->Tick(stage3);

				// Oberheim variations
				output[s] =
					oberheimCoefs[0] * u +
					oberheimCoefs[1] * stage1 +
					oberheimCoefs[2] * stage2 +
					oberheimCoefs[3] * stage3 +
					oberheimCoefs[4] * stage4;
			}
		}
		void Process(size_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				DspFloatType input = samples[s];

				DspFloatType sigma =
					LPF1->GetFeedbackOutput() +
					LPF2->GetFeedbackOutput() +
					LPF3->GetFeedbackOutput() +
					LPF4->GetFeedbackOutput();

				input *= 1.0 + K;

				// calculate input to first filter
				DspFloatType u = (input - K * sigma) * alpha0;

				u = tanh(saturation * u);

				DspFloatType stage1 = LPF1->Tick(u);
				DspFloatType stage2 = LPF2->Tick(stage1);
				DspFloatType stage3 = LPF3->Tick(stage2);
				DspFloatType stage4 = LPF4->Tick(stage3);

				// Oberheim variations
				samples[s] =
					oberheimCoefs[0] * u +
					oberheimCoefs[1] * stage1 +
					oberheimCoefs[2] * stage2 +
					oberheimCoefs[3] * stage3 +
					oberheimCoefs[4] * stage4;
			}
		}
		DspFloatType Tick(DspFloatType input) {
			DspFloatType r;
			Process(1,&input,&r);
			return r;
		}
		virtual void SetResonance(DspFloatType r) override
		{
				// this maps resonance = 1->10 to K = 0 -> 4
				K = (4.0) * (r - 1.0)/(10.0 - 1.0);
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c;

			// prewarp for BZT
			DspFloatType wd = 2.0 * MOOG_PI * cutoff;
			DspFloatType T = 1.0 / sampleRate;
			DspFloatType wa = (2.0 / T) * tan(wd * T / 2.0);
			DspFloatType g = wa * T / 2.0;

			// Feedforward coeff
			DspFloatType G = g / (1.0 + g);

			LPF1->SetAlpha(G);
			LPF2->SetAlpha(G);
			LPF3->SetAlpha(G);
			LPF4->SetAlpha(G);

			LPF1->SetBeta(G*G*G / (1.0 + g));
			LPF2->SetBeta(G*G / (1.0 + g));
			LPF3->SetBeta(G / (1.0 + g));
			LPF4->SetBeta(1.0 / (1.0 + g));

			gamma = G*G*G*G;
			alpha0 = 1.0 / (1.0 + K * gamma);

			// Oberheim variations / LPF4
			oberheimCoefs[0] = 0.0;
			oberheimCoefs[1] = 0.0;
			oberheimCoefs[2] = 0.0;
			oberheimCoefs[3] = 0.0;
			oberheimCoefs[4] = 1.0;
		}

	private:

		VAOnePole * LPF1;
		VAOnePole * LPF2;
		VAOnePole * LPF3;
		VAOnePole * LPF4;

		DspFloatType K;
		DspFloatType gamma;
		DspFloatType alpha0;
		DspFloatType Q;
		DspFloatType saturation;

		DspFloatType oberheimCoefs[5];
	};
}
