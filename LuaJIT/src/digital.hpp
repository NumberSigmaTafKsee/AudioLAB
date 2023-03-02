pragma once

struct ClockMultiplier {
	uint32_t clock = 0;
	uint32_t lastTickSamples = 0;
	double division = 0.f;
	double divisionMult = 0.f;
	uint32_t currentDivision = 0;

	bool process() {
		lastTickSamples++;
		if (division > 0 && currentDivision > divisionMult && currentDivision < clock) {
			currentDivision++;
			divisionMult += division;
			return true;
		}
		currentDivision++;
		return false;
	}

	void tick() {
		clock = lastTickSamples;
		lastTickSamples = 0;
		division = 0.f;
		divisionMult = 0.f;
		currentDivision = 0;
	}

	void trigger(uint32_t div) {
		if (clock == 0) return;
		division = 0.f;
		divisionMult = 0.f;
		if (div == 0) return;
		division = clock / double(div);
	}

	void reset() {
		clock = 0;
		lastTickSamples = 0;
		division = 0.f;
		divisionMult = 0.f;
		currentDivision = 0;
	}
};


struct LinearFade {
	double rise = 1.f;
	double fall = 1.f;
	double currentRise;
	double currentFall;
	double last = 0.f;

	void reset(double last) {
		currentRise = rise;
		currentFall = 0.f;
		this->last = last;
	}

	void triggerFadeIn() {
		currentRise = (fall > 0.f ? (currentFall / fall) : 0.f) * rise;
		currentFall = 0.f;
		last = 1.f;
	}

	void triggerFadeOut() {
		currentFall = (rise > 0.f ? (currentRise / rise) : 0.f) * fall;
		currentRise = rise;
		last = 0.f;
	}

	inline void setRise(double rise) {
		if (currentRise == this->rise) currentRise = rise;
		this->rise = rise;
	}

	inline void setFall(double fall) {
		currentFall = std::min(fall, currentFall);
		this->fall = fall;
	}

	inline void setRiseFall(double rise, double fall) {
		setRise(rise);
		setFall(fall);
	}

	inline double process(double deltaTime) {
		if (currentRise < rise) {
			currentRise += deltaTime;
			return (currentRise / rise);
		}
		else if (currentFall > 0.f) {
			currentFall = std::max(currentFall - deltaTime, 0.f);
			return (currentFall / fall);
		}
		else {
			return last;
		}
	}
};



struct StoermelderSlewLimiter {
	// Minimum and maximum slopes in volts per second
	const double slewMin = 0.1;
	const double slewMax = 10000.f;
	// Amount of extra slew per voltage difference
	const double shapeScale = 1/10.f;

	double shape	= 0.5f;
	double rise = 0.0f;
	double fall = 0.0f;

	double out = 0.0;

	inline void reset() {
		out = 0.f;
	}
	inline void setShape(double shape) {
		this->shape = shape;
	}
	inline void setRise(double rise) {
		this->rise = rise;
	}
	inline void setFall(double fall) {
		this->fall = fall;
	}
	inline void setRiseFall(double rise, double fall) {
		this->rise = rise;
		this->fall = fall;
	}

	double process(double in, double sampleTime) {
		// Rise
		if (in > out) {
			double slew = slewMax * std::pow(slewMin / slewMax, rise);
			out += slew * crossfade(1.f, shapeScale * (in - out), shape) * sampleTime;
			if (out > in)
				out = in;
		}
		// Fall
		else if (in < out) {
			double slew = slewMax * std::pow(slewMin / slewMax, fall);
			out -= slew * crossfade(1.f, shapeScale * (out - in), shape) * sampleTime;
			if (out < in)
				out = in;
		}
		return out;
	}
};