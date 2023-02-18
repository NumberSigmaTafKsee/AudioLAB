// https://github.com/MrJDunn/Compressor
class CompressorBase
{
public:
	CompressorBase();
	virtual ~CompressorBase();

	virtual void setup(double sampleRate, int blockSize);
	virtual void process(AudioBuffer<float>&);

	void setAttack(float);
	void setRelease(float);
	void setRatio(float);
	void setThreshold(float);
		
	float getAttack();
	float getRelease();
	float getRatio();
	float getThreshold();

	float getCompressionAmountDB();

protected:
	/**
		0 < attack < 5
	*/
	float attack;

	/**
		0 < release < 5
	*/
	float release;

	/**
		2 < ratio < 16
	*/
	float ratio;

	/** 
		-96 < threshold < 0
	*/
	float threshold;

private:
	float attackAlpha = 0.99f;
	float releaseAlpha = 0.99f;

	float gainSmoothPrevious = 0.0f;

	double sampleRate = 44100.0;
	int blockSize = 512;
	float reductionAmount = 0.f;
	float minimumDb = -40.f;
	std::atomic<float> lastGreatestReductionAmount = 0.f;
};

CompressorBase::CompressorBase(): threshold(-40), ratio(4), attack(1), release(1)
{
}

CompressorBase::~CompressorBase()
{
}

void CompressorBase::setup(double sampleRate, int blockSize)
{
	this->sampleRate = sampleRate;
	this->blockSize = blockSize;
}

void CompressorBase::process(AudioBuffer<float>& buffer)
{
	int channels = buffer.getNumChannels();
	int bufferSize = buffer.getNumSamples();

	attackAlpha = exp(-1.f / (sampleRate * 0.001f * attack));
	releaseAlpha = exp(-1.f / (sampleRate * 0.001f * release));

	lastGreatestReductionAmount = 0.f;
	
	for (int i = 0; i < channels; ++i)
	{
		float* channelData = buffer.getWritePointer(i);
		for (int j = 0; j < bufferSize; ++j)
		{
			// Calculate db value of sample
			float unipolarInput = abs(channelData[j]);
			float inputDb = 20 * log10(unipolarInput);

			// Avoid negative infinity
			inputDb = std::max(inputDb, minimumDb);

			// Calculate new target gain if greater than threshold
			float outputDb = inputDb > threshold ? threshold + ((inputDb - threshold) / ratio) : inputDb;

			reductionAmount = inputDb - outputDb;

			lastGreatestReductionAmount.store(std::max(reductionAmount, lastGreatestReductionAmount.load()));
		//	DBG(reductionAmount);

			float gainSmooth = 0.f;

			if (reductionAmount < gainSmoothPrevious)
			{
				// Attack
				gainSmooth = (((1.f - attackAlpha) * reductionAmount) + (attackAlpha * gainSmoothPrevious));
			}
			else
			{
				// Release
				gainSmooth = (((1.f - releaseAlpha) * reductionAmount) + (releaseAlpha * gainSmoothPrevious));
			}

			float newSampleValue = powf(10.f, (gainSmooth / 20.f));

			//newSampleValue = std::max(-1.f,std::min(1.f,newSampleValue));
			//jassert(newSampleValue <= 1.f && newSampleValue >= -1.f);

			channelData[j] *= newSampleValue;

			gainSmoothPrevious = gainSmooth;

		}
	}
}

void CompressorBase::setAttack(float val)
{
	attack = std::min(5.f,std::max(1.f, val * 5.f));
}

void CompressorBase::setRelease(float val)
{
	release = std::min(5.f,std::max(1.f, val * 5.f));
}

void CompressorBase::setRatio(float val)
{
	ratio = std::min(16.f,std::max(2.f, val));
}

void CompressorBase::setThreshold(float val)
{
	threshold = std::min(0.f, std::max(minimumDb, val));
}

float CompressorBase::getAttack()
{
	return attack;
}

float CompressorBase::getRelease()
{
	return release;
}

float CompressorBase::getRatio()
{
	return ratio;
}

float CompressorBase::getThreshold()
{
	return threshold;
}

float CompressorBase::getCompressionAmountDB()
{
	return lastGreatestReductionAmount.load();
}
