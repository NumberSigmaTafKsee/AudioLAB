// https://github.com/Hippasusss/Auto
class Helpers
{
public:
	// Returns the absolute average value of samples from an AudioBlock or AudioBuffer
    template<typename SampleType> 
	static SampleType getAverageMagnitude(const dsp::AudioBlock<SampleType>& block)
    {
		SampleType sum {0};
	    const auto numChannels = block.getNumChannels();
	    const auto blockSize = block.getNumSamples();

	    for(auto i = 0; i < numChannels; ++i)
	    {
		    const auto chan = block.getChannelPointer(i);
	        for(auto j = 0; j < blockSize; ++j)
	        {
		        sum += abs(chan[j]);
	        }
	    }

		return sum / (blockSize * numChannels);
    }

	// ew
    template<typename SampleType> 
	static SampleType getAverageMagnitude(const AudioBuffer<SampleType>& block)
    {
		SampleType sum {0};
	    const auto numChannels = block.getNumChannels();
	    const auto blockSize = block.getNumSamples();

	    for(auto i = 0; i < numChannels; ++i)
	    {
		    const auto chan = block.getReadPointer(i);
	        for(auto j = 0; j < blockSize; ++j)
	        {
		        sum += abs(chan[j]);
	        }
	    }

		return sum / (blockSize * numChannels);
    }

	// Returns the highest sample magnitude of all channels
    template<typename SampleType>
	static SampleType getMagnitude(const dsp::AudioBlock<SampleType>& block)
    {
	    
		float value {};
	    const auto numChannels = block.getNumChannels();
	    const auto blockSize = block.getNumSamples();

	    for(auto i = 0; i < numChannels; ++i)
	    {
		    const auto chan = block.getChannelPointer(i);
	        for(auto j = 0; j < blockSize; ++j)
	        {
				const SampleType sample = abs(chan[j]);
		        value = sample > value ? sample : value;
	        }
	    }

		return value / (blockSize * numChannels);
    }

	// Copies the audio from the provided AudioBlock into the provided audio buffer
	template<typename SampleType>
	static void copyAudioBlockIntoBuffer(const dsp::AudioBlock<const SampleType>& sourceBlock,
	                                     AudioBuffer<SampleType>& destinationBuffer,
	                                     const size_t numSamples,
	                                     const size_t sourceStartSample = 0,
	                                     const size_t destStartSample = 0)
	{
	    jassert(sourceBlock.getNumChannels() == destinationBuffer.getNumChannels());
		for(auto channel = 0; channel < destinationBuffer.getNumChannels(); ++channel)
		{
			auto* channelPointer = sourceBlock.getChannelPointer(channel);
			channelPointer += sourceStartSample;
			destinationBuffer.copyFrom(channel, destStartSample, channelPointer, numSamples);
		}
	}

	template<typename SampleType>
	static void sumChannelsToFirstChannel(AudioBuffer<SampleType>& buffer)
    {
		const size_t numChannels = buffer.getNumChannels();
		const size_t numSamples = buffer.getNumSamples();
		auto* const dataChannelWrite = buffer.getWritePointer(0);
		for (size_t i = 1; i < numChannels; ++i)
		{
			const auto* const dataChannelRead = buffer.getReadPointer(i);
			for (size_t j = 0; j < numSamples; ++j)
			{
				dataChannelWrite[j] += dataChannelRead[j];
			}
		}

		// Average in place
		for (size_t i = 0; i < numSamples; ++i)
		{
			dataChannelWrite[i] /= numChannels;
		}
	}

	template<typename ValueType>
	static ValueType getNormalisedDB(ValueType value, ValueType dbMinusInfinity = -100)
    {
	    return jlimit<ValueType>(0.0f, 1.0f, ((ValueType(20.0) * std::log10f(abs(value)) / ValueType(-dbMinusInfinity)) + 1));
    }

private:
	Helpers() = default;
	~Helpers() = default;
};


template <typename ValueType, typename ContainerType>
class RingBuffer
{
public:
	virtual ~RingBuffer() = default;
	virtual void writeValue(ValueType);
	virtual ValueType readValue();
	virtual void readPreviousValues(ContainerType& values);

	virtual size_t getSize() const;
	virtual ValueType operator[](size_t i);
protected:
	RingBuffer();
	RingBuffer(size_t size);
	ContainerType valueArray;
	size_t writeIndex;
	size_t readIndex;
	size_t size;
};

//===============================================================================

template <typename ValueType, typename ContainerType>
RingBuffer<ValueType, ContainerType>::RingBuffer() : writeIndex(0), readIndex(0), size(0)
{
}

template <typename ValueType, typename ContainerType>
RingBuffer<ValueType, ContainerType>::RingBuffer(size_t size) : valueArray(size, 0), writeIndex(0), readIndex(0), size(size)
{
}

template <typename ValueType, typename ContainerType>
void RingBuffer<ValueType, ContainerType>::writeValue(ValueType value)
{
	valueArray[writeIndex] = value;
	++writeIndex %= size;
}

template <typename ValueType, typename ContainerType>
ValueType RingBuffer<ValueType, ContainerType>::readValue()
{
	const ValueType returnValue = valueArray[readIndex];
	++readIndex %= size;
	return returnValue;
}


template <typename ValueType, typename ContainerType>
size_t RingBuffer<ValueType, ContainerType>::getSize() const
{
	return size;
}


template <typename ValueType, typename ContainerType>
ValueType RingBuffer<ValueType, ContainerType>::operator[](size_t i)
{
	return valueArray[i];
}

// copies into provided array/container. copies as many as can into given size
// --------------------|---------------------------|-|----------
//                     ^ size             <i--- -1 ^ ^ write pointer
//                     |------This is copied-------|
//
template <typename ValueType, typename ContainerType>
void RingBuffer<ValueType, ContainerType>::readPreviousValues(ContainerType& values)
{
	const size_t inputSize = values.size();
	const size_t writeIndexLocal = writeIndex; // take local value in case class member is changed in separate thread.

	for (size_t i = 0; i < inputSize; i++)
	{
		// vile one liner takes care of wraparound of index when copying
		const size_t copyIndex = ((((writeIndexLocal - (i + 1)) % size) + size) % size);
		values[inputSize - 1 - i] = valueArray[copyIndex];
	}
}

//===============================================================================
// Vector
//===============================================================================

template <typename ValueType>
class RingBufferVector : public RingBuffer<ValueType, std::vector<ValueType>>
{
public:
	RingBufferVector();
	RingBufferVector(size_t newSize);
	void resize(size_t newSize);

};

//===============================================================================

template <typename ValueType>
RingBufferVector<ValueType>::RingBufferVector()
{
}

template <typename ValueType>
RingBufferVector<ValueType>::RingBufferVector(size_t newSize) : RingBuffer<ValueType, std::vector<ValueType>>(newSize)
{
	this->valueArray.resize(newSize);
}

template <typename ValueType>
void RingBufferVector<ValueType>::resize(size_t newSize)
{
	this->size = newSize;
	this->valueArray.resize(newSize);
}

//===============================================================================
// Array 
//===============================================================================

template <typename ValueType, size_t Size>
class RingBufferArray : public RingBuffer<ValueType, std::array<ValueType, Size>>
{
public:
	RingBufferArray();
};

template <typename ValueType, size_t Size>
RingBufferArray<ValueType, Size>::RingBufferArray()
{
	this->size = Size;
}

//===============================================================================
// Audio Block
//===============================================================================
template <typename SampleType>
class RingBufferAudio
{
public:
	RingBufferAudio() = default;
	RingBufferAudio(size_t channels, size_t length);
	void resize(size_t channels, size_t length);
	void writeBlock(const dsp::AudioBlock<const SampleType>& newBlock);
	void readBlock(AudioBuffer<SampleType>& bufferToFill);
	void getPreviousSamples(AudioBuffer<SampleType>& bufferToFill);
	dsp::AudioBlock<SampleType> getBlock();
	AudioBuffer<SampleType> getBuffer();
private:
	AudioBuffer<SampleType> aggregateBuffer;
	size_t numSamples = 0;
	size_t numChannels = 0;
	size_t writeIndex = 0;
	size_t readIndex = 0;
};

template <typename SampleType>
RingBufferAudio<SampleType>::RingBufferAudio(size_t channels, size_t length) :
	aggregateBuffer(channels, length), numSamples(length), numChannels(channels)
{
}

template <typename SampleType>
void RingBufferAudio<SampleType>::resize(size_t channels, size_t length)
{
	aggregateBuffer.setSize(channels, length);
	numSamples = length;
	numChannels = channels;
}

template <typename SampleType>
void RingBufferAudio<SampleType>::writeBlock(const dsp::AudioBlock<const SampleType>& newBlock)
{
	const size_t numSamplesToCopy = newBlock.getNumSamples();
	const size_t remainingSpace = numSamples - writeIndex;

	if(numSamplesToCopy > remainingSpace)
	{
		Helpers::copyAudioBlockIntoBuffer(newBlock, aggregateBuffer, remainingSpace, 0, writeIndex);
		Helpers::copyAudioBlockIntoBuffer(newBlock, aggregateBuffer, numSamplesToCopy - remainingSpace, remainingSpace, 0);
	}
	else
	{
		Helpers::copyAudioBlockIntoBuffer(newBlock, aggregateBuffer, numSamplesToCopy, 0, writeIndex);
	}

	writeIndex = (writeIndex + numSamplesToCopy) % numSamples;
}

template <typename SampleType>
void RingBufferAudio<SampleType>::readBlock(AudioBuffer<SampleType>& bufferToFill)
{
	const size_t difference = writeIndex - readIndex;
	const size_t remainingSpace = numSamples - readIndex;
	if (difference == 0) return;

	const size_t numSamplesToCopy = writeIndex > readIndex ? difference : remainingSpace + writeIndex;

	bufferToFill.setSize(aggregateBuffer.getNumChannels(), numSamplesToCopy, true, false, true);
	const dsp::AudioBlock<const SampleType> tempBlock = dsp::AudioBlock<SampleType>(aggregateBuffer);

	if(numSamplesToCopy <= remainingSpace)
	{
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, numSamplesToCopy, readIndex);
	}
	else
	{
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, remainingSpace, readIndex);
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, numSamplesToCopy - remainingSpace, 0, remainingSpace);
	}

	readIndex = (readIndex + numSamplesToCopy) % numSamples;
}

template <typename SampleType>
void RingBufferAudio<SampleType>::getPreviousSamples(AudioBuffer<SampleType>& bufferToFill)
{
	const size_t numSamples = bufferToFill.getNumSamples();
	bufferToFill.setSize(numChannels, numSamples);
	const dsp::AudioBlock<const SampleType> tempBlock = dsp::AudioBlock<SampleType>(aggregateBuffer);

	if (writeIndex > numSamples)
	{
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, numSamples, writeIndex - numSamples, 0);
	}
	else
	{
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, writeIndex, 0, 0);
		const size_t samplesLeft = numSamples - writeIndex;
		Helpers::copyAudioBlockIntoBuffer(tempBlock, bufferToFill, samplesLeft, tempBlock.getNumSamples() - samplesLeft , writeIndex);
	}
	
}

template <typename SampleType>
dsp::AudioBlock<SampleType> RingBufferAudio<SampleType>::getBlock()
{
	return dsp::AudioBlock<SampleType>(aggregateBuffer);
}

template <typename SampleType>
AudioBuffer<SampleType> RingBufferAudio<SampleType>::getBuffer()
{
	return aggregateBuffer;
}

class EnvelopeFollower: dsp::ProcessorBase
{

public:
    EnvelopeFollower();
    ~EnvelopeFollower();

    void process(const dsp::ProcessContextReplacing<float>&) override;
    void prepare(const dsp::ProcessSpec&) override;
    void reset() override;

    void setAttack(float milliseconds);
    void setRelease(float milliseconds);
    void setAmount(float newAmount);
    float getValue() const;
    std::function<void(float)> onValueCalculated;

private:
    double sampleRate;
    unsigned int numChannels;
    unsigned int maxBlockSize;

    float amount;
    float attackTime, releaseTime;

    dsp::ProcessorDuplicator<dsp::IIR::Filter<float>, dsp::IIR::Coefficients<float>> filter;
    AudioBuffer<float> copyBuffer;
	RingBufferAudio<float> audioBuffer;


};

EnvelopeFollower::EnvelopeFollower(): 
    sampleRate(44100),
	numChannels(2),
	maxBlockSize(0),
	amount(0),
    attackTime(0.5f),
	releaseTime(1.0f),
    audioBuffer(numChannels, sampleRate),
    copyBuffer(numChannels, 0)
{
}

EnvelopeFollower::~EnvelopeFollower() = default;

void EnvelopeFollower::prepare(const dsp::ProcessSpec& spec) 
{
    numChannels = spec.numChannels;
    sampleRate = spec.sampleRate;
    maxBlockSize = spec.maximumBlockSize;
    audioBuffer.resize(numChannels, spec.sampleRate);

    copyBuffer.setSize(numChannels, maxBlockSize);
    
    filter.state = dsp::IIR::Coefficients<float>::makeLowPass(sampleRate, 10);
    filter.prepare(spec);
}

void EnvelopeFollower::process(const dsp::ProcessContextReplacing<float>& context) 
{
	 // IIR Filter
     // get audio block copy
     // turn it into normal vector (sum channels)
     // filter it
     // put it in ring buffer

     const dsp::AudioBlock<const float>& block = context.getInputBlock();
     dsp::AudioBlock<float> copyBlock(copyBuffer);
     Helpers::copyAudioBlockIntoBuffer(block, copyBuffer, copyBuffer.getNumSamples());
     const dsp::ProcessContextReplacing<float> copyContext(copyBlock);
     filter.process(copyContext);
     auto max = copyContext.getOutputBlock().findMinAndMax().getEnd();
     onValueCalculated(max);

}

void EnvelopeFollower::reset()
{
    filter.reset();
}

void EnvelopeFollower::setAttack(const float milliseconds)
{
    attackTime = milliseconds;
}

void EnvelopeFollower::setRelease(const float milliseconds)
{
    releaseTime = milliseconds;
}

void EnvelopeFollower::setAmount(float newAmount)
{
    amount = newAmount;
}
