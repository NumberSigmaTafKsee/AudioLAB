#include "CombFilter.h"

CombFilter::CombFilter(float sampleRate)
{
	mParamRanges[gain][0] = -1.0f;
	mParamRanges[gain][1] = 1.0f;
	mParamRanges[delayInSec][0] = 0.0f;
	mParamRanges[delayInSec][1] = 10.0f;

	assert(sampleRate > 0);
	mSampleRate = sampleRate;
	mDelayLine.reset(new CRingBuffer<float>(CUtil::float2int<int>(mParamRanges[delayInSec][1] * mSampleRate)));
	mFilterType = FilterType_t::fir;
}

CombFilter::~CombFilter()
{
	mDelayLine.reset();
}

Error_t CombFilter::setParam(CombFilter::Param_t param, float value)
{
	if (!isInParamRange(param, value))
		return Error_t::kFunctionInvalidArgsError;

	if (param == CombFilter::Param_t::delayInSec) {
		mDelayLine->setReadIdx(mDelayLine->getWriteIdx() - CUtil::float2int<int>(value * mSampleRate));
	}
	mParamValues[param] = value;
	return Error_t::kNoError;
}

float CombFilter::getParam(CombFilter::Param_t param) const
{
	return mParamValues[param];
}

Error_t CombFilter::setFilterType(CombFilter::FilterType_t filterType)
{
	if (filterType == numFilterTypes)
		return Error_t::kFunctionInvalidArgsError;

	mFilterType = filterType;
	return Error_t::kNoError;
}

CombFilter::FilterType_t CombFilter::getFiltertType() const
{
	return mFilterType;
}

Error_t CombFilter::process(const float* inputBuffer, float* outputBuffer, int numSamples)
{
	if (!inputBuffer || !outputBuffer)
		return Error_t::kMemError;
	if (numSamples < 0)
		return Error_t::kFunctionInvalidArgsError;

	switch (mFilterType) {
	case fir:
		for (int i = 0; i < numSamples; i++) {
			mDelayLine->putPostInc(inputBuffer[i]);
			outputBuffer[i] = inputBuffer[i] + mParamValues[CombFilter::Param_t::gain] * mDelayLine->getPostInc();
		}
		break;
	case iir:
		for (int i = 0; i < numSamples; i++) {
			outputBuffer[i] = inputBuffer[i] + mParamValues[CombFilter::Param_t::gain] * mDelayLine->getPostInc();
			mDelayLine->putPostInc(outputBuffer[i]);
		}
		break;
	default:
		return Error_t::kUnknownError;
	}
	return Error_t::kNoError;
}

bool CombFilter::isInParamRange(CombFilter::Param_t param, float value) const
{
	return (mParamRanges[param][0] <= value && value <= mParamRanges[param][1]);
}