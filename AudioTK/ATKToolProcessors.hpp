#pragma once

#include "ATK.hpp"

namespace AudioTK
{

	
    struct MonoApplyGainFilter : public ATKFilter
    {
		using Filter = ATK::ApplyGainFilter<DspFloatType>;
		ATK::ApplyGainFilter<DspFloatType> * filter;
		MonoApplyGainFilter() : ATKFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoApplyGainFilter() {
			if(filter) delete filter;
		}
    };
    struct ApplyGainFilter : public ATKFilter
    {
		using Filter = ATK::ApplyGainFilter<DspFloatType>;
		ATK::ApplyGainFilter<DspFloatType> * filter;
		ApplyGainFilter() : ATKFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~ApplyGainFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoBufferFilter : public ATKFilter
    {
		using Filter = ATK::BufferFilter<DspFloatType>;
		Filter * filter;
		MonoBufferFilter() : ATKFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoBufferFilter() {
			if(filter) delete filter;
		}
    };
    struct BufferFilter : public ATKFilter
    {
		using Filter = ATK::BufferFilter<DspFloatType>;
		Filter * filter;
		BufferFilter(size_t nChannels=2) : ATKFilter() {
			filter = new Filter(nChannels);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~BufferFilter() {
			if(filter) delete filter;
		}
    };
    
    struct CachedCosinusGenerator : public ATKMonoFilter
    {
		using Filter = ATK::CachedCosinusGeneratorFilter<DspFloatType>;
		Filter * filter;
		CachedCosinusGenerator(int periods, int seconds=1) : ATKMonoFilter() {
			filter = new Filter(periods,seconds);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~CachedCosinusGenerator() {
			if(filter) delete filter;
		}
			
    };
    struct CachedSinusGenerator : public ATKMonoFilter
    {
		using Filter = ATK::CachedSinusGeneratorFilter<DspFloatType>;
		Filter * filter;
		CachedSinusGenerator(int periods, int seconds=1) : ATKMonoFilter() {
			filter = new Filter(periods,seconds);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~CachedSinusGenerator() {
			if(filter) delete filter;
		}
    };
    struct MonoDecimatorFilter : public ATKMonoFilter
    {
		using Filter = ATK::DecimationFilter<DspFloatType>;
		Filter * filter;
		
		MonoDecimatorFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoDecimatorFilter() {
			if(filter) delete filter;
		}
    };
    struct DecimatorFilter : public ATKFilter
    {
		using Filter = ATK::DecimationFilter<DspFloatType>;
		Filter * filter;
		
		DecimatorFilter(size_t channels=1) {
			filter = new Filter(channels);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~DecimatorFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoDerivativeFilter : public ATKMonoFilter
    {
		using Filter = ATK::DerivativeFilter<DspFloatType>;
		Filter * filter;
		MonoDerivativeFilter() {
			filter = new Filter();
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoDerivativeFilter() {
			if(filter) delete filter;
		}		
    };
    struct DerivativeFilter : public ATKFilter
    {
		using Filter = ATK::DerivativeFilter<DspFloatType>;
		Filter * filter;
		DerivativeFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~DerivativeFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoDryWetFilter : public ATKFilter
    {
		using Filter = ATK::DryWetFilter<DspFloatType>;
		Filter * filter;
		MonoDryWetFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoDryWetFilter() {
			if(filter) delete filter;
		}
    };
    struct DryWetFilter : public ATKFilter
    {
		using Filter = ATK::DryWetFilter<DspFloatType>;
		Filter * filter;
		DryWetFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~DryWetFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoMSFilter : public ATKMonoFilter
    {
		using Filter = ATK::DryWetFilter<DspFloatType>;
		Filter * filter;
		MonoMSFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoMSFilter() {
			if(filter) delete filter;
		}
    };
    struct MSFilter : public ATKFilter
    {
		using Filter = ATK::DryWetFilter<DspFloatType>;
		Filter * filter;
		MSFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MSFilter() {
			if(filter) delete filter;
		}
    };
    struct MaxFilter : public ATKFilter
    {
		using Filter = ATK::MaxFilter<DspFloatType>;
		Filter * filter;
		MaxFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MaxFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoMuteSoloBufferFilter : public ATKFilter
    {
		using Filter = ATK::MuteSoloBufferFilter<DspFloatType>;
		Filter * filter;
		MonoMuteSoloBufferFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoMuteSoloBufferFilter() {
			if(filter) delete filter;
		}
    };
    struct MuteSoloBufferFilter : public ATKFilter
    {
		using Filter = ATK::MuteSoloBufferFilter<DspFloatType>;
		Filter * filter;
		MuteSoloBufferFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MuteSoloBufferFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoOffsetVolumeFilter : public ATKFilter
    {
		using Filter = ATK::OffsetVolumeFilter<DspFloatType>;
		Filter * filter;
		MonoOffsetVolumeFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoOffsetVolumeFilter() {
			if(filter) delete filter;
		}
    };
    struct OffsetVolumeFilter : public ATKFilter
    {
		using Filter = ATK::OffsetVolumeFilter<DspFloatType>;
		Filter * filter;
		OffsetVolumeFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~OffsetVolumeFilter() {
			if(filter) delete filter;
		}
    };
    struct MonoOneMinusFilter : public ATKFilter
    {
		using Filter = ATK::OneMinusFilter<DspFloatType>;
		Filter * filter;
		MonoOneMinusFilter() {
			filter = new Filter(1);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~MonoOneMinusFilter() {
			if(filter) delete filter;
		}
    };
    struct OneMinusFilter : public ATKFilter
    {
		using Filter = ATK::OneMinusFilter<DspFloatType>;
		Filter * filter;
		OneMinusFilter() {
			filter = new Filter(2);
			assert(filter != nullptr);
			setFilter(filter);
		}
		~OneMinusFilter() {
			if(filter) delete filter;
		}
    };
    struct OversamplingFilter : public ATKFilter
    {
	
    };
    struct PanFilter : public ATKFilter
    {

    };    
    struct TanFilter : public ATKFilter
    {

    };
    struct VolumeFilter : public ATKFilter
    {

    };
    struct WhiteNoiseGenerator : public ATKFilter
    {

    };
    
}
