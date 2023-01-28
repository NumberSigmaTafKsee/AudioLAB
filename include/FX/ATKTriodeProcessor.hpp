#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Preamp
{
    struct KorenTriode : public ATKFilter
    {
    protected:
        using Filter = ATK::TriodeFilter<DspFloatType,ATK::KorenTriodeFunction<DspFloatType>>;
        Filter filter;
    public:
        KorenTriode() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~KorenTriode() {
            
        }
    };
    struct EnhancedKorenTriode : public ATKFilter
    {
    protected:
        using Filter = ATK::TriodeFilter<DspFloatType,ATK::EnhancedKorenTriodeFunction<DspFloatType>>;
        Filter filter;
    public:
        EnhancedKorenTriode() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~EnhancedKorenTriode() {
            
        }
    };
    struct LeachTriode : public ATKFilter
    {
    protected:
        using Filter = ATK::TriodeFilter<DspFloatType,ATK::LeachTriodeFunction<DspFloatType>>;
        Filter filter;
    public:
        LeachTriode() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~LeachTriode() {
            
        }
    };    
    struct MunroPiazzaTriode : public ATKFilter
    {
    protected:
        using Filter = ATK::TriodeFilter<DspFloatType,ATK::MunroPiazzaTriodeFunction<DspFloatType>>;
        Filter filter;
    public:
        MunroPiazzaTriode() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~MunroPiazzaTriode() {
            
        }
    };
    struct DempWolfTriode : public ATKFilter
    {
    protected:
        using Filter = ATK::TriodeFilter<DspFloatType,ATK::DempwolfTriodeFunction<DspFloatType>>;
        Filter filter;
    public:
        DempWolfTriode() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~DempWolfTriode() {
            
        }
    };

}
