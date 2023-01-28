#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Preamp
{
    struct KorenTriode2 : public ATKFilter
    {
    protected:
        using Filter = ATK::Triode2Filter<DspFloatType,ATK::KorenTriodeFunction<DspFloatType>>;
        Filter filter;
    public:

        KorenTriode2() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~KorenTriode2() {
            
        }
    };
    struct EnhancedKorenTriode2 : public ATKFilter
    {
    protected:
        using Filter = ATK::Triode2Filter<DspFloatType,ATK::EnhancedKorenTriodeFunction<DspFloatType>>;
        Filter filter;
    public:

        EnhancedKorenTriode2() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~EnhancedKorenTriode2() {
            
        }
    };
    struct LeachTriode2 : public ATKFilter
    {
    protected:
        using Filter = ATK::Triode2Filter<DspFloatType,ATK::LeachTriodeFunction<DspFloatType>>;
        Filter filter;
    public:

        LeachTriode2() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~LeachTriode2() {
            
        }
    };    
    struct MunroPiazzaTriode2 : public ATKFilter
    {
    protected:
        using Filter = ATK::Triode2Filter<DspFloatType,ATK::MunroPiazzaTriodeFunction<DspFloatType>>;
        Filter filter;

    public:
        MunroPiazzaTriode2() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~MunroPiazzaTriode2() {
            
        }
    };
    struct DempWolfTriode2 : public ATKFilter
    {
    protected:
        using Filter = ATK::Triode2Filter<DspFloatType,ATK::DempwolfTriodeFunction<DspFloatType>>;
        Filter filter;

    public:
        DempWolfTriode2() : ATKFilter(),
        filter(Filter::build_standard_filter())
        {            
        }
        ~DempWolfTriode2() {
            
        }
    };

}   
 