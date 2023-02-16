#pragma once

namespace KfrDSP1
{
    template<typename T>
    struct FIRHighpassFilter
    {
        FIRFilter<T> * filter;

        FIRHighpassFilter(size_t num_taps, T x, kfr::expression_pointer<T> & window, bool normalize = true) {
            filter = new FIRFilter<T>(num_taps);
            assert(filter != NULL);
            filter->highpass(x,window,normalize);
        }
        ~FIRHighpassFilter() {
            if(filter) delete filter;
        }
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }
    };
}