#pragma once

namespace KfrDSP1
{
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Chebyshev2
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Chebyshev2Filter
    {
        std::vector<Biquad12DB*> filters;
        
        double fc,fs,W;
        double low,high;
        int order;
        
        FilterType filterType;

        Chebyshev2Filter(FilterType type, int Order, double Fs, double Fc, float w=80)
        {
            order = Order;
            fs    = Fs;
            fc    = Fc;
            filterType = type;
            low   = 0;
            high  = Fs/2;
            W     = w;
            initFilter();
            
        }
        void initFilter()
        {
            switch(filterType)
            {
            case Lowpass:  lowpass(fc,fs); break;
            case Highpass: highpass(fc,fs); break;
            case Bandpass: bandpass(low,high,fs); break;
            case Bandstop: bandstop(low,high,fs); break;
            }

        }
        void setCutoff(double f) {
            fc = f;
            low = f;
            switch(filterType)
            {
                case Lowpass:  doLowpassCutoff(fc,fs); break;
                case Highpass: doHighpassCutoff(fc,fc); break;
                case Bandpass: doBandpassCutoff(low,high,fs); break;
                case Bandstop: doBandstopCutoff(low,high,fs); break;
            }
        }
        void setCutoff(double lo,double hi) {
            low = lo;
            high = hi;
            switch(filterType)
            {
                case Bandpass: doBandpassCutoff(lo,hi,fc); break;
                case Bandstop: doBandstopCutoff(lo,hi,fc); break;
            }
        }

        void doLowpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void lowpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doHighpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void highpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandpassCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandpass(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandstopCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandstop(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            double R = I;
            for(typename std::vector<Biquad12DB*>::reverse_iterator i = filters.rbegin(); i != filters.rend(); i++)
            {
                R = (*i)->Tick(R);
            }
            return R;
        }
    };
}