#include "Kfr/kfrcore.hpp"

using namespace DSP;
using namespace std;
using namespace kfr;

// normalize_frequency(f/sr)
// BiquadParams
// Biquad
// NotchFilter
// LowPassFilter
// HighPassFilter 
// LowShelfFilter
// HighShelfFilter
// BandPassFilter
// PeakFilter
// biquad(params)
// lowpassfilter
// highpassfilter
// notchfilter
// bandpassfilter
// peakfilter
// lowshelffilter
// highshelffilter

int main() 
{
    univector<float> wav = DSP::load_wav<float>("test.wav"),out;
    DSP::Biquad<float> low;
    low.lowpass(100.0/44100.0,0.9);
    low.apply(wav,out);
    DSP::write_wav<float>(out,"low.wav",1,(int)kfr::audio_sample_type::f32,44100);
}