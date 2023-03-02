#pragma once
#include <cmath>

namespace kfr3
{
    ///////////////////////////////////////////////////////////////
    // Utils
    ///////////////////////////////////////////////////////////////
    inline double freq_to_midi(double f) {
        return 12.0*std::log2(f/440.0) + 69;
    }

    inline double midi_to_freq(double m) {
        return std::pow(2.0, (m-69)/12)*440.0;
    }
    inline double cv_to_freq(double cv)
    {
        return std::pow(2,cv);
    }
    inline double semitone(int semi, double f)
    {
        double m = freq_to_midi(f);
        return midi_to_freq(m + semi);
    }
    inline double octave(int octave, double f) {
        double m = freq_to_midi(f);
        return midi_to_freq(m + octave*12);
    }

    ///////////////////////////////////////////////////////////////
    // DC Block
    ///////////////////////////////////////////////////////////////

    template<typename T>
    kfr::univector<T> dcremove(T cutoff, kfr::univector<T> in)
    {        
        return DSP::dcremove(in,cutoff);     
    }

    ///////////////////////////////////////////////////////////////
    // File I/O
    ///////////////////////////////////////////////////////////////

    template<typename T>
    kfr::univector<T> load_wav(const char * filename)
    {
        return DSP::load_wav<T>(filename);        
    }
    template<typename T>
    kfr::univector<T> load_mp3(const char * filename)
    {
        return DSP::load_mp3<T>(filename);        
    }
    template<typename T>
    kfr::univector<T> load_flac(const char * filename)
    {
        return DSP::load_flac<T>(filename);        
    }
    template<typename T>
    void save_wav(kfr::univector<T>  in, const char * filename, size_t channels, int sample_type, double sample_rate, bool use_w64=false)
    {        
        DSP::write_wav(in,filename,channels,sample_type,sample_rate,use_w64);
        
    }

    template<typename T> using WAVFileReader = DSP::WavReader<T>;
    template<typename T> using WAVFileWriter = DSP::WavWriter<T>;
    template<typename T> using MP3FileReader = DSP::MP3Reader<T>;
    template<typename T> using FLACFileReader = DSP::FlacReader<T>;

    ///////////////////////////////////////////////////////////////
    // Plot
    ///////////////////////////////////////////////////////////////
    template<typename T>
    void plot(kfr::univector<T> in, std::string& name = "", std::string& options = "") {
        kfr::plot_show(name,in,options);
    }
}
