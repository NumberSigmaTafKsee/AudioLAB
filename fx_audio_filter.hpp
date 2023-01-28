#include <cstdio>
#include <cstring>
#include <sndfile.hh>
#include <string>
#include <vector>
#include <fftw3.h>
#include <istream>
#include <tuple>
#include <random>
#include <algorithm>
#include <stdexcept>

//Supports only one channel
class IO{
public:
    IO(const std::string &input_path, const std::string &output_path);
    
    //Delete copying
    IO(const IO &io) = delete;
    IO &operator=(const IO &io) = delete;

    //Move
    IO(IO &&io);
    IO &operator=(IO &&io);

    ~IO();
    sf_count_t read_n(double *dest, size_t n) const;
    sf_count_t write_n(double *data, size_t n) const;
    SF_INFO get_info() const{return info;}
private:
    SNDFILE *i_file_, *o_file_;
    SF_INFO info;
};

///Writes to file with overlap
class WriteOverlap{
public:
    ///Default constructor
    /*
     * @param size the size of the window to be overlaped
     * @param max_count the maximum count of datapoints to be written at once
     */
    WriteOverlap(size_t size);
    ///Adds the data to cache and then shifts the cache by shift datapoints. The shifted data is written to a file via io
    void operator()(size_t shift, double *data, const IO &io);
private:
    size_t size_;
    std::vector<double> cache_;
};

///Basic implementation of the Hann window
/*
 * Supposed to be used for reading
 */
class Window{
public:
    Window(size_t size);
    ///Reads data and returns the windowed result, includint the 50% overlap
    /*
     * Returns the number of samples that can be written to file after processing.
     * Output buffer for write should be used in conjuncion.
     */
    size_t operator()(double *output, const IO &io);
    ///Writes the buffered data to file
private:
    size_t size_;
    bool next_end_ = false, first_ = true;
    std::vector<double> older_data_, newer_data_;
    double hann_val(size_t pos) const;
};


/// Basic wrapper for FFT operations form real data
/*
 * Holds two buffers of the needed size. To transform data, copy in into the right buffer
 * and call the corresponding function. Beware that c2r changes input data.
 * 
 * The wrapper remembers the current "state" of the data,
 * so consecutive calls to the same transform transforms the data only once.
 * This also ensures correct rescaling of the data.
 * Rescaling occurs only on c2r transform.
 */
class FFTwrapper{
public:
    ///Default constructor that wraps transforms of size n
    /*
     * Expects the first transform to be from real to complex
     */
    FFTwrapper(int n);
    ///Constructor for both directions of first transform
    /*
     * The is_complex_state parameter indicates wheter the initial data will be complex or real.
     * Beware that the data are always rescaled after the c2r transform, so starting from complex
     * may introduce some scaling error
     */
    FFTwrapper(int n, bool is_complex_state);
    
    //Delete copying
    FFTwrapper(const FFTwrapper &fft) = delete;
    FFTwrapper &operator=(const FFTwrapper &fft) = delete;

    //Move constructors and desctructors
    FFTwrapper(FFTwrapper &&fft);
    FFTwrapper &operator=(FFTwrapper &&fft);

    ~FFTwrapper();
    ///Returns the pointer to the real time-donain data
    double *get_real_ptr();
    ///Returns the length of the real buffer
    size_t get_real_c() const;
    ///Returns the pointer to the complex frequency-domain data
    fftw_complex *get_complex_ptr();
    ///Returns the length of the complex buffer
    size_t get_complex_c() const;
    ///Transforms data from real buffer into complex buffer
    void r2c_transform();
    ///Transforms data from complex buffer into real buffer
    void c2r_transform();

    //Can be changed manually to avoid unnecesarry transforms
    bool is_complex_state_ = false;

private:
    int n_r_, n_c_;
    double *real_;
    fftw_complex *complex_;
    fftw_plan r2c_, c2r_;
};

/*
 *  Frequency is normalized to be independent of sample rate, calculating with sample rate of 44 100 Hz
 *      Frequency of 1 means that there is one period in 44 100 samples
 *      Frequency n means there are n periods in 44 100 samples
 */

class PitchShift{
public:
    ///Ammount must be >1;
    PitchShift(double ammount):ammount_(ammount){};
    void operator()(FFTwrapper &fft);
    ///Pitch shifts the data into the output FFTwrapper complex buffer, overwriting any data stored
    static void shift_complex(const std::vector<std::tuple<double,double>> &in, FFTwrapper &out_fft, double ammount);
private:
    double ammount_;
};

PitchShift make_pitchshift(std::istream &is, const IO &io);


class Lowpass{
public:
    ///Cutoff = what percenatge is cut off
    //Cutoff must be between 0 and 1;
    Lowpass(double cutoff): cutoff_(cutoff){};
    void operator()(FFTwrapper &fft);
private:
    double cutoff_;
};

Lowpass make_lowpass(std::istream &is, const IO &io);

class Highpass{
public:
    ///Cutoff = what percenatge is cut off
    //Cutoff must be between 0 and 1;
    Highpass(double cutoff): cutoff_(cutoff){};
    void operator()(FFTwrapper &fft);
private:
    double cutoff_;
};

Highpass make_highpass(std::istream &is, const IO &io);

///Simple delay of samples, accounting for 50 % overlap
class Delay{
public:
    Delay(
        const std::vector<size_t> &delays,   //delay in samples per one delay "instance"
        const std::vector<double> &gains     //gains of the individual delay "instances"
    );
    void operator()(FFTwrapper &fft);
    //Input and output can be the same pointer
    void operator()(double *input, double* output, size_t size);
    static Delay create_linear(size_t delay, size_t count);
    static Delay create_exp(size_t delay, size_t count, double base);
private:
    std::vector<double> older_data_, newer_data_;
    std::vector<size_t> delays_;
    std::vector<double> gains_;
};

Delay make_delay(std::istream &is, const IO &io);

class Chorus{
public:
    Chorus(
        const std::vector<double> &voice_LFO_pos,
        const std::vector<double> &voice_LFO_delta,
        const std::vector<double> &voice_LFO_amp,
        const std::vector<double> &voice_pitch,
        const std::vector<Delay> &delays
    ):voice_LFO_pos_(voice_LFO_pos),
      voice_LFO_delta_(voice_LFO_delta),
      voice_LFO_amp_(voice_LFO_amp),
      voice_pitch_(voice_pitch),
      delays_(delays)
    {}
    static Chorus create_random(
        size_t voice_count,             //Number of voices to add
        size_t max_delay,               //The maximum delay of the voices
        double LFO_freq_base,           //Base frequency of the pitch LFO
        double LFO_freq_variance,       //The range of the pitch LFO frequencies [base-variance, base+variance]
        double LFO_amplitude_base,      //Base amp for pitch LFO
        double LFO_amplitude_variance,  //amp variance for pitch LFO
        double pitch_variance,          //the variance of the midpoint of the pitch LFO
        double gain                     //The gain of the added voices 1...no volume difference
    );
    static Chorus create_uniform(
        size_t voice_count,             //Number of voices to add
        size_t max_delay,               //The maximum delay of the voices
        double LFO_freq_base,           //Base frequency of the pitch LFO
        double LFO_freq_variance,       //The range of the pitch LFO frequencies [base-variance, base+variance]
        double LFO_amplitude_base,      //Base amp for pitch LFO
        double LFO_amplitude_variance,  //amp variance for pitch LFO
        double pitch_variance,          //the variance of the midpoint of the pitch LFO
        double gain                     //The gain of the added voices 1...no volume difference
    );
    void operator()(FFTwrapper &fft);
private:
    std::vector<double> voice_LFO_pos_;
    std::vector<double> voice_LFO_delta_;
    std::vector<double> voice_LFO_amp_;
    std::vector<double> voice_pitch_;
    std::vector<Delay> delays_;
};

Chorus make_chorus(std::istream &is, const IO &io);


void Lowpass::operator()(FFTwrapper &fft){
    fft.r2c_transform();
    for(size_t i = (size_t)(fft.get_complex_c()*(1-cutoff_)); i<fft.get_complex_c(); ++i){
        fft.get_complex_ptr()[i][0]=0;
        fft.get_complex_ptr()[i][1]=0;
    }
}

void Highpass::operator()(FFTwrapper &fft){
    fft.r2c_transform();
    for(size_t i = 0; i<((size_t)fft.get_complex_c()*cutoff_); ++i){
        fft.get_complex_ptr()[i][0]=0;
        fft.get_complex_ptr()[i][1]=0;
    }
}

const double PI  = 3.141592653589793238463;
const size_t DEFAULT_SAMPLE_RATE = 44100;

Chorus Chorus::create_random(
    size_t voice_count,             //Number of voices to add
    size_t max_delay,               //The maximum delay of the voices
    double LFO_freq_base,           //Base frequency of the pitch LFO
    double LFO_freq_variance,       //The range of the pitch LFO frequencies [base-variance, base+variance]
    double LFO_amp_base,            //Base amplitude for pitch LFO
    double LFO_amp_variance,        //amplitude variance for pitch LFO
    double pitch_variance,          //the variance of the midpoint of the pitch LFO
    double gain                     //The gain of the added voices 1...no volume difference
){
    //Random number distributions
    std::uniform_real_distribution<double> freq_unif(LFO_freq_base-LFO_freq_variance, LFO_freq_base+LFO_freq_variance);
    std::uniform_real_distribution<double> amp_unif(LFO_amp_base-LFO_amp_variance, LFO_amp_base+LFO_amp_variance);
    std::uniform_int_distribution<size_t> delay_unif(0, max_delay);
    std::uniform_real_distribution<double> pos_unif(0, 2 * PI);
    std::uniform_real_distribution<double> ps_unif(1-pitch_variance, 1+pitch_variance);
    std::random_device rd;
    std::default_random_engine re(rd());

    std::vector<double> voice_LFO_pos;
    std::vector<double> voice_LFO_delta;
    std::vector<double> voice_LFO_amp;
    std::vector<double> voice_pitch;
    std::vector<Delay> delays;

    for(size_t i=0; i<voice_count; ++i){
        voice_LFO_delta.push_back((2*PI*freq_unif(re))/DEFAULT_SAMPLE_RATE);
        voice_LFO_amp.push_back(amp_unif(re));
        std::vector<size_t> d={delay_unif(re)};
        std::vector<double> g={gain};
        delays.emplace_back(d,g);
        voice_pitch.push_back(ps_unif(re));
        voice_LFO_pos.push_back(pos_unif(re));
    }

    Chorus ret(voice_LFO_pos,voice_LFO_delta,voice_LFO_amp,voice_pitch,delays);
    return ret;
}
//Uniform = liniar, because its the easiest to implement
Chorus Chorus::create_uniform(
    size_t voice_count,             //Number of voices to add
    size_t max_delay,               //The maximum delay of the voices
    double LFO_freq_base,           //Base frequency of the pitch LFO
    double LFO_freq_variance,       //The range of the pitch LFO frequencies [base-variance, base+variance]
    double LFO_amp_base,            //Base amplitude for pitch LFO
    double LFO_amp_variance,        //amplitude variance for pitch LFO
    double pitch_variance,          //the variance of the midpoint of the pitch LFO
    double gain                     //The gain of the added voices 1...no volume difference
){
    std::vector<double> voice_LFO_pos;
    std::vector<double> voice_LFO_delta;
    std::vector<double> voice_LFO_amp;
    std::vector<double> voice_pitch;
    std::vector<Delay> delays;

    for(size_t i=0; i<voice_count; ++i){
        //Lower frequency ~ higher amp
        double freq = (LFO_freq_base - LFO_freq_variance) + ((2 * LFO_freq_variance)/(voice_count-1))*i;
        voice_LFO_delta.push_back((2*PI*freq)/DEFAULT_SAMPLE_RATE);
        voice_LFO_amp.push_back((LFO_amp_base + LFO_amp_variance) - ((2 * LFO_amp_variance)/(voice_count-1))*i);
        voice_pitch.push_back((1-pitch_variance)+i*(2*pitch_variance)/(voice_count-1));
        std::vector<size_t> d={(max_delay/(voice_count-1))*i};
        std::vector<double> g={gain};
        delays.emplace_back(d,g);
    }
    voice_LFO_pos.resize(voice_count,0);
    Chorus ret(voice_LFO_pos,voice_LFO_delta,voice_LFO_amp,voice_pitch,delays);
    return ret;
}
void Chorus::operator()(FFTwrapper &fft){
    fft.c2r_transform();

    //Copy samples to save them, so we dont overwrite them
    std::vector<double> res(fft.get_real_ptr(),fft.get_real_ptr()+fft.get_real_c());

    fft.r2c_transform();
    
    //Copy spectrum to change it multiple times
    std::vector<std::tuple<double,double>> spectrum;
    for(size_t j=0; j<fft.get_complex_c(); ++j){
        spectrum.push_back(std::make_tuple(fft.get_complex_ptr()[j][0],fft.get_complex_ptr()[j][1]));
    }

    for(size_t i=0; i<voice_LFO_pos_.size();++i){
        double pitch = voice_pitch_[i] * std::exp2(voice_LFO_amp_[i]*std::sin(voice_LFO_pos_[i]));
        PitchShift::shift_complex(spectrum, fft, pitch);
        voice_LFO_pos_[i] += voice_LFO_delta_[i] * fft.get_real_c();
        //Reseting the pos to stay in range [0,2 pi], so we dont lose precision
        if(voice_LFO_pos_[i] > 2*PI)
            voice_LFO_pos_[i] -= 2*PI;
        
        fft.c2r_transform();        

        delays_[i](fft.get_real_ptr(), res.data(), fft.get_real_c());    

    }
    fft.is_complex_state_ = false;
    //overwrite with created data
    for(size_t i=0; i<res.size(); ++i){
        fft.get_real_ptr()[i] = res[i];
    }
}

Delay::Delay(
        const std::vector<size_t> &delays,
        const std::vector<double> &gains
):delays_(delays), gains_(gains){
    size_t max_delay = 0;
    for(auto &&i: delays){
        max_delay = std::max(i,max_delay);
    }
    older_data_.resize(max_delay,0);
    newer_data_.resize(max_delay,0);
}

void Delay::operator()(FFTwrapper &fft){
    fft.c2r_transform();
    operator()(fft.get_real_ptr(), fft.get_real_ptr(), fft.get_real_c());
}

//Only adds the delayed samples to output, does not copy original signal
//Input and output can be the same
void Delay::operator()(double *input, double* output, size_t size){
    //Input and output can be the same pointer, for simplicity, copy the samples
    std::vector<double> input_copy(input, input+size);

    //Add the delay from the second previous window, and zero older_data_
    size_t min = std::min(older_data_.size(), size);
    for(size_t i=0; i < min; ++i){
        output[i] += older_data_[i];
        older_data_[i] = 0;
    }
    if(size < older_data_.size())
        std::rotate(older_data_.begin(), older_data_.begin()+size, older_data_.end());
    
    std::swap(older_data_, newer_data_);

    size_t delaysSize = delays_.size();
    for(size_t i=0; i<size; ++i){
        for(size_t j=0; j< delaysSize; ++j){
            if(i<delays_[j])
                continue;
            output[i] += input_copy[i-delays_[j]] * gains_[j]; 
        }
    }
    
    size_t n_s = newer_data_.size();
    //Save the delayed samples that would lie out of the window
    for(size_t i=size; i<size + n_s; ++i){
        for(size_t j=0; j<delaysSize; ++j){
            if(i-delays_[j]>=size)
                continue;
            newer_data_[i-size] += input_copy[i-delays_[j]] * gains_[j];
        }
    }
}

Delay Delay::create_linear(size_t delay, size_t count){
    std::vector<size_t> delays;
    std::vector<double> gains;
    
    for(size_t i=0; i<count; ++i){
        delays.push_back((i+1)*delay);
    }
    for(size_t i=0; i<count; ++i){
        gains.push_back(1-((i+1.0)/(count+1.0)));
    }
    return Delay(delays, gains);
}
Delay Delay::create_exp(size_t delay, size_t count, double base){
    std::vector<size_t> delays;
    std::vector<double> gains;
    double gain;
    gain = base;
    for(size_t i=0; i<count; ++i){
        delays.push_back((i+1)*delay);
    }
    for(size_t i=0; i<count; ++i){
        gains.push_back(gain);
        gain *= base;
    }
    return Delay(delays, gains);
}

FFTwrapper::FFTwrapper(int n){
    n_r_ = n;
    n_c_ = n/2+1;
    real_ = (double*)fftw_malloc(sizeof(double)*n_r_);
    complex_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n_c_);
    r2c_ = fftw_plan_dft_r2c_1d(n_r_, real_, complex_, FFTW_ESTIMATE);
    c2r_ = fftw_plan_dft_c2r_1d(n_r_, complex_, real_, FFTW_ESTIMATE);
}
FFTwrapper::FFTwrapper(int n, bool is_complex_state): FFTwrapper(n){
    is_complex_state_ = is_complex_state;
}
FFTwrapper::FFTwrapper(FFTwrapper &&fft){
    n_r_ = fft.n_r_;
    n_c_ = fft.n_c_;
    real_ = fft.real_;
    complex_ = fft.complex_;
    r2c_ = fft.r2c_;
    c2r_ = fft.c2r_;

    fft.real_ = nullptr;
    fft.complex_ = nullptr;
    fft.r2c_ = nullptr;
    fft.c2r_ = nullptr;
}
FFTwrapper &FFTwrapper::operator=(FFTwrapper &&fft){
    n_r_ = fft.n_r_;
    n_c_ = fft.n_c_;
    std::swap(real_, fft.real_);
    std::swap(complex_, fft.complex_);
    std::swap(r2c_, fft.r2c_);
    std::swap(c2r_, fft.c2r_);
    return *this;
}

FFTwrapper::~FFTwrapper(){
    //fftw_plan is just a typedefed pointer
    if(r2c_)
        fftw_destroy_plan(r2c_);
    if(c2r_)
        fftw_destroy_plan(c2r_);
    if(real_)
        fftw_free(real_);
    if(complex_)
        fftw_free(complex_);
}
double *FFTwrapper::get_real_ptr(){
    return real_;
}
size_t FFTwrapper::get_real_c() const{
    return n_r_;
}
fftw_complex *FFTwrapper::get_complex_ptr(){
    return complex_;
}
size_t FFTwrapper::get_complex_c() const{
    return n_c_;
}
void FFTwrapper::r2c_transform(){
    if(!is_complex_state_){
        fftw_execute(r2c_);
        is_complex_state_ = true;
    }
}
void FFTwrapper::c2r_transform(){
    if(is_complex_state_){
        fftw_execute(c2r_);
        is_complex_state_ = false;
        //Rescaling
        for(size_t i=0; i < (size_t)n_r_; ++i)
            real_[i] /= n_r_;
    }
}
PitchShift make_pitchshift(std::istream &is, const IO &io){
    double ammount;
    is>>ammount;
    PitchShift ps(ammount);
    return ps;
}
Lowpass make_lowpass(std::istream &is, const IO &io){
    double cutoff;
    is>>cutoff;
    Lowpass lp(cutoff);
    return lp;
}
Highpass make_highpass(std::istream &is, const IO &io){
    double cutoff;
    is>>cutoff;
    Highpass hp(cutoff);
    return hp;
}

Delay make_delay(std::istream &is, const IO &io){
    std::string mode;
    is>>mode;
    std::vector<size_t> delays;
    std::vector<double> gains;
    if(mode=="auto"){
        std::string type;
        size_t delay, count;
        is>>type>>delay>>count;
        //ms -> samples
        delay *= io.get_info().samplerate/1000;
        if(type=="linear"){
            return Delay::create_linear(delay, count);
        }
        else if(type=="exp"){
            double base;
            is>>base;
            return Delay::create_exp(delay, count, base);
        }
        else{
            throw std::invalid_argument("Bad delay type: "+type);
        }
    }
    else if(mode=="manual"){
        size_t delay;
        double gain;
        while(is){
            is>>delay>>gain;
            //ms -> samples
            delays.push_back(delay*io.get_info().samplerate/1000);
            gains.push_back(gain);
        }
        return Delay(delays, gains);
    }
    throw std::invalid_argument("Bad delay mode: "+mode);
}

Chorus make_chorus(std::istream &is, const IO &io){
    size_t voice_count, max_delay;
    double LFO_freq_base, LFO_freq_variance, LFO_amplitude_base,
           LFO_amplitude_variance, pitch_variance, gain;
    std::string mode;
    is>>mode>>voice_count>>max_delay
      >>LFO_freq_base>>LFO_freq_variance
      >>LFO_amplitude_base>>LFO_amplitude_variance
      >>pitch_variance>>gain;
      
    //ms -> samples
    max_delay *= io.get_info().samplerate/1000;

    //from Hz to normalized freq
    LFO_freq_base *= io.get_info().samplerate/44100;

    //from Hz to normalized freq
    LFO_freq_variance *= io.get_info().samplerate/44100;

    if(mode=="random"){
        return Chorus::create_random(
            voice_count,
            max_delay,
            LFO_freq_base,
            LFO_freq_variance,
            LFO_amplitude_base,
            LFO_amplitude_variance,
            pitch_variance,
            gain
        );
    }
    else if(mode=="uniform"){
        return Chorus::create_uniform(
            voice_count,
            max_delay,
            LFO_freq_base,
            LFO_freq_variance,
            LFO_amplitude_base,
            LFO_amplitude_variance,
            pitch_variance,
            gain
        );
    }
    else{
        throw std::invalid_argument("Bad Chorus mode: "+mode);
    }
}

IO::IO(const std::string &input_path, const std::string &output_path){
    info.format = 0;
    i_file_ = sf_open(input_path.c_str(), SFM_READ, &info);
    o_file_ = sf_open(output_path.c_str(), SFM_WRITE, &info);
}

IO::IO(IO &&io){
    i_file_ = io.i_file_;
    o_file_ = io.o_file_;
    info = io.info;
    io.i_file_ = nullptr;
    io.o_file_ = nullptr;
}
IO &IO::operator=(IO &&io){
    info = io.info;
    std::swap(i_file_, io.i_file_);
    std::swap(o_file_, io.o_file_);
    return *this;
}

IO::~IO(){
    if(i_file_)
        sf_close(i_file_);
    if(o_file_)
        sf_close(o_file_);
}
sf_count_t IO::read_n(double *dest, size_t n) const{
   sf_count_t count = sf_read_double(i_file_, dest, n);
   //Zero the data that the read didnt overwrite
   for(size_t i = count; i<n; ++i)
      dest[i] = 0;
   return count;
}
sf_count_t IO::write_n(double *data, size_t n) const{
    return sf_write_double(o_file_, data, n);
}

WriteOverlap::WriteOverlap(size_t size): size_(size){
    cache_.resize(size_,0);
}
//data buffer must be at least size_ long
void WriteOverlap::operator()(size_t shift, double *data, const IO &io){
    for(size_t i=0; i<size_; ++i){
        cache_[i] += data[i];
    }
    io.write_n(cache_.data(),shift);

    std::rotate(cache_.begin(),cache_.begin()+shift, cache_.end());
    for(size_t i=cache_.size()-shift; i<cache_.size(); ++i)
        cache_[i] = 0;
}

//in.size() has to be equal to out_fft.get_complex_c()
void PitchShift::shift_complex(const std::vector<std::tuple<double,double>> &in, FFTwrapper &out_fft, double ammount){
    out_fft.is_complex_state_ = true;
    //Zero out the complex part
    for(size_t i=0; i<out_fft.get_complex_c(); ++i){
        out_fft.get_complex_ptr()[i][0]=0;
        out_fft.get_complex_ptr()[i][1]=0;
    }
    //Do the pitch shift
    for(size_t i=0; i<out_fft.get_complex_c(); ++i){
        double x = i*ammount;
        size_t new_id = (size_t)std::floor(x+0.5);

        if(x<0 || new_id >= out_fft.get_complex_c())
            break;

        out_fft.get_complex_ptr()[new_id][0] += std::get<0>(in[i]);
        out_fft.get_complex_ptr()[new_id][1] += std::get<1>(in[i]);
    }
}

void PitchShift::operator()(FFTwrapper &fft){
    std::vector<std::tuple<double,double>> spectrum;
    fft.r2c_transform();
    for(size_t j=0; j<fft.get_complex_c(); ++j){
        spectrum.push_back(std::make_tuple(fft.get_complex_ptr()[j][0],fft.get_complex_ptr()[j][1]));
    }
    shift_complex(spectrum, fft, ammount_);
}

const double PI = 3.141592653589793238463;

Window::Window(size_t size): size_(size){
    older_data_.resize(size / 2, 0);
    newer_data_.resize(size / 2 + size % 2, 0);
}
size_t Window::operator()(double *output, const IO &io){  
    if(next_end_)
        return 0;
    //Move data from new to old (and from old to new, but we are rewriting those)
    std::swap(older_data_, newer_data_);

    if(first_){
        //Initialize data to old on first read
        io.read_n(older_data_.data(), older_data_.size());
        first_ = false;
    }
    //Rewrite new
    size_t read_n = io.read_n(newer_data_.data(), newer_data_.size());
    
    for(size_t i=0; i<size_; ++i){
        if(i<older_data_.size()){
            output[i] = hann_val(i) * older_data_[i];
        }
        else{
            output[i] = hann_val(i) * newer_data_[i-older_data_.size()];
        }
    }
    
    if(read_n < newer_data_.size()){
        next_end_ = true;
        //Signalize to write everything to file
        return older_data_.size() + read_n;
    }
    return older_data_.size();
}
//Hann window with a bit of +-1 fiddling so we dont multiply the border values by 0
double Window::hann_val(size_t pos) const{
    double val = std::cos((PI/(size_+1))*(pos-((double)(size_-1)/2)));
    return val*val;
}