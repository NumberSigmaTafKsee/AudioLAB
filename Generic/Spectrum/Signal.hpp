#pragma once
#include "samples/sample.hpp"
#include "samples/sample_dsp.hpp"
#include "FunctionGenerator.hpp"
#include "SndFile.hpp"

// Signal   = sample_vector<float>
// SignalD  = sample_vector<double>
// SignalC  = complex_vector<float>
// SignalZ  = complex_vector<double>

// Signal2d = Vector2d
// Signal4d = Vector4d
// Signal8d = Vector8d
// Signal4f = Vector4f
// Signal8f = Vector8f 
// Signal6f = Vector16f 

// SignalArray  = mkl::array<float>
// SignalArrayD = mkl::array<double>
// SignalArrayC = mkl::complex_array<float>
// SignalArrayZ = mkl::complex_array<double>

// SignalVectorF = mkl::vector<float>
// SignalVectorD = mkl::vector<double>
// SignalVectorC = mkl::complex_vector<float>
// SignalVectorZ = mkl::complex_vector<double>

// SignalMatrix   = mkl::matrix<float>
// SignalMatrixD  = mkl::matrix<double>
// SignalMatrixC  = mkl::complex_matrix<float>
// SignalMatrixZ  = mkl::complex_matrix<double>

// ipp
// FFT
// IIR
// FIR
// Convolution
// XCorr
// ACorr
// Polyphase Resampler


struct Signal : public sample_vector<float>
{    
    float sample_rate=44100.0f;
    size_t channels=1;

    Signal() = default;
    Signal(size_t n, size_t channels=1) { resize(n); this->channels = channels; }
    Signal(const char * filename) { loadFile(filename); }
    Signal(const Signal& s) { *this = s; sample_rate = s.sample_rate; channels = s.channels; }

    Signal& operator = (const Signal & s) {
        *this = s; 
        sample_rate = s.sample_rate; 
        channels = s.channels; 
        return *this;
    }
    void loadFile(const char * filename) {
        SndFileReaderFloat r(filename);
        resize(r.size());
        channels = r.channels();
        sample_rate = r.samplerate();
        r >> *this;
    }
    void saveWavFile(const char * filename) {
        SndFileWriterFloat w(filename+".wav",0x10006,channels,sample_rate);
    }

    size_t numChannels() const { return channels; }
    float  sampleRate()  const { return sample_rate; }

    Signal getLeftChannel() const {
        Signal s(size()/2);
        if(channels != 2) s = *this;
        else {
            size_t n = 0;
            for(size_t i = 0; i < size(); i+=2) {
                s[n++] = (*this)[i];
            }
        }
        return s;
    }
    Signal getRightChannel() const {
        Signal s(size()/2);
        if(channels != 2) s = *this;
        else {
            size_t n = 0;
            for(size_t i = 1; i < size(); i+=2) {
                s[n++] = (*this)[i];
            }
        }
        return s;
    }
    void setLeftChannel(const Signal & s) {
        if(channels != 2) return;
        size_t n=0;
        for(size_t i = 0; i < s.size(); i++)
            (*this)[n+=2] = s[i];
    }
    void setRightChannel(const Signal & s) {
        if(channels != 2) return;
        size_t n=1;
        for(size_t i = 0; i < s.size(); i++)
            (*this)[n+=2] = s[i];
    }
    // now can perform the STFT
    // which makes no sense on small real-time buffers (doh)

    Signal& operator + (const Signal& s) {
        for(size_t i = 0; i < n; i++) (*this)[i] += s;
        return *this;
    }
    Signal& operator - (const Signal& s) {
        for(size_t i = 0; i < n; i++) (*this)[i] -= s;
        return *this;
    }
    Signal& operator * (const Signal& s) {
        for(size_t i = 0; i < n; i++) (*this)[i] *= s;
        return *this;
    }
    Signal& operator / (float f) {
        for(size_t i = 0; i < n; i++) (*this)[i] /= f;
        return *this;
    }

    Signal cos() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = cos(s[i]);
        return s;
    }
    Signal sin() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = sin(s[i]);
        return s;
    }
    Signal tan() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = tan(s[i]);
        return s;
    }

    Signal acos() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = acos(s[i]);
        return s;
    }
    Signal asin() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = asin(s[i]);
        return s;
    }
    Signal atan() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = atan(s[i]);
        return s;
    }
    Signal atan2(const Signal& s2) {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = atan2(s[i],s2[i]);
        return s;
    }
    Signal atan2(const float s2) {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = atan2(s[i],s2);
        return s;
    }

    Signal cosh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = cosh(s[i]);
        return s;
    }
    Signal sinh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = sinh(s[i]);
        return s;
    }
    Signal tanh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = tanh(s[i]);
        return s;
    }

    Signal acosh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = acosh(s[i]);
        return s;
    }
    Signal asinh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = asinh(s[i]);
        return s;
    }
    Signal atanh() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = atanh(s[i]);
        return s;
    }

    Signal exp() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = exp(s[i]);
        return s;
    }
    Signal log() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = log(s[i]);
        return s;
    }
    Signal log10() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = log10(s[i]);
        return s;
    }
    Signal pow(const float s1) {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = pow(s[i],s1);
        return s;
    }
    Signal pow(const Signal& s1) {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = pow(s[i],s1[i]);
        return s;
    }

    Signal sqrt() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = sqrt(s[i]);
        return s;
    }
    Signal rsqrt() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = 1.0/sqrt(s[i]);
        return s;
    }
    Signal cbrt() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = cbrt(s[i]);
        return s;
    }
    Signal rcbrt() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = 1.0/cbrt(s[i]);
        return s;
    }
    Signal fabs() {
        Signal s(*this);
        for(size_t i = 0; i < size(); i++) s[i] = fabs(s[i]);
        return s;
    }    

    void copy(size_t n, float * data) {
        resize(n);
        memcpy(this->data(),data,n*sizeof(float));
    }


};

Signal make_sinewave(size_t n, float freq, float sample_rate = 44100.0f) {
    Signal s(n);
    SineGenerator sinwave(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = sinwave.Tick();
    return s;
}
Signal make_coswave(size_t n, float freq, float sample_rate = 44100.0f) {
    Signal s(n);
    CosGenerator coswave(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = coswave.Tick();
    return s;
}
Signal make_phasor(size_t n, float freq, float sample_rate = 44100.0f) {
    Signal s(n);
    PhasorGenerator phasor(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = phasor.Tick();
    return s;
}
Signal make_squarewave(size_t n, float freq, float duty=0.5, float sample_rate = 44100.0f) {
    Signal s(n);
    SquareGenerator sqrwave(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = sqrwave.Tick();
    return s;
}
Signal make_sawwave(size_t n, float freq, float sample_rate = 44100.0f) {
    Signal s(n);
    SawGenerator sawave(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = sawwave.Tick();
    return s;
}
Signal make_trianglewave(size_t n, float freq, float sample_rate = 44100.0f) {
    Signal s(n);
    TriangleGenerator triwave(freq,sample_rate);
    for(size_t i = 0; i < n; i++) s[i] = triwave.Tick();
    return s;
}
Signal make_noise(size_t n) {
    Signal s(n);    
    for(size_t i = 0; i < n; i++) 
    {
        s[i] = std::rand() / (float)RAND_MAX;
        if(std::rand() < RAND_MAX/2) s[i] *= 1.0f;
    }
    return s;
}

Signal loadfile(const char * filename) {
    Signal s(filename);
    return s;
}

void savefile(Signal & s, const char * filename) {
    s.saveFile(filename);
}

Vector2d& operator >> (Signal & s, Vector2d & v) {
    v.resize(s.size());
    for(size_t i = 0; i < s.size(); i++) v[i] = s[i];
    return v;
}
Vector4d& operator >> (Signal & s, Vector4d & v) {
    v.resize(s.size());
    for(size_t i = 0; i < s.size(); i++) v[i] = s[i];
    return v;
}
Vector4f& operator >> (Signal & s, Vector4f & v) {
    v.resize(s.size());
    memcpy(v.data(),s.data(),s.size()*sizeof(float));
    return v;
}
Vector8f& operator >> (Signal & s, Vector8f & v) {
    v.resize(s.size());
    memcpy(v.data(),s.data(),s.size()*sizeof(float));
    return v;
}
