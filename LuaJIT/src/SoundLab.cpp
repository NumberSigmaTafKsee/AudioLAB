#include <vector>
#include "SampleVector.h"
#include "Std.h"
#include "SndFile.h"
#include "Nyquist.h"
#include "Resampler.h"

typedef SampleType float;

enum SampleFormat {
    FORMAT_FLOAT32 = paFloat32,
    FORMAT_INT16   = paint16,
    FORMAT_INT32   = paint32,
    FORMAT_INT24   = paint24,
    FORMAT_INT8    = paint8,
    FORMAT_UINT8   = paUint8,
};

enum WavType {
    WAV_WAVE,
    WAV_MP3,
    WAV_FLAC,
    WAV_OGG,
    WAV_RAW,
    WAV_AU,
    WAV_AIFF
};

enum Quality
{
	QUALITY_BEST = SRC_SINC_BEST_QUALITY,
	QUALITY_MEDIUUM = SRC_SINC_MEDIUM_QUALITY,
	QUALITY_FASTEST = SRC_SINC_FASTEST,
	QUALITY_ZERO_ORDER_HOLD = SRC_ZERO_ORDER_HOLD,
	QUALITY_LINEAR = SRC_LINEAR,
} ;


SampleVector mix(size_t block, size_t channels, const StdVector<SampleVector> & samples) {
    SampleVector r(block,channels);
    for(size_t i = 0; i < samples.size(); i++)
        r += samples[i];
    r /= (float)samples.size();
    return r;
}

SampleVector level(float gain, const SampleVector & samples) {
    return samples * gain;    
}


// stream from disk
struct ChunkSample {

};

// holds entire sample in memory
struct Sample {

    SampleVector samples;
    bool    isLoop;
    int     current_position;

    Sample(SampleVector & s, bool loop = false) {
        isLoop = loop;
        current_position = 0;
        samples = s;
    }

    void Loop() { isLoop = true; } 
    void OneShot() { isLoop = false; }

    void Load(const char * filename, WavType type) 
    {

    }
    void SaveRaw(const char * filename) {

    }
    void SaveWav(const char * filename) {

    }
    void Process(size_t block, SampleVector & output) {                
        if(isLoop == true && current_position >= samples.size()) current_position = 0;
        else if(current_position >= samples.size()) return;
        if(current_position+block > samples.size()) block = samples.size() - current_position;
        output = samples.slice(current_position,block);
        current_position += block;
    }
    SampleBuffer Buffer(size_t block) {
        SampleBuffer r(block,samples.num_channels());
        r.read(block,samples);
    }
};

// libsamplerate
struct Resampler : public Node {
    SRC_STATE * src;

    Resampler(Quality converter_type, int channels) {
        int error = 0;
        src = src_new(converter_type,channels,&error);
        assert(src != NULL);    
    }
    ~Resampler() {
        if(src != NULL) src_delete(src);
    }

    void Process(const SampleVector & samples, SampleVector & output) {
        
    }
};



struct SoundSystem
{
    std::unique_ptr<SoundDevice> input;
    std::unique_ptr<SoundDevice> output;

    SoundSystem() {
        AudioInitialize();        
        input  = unique_ptr<SoundDevice>(NULL);
        output = unique_ptr<SoundDevice>(NULL);
    }
    ~SoundSystm() {
        AudioStop();
    }
    
        
    void SetInput() {
        this->input  = std::unique_ptr<SoundDevice>(input, [](SoundDevice *p){ delete p; });
    }
    void SetOuttput() {
        this->output = std::unique_ptr<SoundDevice>(output, [](SoundDevice *p){ delete p; });
    }
    SoundDevice* GetInputDevice() { return input.get(); }
    SoundDevice* GetOutputDevice() { return output.get(); }    
};


struct Control {
    SampleType *port;

    Control(SampleType * p) : port(p) { }
    
    void        SetValue(SampleType x) { *port = x; }
    SampleType& GetValue() { return *port; }

    operator SampleType() { return *port; } 

};



struct Port {    
    SampleVector input,output;
    SampleMatrix ins,outs;

    Port(size_t n, size_t c = 1) { input.resize(n,c); output.resize(n,c); }    
};


struct Node {

    bool isPlaying;
    SampleType gain;

    std::vector<std::shared_ptr<Port>> port_vector;
    std::vector<std::shared_ptr<Channel>> channel_vector;

    std::map<std::string,Port*>         ports;
    std::map<std::string,Control*>      controls;

    Node() {
        isPlaying = false;
        gain = 1.0;
    }
    virtual void Process(size_t block, SampleVector & output) = 0;
    
    void AddPort(const std::string & name, Port * p) {
        std::shared_ptr<Port> spp = { p, [](Port * p) { delete p; }; }
        port_vector.push_back(spp);
        ports[name] = p;
    }
    void AddControl(const std::string & name, Channel * p) {
        std::shared_ptr<Channel> cpp = { p, [](Channel * c) { delete c; }; }
        channel_vector.push_back(cpp);
        channels[name] = c;
    }
    Control* GetControl(const std::string & name) {
        return controls[name];
    }
    Port* GetPort(const std::string & name) {
        return ports[name];
    }

    void SetVolume(SampleType level) {
        gain = level;
    }
    void Loop() { 
        sampler->Loop();
    }
    void OneShot() {
        sampler->OneShot();
    }
    void Play() { 
        // schedule for playing
        isPlaying = true;
    }
    void Stop() { 
        // stop playing 
        isPlaying = false;
    }
    void Reset() {
        Stop();    
        isPlaying = false;
    }
    SampleVector Level(const SampleVector * in) {
        return in * gain;
    }

};

struct LuaNode : public Node {
    // SWIGLUA_REF func;
};

struct Recorder : public Node {

};


struct Oscillator : public Node {
    // polyblep    
};

struct LFO : public Node {
    // rchoscillators
};

struct NoiseGenerator : public Node {

};


struct Processor : public Node {

    
};

struct Filter : public Node {
    
    
};

struct Mixer : public Node {


};

struct WaveTable : public Generator {
    
};

struct Instrument : public Generator {
    // instrument module made with synthesis libraries etc
};

struct LV2Plugin : public Node {
    // LV2Host
};

struct ModHost : public Node {
    // connection to modhost
};

struct SampleBuffer : public Node {    
    Samples      sampler;
    
    SampleBuffer(size_t s, size_t c) { buffer.resize(s,c);  }
    SampleBuffer(const SampleBuffer & b) { buffer = b.buffer; }
    SampleBuffer(const SampleVector & v) { buffer = v; }

    SampleBuffer& operator = (const SampleBuffer *&b) {
        buffer = b.buffer;        
        return *this;
    }

    SampleType& operator[](size_t i) { return buffer[i]; }

    SampleVector Read(size_t block) { 
        SampleVector r(block,buffer.num_channels());
        r.zero();
        r.copy(buffer,block);        
        return r;
    }
    void Write(size_t block const SampleVector & input) {        
        buffer.copy(input,block);
    }

    void Process(size_t block, SampleVector & out) {
        if(!isPlaying) return;        
        sampler.process(block,out);            
    }

    
};



