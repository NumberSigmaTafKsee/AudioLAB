// Audio Effects Suite
// SimpleResampler (yes)

// IIR
// ZPK/SOS
// Resonant Filters 
// RBJ Coefficients
// Zolzer Coefficients

// FIR
// Kfr standard FIR
// Remez-Exchange 
// Windowed Sinc
// Raised Cosine

// FIR Comb Filter
// IIR Comb Filter
// Universal Comb Filter 

// Cascasder 
// Serializer
// Parallelizer 

// Delay Line
// Fractional Delay Line 

// DFT/IDFT
// DCT/IDCT 
// STFT/ISTFT (libstft)

// QMDSP 
// Aubio 
// Nanosnap
// World 
// libstft 

// Interpolator/Decimator = Resampler

// Convolution Filter
// Singal to Dirac Convolution filter
// Impulse Response Convolution Filter

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>

#include "audiosystem.h"
#include "samples/kfr_sample.hpp"
#include "samples/kfr_sample_dsp.hpp"

#include "Std/StdObject.h"
#include "Std/StdRandom.h"

#include "PolyBLEP.h"

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

//using namespace SoundWave;
using namespace SoundAlchemy;

Std::RandomMersenne noise;

float freq_to_midi(float f) {
    return 12.0*log2(f/440.0) + 69;
}
float midi_to_freq(float m) {
    return pow(2.0, (m-69)/12)*440.0;
}
float cv_to_freq(float cv)
{
    return std::pow(2,cv);
}
float semitone(int semi, float f)
{
    float m = freq_to_midi(f);
    return midi_to_freq(m + semi);
}
float octave(int octave, float f) {
    float m = freq_to_midi(f);
    return midi_to_freq(m + octave*12);
}



DSP::LowPassFilter<float> *biquad;
PolyBLEP osc(44100.0f,PolyBLEP::SAWTOOTH);

int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    LockMidi();    
    float * output = (float*)outputBuffer;
    sample_vector<float> buffer(framesPerBuffer);
    for(size_t i = 0; i < framesPerBuffer; i++)
    {
        buffer[i] = osc.Tick();
    }    
    biquad->apply(buffer);
    memcpy(output, buffer.data(), framesPerBuffer * sizeof(float));
    UnlockMidi();
    return 0;
}            


void midi_msg_print(MidiMsg * msg) {
    printf("%d %d %d\n",msg->msg,msg->data1,msg->data2);
}

void note_on(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    osc.setFrequency(freq);
}
void note_off(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;

}

void control_change(MidiMsg * msg) {
    midi_msg_print(msg);
    if(msg->data1 == 14)
    {
        float fc = midi_to_freq(12)*12*std::log(1 + (float)msg->data2/127.0f)/std::log(2);
        DSP::LowPassFilter<float> new_filter(1,fc,44100.0f);        
        *biquad->filter.filter = *new_filter.filter.filter;
        
    }
}

void repl() {

}


int main()
{
    //set_audio_func(audio_callback);
    Init();
    noise.seed_engine();    
    biquad = new DSP::LowPassFilter<float>(1,0,44100.0f);

    int num_midi = GetNumMidiDevices();
    ITERATE(i,0,num_midi)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    int num_audio = GetNumAudioDevices();
    int pulse = 10;

    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"pulse")) pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }
    
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    
    InitMidiDevice(1,3,3);
    InitAudioDevice(pulse,-1,1,44100,256);
    RunAudio();
    StopAudio();
}
