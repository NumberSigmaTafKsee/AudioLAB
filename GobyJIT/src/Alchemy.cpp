#include "audiosystem.h"
#include "SoundAlchemy.hpp"
#include "TADSR.hpp"
#include "TBasicDelayLine.hpp"
#include "TDistortion.hpp"
#include "TFractionalDelayLine.hpp"
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace SoundAlchemy;


int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{
    
    return 0;
}            


float freq_to_midi(float f) {
    return 12.0*log2(f/440.0) + 69;
}
float midi_to_freq(float m) {
    return pow(2.0, (m-69)/12)*440.0;
}

void note_on(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    
}
void note_off(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;

}

void repl() {

}

TADSR<float> adsr;
BasicDelayLine<float> delay;
Distortion<float> distortion;
FractionalDelayBuffer<float> fracdelay;

int main()
{
    //set_audio_func(audio_callback);
    Init();
    int num_midi = GetNumMidiDevices();
    for(size_t i = 0; i < num_midi; i++)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    int num_audio = GetNumAudioDevices();
    int pulse = 10;
    for(size_t i = 0; i < num_audio; i++)
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }
    
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);

    
    InitMidiDevice(1,3,3);
    InitAudioDevice(pulse,-1,1,44100,256);
    RunAudio();
    StopAudio();    

}