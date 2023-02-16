#include "audiosystem.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <cassert>
#include <signal.h>

using std::string;

struct LFO 
{
    WaveTable * sine;
    WaveTable * saw;
    WaveTable * square;
    WaveTable * triangle;
    float       frequency;
    float       sampleRate;
    enum WaveType {
        SINE,
        SAW,
        SQUARE,
        TRIANGLE
    } current_wave;

    LFO(size_t samples, WaveType type, float freq, float sr) {
        std::vector<float> wav(samples);
        frequency = freq;
        sampleRate= sr;
        current_wave = type;

        MakeSine(wav,1.0,sr);
        sine = new WaveTable();
        sine->addWaveTable(wav,22050.0f);
        
        MakeSaw(wav,1.0,sr);
        saw = new WaveTable();
        saw->addWaveTable(wav,22050.0f);

        MakeSquare(wav,1.0,sr);
        square = new WaveTable();
        square->addWaveTable(wav,22050.0f);

        MakeTriangle(wav,1.0,sr);
        triangle = new WaveTable();
        triangle->addWaveTable(wav,22050.0f);
    }

    void setFrequency(float f) {
        frequency = f;
        sine->setFrequency(f);
        saw->setFrequency(f);
        square->setFrequency(f);
        triangle->setFrequency(f);
    }
    void setPhaseOffset(float offset)
    {
        sine->setPhaseOffset(offset);
        saw->setPhaseOffset(offset);
        square->setPhaseOffset(offset);
        triangle->setPhaseOffset(offset);
    }
    float Tick() {
        float r = 0.0;
        switch(current_wave) {
            case SINE: r = sine.Tick(); break;
            case SAW: r = saw.Tick(); break;
            case SQUARE: r = square.Tick(); break;
            case TRIANGLE: r = triangle.Tick(); break;
        }
        return r;
    }
    void Process(size_t n, float * samples) {
        for(size_t i = 0; i < n; i++) samples[i] = Tick();
    }
};

struct SimpleSynth
{
    BandlimitedOscillator *osc1,*osc2;
    MoogLadderFilter      *filter;
    LFO                   *lfo;
    ADSR                  *env1,*env2;
    PinkNoiseGenerator    *noise;
};

int note_on(MidiMsg * msg) {

    return 0;
}
int note_off(MidiMsg * msg) {

    return 0;
}

int audio_callback(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer,
                    const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData)
{

    return paContinue;
}

int main(int argc, char * argv[]) 
{    
    Init();
    Pm_Initialize();

    set_note_on_callback(note_on);
    set_note_off_callback(note_off);

    int num_midi_device = GetNumMidiDevices();
    for(size_t i = 0; i < num_midi_device; i++)
        printf("Device #%lu : %s\n",i,GetMidiDeviceName(i));

    set_audio_callback(audio_callback);
    int device=14;
    
    Pa_Initialize();
    int num_audio_device = GetNumAudioDevices();
    for(size_t i = 0; i < num_audio_device; i++)
    {
        string s = GetAudioDeviceName(i);
        if(s == "pulse") device = i;
        printf("Audio Device#%lu: %s\n",i,s.c_str());
    }
    InitMidiDevice(1,1,1);
    InitAudioDevice(device,-1,2,44100,256);
    RunAudio();
    StopAudio();

}