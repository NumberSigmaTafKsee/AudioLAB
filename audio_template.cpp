#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <complex>
#include <iostream>
#include <algorithm>

#include "Samples.hpp"
#include "SamplesDSP.hpp"
#include "Octopus.hpp"

#include <Eigen/Core>
//#include "Casino.hpp"
#include "audiosystem.h"
typedef float DspFloatType;
#include "FX/ADSR.hpp"
#include "FX/PolyBLEP.hpp"

#include "audio_iir_filters.hpp"

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

Default noise;
DspFloatType sampleRate=44100.0f;
DspFloatType invSampleRate=1.0/sampleRate;
Octave::Octopus interp;
int blockSize=256;

using namespace Analog::Oscillators::PolyBLEP;
using namespace Envelopes;

PolyBLEP osc(sampleRate,PolyBLEP::SAWTOOTH);
ADSR adsr(0.01,0.1,1.0,0.1,sampleRate);

DspFloatType Freq,Kc,Vel,Fcutoff,Fc=100.0,Qn,Q=0.5,Gain;
Filters::BiquadTransposedTypeII  filter;
Filters::BiquadFilterCascade     cascade;

template<typename Osc>
Eigen::VectorXf osc_tick(Osc & o, size_t n)
{
    Eigen::VectorXf r(n);
    for(size_t i = 0; i < n; i++) r[i] = o.Tick();
    return r;
}
template<typename Envelope>
Eigen::VectorXf env_tick(Envelope & e, size_t n)
{
    Eigen::VectorXf r(n);
    for(size_t i = 0; i < n; i++) r[i] = e.Tick();
    return r;
}
template<typename Envelope>
Eigen::VectorXf env_tick(Envelope & e, Eigen::Map<Eigen::VectorXf>& v, size_t n)
{
    Eigen::VectorXf r(n);
    for(size_t i = 0; i < n; i++) r[i] = v[i]*e.Tick();
    return r;
}
template<typename Filter>
Eigen::VectorXf filter_tick(Filter & f, Eigen::Map<Eigen::VectorXf>& map, size_t n)
{    
    Eigen::VectorXf samples(n);
    for(size_t i = 0; i < n; i++) samples[i] = f.Tick(map[i]);
    return samples;
}



Filters::BiquadSection test(DspFloatType wc, DspFloatType Q)
{

    DspFloatType a0 = tan(M_PI*wc/sampleRate);
    DspFloatType a1 = a0;
    DspFloatType b0 = 1 + a0;
    DspFloatType b1 = -(1 - a0);
    
    Filters::BiquadSection c;
    c.z[0] = a0*a0;
    c.z[1] = 2*(a0*a1);
    c.z[2] = a1*a1;
    c.p[0] = b0*b0;
    c.p[1] = 2*(b0*b1);
    c.p[2] = b1*b1;
    c.z[0] /= c.p[0];    
    c.z[1] /= c.p[0];
    c.z[2] /= c.p[0];
    c.p[1] /= c.p[0];
    c.p[2] /= c.p[0];
    c.p[0]  = c.p[1];
    c.p[1]  = c.p[2];
    
    return c;
}

Filters::BiquadSOS series_filter(Filters::BiquadSection & f1, Filters::BiquadSection & f2)
{
    Filters::BiquadSOS r;
    r.push_back(f1);
    r.push_back(f2);
    return r;
}
Filters::BiquadSOS parallel_filter(Filters::BiquadSection & f1, Filters::BiquadSection & f2)
{
    Filters::BiquadSOS r;
    Filters::BiquadSection c;
    c.z[0] = f1.z[0] + f2.z[0];
    c.z[1] = f1.z[1] + f2.z[1];
    c.z[2] = f1.z[2] + f2.z[2];
    c.p[0] = f1.p[0] * f2.p[0];
    c.p[1] = f1.p[1] * f2.p[1];
    c.p[2] = f1.p[2] * f2.p[2];
    r.push_back(c);
    return r;
}
int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    //LockMidi();
    float * output = (float*)outputBuffer;    
    Eigen::Map<Eigen::VectorXf> ovec(output,framesPerBuffer);    
    Filters::BiquadSection c1 = test(Fc,Q); //butter2(Q);
    Filters::BiquadSection c2 = test(Fc,Q); //butter2(Q);
    Filters::BiquadSOS c  = series_filter(c1,c2);    
    osc.setFrequency(Kc);    
    //filter.setCoefficients(c);    
    cascade.setCoefficients(c);
    ovec = osc_tick(osc,framesPerBuffer);        
    ovec = env_tick(adsr,ovec,framesPerBuffer);
    ovec = filter_tick(cascade,ovec,framesPerBuffer);        
    //UnlockMidi();
    return 0;
}            


float last_freq;
float last_vel;
int   notes_pressed=0;
int   currentNote=69;
int   currentVelocity=0;

void note_on(MidiMsg * msg) {    
    float freq = MusicFunctions::midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    currentNote = msg->data1;
    currentVelocity = msg->data2;
    Freq = MusicFunctions::freq2cv(freq);
    Kc = freq;
    Vel  = velocity;    
    adsr.noteOn();         
    last_freq = Freq;
    last_vel  = velocity;
    notes_pressed++;    
}
void note_off(MidiMsg * msg) {
    float freq = MusicFunctions::midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    notes_pressed--;
    if(notes_pressed <= 0)
    {
        notes_pressed = 0;
        adsr.noteOff();        
    }
}


void midi_msg_print(MidiMsg * msg) {
    printf("%d %d %d\n",msg->msg,msg->data1,msg->data2);
}

void control_change(MidiMsg * msg) {
    midi_msg_print(msg);
    if(msg->data1 == 102)
    {
        double fc = (pow(127.0,((double)msg->data2/127.0f))-1.0)/126.0;
        Fcutoff = 10*fc;        
        Fc = fc*(sampleRate/2);
        printf("Fcutoff=%f Fc=%f\n",Fcutoff,Fc);
    }
    if(msg->data1 == 103)
    {
        double q = (double)msg->data2/127.0f;//(pow(4.0,((double)msg->data2/127.0f))-1.0)/3.0;
        double lg1000 = (log(1000)/log(2));
        Qn = q;                    
        Q = (q*lg1000)+0.5;
        printf("Qn=%f Q=%f\n",Qn,Q);
    }
}


void repl() {
}

int main()
{
    //set_audio_func(audio_callback);
    Init();
    noise.seed_engine();    
    

    int num_midi = GetNumMidiDevices();
    ITERATE(i,0,num_midi)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    int num_audio = GetNumAudioDevices();
    int pulse = 0;
    
    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) { pulse = i; break; }
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }
    
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    
    InitMidiDevice(1,3,3);
    InitAudioDevice(pulse,-1,1,sampleRate,blockSize);
    RunAudio();
    StopAudio();
}
