#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <complex>
#include <iostream>
#include <algorithm>


#include <Eigen/Core>
#include "StdSamples/stdsamples.hpp"
//#include "StdSamples/stdsamples_casino.hpp"
#include "StdSamples/stdsamples_iir_filters.hpp"
#include "audiosystem.h"

typedef float DspFloatType;
#include "FX/ADSR.hpp"
#include "FX/PolyBLEP.hpp"
#include "iir_butterworth_lowpass.hpp"

//#include "Filters/IIRFilters.hpp"
using namespace iir_filters;

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

Default noise;
DspFloatType sampleRate=44100.0f;
DspFloatType invSampleRate=1.0/sampleRate;
//Octave::Octopus interp;
int blockSize=256;

using namespace Analog::Oscillators::PolyBLEP;
using namespace Envelopes;

PolyBLEP osc(sampleRate,PolyBLEP::SAWTOOTH);
ADSR adsr(0.01,0.1,1.0,0.1,sampleRate);

DspFloatType Freq,Kc,Vel,Fcutoff=1,Fc=100.0,Qn,Q=0.5,Gain;
BiquadTransposedTypeII  filter;
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

BiquadSection test(DspFloatType wc, DspFloatType Q)
{

    DspFloatType a0 = tan(M_PI*wc/sampleRate);
    DspFloatType a1 = a0;
    DspFloatType b0 = 1 + a0;
    DspFloatType b1 = -(1 - a0);
    
    BiquadSection c;
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

struct Biquad
{
    DspFloatType z[3];
    DspFloatType p[3];
    DspFloatType x[2];
    DspFloatType y[2];

    Biquad() {
        x[0] = x[1] = 0;
        y[0] = y[1] = 0;    
    }
    void setCoeffs(DspFloatType Z[3], DspFloatType P[3]) {
        memcpy(z,Z,sizeof(z));
        memcpy(p,P,sizeof(p));
    }

    DspFloatType Tick(DspFloatType I)
    {
        DspFloatType r = I*z[0] + x[0]*z[1] + x[1] * z[2] - y[0]*p[0] - y[1]*p[1];
        x[1] = x[0];
        x[0] = I;
        y[1] = y[0];
        y[0] = r;
        return r;
    }
};



struct QuadSection
{
    DspFloatType z[5];
    DspFloatType p[5];
};

struct QuadFilter
{
    DspFloatType z[5];
    DspFloatType p[5];
    DspFloatType x[4];
    DspFloatType y[4];

    QuadFilter() {
        x[0] = x[1] = x[2] = x[3] = 0;
        y[0] = y[1] = y[2] = y[3] = 0;
    }
    void setCoefficients(DspFloatType _z[5], DspFloatType _p[5]) {
        memcpy(z,_z,5*sizeof(DspFloatType));
        memcpy(p,_p,5*sizeof(DspFloatType));
    }
    void setCoefficients(QuadSection & q) {
        memcpy(z,q.z,5*sizeof(DspFloatType));
        memcpy(p,q.p,5*sizeof(DspFloatType));
    }
    DspFloatType Tick(DspFloatType I) {
        DspFloatType out = z[0]*I + z[1]*x[0] + z[2] * x[1] + z[3] * x[2] + z[4]*x[3];
        out = out - p[0]*y[0] - p[1]*y[1] - p[2]*y[2] - p[3]*y[3];
        x[3] = x[2];
        x[2] = x[1];
        x[1] = x[0];
        x[0] = I;
        y[3] = y[2];
        y[2] = y[1];
        y[1] = y[0];
        y[0] = out;
        return out;
    }    
};

// cascade = H1 * H2
// f2(f1(x))

QuadSection series_filter(BiquadSection & f1, BiquadSection & f2)
{
    QuadSection r;
    // z0 + z1 + z2 * z0 + z1 + z2 => z0^2 + 2*z0z1(z)  + 2*z0z2(z2) + z1^2(z2) + 2*z1z2(z3) + z2^2(z4)
    r.z[0] = f1.z[0] * f2.z[0];
    r.z[1] = 2*f1.z[0]*f2.z[1];
    r.z[2] = 2*f1.z[0]*f2.z[2] + f1.z[1]*f2.z[1];
    r.z[3] = 2*f1.z[1]*f2.z[2];
    r.z[4] = f1.z[2] * f2.z[2];
    r.p[0] = f1.p[0] * f2.p[0];
    r.p[1] = 2*f1.p[0]*f2.p[1];
    r.p[2] = 2*f1.p[0]*f2.p[2] + f1.p[1]*f2.p[1];
    r.p[3] = 2*f1.p[1]*f2.p[2];
    r.p[4] = f1.p[2] * f2.p[2];
    r.z[0] /= r.p[0];
    r.z[1] /= r.p[0];
    r.z[2] /= r.p[0];
    r.z[3] /= r.p[0];
    r.z[4] /= r.p[0];
    r.p[1] /= r.p[0];
    r.p[2] /= r.p[0];
    r.p[3] /= r.p[0];
    r.p[4] /= r.p[0];
    r.p[0] = r.p[1];
    r.p[1] = r.p[2];
    r.p[2] = r.p[3];
    r.p[3] = r.p[4];
    r.p[4] = 0;
    return r;
}

// Parallel = H1 + H2
// Lowshelf = G*Hlp + Hhp
// Highshelf= Hlp + G*Hhp
QuadSection parallel_filter(BiquadSection & f1, BiquadSection & f2)
{
    QuadSection r;    
    // z0 + z1 + z2 * p0 + p1 + p2 => z0*p0 + z0*p1(z) + z0*p2(z2) + z1*p0(z) + z1*p1(z2) + z1*p2(z3) + z2*p0(z2) + z2*p1(z3) + z2*p2(z4)
    // H1z * H2p + H2z * H1p
    //      H1s * H2s
    r.z[0] = f1.z[0]*f2.p[0] + f2.z[0]*f1.p[0];
    r.z[1] = f1.z[0]*f2.p[1] + f1.z[1]*f2.p[0]  + f2.z[0]*f1.p[1] + f2.z[1]*f1.p[0];
    r.z[2] = f1.z[0]*f2.p[2] + f1.z[1]*f2.p[1]  + f1.z[2]*f1.p[0] + f2.z[0]*f1.p[2] + f2.z[1]*f1.p[1]  + f2.z[2]*f1.p[0];
    r.z[3] = f1.z[1]*f2.p[2] + f1.z[2]*f2.p[1]  + f2.z[1]*f1.p[2] + f2.z[2]*f1.p[1];
    r.z[4] = f1.z[2]*f2.p[2] + f2.z[2]*f1.p[2];
    r.p[0] = f1.p[0] * f2.p[0];
    r.p[1] = 2*f1.p[0]*f2.p[1];
    r.p[2] = 2*f1.p[0]*f2.p[2] + f1.p[1]*f2.p[1];
    r.p[3] = 2*f1.p[1]*f2.p[2];
    r.p[4] = f1.p[2] * f2.p[2];
    r.z[0] /= r.p[0];
    r.z[1] /= r.p[0];
    r.z[2] /= r.p[0];
    r.z[3] /= r.p[0];
    r.z[4] /= r.p[0];
    r.p[1] /= r.p[0];
    r.p[2] /= r.p[0];
    r.p[3] /= r.p[0];
    r.p[4] /= r.p[0];
    r.p[0] = r.p[1];
    r.p[1] = r.p[2];
    r.p[2] = r.p[3];
    r.p[3] = r.p[4];
    r.p[4] = 0;
    return r;
}

QuadFilter quad;
ButterworthLowpassFilter butters(4,44100.0);


int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    //LockMidi();
    float * output = (float*)outputBuffer;    
    Eigen::Map<Eigen::VectorXf> ovec(output,framesPerBuffer);    

    /*    
    BiquadSection c1 = test(Fc,Q); //butter2(Q);
    BiquadSection c2 = test(Fc,Q); //butter2(Q);
    QuadSection c  = series_filter(c1,c2);    
    quad.setCoefficients(c);    
    */
    osc.setFrequency(Kc);    
    
    butters.setFilter(Fcutoff/10.0*(Kc+Fc),Q);
    //cascade.setCoefficients(c);
    ovec = osc_tick(osc,framesPerBuffer);        
    ovec = env_tick(adsr,ovec,framesPerBuffer);
    //ovec = filter_tick(cascade,ovec,framesPerBuffer);        
    //ovec = filter_tick(filter,ovec,framesPerBuffer);        
    /*
    for(size_t i = 0; i < framesPerBuffer; i++)
        ovec[i] = quad.Tick(ovec[i]);
    */
    butters.ProcessBlock(framesPerBuffer,output,output);
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
    std::cout << "Start\n";

    int num_midi = GetNumMidiDevices();
    ITERATE(i,0,num_midi)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    int num_audio = GetNumAudioDevices();
    int pulse = 0;
    
    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) { pulse = i;  }
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
