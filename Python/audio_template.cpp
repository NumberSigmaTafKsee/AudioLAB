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


#include "SoundObject.hpp"
#include "AudioMidi/audiosystem.h"

#include "FX/ADSR.hpp"
#include "FX/PolyBLEP.hpp"

//#include "audio_iirfilters.hpp"
#include "audio_iir_butterworth.hpp"
#include "audio_iir_bessel.hpp"

#include "Carlo/carlo_casino.hpp"
#include "Filters/IIRFilters.hpp"
#include "octopus.hpp"
#include "octopus_functions.cpp"


#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

Default noise;
DspFloatType sampleRate=44100.0f;
DspFloatType invSampleRate=1.0/sampleRate;
int blockSize=256;

using namespace Analog::Oscillators::PolyBLEP;
using namespace Envelopes;

PolyBLEP osc(sampleRate,PolyBLEP::SAWTOOTH);
ADSR adsr(0.01,0.1,1.0,0.1,sampleRate);

DspFloatType Freq,Kc,Vel,Fcutoff,Fc=100.0,Qn,Q=0.5,Gain;
IIRFilters::BiquadTransposedTypeII  filter;
IIRFilters::BiquadFilterCascade     cascade;
using namespace Octopus;

IIRFilters::BiquadSOS oct_butterlp(int order, DspFloatType Q=1.0)
{
    OctopusValueList values;
    values.vlist(0) = order;
    values.vlist(1) = 1.0;
    values.vlist(2) = 's';
    values = octave_butter(values,3);    
    // octave tf2sos is too slow
    VectorXcf Z = values.vlist(0).float_complex_row_vector_value();
    VectorXcf P = values.vlist(1).float_complex_row_vector_value();
    float    K = values.vlist(2).float_scalar_value();
    IIRFilters::BiquadSOS sos; 
    size_t n = 0;
    if(order %2 != 0) {
        IIRFilters::BiquadSection c;
        std::complex<DspFloatType> p1  = P(n++);
        DspFloatType x1 = abs(p1);
        DspFloatType x2 = 0;
                
        // (s-p1)        
        c.z[0] = K*1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1.0;
        c.p[1] = 1.0/x1;
        c.p[2] = 0.0;
        sos.push_back(c);                        
    }         
    for(size_t i = n; i < order; i += 2)
    {
        std::complex<DspFloatType> p1  = P(n++);
        std::complex<DspFloatType> p2  = P(n++);
        
        DspFloatType x1 = abs(p1*p2);
        DspFloatType x2 = abs(-p1-p2);
        
        // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
        IIRFilters::BiquadSection c;
        c.z[0] = K*1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1;        
        c.p[1] = (1.0/Q)*x2/x1;
        c.p[2] = 1/x1;    
        sos.push_back(c);
    }   
    return sos;
}
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

struct IPPIIRBiquad: public Casino::IPP::IIRBiquad<DspFloatType>
{
    IPPIIRBiquad() = default;
    IPPIIRBiquad(const IIRFilters::BiquadSection &c) {
        setCoefficients(c);
    }
    IPPIIRBiquad(const IIRFilters::BiquadSOS & sos) {
        setCoefficients(sos);
    }
    void setCoefficients(const IIRFilters::BiquadSection & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.z[0];
        buf[1] = c.z[1];
        buf[2] = c.z[2];
        buf[3] = 1.0;
        buf[4] = c.p[0];
        buf[5] = c.p[1];
        this->initCoefficients(blockSize,1,buf);
    }
    void setCoefficients(const IIRFilters::BiquadSOS & sos)
    {        
        DspFloatType buf[6*sos.size()];
        int x = 0;
        for(size_t i = 0; i < sos.size(); i++)
        {    
            buf[x++] = sos[i].z[0];
            buf[x++] = sos[i].z[1];
            buf[x++] = sos[i].z[2];
            buf[x++] = 1.0;
            buf[x++] = sos[i].p[0];
            buf[x++] = sos[i].p[1];
        }
        this->initCoefficients(blockSize,sos.size(),buf);
    }
    void setCoefficients(const Filters::FilterCoefficients & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.b[0];
        buf[1] = c.b[1];
        buf[2] = c.b[2];
        buf[3] = 1.0;
        buf[4] = c.a[0];
        buf[5] = c.a[1];
        this->initCoefficients(blockSize,1,buf);
    }    
    void setCoefficients(const std::vector<Filters::FilterCoefficients> & c)
    {
        DspFloatType buf[6*c.size()];
        int x = 0;
        for(size_t i = 0; i < c.size(); i++)
        {    
            buf[x++] = c[i].b[0];
            buf[x++] = c[i].b[1];
            buf[x++] = c[i].b[2];
            buf[x++] = 1.0;
            buf[x++] = c[i].a[0];
            buf[x++] = c[i].a[1];
        }
        this->initCoefficients(blockSize,c.size(),buf);
    }    
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
    {
        assert(n == this->len);
        this->Execute(in,out);
    }
};

IPPIIRBiquad ipp_filter;

int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{        
    float * output = (float*)outputBuffer;    
    Eigen::Map<Eigen::VectorXf> ovec(output,framesPerBuffer);            
    
    
    //ipp_filter.setCoefficients(x); 

    //IIRFilters::BiquadSOS c = IIRFilters::Bessel::bessellp(5,Q);
    //if(Fc < 30) Fc = 30;
    //if(Fc > sampleRate/2.0) Fc = sampleRate/2.0;
    //auto x = IIRFilters::AnalogBiquadCascade(c,Fc,sampleRate);                
    //cascade.setCoefficients(x);
    
    osc.setFrequency(Kc);
    
    ovec = osc_tick(osc,framesPerBuffer);        
    ovec = env_tick(adsr,ovec,framesPerBuffer);    
    ovec = filter_tick(cascade,ovec,framesPerBuffer);        
    //ovec = filter_tick(filter,ovec,framesPerBuffer);        
    //ipp_filter.ProcessBlock(framesPerBuffer,ovec.data(),ovec.data());
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
        // it needs to be smoothe
        IIRFilters::BiquadSOS sos = oct_butterlp(2,Q);
        auto x = IIRFilters::AnalogBiquadCascade(sos,Fc,sampleRate);        
        cascade.setCoefficients(x);

    }
    if(msg->data1 == 103)
    {
        double q = (double)msg->data2/127.0f;//(pow(4.0,((double)msg->data2/127.0f))-1.0)/3.0;
        double lg1000 = (log(1000)/log(2));
        Qn = q;                    
        Q = (q*lg1000)+0.5;
        printf("Qn=%f Q=%f\n",Qn,Q);
        IIRFilters::BiquadSOS sos = oct_butterlp(2,Q);
        auto x = IIRFilters::AnalogBiquadCascade(sos,Fc,sampleRate);        
        cascade.setCoefficients(x);
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
    int pulse = 6;
    
    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"pulse")) { pulse = i; break; }
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
