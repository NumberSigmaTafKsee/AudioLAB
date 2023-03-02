#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <complex>
#include <iostream>
#include <algorithm>

//#include "Octopus.hpp"

#include <Eigen/Core>
//#include "Casino.hpp"
#include "audiosystem.h"
typedef float DspFloatType;
#include "FX/ADSR.hpp"
#include "FX/PolyBLEP.hpp"

//#include "audio_iir_filters.hpp"
#include "audio_iir_butterworth.hpp"

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

DspFloatType Freq,Kc,Vel,Fcutoff,Fc=100.0,Qn,Q=0.5,Gain;
IIRFilters::BiquadTransposedTypeII  filter;
IIRFilters::BiquadFilterCascade     cascade;

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

IIRFilters::BiquadSection test(DspFloatType wc, DspFloatType Q)
{

    DspFloatType a0 = tan(M_PI*wc/sampleRate);
    DspFloatType a1 = a0;
    DspFloatType b0 = 1 + a0;
    DspFloatType b1 = -(1 - a0);
    
    IIRFilters::BiquadSection c;
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


struct ButterworthLowpassFilter
{
    int order;
    DspFloatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthLowpassFilter(int order, DspFloatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(DspFloatType f, DspFloatType Q) {
        auto x = IIRFilters::Butterworth::butterlp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(DspFloatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(DspFloatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, float * in, float * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<DspFloatType> impulse_response(size_t n)
    {
        std::vector<DspFloatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
};

struct ButterworthHighpassFilter
{
    int order;
    DspFloatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthHighpassFilter(int order, DspFloatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(DspFloatType f, DspFloatType Q) {
        auto x = IIRFilters::Butterworth::butterhp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(DspFloatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(DspFloatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, float * in, float * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<DspFloatType> impulse_response(size_t n)
    {
        std::vector<DspFloatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
};

struct ButterworthBandpassFilter
{
    int order;
    DspFloatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthBandpassFilter(int order, DspFloatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(DspFloatType f, DspFloatType Q) {
        auto x = IIRFilters::Butterworth::butterlp2bp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(DspFloatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(DspFloatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, float * in, float * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<DspFloatType> impulse_response(size_t n)
    {
        std::vector<DspFloatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
};

struct ButterworthBandstopFilter
{
    int order;
    DspFloatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthBandstopFilter(int order, DspFloatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(DspFloatType f, DspFloatType Q) {
        auto x = IIRFilters::Butterworth::butterlp2bs(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(DspFloatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(DspFloatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, float * in, float * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<DspFloatType> impulse_response(size_t n)
    {
        std::vector<DspFloatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
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

QuadSection series_filter(IIRFilters::BiquadSection & f1, IIRFilters::BiquadSection & f2)
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
QuadSection parallel_filter(IIRFilters::BiquadSection & f1, IIRFilters::BiquadSection & f2)
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
ButterworthLowpassFilter butters(3,44100.0);


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
    IIRFilters::BiquadSection c1 = test(Fc,Q); //butter2(Q);
    IIRFilters::BiquadSection c2 = test(Fc,Q); //butter2(Q);
    QuadSection c  = series_filter(c1,c2);    
    quad.setCoefficients(c);    
    */
    osc.setFrequency(Kc);    
    butters.setFilter(Fc,Q);
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
