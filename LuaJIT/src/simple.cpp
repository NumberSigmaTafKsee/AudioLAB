// Audio Effects Suite
// SimpleResampler (yes)

#include "Gamma.hpp"
#include "Stk.hpp"

#include "audiosystem.h"
#include "include/PolyBLEP.h"
#include "include/Moog.h"
#include "include/ADSR.h"
#include "ToneStack.hpp"
#include "include/WaveTable.h"
#include "LFO.h"
#include "SimpleResampler.hpp"
#include "kfr_sample.hpp"
#include "kfr_sample_dsp.hpp"
//#include "cppfilters.hpp"
//#include "Biquad.hpp"
//#include "musicdsp_biquad.h"
//#include "Filter.hpp"
#include "Decimators.hpp"
#include "Resamplers.hpp"
#include "TWaveCycle.hpp"
// not compatible with c++17
//#include "freeverb/nrev.hpp"
#include "MVerb.h"
//#include "BodeShifter.h"
//#include "include/AudioEffectsSuite/DelayEffects/DelayEffects.h"
//#include "include/AudioEffectsSuite/FilterEffects/FilterEffects.h"

#include <random>
#include <chrono>

#include "Std/StdObject.h"
// wtf?
#include "Std/StdRandom.h"


#include <cstring>
#include <cstdlib>
#include <cstdio>

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

using namespace SoundWave;
using namespace SoundAlchemy;

Std::RandomMersenne noise;

float freq_to_midi(float f) {
    return 12.0*log2f(f/440.0) + 69;
}
float midi_to_freq(float m) {
    return 440.0f*powf(2.0, (m-69)/12);
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



struct Simple
{
    PolyBLEP * osc[2];    
    TWaveCycle<float> * cycle;
    MoogLadderFilter * filter;
    ADSR     * env;
    LFO      * lfo;
    ToneStack * tone;
    float      fs;
    float      osc_fine[2];
    float      osc_tune[2];
    float      osc_freq[2];    
    float      filter_cutoff; // this is the knob value
    float      filter_resonance;
    float      lfo_freq;
    float      lfo_osc[2];
    float      lfo_filter_cutoff;
    float      lfo_filter_resonance;
    float      lfo_amp;
    float      input_gain,output_gain;
    float      minClip,maxClip;
    float      freq,velocity;

    Simple(float sr = 44100)
    {
        fs = sr;
        osc[0] = new PolyBLEP(sr,PolyBLEP::SAWTOOTH);
        osc[1] = new PolyBLEP(sr,PolyBLEP::SAWTOOTH);            
        cycle = new TWaveCycle<float>(sr);
        filter = new MoogLadderFilter(FINN_MOOG,sr);
        for(size_t i = 0; i < 2; i++)
        {
            osc_fine[i] = 0.0f;
            osc_tune[i] = 0.0f;
            osc_freq[i] = 440.0f;            
        }
        filter_cutoff = 0;        
        filter_resonance = 0.0f;
        lfo = new LFO(sr);
        env = new ADSR(0.1,0.2,0.75,0.2,sr);
        tone = new ToneStack(880,15000,sr);
        input_gain = 1.0f;
        output_gain = 1.0f;
        minClip = -0.95;
        maxClip = 0.95;
    }
    ~Simple() {
        for(size_t i = 0; i < 2; i++) if(osc[i]) delete osc[i];
        if(filter) delete filter;
        if(lfo) delete lfo;
        if(tone) delete tone;
    }
    void set_osc2_semitone(float semi) {
        osc_tune[1] = semi;
    }
    void set_osc2_fine(float f) {
        osc_fine[1] = f;
    }
    void set_frequency(float f) {
        freq = clamp(f,0,fs/2);        
        for(size_t i = 0; i < 2; i++) 
        {
            osc_freq[i] = f + 12*osc_fine[i];            
            osc[i]->setFrequency(f + semitone(osc_tune[i],f));
        }        
	    cycle->setFrequency(f);                        
    }
    void set_cutoff(float c) {    
        filter_cutoff = clamp(c,0,fs/2);    
        filter->SetCutoff(filter_cutoff);
    }
    void set_resonance(float r) {
        filter_resonance = clamp(r,0,1);
        filter->SetResonance(filter_resonance);
    }
    void set_velocity(float v) {
        velocity = v;
    }
    void note_on() {
        env->gate(1);
    }
    void note_off(float f, float v) {
        if(f == freq)
            env->gate(0);
    }

    float Tick(float I=1,float A=1, float X=0,float Y=0) {        
        float r = osc[0]->Tick() + osc[1]->Tick();
        //float r = cycle->Tick();
        float cut = filter_cutoff;
        float res = filter_resonance;       
        float c   = clamp(freq*12*std::log(1+(filter_cutoff + X))/std::log(2),0,fs/2);        
        set_cutoff(c);
        float t = res + Y*0.5;
        if(t < 0) t = 0.0f;
        if(t > 1) t = 1.0f;
        set_resonance(std::log(1 + t)/std::log(2));
        r = clamp(A*r*input_gain,minClip,maxClip);
        r = filter->Tick(r);
        r *= env->Tick();
        set_cutoff(cut);
        set_resonance(res);
        return std::tanh(output_gain*r);	
    }
};

struct WaveFile
{
    sample_vector<float> buffer;
    int cur_pos;
    bool playing;
    bool loop;

    WaveFile(const char * filename) {
        buffer = load_wav<float>(filename);
        assert(buffer.size() > 0);
        cur_pos = 0;
        playing = false;
        loop = false;
    }
    void process(size_t n, float * buf) {
        if(playing == false) return;
        if(cur_pos >= buffer.size()) 
            if(loop) cur_pos = 0;
            else return;
        for(size_t i = cur_pos; i < cur_pos+n; i++)
        {
            if(cur_pos >= buffer.size()) break;
            buf[i] = 0.5*(buf[i] + buffer[i]);
        }
        cur_pos += n;
    }
    void stop() {
        playing = false;
    }
    void start() {
        playing = true;
    }
    void pause() {
        playing = !playing;
    }
    void reset() {
        cur_pos = 0;
    }
    void toogle_loop(bool f) {
        loop = f;
    }
};



Simple * simple;
RBJLowPassFilter lowpass(15000.0f,88200.0f);
KfrBiquadFilter<float>  filter(2,15000.0f,88200.0f);
Decimateur7 decimator;
std::vector<WaveFile*> wavefiles;

MVerb<float> mverb;

float ** temp1;
float ** temp2;
int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    float * output = (float*)outputBuffer;
    sample_vector<float> temp(framesPerBuffer*2);    
    zeros(temp);
    for(size_t i = 0; i < framesPerBuffer*2; i++)
        temp[i] = simple->Tick();
    //lowpass.Process(temp.data(),framesPerBuffer*2);
    filter.dolowpass(temp);    
    
    //sample_vector<float> out = downsampler_high(2,88200.0f,temp);
    //memcpy(output,out.data(),out.size()*sizeof(float));
    for(size_t i = 0; i < framesPerBuffer; i++) {
       temp1[0][i] = temp1[1][i] = decimator.Calc(temp[i*2],temp[i*2+1]);
    }
    //for(size_t i = 0; i < wavefiles.size(); i++)
    //    wavefiles[i]->process(framesPerBuffer,output);
    
    
    mverb.process(temp1,temp2,framesPerBuffer);
    

    float * p = output;
    for(size_t i = 0; i < framesPerBuffer; i++)
    {
        *p++ = temp2[0][i];
        *p++ = temp2[1][i];
    }
    return 0;
}            



void note_on(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    simple->set_frequency(freq);
    simple->set_velocity(velocity);
    simple->note_on();
    printf("%f %f\n",freq,velocity);
}
void note_off(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = msg->data2/127.0f;
    simple->note_off(freq,velocity);
}
void midi_msg_print(MidiMsg * msg) {
    printf("%d %d %d\n",msg->msg,msg->data1,msg->data2);
}
void control_change(MidiMsg * msg) {
    //midi_msg_print(msg);
    if(msg->data1 == 14) {
        float v = (float)msg->data2/127.0f;
        simple->filter_cutoff = v;        
    }
    if(msg->data1 == 15) {
        float v = (float)msg->data2/127.0f;
        simple->filter_resonance = v;        
    }
}

void repl() {

}


int main()
{
    //set_audio_func(audio_callback);
    Init();
    noise.seed_engine();    
    WaveFile * f = new WaveFile("sine1k.wav");
    //f->start();
    wavefiles.push_back(f);
    
    temp1 = (float**)calloc(2,sizeof(float*));
    temp1[0] = (float*)calloc(1024,sizeof(float));
    temp1[1] = (float*)calloc(1024,sizeof(float));

    temp2 = (float**)calloc(2,sizeof(float*));
    temp2[0] = (float*)calloc(1024,sizeof(float));
    temp2[1] = (float*)calloc(1024,sizeof(float));
    
    mverb.setSampleRate(44100.0f);
    mverb.setParameter(MVerb<float>::MIX,0.5);
    mverb.setParameter(MVerb<float>::SIZE,0.2);
    mverb.setParameter(MVerb<float>::PREDELAY,0.01);

    int num_midi = GetNumMidiDevices();
    for(size_t i = 0; i < num_midi; i++)
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

    simple = new Simple(44100*2);
    simple->set_osc2_semitone(0);
    simple->set_osc2_fine(0.1);
    InitMidiDevice(1,3,3);
    InitAudioDevice(pulse,-1,2,44100,256);
    RunAudio();
    StopAudio();
    delete simple;

}
