///////////////////////////////////////////////////////////////////
// Faust FX
// GTK
// No Gui
///////////////////////////////////////////////////////////////////

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <list>
#include <new>

#include "AudioMidi/audiosystem.h"
#include "Threads.hpp"
#include "Std/StdObject.h"
#include "Std/StdRandom.h"
#include "samples/sample.hpp"
#include "samples/sample_dsp.hpp"

#include "DSP/Decimators.hpp"
#include "MusicFunctions.hpp"

// need to improve these more
//#include "Resamplers.hpp"
//#include "SndFile.hpp"
//#include "WaveFile.hpp"

#include "Faust/FaustFX.hpp"

using namespace std;

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)

//using namespace SoundWave;
//using namespace SoundAlchemy;


float sampleRate = 44100.0f;
float invSampleRate = 1.0f/sampleRate;
float sampleRate2x = 2*sampleRate;

Std::RandomMersenne noise;



// faust
float ** temp_in;
float ** temp_out;


//std::vector<WaveFile*> wavefiles;
//SimpleResampler resampler;
//kfr::samplerate_converter<float> kfr_resampler(kfr::resample_quality::normal,44100,88200);

FaustPolyFX * dx7;
FaustFX * fverb;

float phase = 0.0f;

//BasicDelay basic_delay(44100.0f,2,500,0.5,0.5);
int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    float * input = (float*)inputBuffer;
    float * output = (float*)outputBuffer;
            
    for(size_t i = 0; i < framesPerBuffer; i++)
    {
        float t = 0; //input[2*i];
        temp_in[0][i] = t;
        temp_in[1][i] = t;
    }

    dx7->Run(framesPerBuffer,temp_in,temp_out);    
    fverb->Run(framesPerBuffer,temp_out,temp_out);
    
    float * p = output;
    for(size_t i = 0; i < framesPerBuffer; i++)
    {
        //float x = out[i];
        *p++ = temp_out[0][i];
        *p++ = temp_out[1][i];
    }
            
    return 0;
}            


void midi_msg_print(MidiMsg * msg) {
    printf("%d %d %d\n",msg->msg,msg->data1,msg->data2);
}

float last_freq;
float last_velocity;
int   notes_pressed=0;

float log_scale(float x) {
    return std::log(1+x)/std::log(2);
}

void note_on(MidiMsg * msg) {    
    float freq = midi_to_freq(msg->data1);
    float velocity = log_scale(0.5+0.5*msg->data2/127.0f);

    notes_pressed++;
    dx7->noteOn(msg->data1, msg->data2);
    last_freq = freq;
    last_velocity  = velocity;
}
void note_off(MidiMsg * msg) {
    float freq = midi_to_freq(msg->data1);
    float velocity = log_scale(0.5 + 0.5*msg->data2/127.0f);
    notes_pressed--;
    
    dx7->noteOff(msg->data1,msg->data2);
    // it shouldn't but it does
    if(notes_pressed <= 0) {
        notes_pressed = 0;        
    }    
}

void control_change(MidiMsg * msg) {
    if(msg->data1 == 14) {
        float f = (float)msg->data2/127.0f;
        f = std::log(1 + f)/std::log(2);
        f *= 22050.0f/12.0f;        
        f += last_freq;                
    }
    if(msg->data1 == 15) {
        float f = (float)msg->data2/127.0f;        
        f = std::log(1 + f)/std::log(2);                
        
    }
}

#include "LuaJIT.hpp"
LuaJIT * lua;

int dx7Rand(lua_State * L)
{
    for(std::string& name : dx7->interface->names)
    {
        //std::cout << name << std::endl;
        float r = dx7->interface->controls[name].fMin + 
                    (dx7->interface->controls[name].fMax - dx7->interface->controls[name].fMin)*noise.rand();
        std::cout << "Before= " << dx7->interface->getControl(name) << std::endl;
        dx7->interface->setControl(name,r);
        std::cout << "After=" << dx7->interface->getControl(name) << std::endl;
    }
    return 0;
}

int revRand(lua_State * L)
{
    for(std::string& name : fverb->interface->names)
    {
        //std::cout << name << std::endl;
        float r = fverb->interface->controls[name].fMin + 
                    (fverb->interface->controls[name].fMax - fverb->interface->controls[name].fMin)*noise.rand();
        
        fverb->interface->setControl(name,r);
        
    }
    
    return 0;
};

void connectLua()
{
    lua = new LuaJIT("main.lua");
    lua->CreateCFunction("dx7Rand",dx7Rand);
    lua->CreateCFunction("revRand",revRand);
}

void repl() {
    std::string cmd;
    std::cin >> cmd;
    lua->DoCmd(cmd);
}


int main()
{
  
    Init();
    connectLua();

    dx7 = new FaustPolyFX("Faust/dx7.dsp");
    fverb = new FaustFX("Faust/fverb.dsp");

    temp_in = new float*[2];
    temp_in[0] = new float[1024];
    temp_in[1] = new float[1024];
    
    temp_out = new float*[2];
    temp_out[0] = new float[1024];
    temp_out[1] = new float[1024];
    
    
    noise.seed_engine();    

    //resampler.setup(44100.0f,2);        
        
        
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    
    int num_midi = GetNumMidiDevices();
    ITERATE(i,0,num_midi)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }

    int num_audio = GetNumAudioDevices();
    int pulse = 8;
    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }
    
    InitMidiDevice(0,3,2);
    InitAudioDevice(pulse,pulse,2,44100,256);
    RunAudio();    
    StopAudio();
    delete fverb;
}
