/* 
SimpleSynth - The Simple Synthesizer 
Template for building crap with this junk

AudioDeviceLab.cpp = Testing codes
    PortAudioDevice
    PortMidiDevice 
    RtAudioDevice
    RtMidiDevice 

Audio 
    Play a file 
    Record from input/microphone
    Resampling 
    Mixer
    Panning
*/
/*
MIDI CC - Cutoff, Resonance
LuaJIT  - SetCutoff,SetResonance,SetOsc1Tune,SetOsc2Tune,SetADSR,SetLFORate,SetParameters
*/

#include <iostream>
#include "AudioDevice.h"
//#include "DspGraph.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <chrono>
#include <string>

// lua_register(L, “add2”, add2);
/*
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <luajit.h>
*/



#include "PolyBLEP.h"
#include "Moog.h"
#include "FunctionGenerator.hpp"
#include "ADSR.h"
#include "Amplifier.hpp"
#include "ToneStack.hpp"
#include "MVerb.h"


using namespace std;
//using namespace SoundWave;

struct HardClip
{
	float min,max;
	float gI,gO;

	HardClip() : min(-1),max(1),gI(1),gO(1) {}
	HardClip(float a, float b) : min(a),max(b),gI(1),gO(1) {} 

	float Tick(float I, float A = 1, float X = 1, float Y = 1) {
		float r = A*gI*I;
		X *= min;
		Y *= max;
		if(r < X) r = X;
		if(r > Y) r = Y;
		return gO*r;
	}
};



struct SimpleSynth
{
	PolyBLEP	 	  * osc[2];
	MoogLadderFilter  * filter;
	FunctionGenerator * lfo;
	ADSR			  * envelope;
	MVerb			  * reverb;
	Cubic			  * prefilter;
	Cubic			  * postfiler;
	HardClip		  * amp_clip;
    float osc_mix;
    float osc1_osc2_fm;
    float osc1_osc2_pm;
    float osc2_osc1_fm;
    float osc2_osc1_pm;
	float lfo_osc[2];    
	float lfo_filter_cut;
	float lfo_filter_q;
	float lfo_amplifier;
	float envelope_filter_cut;
	float envelope_filter_q;
	float envelope_amp;
	float input_gain,output_gain;
	float sr;

	SimpleSynth(float sample_rate=44100.0f) {
		sr = sample_rate;
		osc[0] = new PolyBLEP(sr);
		osc[1] = new PolyBLEP(sr);
	}
};

/*
LuaJIT Functions
SimpleSynth synth;

int SetCutoff(lua_State *L) {

}
*/



using std::cout;
using std::cin;
using std::endl;

struct StreamCallbackData {
    const float * input;
    float * output;
    unsigned long frameCount;
    const PaStreamCallbackTimeInfo *timeinfo;
    StreamCallbackFlags flags;
    void * userData;
};

int ProcessFunction(StreamCallbackData * ptr) 
{
    PolyBLEP * wav = (PolyBLEP*)ptr->userData;
    //SampleVector input(ptr->input,ptr->frameCount,1);
    //SampleVector output(ptr->output, ptr->frameCount,1);
    //wav->Run(ptr->frameCount,input, output);   
    //output *= 0.5;  
    //memcpy(ptr->output,output.data(),ptr->frameCount*sizeof(float));
    for(size_t i = 0; i < ptr->frameCount; i++) ptr->output[i] = wav->Tick();
    return 0;
}

int StreamCallback(const void * input, void * output, unsigned long frameCount, const PaStreamCallbackTimeInfo * timeinfo, PaStreamCallbackFlags statusFlags, void * userData) {
    StreamCallbackData data;
    data.input = (float*)input;
    data.output = (float*)output;
    data.frameCount = frameCount;
    data.timeinfo = timeinfo;
    data.flags = (StreamCallbackFlags)statusFlags;
    data.userData = userData;
    return ProcessFunction(&data);
}

int main()
{
    AudioDevice device;
    std::cout << device.GetDeviceCount() << std::endl;
    for(size_t i = 0; i < device.GetHostApiCount(); i++) {
        AudioHostApiInfo api(i);
        cout << "API: " << api.name() << endl;
    }
    int device_id = 14;
    for(size_t i = 0; i < device.GetDeviceCount(); i++)
    {
        AudioDeviceInfo info(i);
        cout << "Device # " << i << ": " << info.name() << endl;
        if(!strcmp(info.name(),"pulse")) device_id = i;
    }
    
    /*
    Node * n = new Node;
    Connection * c = new Connection(256,n);
    WaveTable *wav = new WaveTable(c,44100,256,1);    
    std::vector<float> v(4096);
    SoundWave::MakeSine<float>(v,4096,220,44100);
    wav->AddWaveTable(4096,v,22050);
    wav->SetFrequency(220);
    n->SetDSP(wav);
    */
    PolyBLEP blep(44100);
    AudioStreamParameters input(device_id,1,SAMPLE_FLOAT32,-1.0);
    AudioStreamParameters output(device_id,1,SAMPLE_FLOAT32,-1.0);
    AudioStream stream(&input,&output,44100,256,STREAM_NOFLAG,StreamCallback,&blep);
    stream.StartStream();
    char ch;
    std::cin >> ch;
    stream.StopStream();
    
}