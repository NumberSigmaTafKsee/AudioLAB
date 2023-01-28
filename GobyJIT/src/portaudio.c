#include <portaudio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

// PortAudio
struct AudioDeviceInfo 
{
    PaDeviceInfo input,output;
};

struct AudioStreamParameters 
{
    PaStreamParameters parameters;
};


static int version = 0;
static const char* version_text = NULL;
static int host_api_count = 0;
static int default_host_api = 0;
static int device_count=-1;
static int default_input=-1;
static int default_output=-1;
static bool audio_init = false;
static PaStream * stream;

void InitAudio();
void TerminateAudio();
void Sleep(long msec);

static void (*process_f)(float * samples, size_t n);


void SetProcess(void (*process)(float*,size_t)) {
    process_f = process;
}

void Sleep(long msec) {
    Pa_Sleep(msec);
}
void InitAudio() {
   int err =  Pa_Initialize();
   if(err != paNoError) {
       printf("Can't initialize portaudio\n");
       exit(-1);
   }
   audio_init = true;
   process_f = NULL;
   version = Pa_GetVersion();
   version_text = Pa_GetVersionText();
   host_api_count = Pa_GetHostApiCount();
   default_host_api = Pa_GetDefaultHostApi();
   device_count = Pa_GetDeviceCount();
   default_input = Pa_GetDefaultInputDevice();
   default_output = Pa_GetDefaultOutputDevice();
}
void TerminateAudio() {
    if(!audio_init) return;
    audio_init = false;
    int err = Pa_Terminate();
    if(err != paNoError) {
        printf("Failed to terminate portaudio\n");
    }
}

int GetNumAudioDevices();
void GetAudioDeviceInfo(AudioDeviceInfo * info);
bool IsFormatSupported(AudioDeviceParameters * p, double sample_rate);
void SetInputParameters(AudioDeviceParameters * p, int device_index, int channels, SampleFormat format, double latency);
void SetOutputParameters(AudioDeviceParameters * p, int device_index, int channels, SampleFormat format, double latency);

int GetNumAudioDevices() {
    int numDevices = Pa_GetDeviceCount();
    return numDevices;
}

void GetAudioDeviceInfo(AudioDeviceInfo * info)
{
    assert(info != NULL);
    Pa_GetDeviceInfo(&info->info);
}

bool IsFormatSupported(AudioDeviceParameters * p, double sample_rate) {
    assert(p != NULL);
    PaError err = Pa_IsFormatSupported(&p->input,&p->output, sample_rate)
    return err == paFormatIsSupported;
}

void SetInputParameters(AudioDeviceParameters * p, int device_index, int channels, SampleFormat format, double latency)
{
    assert(p != NULL);
    p->input.channelCount = channels;
    p->input.device       = device_index;    
    p->input.sampleFormat = format;
    p->input.suggestedLatency = latency == 1.0? Pa_GetDeviceInfo(device_index)->defaultLowOutputLatency : latency;
    p->input.hostApiSpecificStreamInfo = NULL;    
}
void SetOutputParameters(AudioDeviceParameters * p, int device_index, int channels, SampleFormat format, double latency) 
{
    assert(p != NULL);
    p->output.channelCount = channels;
    p->output.device       = device_index;    
    p->output.sampleFormat = format;
    p->output.suggestedLatency = latency == 1.0? Pa_GetDeviceInfo(device_index)->defaultLowOutputLatency : latency;
    p->output.hostApiSpecificStreamInfo = NULL;    
}


typedef struct {
    float left_phase;
    float right_phase;
}
paTestData;
static paTestData data;

static int StreamCallback( const void *input,
                            void *output,
                            unsigned long frameCount,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData ) 
{
    /* Cast data passed through stream to our structure. */
    paTestData *data = (paTestData*)userData; 
    float *out = (float*)output;
    unsigned int i;
    (void) input; /* Prevent unused variable warning. */
    
    for( i=0; i<frameCount; i++ )
    {
        *out++ = data->left_phase;  /* left */
        *out++ = data->right_phase;  /* right */
        /* Generate simple sawtooth phaser that ranges between -1.0 and 1.0. */
        data->left_phase += 0.01f;
        /* When signal reaches top, drop back down. */
        if( data->left_phase >= 1.0f ) data->left_phase -= 2.0f;
        /* higher pitch so we can distinguish left and right. */
        data->right_phase += 0.03f;
        if( data->right_phase >= 1.0f ) data->right_phase -= 2.0f;
    }
    return 0;
}                            


void print_info(const PaDeviceInfo * info) {
    printf("Name : %s\n", info->name);
    printf("hostApi: %d\n", info->hostApi);
    printf("maxInputChannels: %d\n", info->maxInputChannels);
    printf("maxOutputChannels: %d\n", info->maxOutputChannels);
    printf("defaultLowInputLatency: %f\n", info->defaultLowInputLatency);
    printf("defaultLowOutputLatency: %f\n", info->defaultLowOutputLatency);
    printf("defaultHighInputLatency: %f\n", info->defaultHighInputLatency);
    printf("defaultHighInputLatency: %f\n", info->defaultHighInputLatency);
    printf("defaultSampleRate: %f\n", info->defaultSampleRate);
}
void print_portaudio() {
    printf("version=%d\n",version);
    printf("version text=%s\n",version_text);
    printf("host api count=%d\n", host_api_count);
    printf("default host api=%d\n", default_host_api);
    printf("device count=%d\n", device_count);
    printf("default input=%d\n", default_input);
    printf("default output=%d\n", default_output);    
}
int main() {
    InitAudio();
    print_portaudio();
    const PaDeviceInfo * info;
    for(size_t i = 0; i < device_count; i++)
    {
        info = Pa_GetDeviceInfo(i);
        print_info(info);
    }
    int err = Pa_OpenDefaultStream( &stream, 0,2,paFloat32,44100,256,StreamCallback,&data);
    if(err != paNoError) {
        printf("Failed to open stream\n");
        exit(-1);
    }
    err = Pa_StartStream(stream);
    if(err != paNoError) {
        printf("Failed to start stream\n");
        exit(-1);
    }
    Sleep(5000);
    err = Pa_StopStream(stream);
    if(err != paNoError) {
        printf("Failed to stop stream\n");
        exit(-1);
    }
    TerminateAudio();
}