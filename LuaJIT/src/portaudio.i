%module luapa
%{ 
#include <portaudio.h>
#include <stdbool.h>
%}
%include "stdint.i"
%include "lua_fnptr.i"
%include <portaudio.h>


int InitAudio(int sample_rate, int frames_per_second);
int StopAudio();
void RunAudio();
float* float_new(size_t size);
void float_free(float * p);
float float_get(float * p, size_t index);
void float_set(float *p, size_t index, float value);

%inline 
%{

float* float_new(size_t size)
{
    return calloc(sizeof(float),size);
}
void float_free(float * p)
{
    free(p);
}
float float_get(float * p, size_t index)
{
    return p[index];
}
void float_set(float *p, size_t index, float value)
{
    p[index] = value;
}


SWIGLUA_REF audio_func={0,0};
SWIGLUA_REF callback_func={0,0};
lua_State * L;
PaStreamParameters outputParameters;
PaStreamParameters inputParameters;
PaStream *stream;
int isAudioRunning=false;
int i;

/*
typedef struct _audio_t
{
SWIGLUA_REF audio_func={0,0};
SWIGLUA_REF callback_func={0,0};
lua_State * L;
PaStreamParameters outputParameters;
PaStreamParameters inputParameters;
PaStream *stream;
int isAudioRunning=false;
}
Audio;
*/
void set_audio_func(SWIGLUA_REF f)
{
    audio_func = f;
    L = lua_newthread(f.L);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    audio_func.L = L;
}

void set_callback_func(SWIGLUA_REF f)
{
    callback_func = f;    
}

void RunAudio()
{
    while(1)    
    {
        if(callback_func.ref != 0)
        {
            swiglua_ref_get(&callback_func);	
            if(lua_pcall(callback_func.L,0,0,0) != 0) 
            {
                printf("%s", lua_tostring(L,-1));
            }    
        }
        Pa_Sleep(10);
    }
}
static int luaCallback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{
    swig_type_info * info = SWIGTYPE_p_float;
    
    if(L == NULL) { return -1; }
    isAudioRunning=true;
    swiglua_ref_get(&audio_func);	
    SWIG_NewPointerObj(L,inputBuffer,info,0);
    SWIG_NewPointerObj(L,outputBuffer,info,0);
    lua_pushinteger(L,framesPerBuffer);
	if(lua_pcall(L,3,0,0) != 0) 
    {
        printf("%s", lua_tostring(L,-1));
    }
    isAudioRunning=false;
    return paContinue;
}   

static void StreamFinished( void* userData )
{

}

int GetNumDevices()
{
    return Pa_GetDeviceCount();
}
const char* GetDeviceName(size_t i)
{
    const PaDeviceInfo* di = Pa_GetDeviceInfo(i);
    return di->name;
}


int InitAudioDevice(int output_device_number, int input_device_number, size_t num_channels, int sample_rate, int frames_per_second)
{
    PaError err;
    err = Pa_Initialize();    
    
    if( err != paNoError ) goto error;

    if(output_device_number > -1)
    {        
        outputParameters.device = output_device_number;
        if (outputParameters.device == paNoDevice) {
            fprintf(stderr,"Error: No default output device.\n");
            goto error;
        }
        outputParameters.channelCount = num_channels;       /* stereo output */
        outputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
        outputParameters.suggestedLatency = 0.001; //Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
        outputParameters.hostApiSpecificStreamInfo = NULL;        
    }
    if(input_device_number > -1)
    {        
        inputParameters.device = input_device_number;
        if (inputParameters.device == paNoDevice) {
            fprintf(stderr,"Error: No default output device.\n");
            goto error;
        }
        inputParameters.channelCount = num_channels;       /* stereo output */
        inputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
        inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowOutputLatency;
        inputParameters.hostApiSpecificStreamInfo = NULL;        
    }

    err = Pa_OpenStream(
              &stream,
              input_device_number > -1? &inputParameters : NULL, /* no input */
              output_device_number > -1? &outputParameters : NULL,
              sample_rate,
              frames_per_second,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              luaCallback,
              NULL );
              
    if( err != paNoError ) goto error;

    printf("Start\n") ;
    //err = Pa_SetStreamFinishedCallback( stream, &StreamFinished );
    //if( err != paNoError ) goto error;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto error;

    
    return err;
error:
    Pa_Terminate();
    fprintf( stderr, "An error occurred while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return err;
}

int InitAudio(int sample_rate, int frames_per_second)
{      
    PaError err;
    err = Pa_Initialize();    
    printf("Init\n");
    if( err != paNoError ) goto error;
   
    
    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto error;
    }
    outputParameters.channelCount = 2;       /* stereo output */
    outputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    
    err = Pa_OpenStream(
              &stream,
              NULL, /* no input */
              &outputParameters,
              sample_rate,
              frames_per_second,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              luaCallback,
              NULL );
    if( err != paNoError ) goto error;

    printf("Start\n") ;
    //err = Pa_SetStreamFinishedCallback( stream, &StreamFinished );
    //if( err != paNoError ) goto error;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto error;

    
    return err;
error:
    Pa_Terminate();
    fprintf( stderr, "An error occurred while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return err;
}


int StopAudio()
{
    PaError err;
    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    Pa_Terminate();
    printf("Test finished.\n");
error:
    Pa_Terminate();
    fprintf( stderr, "An error occurred while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return err;

}


%}