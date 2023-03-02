#include "audiosystem.h"
/////////////////////////////////////////////
// MIDI
/////////////////////////////////////////////

int midi_channel = 0;
MidiMsg * midi_queue = NULL;
PmStream * pm_midi_input = NULL;
PmStream * pm_midi_output = NULL;
pthread_mutex_t Lock;

// most of these will probably never get used
MIDI_CALLBACK note_on  = NULL;
MIDI_CALLBACK note_off = NULL;
MIDI_CALLBACK control_change  = NULL;
MIDI_CALLBACK program_change  = NULL;
MIDI_CALLBACK polyphonic_key_pressure=NULL;
MIDI_CALLBACK channel_pressure  = NULL;
MIDI_CALLBACK pitch_bend  = NULL;
MIDI_CALLBACK realtime_clock = NULL;
MIDI_CALLBACK start_sequence = NULL;
MIDI_CALLBACK continue_sequence = NULL;
MIDI_CALLBACK stop_sequence = NULL;
MIDI_CALLBACK active_sensing = NULL;
MIDI_CALLBACK system_reset = NULL;
MIDI_CALLBACK system_exclusive = NULL;
MIDI_CALLBACK local_control = NULL;
MIDI_CALLBACK all_notes_off = NULL;
MIDI_CALLBACK omni_off = NULL;
MIDI_CALLBACK omni_on = NULL;
MIDI_CALLBACK mono_mode = NULL;
MIDI_CALLBACK poly_mode = NULL;
MIDI_CALLBACK midi_clock = NULL;
MIDI_CALLBACK midi_timing_code = NULL;
MIDI_CALLBACK reset_all_controllers = NULL;
MIDI_CALLBACK song_position = NULL;
MIDI_CALLBACK song_select = NULL;
MIDI_CALLBACK tuning_request = NULL;


/////////////////////////////////////////////
// Audio
/////////////////////////////////////////////


AUDIO_CALLBACK audio_func;
AUDIO_CALLBACK callback_func;
PaStreamParameters outputParameters;
PaStreamParameters inputParameters;
PaStream *stream;
int isAudioRunning=false;
pthread_t repl_thread;


void dummy_repl() {

}
void (*repl_func)() = dummy_repl;

void set_repl_func(void (*r)()) { repl_func = r; }

void set_note_on_func(MIDI_CALLBACK f)
{
    note_on = f;         
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    note_on = L;        
    */
}
void set_note_off_func(MIDI_CALLBACK f)
{
    note_off = f;            
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    note_off = L;        
    */
}
void set_control_change_func(MIDI_CALLBACK f)
{
    control_change = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    control_change = L;        
    */
}
void set_program_change_func(MIDI_CALLBACK f)
{
    program_change = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    program_change = L;        
    */
}
void set_polyphonic_key_pressure_func(MIDI_CALLBACK f)
{
    polyphonic_key_pressure = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    polyphonic_key_pressure = L;        
    */
}

void set_channel_pressure_func(MIDI_CALLBACK f)
{
    channel_pressure = f;   
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    channel_pressure = L;         
    */
}

void set_pitch_bend_func(MIDI_CALLBACK f)
{
    pitch_bend = f; 
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    pitch_bend = L;           
    */
}

void set_realtime_clock_func(MIDI_CALLBACK f)
{
    realtime_clock = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    realtime_clock = L;        
    */
}
void set_start_sequence_func(MIDI_CALLBACK f)
{
    start_sequence = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    start_sequence = L;        
    */
}
void set_continue_sequence_func(MIDI_CALLBACK f)
{
    continue_sequence = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    continue_sequence = L;        
    */
}
void set_stop_sequence_func(MIDI_CALLBACK f)
{
    stop_sequence = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    stop_sequence = L;        
    */
}
void set_active_sensing_func(MIDI_CALLBACK f)
{
    active_sensing = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    active_sensing = L;        
    */
}
void set_system_reset_func(MIDI_CALLBACK f)
{
    system_reset = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    system_reset = L;        
    */
}
void set_system_exclusive_func(MIDI_CALLBACK f)
{
    system_exclusive = f;    
    /*
    lua_State *L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    system_exclusive = L;        
    */
}
void set_local_control_func(MIDI_CALLBACK f)
{
    local_control = f;    
}
void set_all_notes_off_func(MIDI_CALLBACK f)
{
    all_notes_off = f;    
}
void set_omni_off_func(MIDI_CALLBACK f)
{
    omni_off = f;    
}
void set_omni_on_func(MIDI_CALLBACK f)
{
    omni_on = f;    
}
void set_mono_mode_func(MIDI_CALLBACK f)
{
    mono_mode = f;    
}
void set_poly_mode_func(MIDI_CALLBACK f)
{
    poly_mode = f;    
}
void set_clock_func(MIDI_CALLBACK f)
{
    midi_clock = f;    
}
void set_midi_timing_code_func(MIDI_CALLBACK f)
{
    midi_timing_code = f;    
}
void set_reset_all_controllers_func(MIDI_CALLBACK f)
{
    reset_all_controllers = f;    
}
void set_song_position_func(MIDI_CALLBACK f)
{
    song_position = f;    
}
void set_select_func(MIDI_CALLBACK f)
{
    song_select = f;    
}
void set_tuning_request_func(MIDI_CALLBACK f)
{
    tuning_request = f;    
}


size_t GetNumMidiDevices()
{
    return Pm_CountDevices();
}
const char* GetMidiDeviceName(size_t i)
{
    const PmDeviceInfo * pm = Pm_GetDeviceInfo(i);
    return pm->name;
    
}
void LockMidi()
{
    while( pthread_mutex_lock(&Lock) != 0);
}
void UnlockMidi()
{
    while(pthread_mutex_unlock(&Lock) != 0);
}
    
MidiMsg* NewMessage(int status, int data1, int data2, int msg, int channel) {
    MidiMsg * p = (MidiMsg*)calloc(1,sizeof(MidiMsg));
    p->status = status;
    p->data1  = data1;
    p->data2  = data2;
    p->msg    = msg;
    p->channel = channel;
    p->next   = NULL;
    return p;
}
void AddMessage(MidiMsg * head, MidiMsg * last) {
    MidiMsg * p = head;
    if(p == NULL) return;
    while(p->next != NULL) {
        p = p->next;
    }
    p->next = last;
    last->next = NULL;
}


void Callback(MIDI_CALLBACK f, MidiMsg * msg ) {
    if(f == NULL) return;
    f(msg);
}
void ExecQueue(MidiMsg * msgs) 
{
    MidiMsg * p = msgs, *t;
    while(p != NULL) 
    {
        int status = p->status;
        int data1  = p->data1;
        int data2  = p->data2;
        int msg    = p->msg & 0xF0;
        int channel= p->msg & 0x0F;
        
        if( msg == 0x90 && note_on != 0)
        {                
            // note on
            Callback(note_on,p);
        }
        else if( msg == 0x80 && note_off != 0)
        {
            // note off                
            Callback(note_off,p);
        }
        else if(msg == 0xA0)
        {
            // polyphonic pressure
            Callback(polyphonic_key_pressure,p);
        }
        else if(msg == 0xB0)
        {
            // control change
            Callback(control_change,p);
        }
        else if(msg == 0xC0)
        {
            // program change        
            Callback(program_change,p);
        }
        else if(msg == 0xD0)
        {
            // channel pressure
            Callback(channel_pressure,p);
            
        }
        else if(msg == 0xE0)
        {
            // pitchbend
            Callback(pitch_bend,p);
        }
        else if(status == 0x79)
        {
            // reset all conrollers
            Callback(reset_all_controllers,p);
        }
        else if(status == 0x7A)
        {
            // local control
            Callback(local_control,p);
        }
        else if(status == 0x7B)
        {
            // all notes off
            Callback(all_notes_off,p);
        }
        else if(status == 0x7C)
        {
            // omni off
            Callback(omni_off,p);
        }
        else if(status == 0x7D)
        {
            // omni on
            Callback(omni_on,p);
        }
        else if(status == 0x7E)
        {
            // mono mode
            Callback(mono_mode,p);
        }
        else if(status == 0x7F)
        {
            // poly mode
            Callback(poly_mode,p);
        }
        else if(status == 0xF8)
        {
            // clock
            Callback(midi_clock,p);
        }
        else if(status == 0xFA)
        {
            // start sequences
            Callback(start_sequence,p);
        }
        else if(status == 0xFB)
        {
            // continue sequence
            Callback(continue_sequence,p);
        }
        else if(status == 0xFC)
        {
            // stop sequence
            Callback(stop_sequence,p);
        }
        else if(status == 0xFE)
        {
            // active sense
            Callback(active_sensing,p);
        }
        else if(status == 0xFF)
        {
            // system reset
            Callback(system_reset,p);
        }
        else if(status == 0xF1)
        {
            // midi timing code
            Callback(midi_timing_code,p);
        }
        else if(status == 0xF2)
        {
            // song position
            Callback(song_position,p);
        }
        else if(status == 0xF3)
        {
            // song select
            Callback(song_select,p);
        }
        else if(status == 0xF6)
        {
            // tune request
            Callback(tuning_request,p);
        }
        else if(status == 0xF0)
        {
            // system exclusive
            Callback(system_exclusive,p);
        }
        t = p->next;
        free(p);
        p = t;
    }
} 


void RunQueue() {    
    MidiMsg * msgs = midi_queue;
    midi_queue = NULL;
    ExecQueue(msgs);        
}

static void process_midi(PtTimestamp timestamp, void * userData)
{
    PmError result;
    PmEvent buffer;
    int channel;
   
    LockMidi();
    do
    {
        result = Pm_Poll(pm_midi_input);        
        if(result)
        {
            int status,data1,data2,msg;
            if(Pm_Read(pm_midi_input, &buffer, 1) == pmBufferOverflow)
                continue;
            status = Pm_MessageStatus(buffer.message);
            data1  = Pm_MessageData1(buffer.message);
            data2  = Pm_MessageData2(buffer.message);
            channel = status & 0x0F;
            msg = status & 0xF0;   
            MidiMsg * pMsg = NewMessage(status,data1,data2,msg,channel);
            if(midi_queue == NULL) midi_queue = pMsg;
            else AddMessage(midi_queue,pMsg);
        }
    }        
    while(result);   
    UnlockMidi();
    
}



void InitMidi() {
    pthread_mutex_init(&Lock,NULL);
    LockMidi();
}
void InitMidiDevice(int channel, int input_id, int output_id)
{
    const PmDeviceInfo *in_info,*out_info;
    midi_channel = channel;
    int output_buffer_size = 100;
    int latency = 0;

    
    Pt_Start(10,&process_midi,0);
    Pm_Initialize();
    
    if(input_id > -1)
    {
        in_info  = Pm_GetDeviceInfo(input_id);
        Pm_OpenInput(&pm_midi_input, input_id, NULL, output_buffer_size, NULL, NULL);
    }
    if(output_id > -1)
    {
        out_info  = Pm_GetDeviceInfo(output_id);
        Pm_OpenOutput(&pm_midi_output, output_id, NULL, output_buffer_size, NULL, NULL, latency);
    }
}

void StopMidi()
{    
    Pt_Stop();
    if(pm_midi_input)  Pm_Close(pm_midi_input);
    if(pm_midi_output) Pm_Close(pm_midi_output);    
}

void set_audio_func(AUDIO_CALLBACK f)
{   
    /*     
    audio_func = f;    
    Lua = f;
    L = lua_newthread(f);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    audio_func = L;    
    */
    audio_func = f;
}

void set_callback_func(AUDIO_CALLBACK f)
{
    callback_func = f;      
}


void *REPL(void * args) {    
    while(1)
    {
        repl_func();
    }
    return args;    
}


void RunAudio()
{
    int r = pthread_create(&repl_thread, NULL, REPL,NULL);
    UnlockMidi();
    while(1)    
    {               
        LockMidi();     
        RunQueue();
        UnlockMidi();            
        Pa_Sleep(10);
    }
}

void StreamFinished( void* userData )
{

}


int GetNumAudioDevices()
{
    return Pa_GetDeviceCount();
}

const char* GetAudioDeviceName(size_t i)
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
              audio_func,
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
    exit(-1);
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
              audio_func,
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
    exit(-1);
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
    return 0;
error:
    Pa_Terminate();
    fprintf( stderr, "An error occurred while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    exit(-1);
    return err;

}


void Stop()
{
    LockMidi();
    StopMidi();
    StopAudio();    
    pthread_mutex_destroy(&Lock);
}
void Init()
{
    InitMidi();
    Pa_Initialize();
}