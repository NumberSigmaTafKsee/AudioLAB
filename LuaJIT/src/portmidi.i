%module luapm 
%{
#include <portmidi.h>
#include <porttime.h>
#include <pthread.h>

%}

%include <portmidi.h>
%include <porttime.h>
%include "stdint.i"
%include "lua_fnptr.i"

%inline
%{
SWIGLUA_REF note_on  = {0,0};
SWIGLUA_REF note_off = {0,0};
SWIGLUA_REF control_change  = {0,0};
SWIGLUA_REF program_change  = {0,0};
SWIGLUA_REF polyphonic_key_pressure={0,0};
SWIGLUA_REF channel_pressure  = {0,0};
SWIGLUA_REF pitch_bend  = {0,0};
SWIGLUA_REF realtime_clock = {0,0};
SWIGLUA_REF start_sequence = {0,0};
SWIGLUA_REF continue_sequence = {0,0};
SWIGLUA_REF stop_sequence = {0,0};
SWIGLUA_REF active_sensing = {0,0};
SWIGLUA_REF system_reset = {0,0};
SWIGLUA_REF system_exclusive = {0,0};

void set_note_on_func(SWIGLUA_REF f)
{
    note_on = f; 
    /*
    lua_State *L = lua_newthread(f.L);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    note_on.L = L;
    */
}
void set_note_off_func(SWIGLUA_REF f)
{
    note_off = f;    
    /*
    lua_State *L = lua_newthread(f.L);        
    luaL_ref(L,LUA_REGISTRYINDEX);
    note_off.L = L;
    */
}
void set_control_change_func(SWIGLUA_REF f)
{
    control_change = f;    
}
void set_program_change_func(SWIGLUA_REF f)
{
    program_change = f;    
}
void set_polyphonic_key_pressure_func(SWIGLUA_REF f)
{
    polyphonic_key_pressure = f;    
}

void set_channel_pressure_func(SWIGLUA_REF f)
{
    channel_pressure = f;    
}

void set_pitch_bend_func(SWIGLUA_REF f)
{
    pitch_bend = f;    
}

void set_realtime_clock_func(SWIGLUA_REF f)
{
    realtime_clock = f;    
}

void set_start_sequence_func(SWIGLUA_REF f)
{
    start_sequence = f;    
}
void set_continue_sequence_func(SWIGLUA_REF f)
{
    continue_sequence = f;    
}
void set_stop_sequence_func(SWIGLUA_REF f)
{
    stop_sequence = f;    
}
void set_active_sensing_func(SWIGLUA_REF f)
{
    active_sensing = f;    
}
void set_system_reset_func(SWIGLUA_REF f)
{
    system_reset = f;    
}
void set_system_exclusive_func(SWIGLUA_REF f)
{
    system_exclusive = f;    
}

int midi_channel = 0;

PmStream * pm_midi_input = NULL;
PmStream * pm_midi_output = NULL;
pthread_mutex_t Lock;

size_t GetNumDevices()
{
    return Pm_CountDevices();
}

const char* GetDeviceName(size_t i)
{
    const PmDeviceInfo * pm = Pm_GetDeviceInfo(i);
    return pm->name;
    
}

static void process_midi(PtTimestamp timestamp, void * userData)
{
    PmError result;
    PmEvent buffer;
    int channel;
    lua_State *L;

    pthread_mutex_lock(&Lock);
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
                        
            if( msg == 0x90 && note_on.L != 0)
            {                
                // note on
                swiglua_ref_get(&note_on);	
                lua_pushnumber(note_on.L,channel);
                lua_pushnumber(note_on.L,data1);
                lua_pushnumber(note_on.L,data2);
                if(lua_pcall(note_on.L,3,0,0) != 0)
                {
                    printf("%s", lua_tostring(note_on.L,-1));
                }
            }
            else if( msg == 0x80 && note_off.L != 0)
            {
                // note off                
                swiglua_ref_get(&note_off);	
                lua_pushnumber(note_off.L,channel);
                lua_pushnumber(note_off.L,data1);
                lua_pushnumber(note_off.L,data2);
                if(lua_pcall(note_off.L,3,0,0) != 0)
                {
                    printf("%s", lua_tostring(note_off.L,-1));
                }

            }
            else if(msg == 0xA0)
            {
                // polyphonic pressure
            }
            else if(msg == 0xB0)
            {
                // control change
            }
            else if(msg == 0xC0)
            {
                // program change
            }
            else if(msg == 0xD0)
            {
                // channel pressure
            }
            else if(msg == 0xE0)
            {
                // pitchbend
            }
            else if(status == 0x79)
            {
                // reset all conrollers
            }
            else if(status == 0x7A)
            {
                // local control
            }
            else if(status == 0x7B)
            {
                // all notes off
            }
            else if(status == 0x7C)
            {
                // omni off
            }
            else if(status == 0x7D)
            {
                // omni on
            }
            else if(status == 0x7E)
            {
                // mono mode
            }
            else if(status == 0x7F)
            {
                // poly mode
            }
            else if(status == 0xF8)
            {
                // clock
            }
            else if(status == 0xFA)
            {
                // start sequences
            }
            else if(status == 0xFB)
            {
                // continue sequence
            }
            else if(status == 0xFC)
            {
                // stop sequence
            }
            else if(status == 0xFE)
            {
                // active sense
            }
            else if(status == 0xFF)
            {
                // system reset
            }
            else if(status == 0xF1)
            {
                // midi timing code
            }
            else if(status == 0xF2)
            {
                // song position
            }
            else if(status == 0xF3)
            {
                // song select
            }
            else if(status == 0xF6)
            {
                // tune request
            }
            else if(status == 0xF0)
            {
                // system exclusive
            }
        }
    } while(result);
    pthread_mutex_unlock(&Lock);
}

void InitMidi(int channel, int input_id, int output_id)
{
    const PmDeviceInfo *in_info,*out_info;
    midi_channel = channel;
    int output_buffer_size = 100;
    int latency = 0;

    Pt_Start(10,&process_midi,0);
    Pm_Initialize();
    pthread_mutex_init(&Lock,NULL);

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
    pthread_mutex_destroy(&Lock);
}

void LockMidi()
{
    pthread_mutex_lock(&Lock);
}
void UnlockMidi()
{
    pthread_mutex_unlock(&Lock);
}

%}