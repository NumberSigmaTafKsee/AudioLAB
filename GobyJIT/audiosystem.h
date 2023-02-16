
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <portaudio.h>
#include <stdbool.h>
#include <portmidi.h>
#include <porttime.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct midimsg
{
    int status,data1,data2,msg,channel;    
    struct midimsg * next;
} 
MidiMsg;


typedef void (*MIDI_CALLBACK)(MidiMsg * msg);

typedef int (*AUDIO_CALLBACK)( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData );

void set_note_on_func(MIDI_CALLBACK f);
void set_note_off_func(MIDI_CALLBACK f);
void set_control_change_func(MIDI_CALLBACK f);
void set_program_change_func(MIDI_CALLBACK f);
void set_polyphonic_key_pressure_func(MIDI_CALLBACK f);
void set_channel_pressure_func(MIDI_CALLBACK f);
void set_pitch_bend_func(MIDI_CALLBACK f);
void set_realtime_clock_func(MIDI_CALLBACK f);
void set_start_sequence_func(MIDI_CALLBACK f);
void set_continue_sequence_func(MIDI_CALLBACK f);
void set_stop_sequence_func(MIDI_CALLBACK f);
void set_active_sensing_func(MIDI_CALLBACK f);
void set_system_reset_func(MIDI_CALLBACK f);
void set_system_exclusive_func(MIDI_CALLBACK f);
void set_local_control_func(MIDI_CALLBACK f);
void set_all_notes_off_func(MIDI_CALLBACK f);
void set_omni_off_func(MIDI_CALLBACK f);
void set_omni_on_func(MIDI_CALLBACK f);
void set_mono_mode_func(MIDI_CALLBACK f);
void set_poly_mode_func(MIDI_CALLBACK f);
void set_clock_func(MIDI_CALLBACK f);
void set_midi_timing_code_func(MIDI_CALLBACK f);
void set_reset_all_controllers_func(MIDI_CALLBACK f);
void set_song_position_func(MIDI_CALLBACK f);
void set_select_func(MIDI_CALLBACK f);
void set_tuning_request_func(MIDI_CALLBACK f);


void Init();
void InitMidi();
size_t GetNumMidiDevices();
const char* GetMidiDeviceName(size_t i);
void LockMidi();
void UnlockMidi();
void InitMidiDevice(int channel, int input_id, int output_id);
void StopMidi();

void set_repl_func(void (*r)());
void set_audio_func(AUDIO_CALLBACK f);
void set_callback_func(AUDIO_CALLBACK f);
int  InitAudio(int sample_rate, int frames_per_second);
int  StopAudio();
void RunAudio();
int  GetNumAudioDevices();
const char* GetAudioDeviceName(size_t i);
int InitAudioDevice(int output_device_number, int input_device_number, size_t num_channels, int sample_rate, int frames_per_second);
int InitAudio(int sample_rate, int frames_per_second);
int StopAudio();
void Stop();

void SendMidiMessage(int msg, int p1, int p2);
void SendSysex(char * msg);

#ifdef __cplusplus
}
#endif