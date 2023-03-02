     
#include "rtaudio_c.h"
#include "rtmidi_c.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char*  GetAudioVersion();
void GetCompiledApis();
unsigned int GetNumCompiledApis();
const char* GetApiName(size_t i);
const char* GetApiDisplayName(size_t i);


static const rtaudio_api_t *compiled_apis;
static rtaudio_t audio;
unsigned int num_apis = 0;

const char* GetAudioVersion() {
    return rtaudio_version();
}

unsigned int GetNumCompiledApis() {
    return rtaudio_get_num_compiled_apis();
}

void GetCompiledApis() {
    num_apis = GetNumCompiledApis();
    compiled_apis = rtaudio_compiled_api();
}
const char* GetApiName(size_t i) {
    return rtaudio_api_name(compiled_apis[i]);
}
const char* GetApiDisplayName(size_t i) {
    return rtaudio_api_display_name(compiled_apis[i]);
}


int DeviceCount() {
    return rtaudio_device_count(audio);
}
void InitAudio(size_t api) {
    audio =  rtaudio_create(compiled_apis[api]);
}
void ShutdownAudio() {
    rtaudio_destroy(audio);
}

int main() {
    GetCompiledApis();    
    printf("num apis=%d\n", num_apis);
    for(size_t i = 0; i < num_apis; i++) {
        printf("api : %s\n", GetApiName(i));
        printf("name: %s\n", GetApiDisplayName(i));
    }
}