#pragma once

#include <iostream>
#include <future>

#include <SDL2/SDL.h>

#include "WAVReader.h"
#include "../window/Animator.h"

#define RUN(CMD) if(!CMD) { \
                    printf("|\n| [%d] PortAudio error: %s\n|\n", __LINE__, SDL_GetError()); \
                    return; \
                 }

struct AudioData{
    SDL_AudioSpec* spec;

    Uint8* dataPos;
    Uint8* currentPos;

    Uint32 dataLength;
    Uint32 currentLength;
};

class Player{
public:
    Player(Animator* animator);

    ~Player();

    void playSound(const char* filePath);
    void stopSound();
    void restartSound();
    void continueSound();

    static double value;
    static double freq;

    static Animator* animator;
    static void handleSound(bool* playing, AudioData* sound);
    
    void loadData();
private:
    AudioData* soundData;
    const char* filePath;
    bool playing;

    std::future<void> futureCall;
};
