#include "IO/SndFile.hpp"
#include "kissfft/kissfft.hh"
#include <iostream>

int main()
{
    SoundWave::SndFileReader file("test.wav");
    float buffer[1024];
    std::cout << file.size() << std::endl;
    file.read(1024,buffer);
}