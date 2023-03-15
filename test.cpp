#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "SoundObject.hpp"
#include "FX/LiquidNeuron.hpp"
#include "FX/Amplifiers.hpp"
#include "FX/LiquidMoog.hpp"

using namespace Liquid;
int main()
{
	LiquidNeuron * neuron = new LiquidNeuron(44100,0.01);
	LiquidMoog * moog = new LiquidMoog(44100.0,1000.0,0.5);
	std::cout << "ok\n";
}
