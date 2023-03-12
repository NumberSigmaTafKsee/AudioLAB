# AudioLAB
 Generic LV2 Plugins

# I have the theory
* I have the power to make it take my course

# Current Plugin
* AnalogSVF
* MoogLadderFilters
* AmpClipper

# AnalogSVF
* Just your kind of svf
* From the Will Pirkin Theory Society
* Has a built in ADSR envelope that can be midi trigger
* Has a built in LFO
* It will work in frequency or it can be set to control-voltage
* It will have a built in LuaJIT interpreter for control-midi-matrix theory and tap alchemy
* It should be possible to modulate it per sample by the FIAXY equations
* It can be oversampled 2x,4x,6x,8x

# MoogLadderFilters
* All of them are from github : houvilainen,stilson,musicdsp,krajeski,microtrack,don phillips runge kut, more
* Has built in ADSR and LFO
* It will have a built in LuaJIT interpreter for control-midi-matrix theory and tap alchemy
* It will work in frequency or it can be set to control-voltage
* It can be oversampled 2x,4x,6x,8x

# AmpClipperCircuit
* It is something I found which doesn't have alot of practical purpose at this time for anything
* It is supposed to be a diode clipper but it is like a diode filter instead
* It has significance as it demonstrate the theory of a non-linear memory algorithm despite it wasn't supposed to

# LiquidSVF
* It is the AnalogSVF but it uses liquid neurons for the coefficients

# AdalineSVF
* It uses adaptive adaline neuron coefficients

# FIRSvf
* It uses frequency sampling of the state variable filter

# FIRMoog
* It uses frequency sampling of the moog filters

# IIRFilters
* general purpose IIR filters

# DSPFilters
* general purpose DSPfilters

# SpuceFilters
* general purpose Filters

# CppFilters
* Generatel purpose biquad filters

# MultiMassberg
* General purpose multibands massberg

# ZolzerFilter
* Biquads from DAFX books

# RbjFilter
* The standard RBJ equations

# PolyRoot
* It will sovle any polynomials in the frequency domain and create a SoS Biquad cascade for you
* You can inject transformation into the S or Z domain
* It can work in S or Z

# IPPTransferFunc
* It uses the IIR filter in Intel IPP for any arbtrrary transfer function

# IPPBiquad
* It uses the IIR filter for coefficients in biquad sos cascade form

# OctopusIIR
* It will let you run GNU Octave script to generate IIR filters for ZPK or as A,B
* Can use the functions in signal or can design anything in ZPK or however you want
* It must fit in Biquad or it must fit in TransferFunction
* It can not use Octave Tf2Sos in realtime it is not fast enough to do that
* Generally you can just do this offline but habing it in a plugin can be interesting thing to do

# OctopusFIR
* It will let you run GNU octave script to generate FIR filter coefficients if you want to
* If the GNU octave function is to slow you can not do this inreal time
* firpm may not be fast enough for realtime so you put it into a file instead

# Plotonic
* when you push the button it will draw the sample in gnuplot 

# BlitOscillators
* It is the same thing adapted from the Stk algorithms
* You can control it with MIDI
* You can let it run free as frequency 
* It has a built in ADSR that can be bipass
* It can be phase modulated
* It can be oversamples 2x,4x,6x,8x
