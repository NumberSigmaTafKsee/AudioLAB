# AudioLAB
 Generic LV2 Plugins

# There are a tremenedous amount of them
* It is going to take a while to connect the LV2 frames to them all
* They will be highly detailed to do any theory possible

# MIDI CC
* any midi cc use the same format as the DSI Prophet 08/mopho/tetra 

# AnalogSVF
* Will pirkins equations
* mono audio channel
* midi input , cc103 = cutoff, cc104 = resonance
* midi output (for control)
* built in adsr
* built in lfo
* cv control ports for all parameters

# Moog Ladder Filter
* From Github repo
* mono audio
* midi input (note,controls), cc103=cutoff, cc104=resonance
* midi output (to cc etc)
* built in adsr
* built in lfo to cutoff, q
* cv control ports for all parameters

# RBJFilter
* all equations of RBJ
* mono and stereo
* midi input can set the keyboard frequency, cc103=cutoff,cc104=q,cc105=bandwidth,cc106=slope
* midi output to control synthesizer
* cv ports for everybody

# ZolzerFilter
* Equations from DAFX book

# AR
* Attack Release envelope
* Midi note on
* Midi output to control etc
* CV Envelope output
* Audio in * envelope

# ASR
* Attack Sustain Release envelope
* Same as AR but sustains

# ADSR
* Attack Decay Sustain Release Envelope
* audio input * envelope
* cv envelope output
* cv envelope modulation input
* logic trigger on/off

# StereoAnalogSVF
* 2 channel SVF
* Does not have all outputs like the mono

# StereoMoogLadder
* 2 channel Moog Ladder Filter

# BlitOscillator
* Some blit saw and square tooth
* mix saw, square and impulse
* integrator input cv
* differentiator output cv
* comparator wave output cv
* impulse train cv

# QuadraphonicBlit
* 4 x Blit oscillator
* mixed or paralle

# OctaphonicBlit
* 8 x Blit Oscillator
* mixed or parallel

# Blit2Oscillator
* a better blit
* Saw, Triangle, Square

# DPW
* Saw
* Square
* Triangle

# minBlip
* MinBLEP Saw,Square,Triangle

# VCO
* Vectorized PolyBlep
* many waveforms

# QuadVCO
* 4x VCO mixed or parallel

# OctaVCO
* 8x VCO mixed or parallel

# VCO9000
* 8x8 VCO terrain
* Matrix sequencing equation theories

# LFO9000
* Deep Tesla LFO Frequency Terrain

# IIR Filters
* Butterworths
* Bessels
* Chebyshevs
* Ellipticals
* DspFilters
* Spuce Filters
* IPP Filters
