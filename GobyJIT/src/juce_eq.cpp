/*
  ==============================================================================
    This file was auto-generated!
    It contains the basic framework code for a JUCE plugin processor.
  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "Filter.hpp"

using namespace std;

//==============================================================================
GraphicEqAudioProcessor::GraphicEqAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       ), state(*this, nullptr, Identifier("GraphicEQParameters"), createParameterLayout())
#endif
{
    //Add the gain parameters here in the Constructor so they can be automated in the DAW. Sets the name of the variable and how it appears in the DAW automation, as well as the range of automatable numbers (min/max gain value). Add for each gain slider
    
    //LSF
    /*
    addParameter(LSFgain = new AudioParameterFloat("LSFgain",        //Parameter ID
                                                   "Low Shelf Gain", //DAW parameter name
                                                   -24.0,            //Minimum value
                                                   12.0,             //Maximum value
                                                   0.0));            //Default value
     
    
    //PKF1
    addParameter(PKF1gain = new AudioParameterFloat("PKF1gain",      //Parameter ID
                                                   "250 Hz Gain",    //DAW parameter name
                                                   -24.0,            //Minimum value
                                                   12.0,             //Maximum value
                                                   0.0));            //Default value
     
    
    //PKF2
    addParameter(PKF2gain = new AudioParameterFloat("PKF2gain",       //Parameter ID
                                                    "500 Hz Gain",    //DAW parameter name
                                                    -24.0,            //Minimum value
                                                    12.0,             //Maximum value
                                                    0.0));            //Default value
    
    //PKF3
    addParameter(PKF3gain = new AudioParameterFloat("PKF3gain",       //Parameter ID
                                                    "1 kHz Gain",     //DAW parameter name
                                                    -24.0,            //Minimum value
                                                    12.0,             //Maximum value
                                                    0.0));            //Default value
    
    //PKF4
    addParameter(PKF4gain = new AudioParameterFloat("PKF4gain",       //Parameter ID
                                                    "2 kHz Gain",     //DAW parameter name
                                                    -24.0,            //Minimum value
                                                    12.0,             //Maximum value
                                                    0.0));            //Default value
    
    //PKF5
    addParameter(PKF5gain = new AudioParameterFloat("PKF5gain",       //Parameter ID
                                                    "4 kHz Gain",     //DAW parameter name
                                                    -24.0,            //Minimum value
                                                    12.0,             //Maximum value
                                                    0.0));            //Default value
    
    //PKF6
    addParameter(PKF6gain = new AudioParameterFloat("PKF6gain",       //Parameter ID
                                                    "8 kHz Gain",     //DAW parameter name
                                                    -24.0,            //Minimum value
                                                    12.0,             //Maximum value
                                                    0.0));            //Default value
    
    //HSF
    addParameter(HSFgain = new AudioParameterFloat("HSFgain",             //Parameter ID
                                                    "High Shelf Gain",    //DAW parameter name
                                                    -24.0,                //Minimum value
                                                    12.0,                 //Maximum value
                                                    0.0));                //Default value
    
    //Creates an instance of each filter object
    lowShelf = new Filter; peak250 = new Filter; peak500 = new Filter; peak1000 = new Filter;
    peak2000 = new Filter; peak4000 = new Filter; peak8000 = new Filter; highShelf = new Filter;
    
    //Create Instances of Smoothing Objects
    lowShelfSmooth = new LinearSmoothedValue<float>; peakOneSmooth = new LinearSmoothedValue<float>;
    peakTwoSmooth = new LinearSmoothedValue<float>; peakThreeSmooth = new LinearSmoothedValue<float>;
    peakFourSmooth = new LinearSmoothedValue<float>; peakFiveSmooth = new LinearSmoothedValue<float>;
    peakSixSmooth = new LinearSmoothedValue<float>; highShelfSmooth = new LinearSmoothedValue<float>;
     */
}

GraphicEqAudioProcessor::~GraphicEqAudioProcessor()
{
    //Delete Filter Objects from memory
    /*
    delete lowShelf; delete peak250; delete peak500; delete peak1000; delete peak2000; delete peak4000; delete peak8000; delete highShelf;
    
    //Delete smoothing objects from memory
    delete lowShelfSmooth; delete peakOneSmooth; delete peakTwoSmooth; delete peakThreeSmooth; delete peakFourSmooth; delete peakFiveSmooth; delete peakSixSmooth; delete highShelfSmooth;
     */
}

AudioProcessorValueTreeState::ParameterLayout GraphicEqAudioProcessor::createParameterLayout(){
    
    std::vector<std::unique_ptr<RangedAudioParameter>> params;
    
    //Identifier  DAW Display  [  Range  ]  IC
    params.push_back(std::make_unique<AudioParameterFloat>("LSFgain", "Low Shelf Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("250gain", "250Hz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("500gain", "500Hz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("1kgain", "1 kHz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("2kgain", "2 kHz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("4kgain", "4 kHz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("8kgain", "8 kHz Gain", -24.f, 12.f, 0.f));
    params.push_back(std::make_unique<AudioParameterFloat>("HSFgain", "High Shelf Gain", -24.f, 12.f, 0.f));

    return {params.begin(), params.end()};
    
}

//==============================================================================
const String GraphicEqAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool GraphicEqAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool GraphicEqAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool GraphicEqAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double GraphicEqAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int GraphicEqAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int GraphicEqAudioProcessor::getCurrentProgram()
{
    return 0;
}

void GraphicEqAudioProcessor::setCurrentProgram (int index)
{
}

const String GraphicEqAudioProcessor::getProgramName (int index)
{
    return {};
}

void GraphicEqAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
//prepareToPlay is called whenever major DAW preferences are changed, such as changing the buffer size, sampling frequency, or other preferences.
//Is is NOT called everytime you press the play button in the DAW
void GraphicEqAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
    
    //Assign Sampling Rate
    Fs = getSampleRate();
    
    //Set the Q and f values for each filter object
    lowShelf.setQ(Q);
    lowShelf.setFreq(125);
    lowShelf.setOmega(125, Fs);
    lowShelf.setAlpha(lowShelf.getOmega(), Q);
    
    peak250.setQ(Q);
    peak250.setFreq(250);
    peak250.setOmega(250, Fs);
    peak250.setAlpha(peak250.getOmega(), Q);
    
    peak500.setQ(Q);
    peak500.setFreq(500);
    peak500.setOmega(500, Fs);
    peak500.setAlpha(peak500.getOmega(), Q);
    
    peak1000.setQ(Q);
    peak1000.setFreq(1000);
    peak1000.setOmega(1000, Fs);
    peak1000.setAlpha(peak1000.getOmega(), Q);
    
    peak2000.setQ(Q);
    peak2000.setFreq(2000);
    peak2000.setOmega(2000, Fs);
    peak2000.setAlpha(peak2000.getOmega(), Q);
    
    peak4000.setQ(Q);
    peak4000.setFreq(4000);
    peak4000.setOmega(4000, Fs);
    peak4000.setAlpha(peak4000.getOmega(), Q);
    
    peak8000.setQ(Q);
    peak8000.setFreq(8000);
    peak8000.setOmega(8000, Fs);
    peak8000.setAlpha(peak8000.getOmega(), Q);
    
    highShelf.setQ(Q);
    highShelf.setFreq(12000);
    highShelf.setOmega(12000, Fs);
    highShelf.setAlpha(highShelf.getOmega(), Q);
    
    //Reset to new Sample Rate and Ramp Length, whenever any major DAW preferences are changed
    /*
    lowShelfSmooth->reset(Fs, 0.02);
    peakOneSmooth->reset(Fs, 0.02);
    peakTwoSmooth->reset(Fs, 0.02);
    peakThreeSmooth->reset(Fs, 0.02);
    peakFourSmooth->reset(Fs, 0.02);
    peakFiveSmooth->reset(Fs, 0.02);
    peakSixSmooth->reset(Fs, 0.02);
    highShelfSmooth->reset(Fs, 0.02);
     */
}

void GraphicEqAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool GraphicEqAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

//Passes in a reference to buffer and midiMessages values. Any updates to them within the process block will also update the buffer and midiMessages variables elsewhere in the program.
void GraphicEqAudioProcessor::processBlock (AudioBuffer<float> &buffer, MidiBuffer &midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (int i = totalNumInputChannels; i < totalNumOutputChannels; ++i){
        buffer.clear (i, 0, buffer.getNumSamples());
    }
    
    
    
    for (int channel = 0; channel < totalNumInputChannels; ++channel)
    {
        //Looping through all the samples in the input buffer
        for (int sample = 0; sample < buffer.getNumSamples(); ++sample){
            //Current sample in our buffer
            input = buffer.getReadPointer(channel)[sample];
            
            // Filter instantaneous gain values
            float lsfGain = *state.getRawParameterValue("LSFgain");
            float gain250 = *state.getRawParameterValue("250gain");
            float gain500 = *state.getRawParameterValue("500gain");
            float gain1k = *state.getRawParameterValue("1kgain");
            float gain2k = *state.getRawParameterValue("2kgain");
            float gain4k = *state.getRawParameterValue("4kgain");
            float gain8k = *state.getRawParameterValue("8kgain");
            float hsfGain = *state.getRawParameterValue("HSFgain");
            
            // Perform gain smoothing
            lowShelfSmooth = (1.f - alpha)*lsfGain + alpha*lowShelfSmooth;
            peak250Smooth = (1.f - alpha)*gain250 + alpha*peak250Smooth;
            peak500Smooth = (1.f - alpha)*gain500 + alpha*peak500Smooth;
            peak1000Smooth = (1.f - alpha)*gain1k + alpha*peak1000Smooth;
            peak2000Smooth = (1.f - alpha)*gain2k + alpha*peak2000Smooth;
            peak4000Smooth = (1.f - alpha)*gain4k + alpha*peak4000Smooth;
            peak8000Smooth = (1.f - alpha)*gain8k + alpha*peak8000Smooth;
            highShelfSmooth = (1.f - alpha)*hsfGain + alpha*highShelfSmooth;
            
            // Filtering methods
            u1 = lowShelf.lowShelfFilter(input, lowShelfSmooth, channel);
            u2 = peak250.peakFilter(u1, peak250Smooth, channel);
            u3 = peak500.peakFilter(u2, peak500Smooth, channel);
            u4 = peak1000.peakFilter(u3, peak1000Smooth, channel);
            u5 = peak2000.peakFilter(u4, peak2000Smooth, channel);
            u6 = peak4000.peakFilter(u5, peak4000Smooth, channel);
            u7 = peak8000.peakFilter(u6, peak8000Smooth, channel);
            output = highShelf.highShelfFilter(u7, highShelfSmooth, channel);
            
            //Output the final processed sound after all filters have been applied
            buffer.getWritePointer(channel)[sample] = output;
        }
    }
}

//==============================================================================
bool GraphicEqAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* GraphicEqAudioProcessor::createEditor()
{
    return new GraphicEqAudioProcessorEditor (*this);
}

//==============================================================================
void GraphicEqAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
    
    /*
    ScopedPointer<XmlElement> xml(new XmlElement ("GraphicEQAudioProcessorParams"));
    xml -> setAttribute("LSFgain", (double) *LSFgain);
    xml -> setAttribute("PKF1gain", (double) *PKF1gain);
    xml -> setAttribute("PKF2gain", (double) *PKF2gain);
    xml -> setAttribute("PKF3gain", (double) *PKF3gain);
    xml -> setAttribute("PKF4gain", (double) *PKF4gain);
    xml -> setAttribute("PKF5gain", (double) *PKF5gain);
    xml -> setAttribute("PKF6gain", (double) *PKF6gain);
    xml -> setAttribute("HSFgain", (double) *HSFgain);
    
    copyXmlToBinary(*xml, destData);
     */
    
    auto currentState = state.copyState();
    std::unique_ptr<XmlElement> xml(currentState.createXml());
    copyXmlToBinary(*xml, destData);
    
}

void GraphicEqAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // this method will restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
    /*
    ScopedPointer<XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
    
    //This if statement accounts for the condition when we do not have anything saved in memory.
    // if(xmlState != nullptr) returns true if the xmlState memory is not empty, i.e. we have something saved there
    if(xmlState != nullptr){
        
        //This if statement is looking for the string reference "GraphicEQAudioProcessorParams". This MUST equal the same string as the getStateInformation method in order for it to properly find the correct data in memory
        
        if(xmlState->hasTagName("GraphicEQAudioProcessorParams")){
            *LSFgain = xmlState -> getDoubleAttribute("LSFgain", 0.0);
            *PKF1gain = xmlState -> getDoubleAttribute("PKF1gain", 0.0);
            *PKF2gain = xmlState -> getDoubleAttribute("PKF2gain", 0.0);
            *PKF3gain = xmlState -> getDoubleAttribute("PKF3gain", 0.0);
            *PKF4gain = xmlState -> getDoubleAttribute("PKF4gain", 0.0);
            *PKF5gain = xmlState -> getDoubleAttribute("PKF5gain", 0.0);
            *PKF6gain = xmlState -> getDoubleAttribute("PKF6gain", 0.0);
            *HSFgain = xmlState -> getDoubleAttribute("HSFgain", 0.0);
        }
    }
     */
    
    std::unique_ptr<XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
    if (xmlState && xmlState->hasTagName(state.state.getType())){
        state.replaceState(ValueTree::fromXml(*xmlState));
    }
    
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new GraphicEqAudioProcessor();
}