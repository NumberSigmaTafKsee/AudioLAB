#pragma once


#include <iostream>
#include "algorithmfactory.h"
#include "essentiamath.h"
#include "pool.h"

using namespace std;
using namespace essentia;
using namespace standard;

enum ofxAAAlgorithm{
    RMS,
    ENERGY,
    POWER,
    PITCH_FREQ,
    PITCH_CONFIDENCE,
    PITCH_SALIENCE,
    INHARMONICITY,
    HFC,
    CENTROID,
    SPECTRAL_COMPLEXITY,
    DISSONANCE,
    ROLL_OFF,
    ODD_TO_EVEN,
    STRONG_PEAK,
    STRONG_DECAY,
    SPECTRUM,
    MEL_BANDS,
    MFCC,
    HPCP,
    MULTI_PITCHES,
    PITCH_SALIENCE_FUNC_PEAKS,
    TRISTIMULUS,
    ONSETS
    
};

//OSC INDEXES
#define POWER_OSC_IDX 0
#define PITCH_FREQ_OSC_IDX 1
#define PITCH_CONF_OSC_IDX 2
#define PITCH_SALIENCE_OSC_IDX 3
#define HFC_OSC_IDX 4
#define CENTROID_OSC_IDX 5
#define SPEC_COMP_OSC_IDX 6
#define INHARMONICITY_OSC_IDX 7
#define DISSONANCE_OSC_IDX 8
#define ROLLOFF_OSC_IDX 9
#define ODDTOEVEN_OSC_IDX 10
#define ONSETS_OSC_IDX 11



class ofxAABaseAlgorithm{
    
public:
    
    
    void init();
    
    float getValue();
    float getValueDb();
    
    ///Gets the value normalized from 0 to maxEstimatedValue with Clamping
    float getValueNormalized();
    
    float getValueNormalized(float min, float max, bool doClamp=true);
    float getValueDbNormalized(float min, float max, bool doClamp=true);
    float getSmoothedValue(float smthAmnt);
    
    ///Gets the value normalized and smoothed from 0 to maxEstimatedValue with Clamping
    float getSmoothedValueNormalized(float smthAmnt);
    
    float getSmoothedValueNormalized(float smthAmnt, float min, float max, bool doClamp=true);
    float getSmoothedValueDbNormalized(float smthAmnt, float min, float max, bool doClamp=true);
    
    float getMaxEstimatedValue();
    
    bool getIsActive();
    
    void setActive(bool state);
    void setValueZero();
    void setMaxEstimatedValue(float value);
    
    void compute();
    
    void castValueToFloat();
    
    void deleteAlgorithm();
    
    Algorithm* algorithm;
    Real realValue;
    
    
protected:
    
    bool isActivated;
    
private:
    
    float floatValue;
    float smoothedFloatValue;
    float smoothedNormFloatValue;
    float smoothedNormFloatValueDb;
    
    float maxEstimatedValue;
    
};
//---------------------------------------------------------------------
class ofxAAOneVectorOutputAlgorithm : public ofxAABaseAlgorithm{

public:
    
    void init();
    
    void initAndAssignSize(int size, int initValues);
    
    void assignFloatValuesSize(int size, int val);
    
    void castValuesToFloat(bool logarithmic);
    
    void updateLogRealValues();
    
    int getBinsNum();
    vector<float>& getValues();
    vector<float>& getSmoothedValues(float smthAmnt);
    
    vector<Real> realValues;
    vector<Real> logRealValues;
    
private:
    
    vector<float> floatValues;
    vector<float> smoothedFloatValues;
    
};


//---------------------------------------------------------------------
struct SalienceFunctionPeak{
    float bin;//cents
    float value;
    
    SalienceFunctionPeak(){
        bin = 0.0;
        value = 0.0;
    }
};

class ofxAAPitchSalienceFunctionPeaksAlgorithm : public ofxAABaseAlgorithm{
public:
    
    void init();
    
    void castValuesToFloat();
    
    vector<SalienceFunctionPeak>& getPeaks();
    vector<SalienceFunctionPeak>& getSmoothedPeaks(float smthAmnt);
    
    ///change to maxPeaks!!!!
    void setMaxPeaksNum(int maxBins){maxPeaksNum = maxBins;}
    
    vector<Real> realSalienceBins;
    vector<Real> realSalienceValues;
    
private:
    vector<SalienceFunctionPeak> peaks;
    vector<SalienceFunctionPeak> smoothedPeaks;
    
    bool limitPeaksNum;
    
    int maxPeaksNum;
};


//---------------------------------------------------------------------
class ofxAACartToPolAlgorithm : public ofxAABaseAlgorithm{
public:
    vector<Real> magnitudes;
    vector<Real> phases;
    
};
//---------------------------------------------------------------------
//class used for SpectralPeaks & HarmonicPeaks
class ofxAAPeaksAlgorithm : public ofxAABaseAlgorithm{
public:
    vector<Real> frequencies;
    vector<Real> magnitudes;
    
};

//---------------------------------------------------------------------
class ofxAAFftAlgorithm : public ofxAAOneVectorOutputAlgorithm{
    public:
    vector<complex<Real> > fftRealValues;

};
//---------------------------------------------------------------------
class ofxAAPitchDetectAlgorithm : public ofxAABaseAlgorithm{

public:
    
    void init();
    
    void castValuesToFloat();
    
    float getPitchValue();
    float getPitchValueNormalized();
    float getConfidenceValue();
    
    float getSmoothedPitchValue(float smthAmnt);
    float getSmoothedPitchValueNormalized(float smthAmnt);
    float getSmoothedConfidenceValue(float smthAmnt);
    
    float getMaxPitchEstimatedValue(){return pitchMaxEstimatedValue;}
    
    void setMaxPitchEstimatedValue(float value);
    
    Real pitchRealVal, confidenceRealVal;

private:
    
    float pitchFloatVal, confidenceFloatVal;
    
    float smoothedPitchFloatValue;
    float smoothedNormPitchFloatValue;
    float smoothedConfidenceFloatValue;
    
    float pitchMaxEstimatedValue;

};

//---------------------------------------------------------------------
class ofxAATuningFrequencyAlgorithm : public ofxAABaseAlgorithm{
    
public:
    
    void castValuesToFloat();
    
    float getFreqValue();
    float getCentsValue();
    
    Real freqRealVal,   centsRealVal;
    
    private:
    float freqFloatVal, centsFloatVal;
    
};

class ofxAAMultiPitchKlapuriAlgorithm : public ofxAABaseAlgorithm{
    
public:
    
    void setup(ofxAAPitchSalienceFunctionPeaksAlgorithm* saliencePeaksPtr,
               ofxAAOneVectorOutputAlgorithm* spectrumPtr,
               int sampleRate);
    void compute();
    
    vector<float>& getPitches();
    
    int frequencyToCentBin(Real frequency);
    Real getWeight(int centBin, int harmonicNumber);
    
private:
    
    ofxAAPitchSalienceFunctionPeaksAlgorithm* _saliencePeaksAlgthm;
    ofxAAOneVectorOutputAlgorithm* _spectrumAlgthm;
    
    Real _binResolution;
    int _binsInSemitone;
    std::vector<Real> _centSpectrum;
    Real _sampleRate;
    int _binsInOctave;
    Real _referenceTerm;
    Real _referenceFrequency;
    int _numberBins;
    Real _centToHertzBase;
    int _numberHarmonicsMax;
    vector<Real> nearestBinWeights;
    
    vector<float> pitches;
    
};

enum OnsetsTimeThresholdMode{
    TIME_BASED,
    BUFFER_NUM_BASED
};

class ofxAAOnsetsAlgorithm : public ofxAABaseAlgorithm{

public:
    
    void setup(int bufferSize);
    void compute();
    void castValuesToFloat();
    void evaluate();
    
    void reset();
    
    bool getValue(){return _value;}
    float getOnsetSilenceThreshold(){return silenceThreshold;}
    float getOnsetTimeThreshold(){return timeThreshold;}
    float getOnsetAlpha(){return alpha;}
    
    void setOnsetSilenceThreshold(float val){silenceThreshold=val;}
    void setOnsetAlpha(float val){alpha=val;}
    void setOnsetTimeThreshold(float ms){timeThreshold = ms;}
    void setOnsetBufferNumThreshold(int buffersNum){bufferNumThreshold = buffersNum;}
    void setUseTimeThreshold(bool doUse){usingTimeThreshold = doUse;}
    void setOnsetTimeThresholdsMode(OnsetsTimeThresholdMode mode){onsetsMode = mode;}

    ofxAABaseAlgorithm onsetHfc;
    ofxAABaseAlgorithm onsetComplex;
    ofxAABaseAlgorithm onsetFlux;
    
private:
    
    bool _value;//isOnset
    
    bool onsetBufferEvaluation (Real iDetectHfc, Real iDetectComplex, Real iDetectFlux);
    bool onsetTimeThresholdEvaluation();
    bool onsetBufferNumThresholdEvaluation();//framebased threshold eval.
    
    int detecBufferSize;
    vector<vector<Real> > detections;
    vector<Real> detection_sum;
    Real hfc_max, complex_max, flux_max;

    Real silenceThreshold, alpha;
    bool addHfc, addComplex, addFlux;
    
    bool usingTimeThreshold;
    float timeThreshold;
    float lastOnsetTime;
    int bufferNumThreshold;
    int lastOnsetBufferNum;
    
    OnsetsTimeThresholdMode onsetsMode;
    int bufferCounter;
    

};

//----------------------------------
//for scaling values:
#define DB_MIN -6
#define DB_MAX 0
#define MFCC_MAX_ESTIMATED_VALUE 300.0
//----------------------------------
//Vars for algorithm creation:
#define MELBANDS_BANDS_NUM 24
#define DCT_COEFF_NUM 13
#define PITCH_SALIENCE_FUNC_NUM 10
#define TRISTIMULUS_BANDS_NUM 3

#define HPCP_SIZE 12
#define HPCP_MIN_FREQ 40.0//hz
#define HPCP_MAX_FREQ 5000.0//hz

#define PEAKS_MAXPEAKS_NUM 10000
#define PEAKS_MIN_FREQ 40.0//hz
#define PEAKS_MAX_FREQ 5000.0//hz
//----------------------------------

class ofxAudioAnalyzerUnit
{

public:
    
    ofxAudioAnalyzerUnit(int sampleRate, int bufferSize){
        setup(sampleRate, bufferSize);
    }
    ~ofxAudioAnalyzerUnit(){
        exit();
    }
    
    void setup(int sampleRate, int bufferSize);
    void analyze(const vector<float> &  inBuffer);
    void exit();
    
    void resetOnsets();
    
    int getSampleRate() {return _samplerate;}
    int getBufferSize() {return _framesize;}
    
    void setActive(ofxAAAlgorithm algorithm, bool state);
    void setMaxEstimatedValue(ofxAAAlgorithm algorithm, float value);
    void setOnsetsParameters(float alpha, float silenceTresh, float timeTresh, bool useTimeTresh = true);
    void setSalienceFunctionPeaksParameters(int maxPeaks);
    
    float getValue(ofxAAAlgorithm algorithm, float smooth=0.0, bool normalized=false);
    vector<float>& getValues(ofxAAAlgorithm algorithm, float smooth=0.0);
    vector<SalienceFunctionPeak>& getPitchSaliencePeaksRef(float smooth=0.0);
    bool getIsActive(ofxAAAlgorithm algorithm);
    bool getOnsetValue();
    
    int getBinsNum(ofxAAAlgorithm algorithm);
    float getMaxEstimatedValue(ofxAAAlgorithm algorithm);
    
    ofxAAOnsetsAlgorithm* getOnsetsAlgorithmPtr(){return &onsets;}
    
    int   getPitchFreqAsMidiNote(float smooth=0.0);
    string getPitchFreqAsNoteName(float smooth=0.0);
    
   
private:
    //Utils:
    int pitchToMidi(float pitch);
    string midiToNoteName(int midiNote);
    
    //Algorithms:
    vector<Real> audioBuffer;

    //algorithms with return value func
    ofxAABaseAlgorithm rms;
    ofxAABaseAlgorithm energy;
    ofxAABaseAlgorithm power;
    ofxAAPitchDetectAlgorithm pitchDetect;
    ofxAABaseAlgorithm pitchSalience;
    //ofxAATuningFrequencyAlgorithm tuning;
    ofxAABaseAlgorithm inharmonicity;
    ofxAABaseAlgorithm hfc;
    ofxAABaseAlgorithm centroid;
    ofxAABaseAlgorithm spectralComplex;
    ofxAABaseAlgorithm dissonance;
    ofxAABaseAlgorithm rollOff;
    ofxAABaseAlgorithm oddToEven;
    ofxAABaseAlgorithm strongPeak;
    ofxAABaseAlgorithm strongDecay;
    
    ofxAAOneVectorOutputAlgorithm spectrum;
    ofxAAOneVectorOutputAlgorithm melBands;
    ofxAAOneVectorOutputAlgorithm dct;//MFCC
    ofxAAOneVectorOutputAlgorithm hpcp;
    ofxAAOneVectorOutputAlgorithm tristimulus;
    //MultiPitch
    ofxAAPitchSalienceFunctionPeaksAlgorithm pitchSalienceFunctionPeaks;
    ofxAAMultiPitchKlapuriAlgorithm multiPitchKlapuri;
    
    //Algorithms for internal functionality:
    ofxAAOneVectorOutputAlgorithm dcremoval;
    ofxAAOneVectorOutputAlgorithm window;
    ofxAAFftAlgorithm fft;
    ofxAACartToPolAlgorithm cartesian2polar;
    ofxAAPeaksAlgorithm spectralPeaks;
    ofxAAPeaksAlgorithm harmonicPeaks;
    ofxAAOneVectorOutputAlgorithm pitchSalienceFunction;
    ofxAAOnsetsAlgorithm onsets;
    
    //--------
    int _samplerate;
    int _framesize;
    
    int hopsize;
    int zeropadding;
    Real framerate;
    

    
};


class ofxAudioAnalyzer{
 
 public:
    
    void setup(int sampleRate, int bufferSize, int channels);
    void reset(int sampleRate, int bufferSize, int channels);
    void analyze(const std::vector<float> & inBuffer);
    void exit();
    
    int getSampleRate() {return _samplerate;}
    int getBufferSize() {return _buffersize;}
    int getChannelsNum(){return _channels;}
    
    ///Gets value of single output  Algorithms.
    ///\param algorithm
    ///\param channel: starting from 0 (for stereo setup, 0 and 1)
    ///\param smooth: smoothing amount. 0.0=non smoothing, 1.0=fixed value
    float getValue(ofxAAAlgorithm algorithm, int channel, float smooth=0.0, bool normalized=false);
    
    ///Gets values of vector output Algorithms.
    ///\param algorithm
    ///\param channel: starting from 0 (for stereo setup, 0 and 1)
    ///\param smooth: smoothing amount. 0.0=non smoothing, 1.0=fixed value
    vector<float>& getValues(ofxAAAlgorithm algorithm, int channel, float smooth=0.0);
    
    ///Gets the array of pitch salience function peaks: bin/cents & value
    vector<SalienceFunctionPeak>& getSalienceFunctionPeaks(int channel, float smooth=0.0);
    
    ///Returns if there is an onset in the speciefied channel.
    bool getOnsetValue(int channel);
    
    ///Gets the state of an algorithm
    bool getIsActive(int channel, ofxAAAlgorithm algorithm);

    ///Pointers for the audio analyzing units.
    ///Use very carefully!
    vector<ofxAudioAnalyzerUnit*>& getChannelAnalyzersPtrs(){return channelAnalyzerUnits;}
    
    ///Resets onsetsr detections buffer
    void resetOnsets(int channel);
    
    ///Activates and deactives algorithms.
    void setActive(int channel, ofxAAAlgorithm algorithm, bool state);
    
    ///Set max estimated values for algorithms that are not normalized
    void setMaxEstimatedValue(int channel, ofxAAAlgorithm algorithm, float value);
    
    ///Sets onsets detection parameters
    ///\param channel: starting from 0 (for stereo setup, 0 and 1)
    ///\param alpha: the proportion of the mean included to reject smaller peaks--filters very short onsets
    ///\param silenceThreshold: the thhreshold for silence
    ///\param timeThreshold: time threshold in ms.
    ///\param useTimeThreshold: use or note the time threshold.
    void setOnsetsParameters(int channel, float alpha, float silenceTresh, float timeTresh, bool useTimeTresh = true);
    
    void setSalienceFunctionPeaksParameters(int channel, int maxPeaks);
    
    

 private:
    
    int _samplerate;
    int _buffersize;
    int _channels;
    
    vector<ofxAudioAnalyzerUnit*> channelAnalyzerUnits;
    
    
};

#pragma mark - Main funcs
//--------------------------------------------------------------
void ofxAudioAnalyzerUnit::setup(int sampleRate, int bufferSize){
    
    #pragma mark -Init variables:
    
    _framesize = bufferSize;
    hopsize = _framesize/2;
    _samplerate = sampleRate;
    zeropadding = 0;
    framerate = (Real) _samplerate / hopsize;
    
    audioBuffer.resize(bufferSize);

    
    #pragma mark -Init algorithms:
    
    onsets.setup(bufferSize);
    
    rms.init();
    energy.init();
    power.init();
    pitchDetect.init();
    pitchSalience.init();
    inharmonicity.init();
    hfc.init();
    centroid.init();
    spectralComplex.init();
    dissonance.init();
    rollOff.init();
    oddToEven.init();
    strongPeak.init();
    strongDecay.init();
    
    spectrum.initAndAssignSize((bufferSize/2)+1, 0.0);
    melBands.initAndAssignSize(MELBANDS_BANDS_NUM,0);
    dct.initAndAssignSize(DCT_COEFF_NUM, 0);
    hpcp.initAndAssignSize(HPCP_SIZE, 0);
    pitchSalienceFunction.initAndAssignSize(PITCH_SALIENCE_FUNC_NUM, 0.0);
    tristimulus.initAndAssignSize(TRISTIMULUS_BANDS_NUM, 0.0);
    
    dcremoval.init();
    window.init();
    fft.init();
    cartesian2polar.init();
    spectralPeaks.init();
    harmonicPeaks.init();
    
    //-:Set Max Estimated Values for Non Normalized Algorithms
    //default values set from testing with white noise.
    energy.setMaxEstimatedValue(100.0);
    pitchDetect.setMaxPitchEstimatedValue(4186.0);//C8
    hfc.setMaxEstimatedValue(2000.0);
    spectralComplex.setMaxEstimatedValue(20.0);
    centroid.setMaxEstimatedValue(11000.0);
    rollOff.setMaxEstimatedValue(sampleRate/2);
    oddToEven.setMaxEstimatedValue(10.0);
    strongPeak.setMaxEstimatedValue(20.0);
    strongDecay.setMaxEstimatedValue(100.0);
    
    //------Not very useful...
    pitchSalienceFunctionPeaks.init();
    setActive(PITCH_SALIENCE_FUNC_PEAKS, false);
    multiPitchKlapuri.init();
    setActive(MULTI_PITCHES, false);
    //------------------
    
    essentia::init();

    #pragma mark -Create algorithms
    AlgorithmFactory& factory = AlgorithmFactory::instance();
    
    rms.algorithm = factory.create("RMS");
    energy.algorithm = factory.create("Energy");
    power.algorithm = factory.create("InstantPower");

    dcremoval.algorithm = factory.create("DCRemoval", "sampleRate", _samplerate);

    window.algorithm = factory.create("Windowing",
                                 "type", "hann",
                                 "zeroPadding", zeropadding);

    fft.algorithm = factory.create("FFT", "size", _framesize);
    cartesian2polar.algorithm = factory.create("CartesianToPolar");

    onsets.onsetHfc.algorithm     = factory.create("OnsetDetection", "method", "hfc", "sampleRate", _samplerate);
    onsets.onsetComplex.algorithm = factory.create("OnsetDetection", "method", "complex", "sampleRate", _samplerate);
    onsets.onsetFlux.algorithm    = factory.create("OnsetDetection", "method", "flux", "sampleRate", _samplerate);

    spectrum.algorithm = factory.create("Spectrum",
                                   "size", _framesize);

    hfc.algorithm = factory.create("HFC", "sampleRate", _samplerate);

    pitchSalience.algorithm = factory.create("PitchSalience", "sampleRate",_samplerate);

    centroid.algorithm = factory.create("Centroid", "range", _samplerate/2);

    spectralComplex.algorithm = factory.create("SpectralComplexity", "sampleRate", _samplerate);
    
    dissonance.algorithm = factory.create("Dissonance");
    
    rollOff.algorithm = factory.create("RollOff",
                                       "sampleRate", _samplerate);
    oddToEven.algorithm = factory.create("OddToEvenHarmonicEnergyRatio");
    strongPeak.algorithm = factory.create("StrongPeak");
    strongDecay.algorithm = factory.create("StrongDecay",
                                           "sampleRate", _samplerate);
    
    spectralPeaks.algorithm = factory.create("SpectralPeaks",
                                "maxPeaks", PEAKS_MAXPEAKS_NUM,
                                "magnitudeThreshold", 0.00001,
                                "minFrequency", PEAKS_MIN_FREQ,
                                "maxFrequency", PEAKS_MAX_FREQ,
                                "orderBy", "frequency");


    melBands.algorithm = factory.create("MelBands",
                                        "inputSize", (_framesize/2)+1,
                                        "sampleRate", _samplerate,
                                        "highFrequencyBound", _samplerate/2,
                                        "numberBands", MELBANDS_BANDS_NUM);
    
    dct.algorithm = factory.create("DCT",
                                   "inputSize", MELBANDS_BANDS_NUM,
                                   "outputSize", DCT_COEFF_NUM);

    inharmonicity.algorithm = factory.create("Inharmonicity");

    pitchDetect.algorithm = factory.create("PitchYinFFT",
                                      "frameSize", _framesize,
                                      "sampleRate", _samplerate);

    harmonicPeaks.algorithm = factory.create("HarmonicPeaks");

    hpcp.algorithm = factory.create("HPCP",
                                   "size", HPCP_SIZE,
                                   "referenceFrequency", 440,
                                   "bandPreset", false,
                                   "minFrequency", HPCP_MIN_FREQ,
                                   "maxFrequency", HPCP_MAX_FREQ,
                                   "weightType", "squaredCosine",
                                   "nonLinear", false,
                                   "windowSize", 4.0/3.0);

    
    pitchSalienceFunction.algorithm = factory.create("PitchSalienceFunction");
    pitchSalienceFunctionPeaks.algorithm = factory.create("PitchSalienceFunctionPeaks");
    
    tristimulus.algorithm = factory.create("Tristimulus");
    
    
    #pragma mark -Connect algorithms
    
    //DCRemoval
    dcremoval.algorithm->input("signal").set(audioBuffer);
    dcremoval.algorithm->output("signal").set(dcremoval.realValues);
    //RMS
    rms.algorithm->input("array").set(dcremoval.realValues);
    rms.algorithm->output("rms").set(rms.realValue);
    //Energy
    energy.algorithm->input("array").set(dcremoval.realValues);
    energy.algorithm->output("energy").set(energy.realValue);
    //Power
    power.algorithm->input("array").set(dcremoval.realValues);
    power.algorithm->output("power").set(power.realValue);
    //Window
    window.algorithm->input("frame").set(dcremoval.realValues);
    window.algorithm->output("frame").set(window.realValues);
    //Onsets
    fft.algorithm->input("frame").set(window.realValues);
    fft.algorithm->output("fft").set(fft.fftRealValues);
    cartesian2polar.algorithm->input("complex").set(fft.fftRealValues);
    cartesian2polar.algorithm->output("magnitude").set(cartesian2polar.magnitudes);
    cartesian2polar.algorithm->output("phase").set(cartesian2polar.phases);
    //-
    onsets.onsetHfc.algorithm->input("spectrum").set(cartesian2polar.magnitudes);
    onsets.onsetHfc.algorithm->input("phase").set(cartesian2polar.phases);
    onsets.onsetHfc.algorithm->output("onsetDetection").set(onsets.onsetHfc.realValue);
    //-
    onsets.onsetComplex.algorithm->input("spectrum").set(cartesian2polar.magnitudes);
    onsets.onsetComplex.algorithm->input("phase").set(cartesian2polar.phases);
    onsets.onsetComplex.algorithm->output("onsetDetection").set(onsets.onsetComplex.realValue);
    //-
    onsets.onsetFlux.algorithm->input("spectrum").set(cartesian2polar.magnitudes);
    onsets.onsetFlux.algorithm->input("phase").set(cartesian2polar.phases);
    onsets.onsetFlux.algorithm->output("onsetDetection").set(onsets.onsetFlux.realValue);
    //Spectrum
    spectrum.algorithm->input("frame").set(window.realValues);
    spectrum.algorithm->output("spectrum").set(spectrum.realValues);
    //HFC
    hfc.algorithm->input("spectrum").set(spectrum.realValues);
    hfc.algorithm->output("hfc").set(hfc.realValue);
    //Pitch Salience
    pitchSalience.algorithm->input("spectrum").set(spectrum.realValues);
    pitchSalience.algorithm->output("pitchSalience").set(pitchSalience.realValue);
    //Centroid
    centroid.algorithm->input("array").set(spectrum.realValues);
    centroid.algorithm->output("centroid").set(centroid.realValue);
    //Spectral Complexity
    spectralComplex.algorithm->input("spectrum").set(spectrum.realValues);
    spectralComplex.algorithm->output("spectralComplexity").set(spectralComplex.realValue);
    //Peak detection
    spectralPeaks.algorithm->input("spectrum").set(spectrum.realValues);
    spectralPeaks.algorithm->output("frequencies").set(spectralPeaks.frequencies);
    spectralPeaks.algorithm->output("magnitudes").set(spectralPeaks.magnitudes);
    //HPCP
    hpcp.algorithm->input("frequencies").set(spectralPeaks.frequencies);
    hpcp.algorithm->input("magnitudes").set(spectralPeaks.magnitudes);
    hpcp.algorithm->output("hpcp").set(hpcp.realValues);
    //MelBands
    melBands.algorithm->input("spectrum").set(spectrum.realValues);
    melBands.algorithm->output("bands").set(melBands.realValues);
    //DCT
    dct.algorithm->input("array").set(melBands.logRealValues);
    dct.algorithm->output("dct").set(dct.realValues);
    //PitchDetection
    pitchDetect.algorithm->input("spectrum").set(spectrum.realValues);
    pitchDetect.algorithm->output("pitch").set(pitchDetect.pitchRealVal);
    pitchDetect.algorithm->output("pitchConfidence").set(pitchDetect.confidenceRealVal);
    //HarmonicPeaks
    harmonicPeaks.algorithm->input("frequencies").set(spectralPeaks.frequencies);
    harmonicPeaks.algorithm->input("magnitudes").set(spectralPeaks.magnitudes);
    harmonicPeaks.algorithm->input("pitch").set(pitchDetect.pitchRealVal);
    harmonicPeaks.algorithm->output("harmonicFrequencies").set(harmonicPeaks.frequencies);
    harmonicPeaks.algorithm->output("harmonicMagnitudes").set(harmonicPeaks.magnitudes);
    //Inharmonicity
    inharmonicity.algorithm->input("frequencies").set(harmonicPeaks.frequencies);
    inharmonicity.algorithm->input("magnitudes").set(harmonicPeaks.magnitudes);
    inharmonicity.algorithm->output("inharmonicity").set(inharmonicity.realValue);
    //Dissonance
    dissonance.algorithm->input("frequencies").set(spectralPeaks.frequencies);
    dissonance.algorithm->input("magnitudes").set(spectralPeaks.magnitudes);
    dissonance.algorithm->output("dissonance").set(dissonance.realValue);
    //Pitch Salience Function
    pitchSalienceFunction.algorithm->input("frequencies").set(spectralPeaks.frequencies);
    pitchSalienceFunction.algorithm->input("magnitudes").set(spectralPeaks.magnitudes);
    pitchSalienceFunction.algorithm->output("salienceFunction").set(pitchSalienceFunction.realValues);
    //Pitch Salience Function Peaks
    pitchSalienceFunctionPeaks.algorithm->input("salienceFunction").set(pitchSalienceFunction.realValues);
    pitchSalienceFunctionPeaks.algorithm->output("salienceBins").set(pitchSalienceFunctionPeaks.realSalienceBins);
    pitchSalienceFunctionPeaks.algorithm->output("salienceValues").set(pitchSalienceFunctionPeaks.realSalienceValues);

    //RollOff
    rollOff.algorithm->input("spectrum").set(spectrum.realValues);
    rollOff.algorithm->output("rollOff").set(rollOff.realValue);
    //StrongPeak
    strongPeak.algorithm->input("spectrum").set(spectrum.realValues);
    strongPeak.algorithm->output("strongPeak").set(strongPeak.realValue);
    //StrongDecay
    strongDecay.algorithm->input("signal").set(dcremoval.realValues);
    strongDecay.algorithm->output("strongDecay").set(strongDecay.realValue);
    //OddToEven
    oddToEven.algorithm->input("frequencies").set(harmonicPeaks.frequencies);
    oddToEven.algorithm->input("magnitudes").set(harmonicPeaks.magnitudes);
    oddToEven.algorithm->output("oddToEvenHarmonicEnergyRatio").set(oddToEven.realValue);
    //Tristimulus
    tristimulus.algorithm->input("frequencies").set(harmonicPeaks.frequencies);
    tristimulus.algorithm->input("magnitudes").set(harmonicPeaks.magnitudes);
    tristimulus.algorithm->output("tristimulus").set(tristimulus.realValues);
    
    //MultiPitch Kalpuri:
    multiPitchKlapuri.setup(&pitchSalienceFunctionPeaks, &spectrum, _samplerate);
    
    //-Shutdown factory:
    factory.shutdown();
    
    
}


//--------------------------------------------------------------
void ofxAudioAnalyzerUnit::analyze(const vector<float> & inBuffer){
    
    if(inBuffer.size() != _framesize){
        //ofLogWarning()<<"ofxAudioAnalyzerUnit: buffer requested to analyze size(" <<inBuffer.size()<<")doesnt match the buffer size already set: "<<_framesize;
    }
    
    //Cast of incoming audio buffer to Real
    for (int i=0; i<inBuffer.size();i++){
        audioBuffer[i] = (Real) inBuffer[i];
    }
    
    #pragma mark -Compute Algorithms
    
    dcremoval.compute();
    rms.compute();
    energy.compute();
    power.compute();
    window.compute();
    
    if(onsets.getIsActive()){
        fft.compute();
        cartesian2polar.compute();
        onsets.compute();
    }
    
    //spectrum must always be computed as it is neede for other algorithms
    spectrum.algorithm->compute();
    
    hfc.compute();
    pitchSalience.compute();
    pitchDetect.compute();
    centroid.compute();
    spectralComplex.compute();
    if(melBands.getIsActive()){
        melBands.algorithm->compute();
        if(dct.getIsActive()){
            melBands.updateLogRealValues();
            dct.compute();
        }
    }else{
        dct.setActive(false);//dct needs melBands to be active
    }
    spectralPeaks.compute();
    hpcp.compute();
    
    if (inharmonicity.getIsActive()){
        harmonicPeaks.compute();
        inharmonicity.algorithm->compute();
    }

    dissonance.compute();
    pitchSalienceFunction.compute();
    pitchSalienceFunctionPeaks.compute();
    
    multiPitchKlapuri.compute();
    
    rollOff.compute();
    oddToEven.compute();
    strongPeak.compute();
    
    tristimulus.compute();
    if(dcremoval.realValues[0] != 0.0){
        //the strong decay is not defined for a zero signal
        strongDecay.compute();
    }
    
    
    #pragma mark -Cast results to float
    
    spectrum.castValuesToFloat(true);
    
    rms.castValueToFloat();
    energy.castValueToFloat();
    power.castValueToFloat();
    pitchDetect.castValuesToFloat();
    pitchSalience.castValueToFloat();
    
    melBands.castValuesToFloat(true);
    dct.castValuesToFloat(false);
    hpcp.castValuesToFloat(false);
    hfc.castValueToFloat();
    centroid.castValueToFloat();
    spectralComplex.castValueToFloat();
    inharmonicity.castValueToFloat();
    dissonance.castValueToFloat();
    rollOff.castValueToFloat();
    oddToEven.castValueToFloat();
    strongPeak.castValueToFloat();
    strongDecay.castValueToFloat();
    
    pitchSalienceFunctionPeaks.castValuesToFloat();
    
    tristimulus.castValuesToFloat(false);
    
    onsets.castValuesToFloat();
    onsets.evaluate();
    
 
}

//--------------------------------------------------------------
void ofxAudioAnalyzerUnit::exit(){

    rollOff.deleteAlgorithm();
    oddToEven.deleteAlgorithm();
    strongPeak.deleteAlgorithm();
    strongDecay.deleteAlgorithm();
    tristimulus.deleteAlgorithm();
    
    dcremoval.deleteAlgorithm();
    rms.deleteAlgorithm();
    energy.deleteAlgorithm();
    power.deleteAlgorithm();
    window.deleteAlgorithm();;
    spectrum.deleteAlgorithm();
    spectralPeaks.deleteAlgorithm();;
    pitchDetect.deleteAlgorithm();
    pitchSalience.deleteAlgorithm();
    dissonance.deleteAlgorithm();
    
    melBands.deleteAlgorithm();
    dct.deleteAlgorithm();
    harmonicPeaks.deleteAlgorithm();;
    hpcp.deleteAlgorithm();
    hfc.deleteAlgorithm();
    inharmonicity.deleteAlgorithm();
    centroid.deleteAlgorithm();
    spectralComplex.deleteAlgorithm();
    fft.deleteAlgorithm();;
    cartesian2polar.deleteAlgorithm();;
    onsets.onsetComplex.deleteAlgorithm();;
    onsets.onsetHfc.deleteAlgorithm();;
    onsets.onsetFlux.deleteAlgorithm();;
    pitchSalienceFunction.deleteAlgorithm();
    pitchSalienceFunctionPeaks.deleteAlgorithm();
    
    essentia::shutdown();
    
    //ofLogVerbose()<<"AudioAnalyzer exit";

}

//--------------------------------------------------------------
#pragma mark - Activates
//----------------------------------------------
void ofxAudioAnalyzerUnit::setActive(ofxAAAlgorithm algorithm, bool state){
    
    switch (algorithm) {
        case RMS:
            rms.setActive(state);
            break;
        case ENERGY:
            energy.setActive(state);
            break;
        case POWER:
            power.setActive(state);
            break;
        case PITCH_FREQ:
            pitchDetect.setActive(state);
            break;
        case PITCH_CONFIDENCE:
            pitchDetect.setActive(state);
            break;
        case PITCH_SALIENCE:
            pitchSalience.setActive(state);
            break;
        case INHARMONICITY:
            inharmonicity.setActive(state);
            break;
        case HFC:
             hfc.setActive(state);
            break;
        case CENTROID:
            centroid.setActive(state);
            break;
        case SPECTRAL_COMPLEXITY:
            spectralComplex.setActive(state);
            break;
        case DISSONANCE:
            dissonance.setActive(state);
            break;
        case SPECTRUM:
            //ofLogWarning()<<"ofxAudioAnalyzerUnit: Spectrum Algorithm cant be turned off.";
            break;
        case MEL_BANDS:
            melBands.setActive(state);
            if(state==false)dct.setActive(state);//dct needs melBands to be active.
            break;
        case MFCC:
            //dct needs melBands to be active.
            melBands.setActive(state);
            dct.setActive(state);
            break;
        case HPCP:
            hpcp.setActive(state);
            break;
        case MULTI_PITCHES:
            multiPitchKlapuri.setActive(state);
            break;
        case PITCH_SALIENCE_FUNC_PEAKS:
            pitchSalienceFunction.setActive(state);
            pitchSalienceFunctionPeaks.setActive(state);
            break;
        case ONSETS:
            onsets.setActive(state);
            break;
        case ROLL_OFF:
            rollOff.setActive(state);
            break;
        case ODD_TO_EVEN:
            oddToEven.setActive(state);
            break;
        case STRONG_PEAK:
            strongPeak.setActive(state);
            break;
        case STRONG_DECAY:
            strongDecay.setActive(state);
            break;
        case TRISTIMULUS:
            tristimulus.setActive(state);
            break;
            
        default:
            //ofLogWarning()<<"ofxAudioAnalyzerUnit: wrong algorithm to set active.";
            break;
    }
}



//----------------------------------------------
#pragma mark - Get values
//----------------------------------------------
bool ofxAudioAnalyzerUnit::getIsActive(ofxAAAlgorithm algorithm){
    
    switch (algorithm) {
        case RMS:
            return rms.getIsActive();
            break;
        case ENERGY:
            return energy.getIsActive();
            break;
        case POWER:
            return power.getIsActive();
            break;
        case PITCH_FREQ:
            return pitchDetect.getIsActive();
            break;
        case PITCH_CONFIDENCE:
            return pitchDetect.getIsActive();
            break;
        case PITCH_SALIENCE:
            return pitchSalience.getIsActive();
            break;
        case INHARMONICITY:
            return inharmonicity.getIsActive();
            break;
        case HFC:
            return hfc.getIsActive();
            break;
        case CENTROID:
            return centroid.getIsActive();
            break;
        case SPECTRAL_COMPLEXITY:
            return spectralComplex.getIsActive();
            break;
        case DISSONANCE:
            return dissonance.getIsActive();
            break;
        case SPECTRUM:
            //ofLogWarning()<<"ofxAudioAnalyzerUnit: Spectrum Algorithm cant be turned off.";
            break;
        case MEL_BANDS:
            return melBands.getIsActive();
            break;
        case MFCC:
            return dct.getIsActive();
            break;
        case HPCP:
            return hpcp.getIsActive();
            break;
        case MULTI_PITCHES:
            return multiPitchKlapuri.getIsActive();
            break;
        case PITCH_SALIENCE_FUNC_PEAKS:
            return pitchSalienceFunctionPeaks.getIsActive();
            break;
        case ONSETS:
            return onsets.getIsActive();
            break;
        case ROLL_OFF:
            return rollOff.getIsActive();
            break;
        case ODD_TO_EVEN:
            return oddToEven.getIsActive();
            break;
        case STRONG_PEAK:
            return strongPeak.getIsActive();
            break;
        case STRONG_DECAY:
            return strongDecay.getIsActive();
            break;
        case TRISTIMULUS:
            return tristimulus.getIsActive();
            break;
            
        default:
            //ofLogWarning()<<"ofxAudioAnalyzerUnit: wrong algorithm to get if is active.";
            break;
    }
  return true;
}
//----------------------------------------------
float ofxAudioAnalyzerUnit::getValue(ofxAAAlgorithm algorithm, float smooth, bool normalized){
    
    float r = 0.0;
    
    switch (algorithm) {
        
        case RMS:
            r = smooth ?
                rms.getSmoothedValueDbNormalized(smooth, DB_MIN, DB_MAX):
                rms.getValueDbNormalized(DB_MIN, DB_MAX);
            break;
            
        case ENERGY:
            r = smooth ?
                energy.getSmoothedValueNormalized(smooth):
                energy.getValueNormalized();
            break;
            
        case POWER:
            r = smooth ?
                power.getSmoothedValueDbNormalized(smooth, DB_MIN, DB_MAX):
                power.getValueDbNormalized(DB_MIN, DB_MAX);
            break;
            
        case PITCH_FREQ:
            if (normalized){
                r = smooth ?
                pitchDetect.getSmoothedPitchValueNormalized(smooth):
                pitchDetect.getPitchValueNormalized();
            }else{
                r = smooth ?
                pitchDetect.getSmoothedPitchValue(smooth):
                pitchDetect.getPitchValue();
            }
            break;
            
        case PITCH_CONFIDENCE:
            r = smooth ?
                pitchDetect.getSmoothedConfidenceValue(smooth):
                pitchDetect.getConfidenceValue();
            break;
            
        case PITCH_SALIENCE:
            r = smooth ?
                pitchSalience.getSmoothedValue(smooth):
                pitchSalience.getValue();
            break;
            
        case INHARMONICITY:
            r =  smooth ?
                inharmonicity.getSmoothedValue(smooth):
                inharmonicity.getValue();
            break;
            
        case HFC:
            if (normalized){
                r = smooth ?
                hfc.getSmoothedValueNormalized(smooth):
                hfc.getValueNormalized();
            }else{
                r = smooth ?
                hfc.getSmoothedValue(smooth):
                hfc.getValue();
            }
            break;
            
        case SPECTRAL_COMPLEXITY:
            if (normalized){
                r = smooth ?
                spectralComplex.getSmoothedValueNormalized(smooth):
                spectralComplex.getValueNormalized();
            }else{
                r = smooth ?
                spectralComplex.getSmoothedValue(smooth):
                spectralComplex.getValue();
            }
            break;
            
        case CENTROID:
            if (normalized){
                r = smooth ?
                centroid.getSmoothedValueNormalized(smooth):
                centroid.getValueNormalized();
            }else{
                r = smooth ?
                centroid.getSmoothedValue(smooth):
                centroid.getValue();
            }
            break;
            
        case DISSONANCE:
            r = smooth ?
                dissonance.getSmoothedValue(smooth):
                dissonance.getValue();
            break;

        case ROLL_OFF:
            if (normalized){
                r = smooth ?
                rollOff.getSmoothedValueNormalized(smooth):
                rollOff.getValueNormalized();
            }else{
                r = smooth ?
                rollOff.getSmoothedValue(smooth):
                rollOff.getValue();
            }
            break;
        case ODD_TO_EVEN:
            if (normalized){
                r = smooth ?
                oddToEven.getSmoothedValueNormalized(smooth):
                oddToEven.getValueNormalized();
            }else{
                r = smooth ?
                oddToEven.getSmoothedValue(smooth):
                oddToEven.getValue();
                //limit value, because this algorithm reaches huge values (eg: 3.40282e+38)
                r = std::clamp(r, 0.0f, oddToEven.getMaxEstimatedValue());
            }
            break;
        case STRONG_PEAK:
            if (normalized){
                r = smooth ?
                strongPeak.getSmoothedValueNormalized(smooth):
                strongPeak.getValueNormalized();
            }else{
                r = smooth ?
                strongPeak.getSmoothedValue(smooth):
                strongPeak.getValue();
            }
            break;
        case STRONG_DECAY:
            if (normalized){
                r = smooth ?
                strongDecay.getSmoothedValueNormalized(smooth):
                strongDecay.getValueNormalized();
            }else{
                r = smooth ?
                strongDecay.getSmoothedValue(smooth):
                strongDecay.getValue();
            }
            break;
            
            
        default:
            //ofLogWarning()<<"ofxAudioAnalyzerUnit: wrong algorithm for getting value.";
            break;
    }
    
    return r;
}
//----------------------------------------------
bool ofxAudioAnalyzerUnit::getOnsetValue(){
    return onsets.getValue();
}
//----------------------------------------------
vector<float>& ofxAudioAnalyzerUnit::getValues(ofxAAAlgorithm algorithm, float smooth){
    
    switch (algorithm) {
        
        case SPECTRUM:
            return smooth ? spectrum.getSmoothedValues(smooth) : spectrum.getValues();
            break;
            
        case MEL_BANDS:
            return smooth ? melBands.getSmoothedValues(smooth) : melBands.getValues();
            break;
            
        case MFCC:
            return smooth ? dct.getSmoothedValues(smooth) : dct.getValues();
            break;
            
        case HPCP:
            return smooth ? hpcp.getSmoothedValues(smooth) : hpcp.getValues();
            break;
            
        case MULTI_PITCHES:
            return multiPitchKlapuri.getPitches();
            break;
            
        case TRISTIMULUS:
            return smooth ? tristimulus.getSmoothedValues(smooth) : tristimulus.getValues();
            break;
            
        default:
            //ofLogError()<<"ofxAudioAnalyzerUnit: wrong algorithm for getting values.";
            break;
    }
    return multiPitchKlapuri.getPitches(); //Just trying to avoid werror.
}
//----------------------------------------------
vector<SalienceFunctionPeak>& ofxAudioAnalyzerUnit::getPitchSaliencePeaksRef(float smooth){
    
    return smooth ? pitchSalienceFunctionPeaks.getSmoothedPeaks(smooth) : pitchSalienceFunctionPeaks.getPeaks();
    
//    return pitchSalienceFunctionPeaks.getPeaks();
}
//----------------------------------------------
int ofxAudioAnalyzerUnit::getBinsNum(ofxAAAlgorithm algorithm){
    
    switch (algorithm) {
        
        case SPECTRUM:
            return spectrum.getBinsNum();
            break;
        case MEL_BANDS:
            return melBands.getBinsNum();
            break;
        case MFCC:
            return dct.getBinsNum();
            break;
        case HPCP:
            return hpcp.getBinsNum();
            break;
            
        default:
            //ofLogError()<<"ofxAudioAnalyzerUnit: wrong algorithm for getting bins number.";
            break;
    }
    return 0;
}
//----------------------------------------------
float ofxAudioAnalyzerUnit::getMaxEstimatedValue(ofxAAAlgorithm algorithm){
    
    float r = 0.0;
    
    switch (algorithm) {
            
        case ENERGY:
            r = energy.getMaxEstimatedValue();
            break;
        case PITCH_FREQ:
            r = pitchDetect.getMaxPitchEstimatedValue();
            break;
        case HFC:
            r = hfc.getMaxEstimatedValue();
            break;
        case SPECTRAL_COMPLEXITY:
            r = spectralComplex.getMaxEstimatedValue();
            break;
        case CENTROID:
            r = centroid.getMaxEstimatedValue();
            break;
        case ROLL_OFF:
            r = rollOff.getMaxEstimatedValue();
            break;
        case ODD_TO_EVEN:
            r = oddToEven.getMaxEstimatedValue();
            break;
        case STRONG_PEAK:
            r = strongPeak.getMaxEstimatedValue();
            break;
        case STRONG_DECAY:
            r = strongDecay.getMaxEstimatedValue();
            break;
            
        default:
            //ofLogError()<<"ofxAudioAnalyzerUnit: wrong algorithm for getting max estimated value.";
            break;
    }
    
    return r;
}
//----------------------------------------------
void ofxAudioAnalyzerUnit::setMaxEstimatedValue(ofxAAAlgorithm algorithm, float value){
    
    switch (algorithm) {
            
        case ENERGY:
            energy.setMaxEstimatedValue(value);
            break;
        case PITCH_FREQ:
            pitchDetect.setMaxPitchEstimatedValue(value);
            break;
        case HFC:
            hfc.setMaxEstimatedValue(value);
            break;
        case SPECTRAL_COMPLEXITY:
            spectralComplex.setMaxEstimatedValue(value);
            break;
        case CENTROID:
            centroid.setMaxEstimatedValue(value);
            break;
        case ROLL_OFF:
            rollOff.setMaxEstimatedValue(value);
            break;
        case ODD_TO_EVEN:
            oddToEven.setMaxEstimatedValue(value);
            break;
        case STRONG_PEAK:
            strongPeak.setMaxEstimatedValue(value);
            break;
        case STRONG_DECAY:
            strongDecay.setMaxEstimatedValue(value);
            break;
            
        default:
             //ofLogError()<<"ofxAudioAnalyzerUnit: wrong algorithm for setting max estimated value.";
            break;
    }
}
//----------------------------------------------
void ofxAudioAnalyzerUnit::resetOnsets(){
    onsets.reset();
}
//----------------------------------------------
void ofxAudioAnalyzerUnit::setOnsetsParameters(float alpha, float silenceTresh, float timeTresh, bool useTimeTresh){
    
    onsets.setOnsetAlpha(alpha);
    onsets.setOnsetSilenceThreshold(silenceTresh);
    onsets.setOnsetTimeThreshold(timeTresh);
    onsets.setUseTimeThreshold(useTimeTresh);
}
//----------------------------------------------
void ofxAudioAnalyzerUnit::setSalienceFunctionPeaksParameters(int maxPeaks){
    pitchSalienceFunctionPeaks.setMaxPeaksNum(maxPeaks);
}
//----------------------------------------------
#pragma mark - Utils
//----------------------------------------------
int ofxAudioAnalyzerUnit::getPitchFreqAsMidiNote(float smooth){
    return pitchToMidi(getValue(PITCH_FREQ, smooth));
}
//----------------------------------------------
string ofxAudioAnalyzerUnit::getPitchFreqAsNoteName(float smooth){
    return midiToNoteName(getValue(PITCH_FREQ, smooth));
}
//----------------------------------------------
int ofxAudioAnalyzerUnit::pitchToMidi(float pitch){
    return round (12*log2(pitch/440) + 69);
}
//--------------------------------------------------------------
string ofxAudioAnalyzerUnit::midiToNoteName(int midiNote){
    
    string noteName;
    int mod = midiNote%12;
    switch (mod){
        case 0: noteName = "C";
            break;
        case 1: noteName = "C#";
            break;
        case 2: noteName = "D";
            break;
        case 3: noteName = "D#";
            break;
        case 4: noteName = "E";
            break;
        case 5: noteName = "F";
            break;
        case 6: noteName = "F#";
            break;
        case 7: noteName = "G";
            break;
        case 8: noteName = "G#";
            break;
        case 9: noteName = "A";
            break;
        case 10: noteName = "Bb";
            break;
        case 11: noteName = "B";
            break;
        default:
            break;
            
    }
    return (noteName);
    
}

//-------------------------------------------
#pragma mark - ofxAABaseAlgorithm
//-------------------------------------------
void ofxAABaseAlgorithm::init(){
    
    isActivated = true;
    floatValue = 0.0;
    smoothedFloatValue = 0.0;
    smoothedNormFloatValue = 0.0;
    smoothedNormFloatValueDb = 0.0;

}
//-------------------------------------------
float ofxAABaseAlgorithm::getValue(){
    return floatValue;
}
//-------------------------------------------
float ofxAABaseAlgorithm::getValueDb(){
    //returns floatValue in a logaritmic scale
    //0.000001 to 1 -->  -6 to 0
    return log10(floatValue);
}
//-------------------------------------------
float ofxAABaseAlgorithm::getValueNormalized(){
    //return ofMap(floatValue, 0.0, maxEstimatedValue, 0.0, 1.0, true);
    return std::clamp(floatValue(((1.0f/maxEstimatedValue)),0.0f,1.0f);
}
//-------------------------------------------
float ofxAABaseAlgorithm::getValueNormalized(float min, float max, bool doClamp){
    //return ofMap(floatValue, min, max, 0.0, 1.0, doClamp);
    return( std::clamp(floatValue*(1.0f/(max-min))),0.0f,1.0f);
}
//-------------------------------------------
float ofxAABaseAlgorithm::getValueDbNormalized(float min, float max, bool doClamp){
    //return ofMap(getValueDb(), min, max, 0.0, 1.0, doClamp);
    return std::clamp( getValueDb()*(1.0f/(max-min)),0.0f,1.0f);
}
//-------------------------------------------
float ofxAABaseAlgorithm::getSmoothedValue(float smthAmnt){
    smoothedFloatValue =  smoothedFloatValue * smthAmnt + (1-smthAmnt) * floatValue;
    return smoothedFloatValue;
}
//-------------------------------------------
float ofxAABaseAlgorithm::getSmoothedValueNormalized(float smthAmnt){
    //float normVal = ofMap(floatValue, 0.0, maxEstimatedValue, 0.0, 1.0, true);
    float normVal = std::clamp(floatValue/maxEstimatedValue,0.0f,1.0f);
    smoothedNormFloatValue =  smoothedNormFloatValue * smthAmnt + (1-smthAmnt) * normVal;
    return smoothedNormFloatValue;
}
//-------------------------------------------
float ofxAABaseAlgorithm::getSmoothedValueNormalized(float smthAmnt, float min, float max, bool doClamp){
    //float normVal = ofMap(floatValue, min, max, 0.0, 1.0, doClamp);
    float normVal = std::clamp(floatValue*(1.0f/(max-min)),0.0f,1.0f);
    smoothedNormFloatValue =  smoothedNormFloatValue * smthAmnt + (1-smthAmnt) * normVal;
    return smoothedNormFloatValue;
}
//-------------------------------------------
float ofxAABaseAlgorithm::getSmoothedValueDbNormalized(float smthAmnt, float min, float max, bool doClamp){
    
    //float normVal = ofMap(getValueDb(), min, max, 0.0, 1.0, doClamp);
    float normVal = std::clamp(getValueDb()*(1.0f/(max-min)),0.0f,1.0f);
    smoothedNormFloatValueDb = smoothedNormFloatValueDb * smthAmnt + (1-smthAmnt) * normVal;
    return smoothedNormFloatValueDb;
    
}
//-------------------------------------------
bool ofxAABaseAlgorithm::getIsActive(){
    return isActivated;
}
//-------------------------------------------
float ofxAABaseAlgorithm::getMaxEstimatedValue(){
    return maxEstimatedValue;
}
//-------------------------------------------
void ofxAABaseAlgorithm::setActive(bool state){
    isActivated = state;
}
//-------------------------------------------
void ofxAABaseAlgorithm::setValueZero(){
    floatValue = 0.0;
}
//-------------------------------------------
void ofxAABaseAlgorithm::setMaxEstimatedValue(float value){
    maxEstimatedValue = value;
}
//-------------------------------------------
void ofxAABaseAlgorithm::compute(){
    if(isActivated){
        algorithm->compute();
    }
}
//-------------------------------------------
void ofxAABaseAlgorithm::castValueToFloat(){
    if(isActivated)
        floatValue = (float) realValue;
    else
        floatValue = 0.0;
}
//-------------------------------------------
void ofxAABaseAlgorithm::deleteAlgorithm(){
    delete algorithm;
}
//-------------------------------------------
#pragma mark - ofxAAOneVectorOutputAlgorithm
//-------------------------------------------
void ofxAAOneVectorOutputAlgorithm::init(){
    isActivated = true;
}
//-------------------------------------------
void ofxAAOneVectorOutputAlgorithm::initAndAssignSize(int size, int initValues){
    isActivated = true;
    assignFloatValuesSize(size, initValues);
}
//-------------------------------------------
void ofxAAOneVectorOutputAlgorithm::assignFloatValuesSize(int size, int val){
    floatValues.assign(size, val);
    smoothedFloatValues.assign(size, val);
}
//-------------------------------------------
void ofxAAOneVectorOutputAlgorithm::castValuesToFloat(bool logarithmic){
    
    for (int i=0; i<realValues.size(); i++){
        if(getIsActive()){
            if(logarithmic){
                
                if(realValues[i] == 0.0){
                    floatValues[i] = log10(0.000001);//DB_MIN
                }else{
                    floatValues[i] = log10((float) realValues[i]);
                }
                
            }else{
                floatValues[i] = (float) realValues[i];
            }
        }else{
            if(logarithmic){
                floatValues[i] = log10(0.000001);//DB_MIN
            }else{
                floatValues[i] = 0.0;
            }
        }
        
    }
    
}
//-------------------------------------------
void ofxAAOneVectorOutputAlgorithm::updateLogRealValues(){
    logRealValues.resize(realValues.size());
    for (int i=0; i<realValues.size(); ++i)
        logRealValues[i] = amp2db(realValues[i]);
    
}
//-------------------------------------------
int ofxAAOneVectorOutputAlgorithm::getBinsNum(){
    return floatValues.size();
}
//-------------------------------------------
vector<float>& ofxAAOneVectorOutputAlgorithm::getValues(){
    return floatValues;
}
//-------------------------------------------
vector<float>& ofxAAOneVectorOutputAlgorithm::getSmoothedValues(float smthAmnt){
    
    for(int i=0; i<smoothedFloatValues.size(); i++){
        smoothedFloatValues[i] = smoothedFloatValues[i] * smthAmnt + (1-smthAmnt) * floatValues[i];
    }
    
    return smoothedFloatValues;
}
//-------------------------------------------
#pragma mark - ofxAAPitchSalienceFuntionPeaksAlgorithm
//-------------------------------------------
void ofxAAPitchSalienceFunctionPeaksAlgorithm::init(){
    
    limitPeaksNum = true;
    maxPeaksNum = 4;
    
    isActivated = true;

}
//-------------------------------------------
void ofxAAPitchSalienceFunctionPeaksAlgorithm::castValuesToFloat(){
    
    peaks.clear();
    
    peaks.resize(realSalienceBins.size());

    for (int i=0; i<realSalienceBins.size(); i++){
        peaks[i].bin = (float) realSalienceBins[i];
        peaks[i].value = realSalienceValues[i];
    }

}
//-------------------------------------------
vector<SalienceFunctionPeak>& ofxAAPitchSalienceFunctionPeaksAlgorithm::getPeaks(){
    
    if (limitPeaksNum && peaks.size() > maxPeaksNum){
        peaks.resize(maxPeaksNum);
    }
    
    return peaks;
}

//-------------------------------------------
vector<SalienceFunctionPeak>& ofxAAPitchSalienceFunctionPeaksAlgorithm::getSmoothedPeaks(float smthAmnt){
    
    if (limitPeaksNum && peaks.size() > maxPeaksNum){
        peaks.resize(maxPeaksNum);
    }
    
    smoothedPeaks.resize(peaks.size(), SalienceFunctionPeak());
    
    for(int i=0; i<smoothedPeaks.size(); i++){
        smoothedPeaks[i].bin = smoothedPeaks[i].bin * smthAmnt + (1-smthAmnt) * peaks[i].bin;
//        smoothedPeaks[i].bin = peaks[i].bin;
        smoothedPeaks[i].value = smoothedPeaks[i].value * smthAmnt + (1-smthAmnt) * peaks[i].value;
    }
    
    return smoothedPeaks;
}
//-------------------------------------------
#pragma mark - ofxAAPitchDetectAlgorithm
//-------------------------------------------
void ofxAAPitchDetectAlgorithm::init(){
    
    isActivated = true;
    
    pitchFloatVal = 0.0;
    confidenceFloatVal = 0.0;
    
    smoothedPitchFloatValue = 0.0;
    smoothedNormPitchFloatValue = 0.0;
    smoothedConfidenceFloatValue = 0.0;
    
}
//-------------------------------------------
void ofxAAPitchDetectAlgorithm::castValuesToFloat(){
    if(getIsActive()){
        pitchFloatVal = (float) pitchRealVal;
        confidenceFloatVal = (float)confidenceRealVal;
        //avoid negative values...
        if(confidenceFloatVal<0.0) confidenceFloatVal = 0.0;
    }
    else{
        pitchFloatVal = confidenceFloatVal = 0.0;
    }
    
}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getPitchValue(){
    return pitchFloatVal;
}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getPitchValueNormalized(){
    return ofMap(pitchFloatVal, 0.0, pitchMaxEstimatedValue, 0.0, 1.0, true);
}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getConfidenceValue(){
    return confidenceFloatVal;
}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getSmoothedPitchValueNormalized(float smthAmnt){
    float normVal = ofMap(pitchFloatVal, 0.0, pitchMaxEstimatedValue, 0.0, 1.0, true);
    smoothedNormPitchFloatValue =  smoothedNormPitchFloatValue * smthAmnt + (1-smthAmnt) * normVal;
    return smoothedNormPitchFloatValue;

}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getSmoothedPitchValue(float smthAmnt){
    smoothedPitchFloatValue =  smoothedPitchFloatValue * smthAmnt + (1-smthAmnt) * pitchFloatVal;
    return smoothedPitchFloatValue;
}
//-------------------------------------------
float ofxAAPitchDetectAlgorithm::getSmoothedConfidenceValue(float smthAmnt){
    smoothedConfidenceFloatValue =  smoothedConfidenceFloatValue * smthAmnt + (1-smthAmnt) * confidenceFloatVal;
    return smoothedConfidenceFloatValue;
}
//-------------------------------------------
void ofxAAPitchDetectAlgorithm::setMaxPitchEstimatedValue(float value){
    pitchMaxEstimatedValue = value;
}


//-------------------------------------------
#pragma mark - ofxAATuningFrequencyAlgorithm
//-------------------------------------------
void ofxAATuningFrequencyAlgorithm::castValuesToFloat(){
    if(getIsActive()){
        freqFloatVal = (float) freqRealVal;
        centsFloatVal = (float) centsRealVal;
    }
    else{
        freqFloatVal = centsFloatVal = 0.0;
    }
    
}
//-------------------------------------------
float ofxAATuningFrequencyAlgorithm::getFreqValue(){
    return freqFloatVal;
}
//-------------------------------------------
float ofxAATuningFrequencyAlgorithm::getCentsValue(){
    return centsFloatVal;
}
//-------------------------------------------

//-------------------------------------------------------
void ofxAudioAnalyzer::setup(int sampleRate, int bufferSize, int channels){
    
    _samplerate = sampleRate;
    _buffersize = bufferSize;
    _channels = channels;
    
    if(_channels <= 0){
        ofLogWarning()<<"ofxAudioAnalyzer: channels cant be set to none. Setting 1 channel";
        _channels = 1;
    }
    
    for(int i=0; i<_channels; i++){
        ofxAudioAnalyzerUnit * aaUnit = new ofxAudioAnalyzerUnit(_samplerate, _buffersize);
        channelAnalyzerUnits.push_back(aaUnit);
    }
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::reset(int sampleRate, int bufferSize, int channels){
    
    _samplerate = sampleRate;
    _buffersize = bufferSize;
    _channels = channels;
    
    if(_channels <= 0){
        ofLogWarning()<<"ofxAudioAnalyzer: channels cant be set to none. Setting 1 channel";
        _channels = 1;
    }
    
    for (int i=0; i<channelAnalyzerUnits.size(); i++){
        channelAnalyzerUnits[i]->exit();
    }
    channelAnalyzerUnits.clear();
    
    for(int i=0; i<_channels; i++){
        ofxAudioAnalyzerUnit * aaUnit = new ofxAudioAnalyzerUnit(_samplerate, _buffersize);
        channelAnalyzerUnits.push_back(aaUnit);
    }
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::analyze(const ofSoundBuffer & inBuffer){
   
    if(inBuffer.getNumChannels() != _channels){
        ofLogError()<<"ofxAudioAnalyzer: inBuffer channels number incorrect.";
        return;
    }
    
    if(channelAnalyzerUnits.size()!= _channels){
        ofLogError()<<"ofxAudioAnalyzer: wrong number of audioAnalyzerUnits";
        return;
    }
    
    if(inBuffer.getSampleRate() != _samplerate){
        ofLogWarning()<<"ofxAudioAnalyzer: inBuffer sample rate not matching.";
    }
    
    for (int i=0; i<_channels; i++){
        ofSoundBuffer chBuff;
        inBuffer.getChannel(chBuff, i);//copy channel from inBuffer
        if(channelAnalyzerUnits[i]!=nullptr){
            channelAnalyzerUnits[i]->analyze(chBuff.getBuffer());
            //cout<<"analyzer-"<<i<<"rms"<<channelAnalyzerUnits[i]->getRms()<<endl;
        }else{
            ofLogError()<<"ofxAudioAnalyzer: channelAnalyzer NULL pointer";
        }
        
    }
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::exit(){
    
    for(int i=0; i<channelAnalyzerUnits.size();i++){
        channelAnalyzerUnits[i]->exit();
    }
    
}
//-------------------------------------------------------
float ofxAudioAnalyzer::getValue(ofxAAAlgorithm algorithm, int channel, float smooth, bool normalized){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        return 0.0;
    }
    
    return channelAnalyzerUnits[channel]->getValue(algorithm, smooth, normalized);

}
//-------------------------------------------------------
vector<float>& ofxAudioAnalyzer::getValues(ofxAAAlgorithm algorithm, int channel, float smooth){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        vector<float>r (1, 0.0);
        return r;
    }
    
    return channelAnalyzerUnits[channel]->getValues(algorithm, smooth);
    
}
//-------------------------------------------------------
vector<SalienceFunctionPeak>& ofxAudioAnalyzer::getSalienceFunctionPeaks(int channel, float smooth){
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        //SalienceFunctionPeak peak = SalienceFunctionPeak();
        vector<SalienceFunctionPeak> r(1, SalienceFunctionPeak());
        return r;
    }
    
     return channelAnalyzerUnits[channel]->getPitchSaliencePeaksRef(smooth);

}
//-------------------------------------------------------
bool ofxAudioAnalyzer::getOnsetValue(int channel){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        return false;
    }
    
    return channelAnalyzerUnits[channel]->getOnsetValue();
    
}
//-------------------------------------------------------
bool ofxAudioAnalyzer::getIsActive(int channel, ofxAAAlgorithm algorithm){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting if its active is incorrect.";
        return false;
    }
    
    return channelAnalyzerUnits[channel]->getIsActive(algorithm);
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::resetOnsets(int channel){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        return;
    }
    
    channelAnalyzerUnits[channel]->resetOnsets();
}
//-------------------------------------------------------
void ofxAudioAnalyzer::setActive(int channel, ofxAAAlgorithm algorithm, bool state){

    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for setting active is incorrect.";
        return;
    }
    
    channelAnalyzerUnits[channel]->setActive(algorithm, state);
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::setMaxEstimatedValue(int channel, ofxAAAlgorithm algorithm, float value){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for setting max estimated value is incorrect.";
        return;
    }
    
    channelAnalyzerUnits[channel]->setMaxEstimatedValue(algorithm, value);
    
}
//-------------------------------------------------------
void ofxAudioAnalyzer::setOnsetsParameters(int channel, float alpha, float silenceTresh, float timeTresh, bool useTimeTresh){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        return;
    }
    channelAnalyzerUnits[channel]->setOnsetsParameters(alpha, silenceTresh, timeTresh, useTimeTresh);

}
//-------------------------------------------------------
void ofxAudioAnalyzer::setSalienceFunctionPeaksParameters(int channel, int maxPeaks){
    
    if (channel >= _channels){
        ofLogError()<<"ofxAudioAnalyzer: channel for getting value is incorrect.";
        return;
    }
    
    channelAnalyzerUnits[channel]->setSalienceFunctionPeaksParameters(maxPeaks);

}


//-------------------------------------------
void ofxAAOnsetsAlgorithm::setup(int bufferSize){
    /*
     at 44100:
     512 samples = 11.6 ms
     
     detetctBufferSize: 1 detection x buffers.
     512 x 64 = 742.4 ms
     */
    
//    detecBufferSize = bufferSize;
    detecBufferSize = 64;
    
    detection_sum.assign(detecBufferSize, 0.0);
    detections.assign(3, vector<Real> (detecBufferSize));
    silenceThreshold = 0.02;
    alpha = 0.1;
    timeThreshold = 100.0;
    bufferNumThreshold = 7; //116 ms at 60 fps
    lastOnsetTime = 0.0;
    lastOnsetBufferNum = 0;
    addHfc = addComplex = addFlux = true;
    hfc_max = complex_max = flux_max = 0.0;
    usingTimeThreshold = true;
    onsetsMode = TIME_BASED;
    bufferCounter = 0;
    
    _value = false;
    
    init();
    
    onsetHfc.init();
    onsetComplex.init();
    onsetFlux.init();
}
//-------------------------------------------
void ofxAAOnsetsAlgorithm::compute(){
    
    onsetHfc.compute();
    onsetComplex.compute();
    onsetFlux.compute();

}
//-------------------------------------------
void ofxAAOnsetsAlgorithm::castValuesToFloat(){
    
    if(isActivated){
        onsetHfc.castValueToFloat();
        onsetComplex.castValueToFloat();
        onsetFlux.castValueToFloat();
    }else{
        onsetHfc.setValueZero();
        onsetFlux.setValueZero();
        onsetComplex.setValueZero();
    }

}
//-------------------------------------------
void ofxAAOnsetsAlgorithm::evaluate(){
    
    //is current buffer an Onset?
    bool isCurrentBufferOnset = onsetBufferEvaluation(onsetHfc.getValue(), onsetComplex.getValue(), onsetFlux.getValue());
    
    //if current buffer is onset, check for timeThreshold evaluation
    if (usingTimeThreshold && isCurrentBufferOnset){
        
        switch (onsetsMode) {
            case TIME_BASED:
                _value = onsetTimeThresholdEvaluation();
                break;
            case BUFFER_NUM_BASED:
                _value = onsetBufferNumThresholdEvaluation();
                break;
            default:
                break;
        }
    
    }
    else{
        _value = isCurrentBufferOnset;
    }
    
    //update bufferCounter for frameBased timeThreshold evaluation:
    if (onsetsMode == BUFFER_NUM_BASED) bufferCounter++;
    
}

//--------------------------------------------------------------
bool ofxAAOnsetsAlgorithm::onsetBufferEvaluation (Real iDetectHfc, Real iDetectComplex, Real iDetectFlux){
    
    Real prop_hfc, prop_complex, prop_flux;
    
    if (iDetectHfc > hfc_max) {
        prop_hfc = iDetectHfc/hfc_max;
        hfc_max = iDetectHfc;
        for (int i=0; i<detections[0].size(); i++)
            detections[0][i] /= prop_hfc;
    }
    if (iDetectComplex > complex_max){
        prop_complex = iDetectComplex/complex_max;
        complex_max = iDetectComplex;
        for (int i=0; i<detections[1].size(); i++)
            detections[1][i] /= prop_complex;
    }
    if (iDetectFlux > flux_max){
        prop_flux = iDetectFlux/flux_max;
        flux_max = iDetectFlux;
        for (int i=0; i<detections[2].size(); i++)
            detections[2][i] /= prop_flux;
    }
    
    Real hfc_norm = iDetectHfc / hfc_max;
    Real complex_norm = iDetectComplex / complex_max;
    Real flux_norm = iDetectFlux / flux_max;
    
    detections[0].push_back(hfc_norm);
    detections[0].erase(detections[0].begin());
    
    detections[1].push_back(complex_norm);
    detections[1].erase(detections[1].begin());
    
    detections[2].push_back(flux_norm);
    detections[2].erase(detections[2].begin());
    
    for (int i=0; i<detection_sum.size(); i++){
        int n=0;
        detection_sum[i]=0;
        if(addHfc){
            detection_sum[i] += detections[0][i];
            n++;
        }
        if(addComplex){
            detection_sum[i] += detections[1][i];
            n++;
        }
        if(addFlux){
            detection_sum[i] += detections[2][i];
            n++;
        }
        if(n>0) detection_sum[i] /= n;
        if(detection_sum[i] < silenceThreshold) detection_sum[i] = 0.0;
    }
    
    Real buffer_median = median (detection_sum);
    Real buffer_mean = mean (detection_sum);
    Real onset_thhreshold = buffer_median + alpha * buffer_mean;
    
    bool onsetDetection = detection_sum[detection_sum.size()-1] > onset_thhreshold;
    
    return onsetDetection;
    
}
//----------------------------------------------
bool ofxAAOnsetsAlgorithm::onsetTimeThresholdEvaluation(){
    
    bool onsetTimeEvaluation = false;
    
    float currentTime = ofGetElapsedTimeMillis();
    
    //elapsed time since last onset:
    float elapsed = currentTime - lastOnsetTime;
    
    if (elapsed>timeThreshold){
        onsetTimeEvaluation = true;
        lastOnsetTime = currentTime;
    }
    
    return onsetTimeEvaluation;
}
//----------------------------------------------
bool ofxAAOnsetsAlgorithm::onsetBufferNumThresholdEvaluation(){
    
    bool onsetBufferNumEvaluation = false;
    
    //elapsed frames/buffers since last onset:
    int elapsedBuffers = bufferCounter - lastOnsetBufferNum;
    
    if (elapsedBuffers > bufferNumThreshold){
        onsetBufferNumEvaluation = true;
        lastOnsetBufferNum = bufferCounter;
    }
    
    return onsetBufferNumEvaluation;
    
    
}
//----------------------------------------------
void ofxAAOnsetsAlgorithm::reset(){
    
    hfc_max = complex_max = flux_max = 0.0;
    for (int i=0; i<detection_sum.size(); i++){
        detection_sum[i] = 0.0;
    }
    
    //necessary?
    onsetHfc.algorithm->reset();
    onsetComplex.algorithm->reset();
    onsetFlux.algorithm->reset();
    
    
    bufferCounter = 0;
}



