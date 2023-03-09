#pragma once

namespace Loudia
{
        /**
    * WARNING: the process() method is NOT re-entrant.
    */
    class VoiceActivityDetection {

    protected:
        int _lowBand, _highBand;
        Real _sampleRate;
        int _fftSize;
        MatrixXR _memory;
        int _currentMemoryPos;
        int _memorySize;

        BarkBands _barkBands;
        MatrixXR _bands;

        int _halfSize;
            

    public:
        VoiceActivityDetection(int lowBand = 4, int highBand = 16,
                                Real sampleRate = 44100,
                                int fftSize = 1024,
                                int memorySize = 12);
        ~VoiceActivityDetection();
        
        void process(const MatrixXR& frames, MatrixXR* vad);

        void setup();
        void reset();


        // Parameters
        int lowBand() const { return _lowBand; }
        void setLowBand(int lowBand, bool callSetup = true ) { _lowBand = lowBand; if ( callSetup ) setup(); }

        int highBand() const { return _highBand; }
        void setHighBand(int highBand, bool callSetup = true ) { _highBand = highBand; if ( callSetup ) setup(); }

        Real sampleRate() const { return _sampleRate; }
        void setSampleRate(Real sampleRate, bool callSetup = true ) { _sampleRate = sampleRate; if ( callSetup ) setup(); }

        int memorySize() const { return _memorySize; }
        void setMemorySize(int memorySize, bool callSetup = true ) { _memorySize = memorySize; if ( callSetup ) setup(); }

        /**
            Returns the size of the FFT to be performed.  The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const { return _fftSize; };

        /**
            Specifies the @a size of the FFT to be performed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true ) { _fftSize = size; if ( callSetup ) setup(); }

    };    

VoiceActivityDetection::VoiceActivityDetection(int lowBand, int highBand,
                                               Real sampleRate,
                                               int fftSize,
                                               int memorySize) : 
  _lowBand(lowBand),
  _highBand(highBand),
  _sampleRate(sampleRate),
  _fftSize(fftSize),
  _memorySize(memorySize),
  _barkBands(lowBand, highBand, sampleRate, fftSize)
{

  LOUDIA_DEBUG("VoiceActivityDetection: Constructor");

  setup();
  
  LOUDIA_DEBUG("VoiceActivityDetection: Constructed");
}

VoiceActivityDetection::~VoiceActivityDetection(){
  LOUDIA_DEBUG("VoiceActivityDetection: Destroying...");
  LOUDIA_DEBUG("VoiceActivityDetection: Destroyed out");
}

void VoiceActivityDetection::setup(){
  LOUDIA_DEBUG("VoiceActivityDetection: Setting up...");

  _memory = MatrixXR::Zero(_memorySize, _highBand - _lowBand + 1);
  _currentMemoryPos = 0;
  
  _barkBands.setSampleRate(_sampleRate, false);
  _barkBands.setLowBand(_lowBand, false);
  _barkBands.setHighBand(_highBand, false);
  _barkBands.setup();
  
  LOUDIA_DEBUG("VoiceActivityDetection: Finished set up...");
}

void VoiceActivityDetection::process(const MatrixXR& frames, MatrixXR* vad){
  const int rows = frames.rows();
  
  vad->resize(rows, 1);

  for (int i=0; i < rows; i++){
    // compute barkbands
    _barkBands.process(frames.row(0), &_bands);

    // copy frame into memory
    _memory.row(_currentMemoryPos) = _bands.row(0);

    _currentMemoryPos = (_currentMemoryPos + 1) % _memorySize;

    // compute the VAD
    RowXR LTSE = _memory.colwise().maxCoeff();
    RowXR noise = _memory.colwise().sum() / _memorySize;

    (*vad)(i,0) = log10((LTSE.array().square() / noise.array().square()).sum());
  }
}

void VoiceActivityDetection::reset(){
  _barkBands.reset();
}

}