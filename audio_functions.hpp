namespace Oscillators::Functions
{
    struct FunctionGenerator : public GeneratorProcessor
    {
        enum Type
        {
            NOISE,
            SINEWAVE,   
            COSWAVE,
            PHASOR,
            SQUARE,
            PULSE,
            SAW,
            REVERSE_SAW,
            TRIANGLE,
        } 
        type = SAW;

        enum Polarity
        {
            BIPOLAR,
            POSITIVE,
            NEGATIVE,
        }
        polarity = POSITIVE;

        FunctionGenerator(Type t) 
        {
            inc = 440.0/sampleRate;
            type = t;
        }
        FunctionGenerator(DspFloatType f, DspFloatType sr=44100) : GeneratorProcessor(), frequency(f),sampleRate(sr),startphase(0),phase(0.0) {
            inc = f/sr;
        }
        void setFrequency(DspFloatType f)
        {
            frequency = f;
            inc = f/sampleRate;
        }
        void phaseReset(DspFloatType phaseIn) {
            // This allows you to set the phase of the oscillator to anything you like.
            phase = phaseIn;
        }

        inline void phaseIncrement() {        
            phase = fmod(phase+inc,1);
        } 

        DspFloatType phasor(DspFloatType frequency, DspFloatType startphase=0, DspFloatType endphase=1) {
            
            output = phase;
            
            if (phase < startphase) {
                phase = startphase;
            }
            
            if (phase >= endphase)
                phase = startphase;
            
            phase += ((endphase - startphase) / (sampleRate / frequency));
            
            return (output);
        }
            
        DspFloatType noise() {
            // White Noise
            // always the same unless you seed it.
            DspFloatType r = std::rand() / (DspFloatType)RAND_MAX;
            return r;
        }

        DspFloatType sinewave() {
            return 0.5 + 0.5*sin(phase * TWOPI);            
        }

        DspFloatType coswave() {
            return 0.5 + 0.5*cos(phase * TWOPI);            
        }

        DspFloatType phasor() {
            return phase;
        }

        DspFloatType square() {
            return phase < 0.5f ? 0 : 1.0f;            
        }

        DspFloatType pulse() {
            DspFloatType output;
            if (phase < duty)
                output = 0;
            if (phase > duty)
                output = 1.f;
            return output;
        }
        DspFloatType saw() 
        { 
            return phase;
        }

        DspFloatType triangle() {     
            DspFloatType output;   
            if (phase <= 0.5f) {
                output = (phase - 0.25f) * 4;
            } else {
                output = ((1.0f - phase) - 0.25f) * 4;
            }
            return (output);
        }

        DspFloatType function() {
            switch(type)
            {
            case SINEWAVE: return sinewave();
            case COSWAVE: return coswave();
            case PHASOR: return phasor();
            case SQUARE: return square();
            case PULSE: return pulse();
            case SAW: return saw();            
            case REVERSE_SAW: return 1-saw();
            }
            return triangle();
        }
    
        DspFloatType Tick(DspFloatType I=0, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0) {        
            //DspFloatType p = phase;
            //phase = phasor(frequency + 0.5*X*frequency) + Y + phase;
            //phase = std::fmod(phase,1);
            DspFloatType r = function();
            //phase = p;        
            phaseIncrement();      
            if(polarity == POSITIVE) return r;
            else if(polarity == NEGATIVE) return -r;
            return 2*r-1;
        }
        DspFloatType operator()() {
            return Tick();
        }
        DspFloatType sampleRate;
        DspFloatType frequency;
        DspFloatType phase;
        DspFloatType duty = 0.5f;
        DspFloatType startphase;
        DspFloatType endphase;
        DspFloatType output;
        DspFloatType tri;
        DspFloatType inc;
    };

    struct NoiseFunction : public FunctionGenerator
    {
        NoiseFunction() : FunctionGenerator(FunctionGenerator::Type::NOISE)
        {

        }
    };
    struct SinFunction : public FunctionGenerator
    {
        SinFunction() : FunctionGenerator(FunctionGenerator::Type::SINEWAVE)
        {

        }
    };
    struct CosFunction : public FunctionGenerator
    {
        CosFunction() : FunctionGenerator(FunctionGenerator::Type::COSWAVE)
        {

        }
    };
    struct PhasorFunction : public FunctionGenerator
    {
        PhasorFunction() : FunctionGenerator(FunctionGenerator::Type::PHASOR)
        {

        }
    };
    struct SquareFunction : public FunctionGenerator
    {
        SquareFunction() : FunctionGenerator(FunctionGenerator::Type::SQUARE)
        {

        }
    };
    struct PulseFunction : public FunctionGenerator
    {
        PulseFunction() : FunctionGenerator(FunctionGenerator::Type::PULSE)
        {

        }
    };
    struct SawFunction : public FunctionGenerator
    {
        SawFunction() : FunctionGenerator(FunctionGenerator::Type::SAW)
        {

        }
    };    
    struct ReverseSawFunction : public FunctionGenerator
    {
        ReverseSawFunction() : FunctionGenerator(FunctionGenerator::Type::REVERSE_SAW)
        {

        }
    };    
    struct TriangleFunction : public FunctionGenerator
    {
        TriangleFunction() : FunctionGenerator(FunctionGenerator::Type::TRIANGLE)
        {

        }
    };
}
