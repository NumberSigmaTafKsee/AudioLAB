//saul

class IIRFilterFirstOrder;

class IIRCoefficientsSinglePole
{
public:
    //==============================================================================
    IIRCoefficientsSinglePole() noexcept;
    ~IIRCoefficientsSinglePole() noexcept;


    IIRCoefficientsSinglePole (double c0, double c1, double c2, double c3) noexcept;
    
    //==============================================================================
    static IIRCoefficientsSinglePole makeLowPass (double sampleRate,
                                        double frequency) noexcept;
    
    static IIRCoefficientsSinglePole makeLowPass (double sampleRate,
                                        double frequency,
                                        double Q) noexcept;
    
    //==============================================================================
    static IIRCoefficientsSinglePole makeHighPass (double sampleRate,
                                         double frequency) noexcept;
    
    static IIRCoefficientsSinglePole makeHighPass (double sampleRate,
                                         double frequency,
                                         double Q) noexcept;
    
    //==============================================================================
    float coefficients[3];
};


class IIRFilterFirstOrder
{
public:
    
    IIRFilterFirstOrder() noexcept;
    ~IIRFilterFirstOrder() noexcept;
    
    //==============================================================================
    void setCoefficients (const IIRCoefficientsSinglePole& newCoefficients) noexcept;
    
    IIRCoefficientsSinglePole getCoefficients() const noexcept    { return coefficients; }
    //==============================================================================
    void reset() noexcept;

    //NOT thread-safe - no locking, etc..
    float processSingleSampleRaw (float sample) noexcept;
    
protected:
    //==============================================================================
    IIRCoefficientsSinglePole coefficients;
    float v1 = 0.f;
};

