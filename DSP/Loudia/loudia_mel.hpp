#pragma once

namespace Loudia
{
        /**
    * @class MelBands
    *
    * @brief Algorithm to get the band values of Mel-scale frequency warpped magnitude spectrums.
    *
    * This class represents an object to Mel bands on vectors
    * of Real values.  This method is a special case of the Bands algorithm and is used
    * in many other spectral algorithms such as MFCC and SpectralWhitening.
    *
    * The method consists in a set triangular 50% overlapping windows spaced evenly on a
    * Mel-frequency scale.
    *
    * The sampleRate and FFT size of the input spectrum are specified using setSampleRate() and
    * setFftSize().
    *
    * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
    * setHighFrequency().
    *
    * The number of bands is specified using setBandCount().
    *
    * @author Ricard Marxer
    *
    * @sa Bands, MFCC, SpectralWhitening
    */
    class MelBands {
    public:
    /**
        @enum ScaleType
        @brief Specifies the type of the scale used.
        @sa scaleType
    */
    enum ScaleType {
        STEVENS = 0 /**< Mel scales computed using the original formula proposed by:
                    *
                    * Stevens, Stanley Smith; Volkman; John; & Newman, Edwin. (1937).
                    * A scale for the measurement of the psychological magnitude of pitch.
                    * Journal of the Acoustical Society of America, 8 (3), 185-190.
                    *
                    */,
        FANT = 1 /**< Mel scales computed using the formula proposed by:
                *
                * Fant, Gunnar. (1968).
                * Analysis and synthesis of speech processes.
                * In B. Malmberg (Ed.), Manual of phonetics (pp. 173-177). Amsterdam: North-Holland.
                *
                */,
        GREENWOOD = 2 /**< Mel scales computed using the Greenwood function:
                    *
                    * Greenwood, DD. (1990)
                    * A cochlear frequency-position function for several species - 29 years later,
                    * Journal of the Acoustical Society of America, vol. 87, pp. 2592-2605.
                    *
                    */
    };

    protected:
        Real _lowFrequency;
        Real _highFrequency;
        int _bandCount;
        Real _sampleRate;
        int _fftSize;
        ScaleType _scaleType;

        Bands _bands;
        MatrixXR _centersLinear;

        Real (*_linearToMel)(Real linearFreq);

        Real (*_melToLinear)(Real melFreq);

        void (*_linearToMelMatrix)(const MatrixXR& linearFreq, MatrixXR* melFreq);

        void (*_melToLinearMatrix)(const MatrixXR& melFreq, MatrixXR* linearFreq);

        void triangleWindow(MatrixXR* window, Real start, Real stop, Real center = -1, Real height = Real(1.0));

    public:
        /**
            Constructs a Mel bands object with the specified @a lowFrequency, @a highFrequency,
            @a bandCount, @a sampleRate, @a fftSize and @a scaleType settings.
            @param lowFrequency frequency of the lowest Mel band,
            must be greater than zero 0 and lower than half the sampleRate.
            @param highFrequency frequency of the highest Mel band,
            must be greater than zero 0 and lower than half the sampleRate.
            @param bandCount number of Mel bands.
            @param sampleRate sampleRate frequency of the input signal.
            @param fftSize size of the FFT.
            @param scaleType scale used for the frequency warping.
        */
        MelBands(Real lowFrequency = 50.0, Real highFrequency = 6000.0, int bandCount = 40, Real sampleRate = 44100.0, int fftSize = 1024, ScaleType scaleType = GREENWOOD);

        void setup();
        void reset();

        /**
            Calculates the bands of @a spectrums.
            @param spectrums matrix of Real values.
            @param bands pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a spectrums and as many columns as the number of bands
            as specified by bandCount.
            Note that if the output matrix is not of the required size it will be resized,
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& spectrums, MatrixXR* bands);

        /**
            Return the vector of weights.
            @sa starts, bandWeights, setStartsWeights
        */
        std::vector<MatrixXR> weights() const;

        /**
            Return in @a bandWeights the weights of the band given by the index @a band.
            @sa weights
        */
        void bandWeights(int band, MatrixXR* bandWeights) const;

        /**
            Return in @a result the single column matrix of start indices of the bands.
        */
        void starts(MatrixXI* result) const;

        /**
            Return number of bands.
        */
        int bands() const;

        /**
            Return in @a result the single column matrix of center fractional indices of the bands.
        */
        void centers(MatrixXR* result) const;

        /**
            Returns the number of bands to be performed.
            The default is 40.
            @sa setBandCount()
        */
        int bandCount() const;

        /**
            Specifies the @a count of bands to be performed.
            @sa bandCount()
        */
        void setBandCount( int count, bool callSetup = true );

        /**
            Return the low frequency.
            The default is 50.0.
            @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
        */
        Real lowFrequency() const;

        /**
            Specifies the low @a frequency.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            @sa lowFrequency, highFrequency, setHighFrequency
        */
        void setLowFrequency( Real frequency, bool callSetup = true );

        /**
            Return the high frequency.
            The default is 6000.0.
            @sa lowFrequency, setLowFrequency, setHighFrequency
        */
        Real highFrequency() const;

        /**
            Specifies the high @a frequency.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            @sa lowFrequency, highFrequency, setLowFrequency
        */
        void setHighFrequency( Real frequency, bool callSetup = true );

        /**
            Return the sampleRate frequency of the input signal.
            The default is 44100.0.
            @sa setSampleRate
        */
        Real sampleRate() const;

        /**
            Specifies the sampleRate @a frequency of the input signal.
            @sa sampleRate
        */
        void setSampleRate( Real frequency, bool callSetup = true );

        /**
            Returns the size of the FFT that has been performed for the input.
            The default is 1024.
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT that has been performed for the input.
            The given @a size must be higher than 0.
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Return the type of the frequency warping scale.
            By default it is GREENWOOD.
            @sa setScaleType()
        */
        MelBands::ScaleType scaleType() const;

        /**
            Specify the type of the frequency warping scale.
            @sa scaleType()
        */
        void setScaleType( MelBands::ScaleType type, bool callSetup = true );
    }
        /**
    *
    * Mel scales computed using the Greenwood function:
    *
    * Greenwood, DD. (1990)
    * A cochlear frequency-position function for several species - 29 years later,
    * Journal of the Acoustical Society of America, vol. 87, pp. 2592-2605.
    *
    */
    Real linearToMelGreenwood1990(Real linearFreq);

    Real melToLinearGreenwood1990(Real melFreq);

    void linearToMelMatrixGreenwood1990(const MatrixXR& linearFreq, MatrixXR* melFreq);

    void melToLinearMatrixGreenwood1990(const MatrixXR& melFreq, MatrixXR* linearFreq);


    /**
    *
    * Mel scales computed using the original formula proposed by:
    *
    * Stevens, Stanley Smith; Volkman; John; & Newman, Edwin. (1937). 
    * A scale for the measurement of the psychological magnitude of pitch.
    * Journal of the Acoustical Society of America, 8 (3), 185-190.
    *
    */
    Real linearToMelStevens1937(Real linearFreq);

    Real melToLinearStevens1937(Real melFreq);

    void linearToMelMatrixStevens1937(const MatrixXR& linearFreq, MatrixXR* melFreq);

    void melToLinearMatrixStevens1937(const MatrixXR& melFreq, MatrixXR* linearFreq);


    /**
    *
    * Mel scales computed using the formula proposed by:
    *  
    * Fant, Gunnar. (1968).
    * Analysis and synthesis of speech processes.
    * In B. Malmberg (Ed.), Manual of phonetics (pp. 173-177). Amsterdam: North-Holland.
    *
    */
    Real linearToMelFant1968(Real linearFreq);

    Real melToLinearFant1968(Real melFreq);

    void linearToMelMatrixFant1968(const MatrixXR& linearFreq, MatrixXR* melFreq);

    void melToLinearMatrixFant1968(const MatrixXR& melFreq, MatrixXR* linearFreq);


    MelBands::MelBands(Real lowFrequency, Real highFrequency, int bandCount, Real sampleRate, int fftSize, ScaleType scaleType) 
    {
    
    LOUDIA_DEBUG("MELBANDS: Constructor lowFrequency: " << _lowFrequency << 
            ", highFrequency: " << _highFrequency << 
            ", bandCount: " << _bandCount << 
            ", sampleRate: " << _sampleRate << 
            ", fftSize: " << _fftSize << 
            ", scaleType:" << _scaleType);

    if ( lowFrequency >= highFrequency ) {
        // Throw an exception, highFrequency must be higher than lowFrequency
    }

    if ( bandCount <= 0 ) {
        // Throw an exception, bandCount must be higher than 0
    }
    
    setLowFrequency( lowFrequency, false );
    setHighFrequency( highFrequency, false );
    setBandCount( bandCount, false );
    setSampleRate( sampleRate, false );
    setFftSize( fftSize, false );
    setScaleType( scaleType, false );
    
    setup();
    
    LOUDIA_DEBUG("MELBANDS: Constructed");
    
    }

    void MelBands::setup(){
    LOUDIA_DEBUG("MELBANDS: Setting up...");
    
    // Set the linearToMel and melToLinear functions
    switch(_scaleType) {
    case STEVENS:
        _linearToMel = &linearToMelStevens1937;
        _melToLinear = &melToLinearStevens1937;
    
        _linearToMelMatrix = &linearToMelMatrixStevens1937;
        _melToLinearMatrix = &melToLinearMatrixStevens1937;

        break;

    case FANT:
        _linearToMel = &linearToMelFant1968;
        _melToLinear = &melToLinearFant1968;
    
        _linearToMelMatrix = &linearToMelMatrixFant1968;
        _melToLinearMatrix = &melToLinearMatrixFant1968;

        break;

    case GREENWOOD:
        _linearToMel = &linearToMelGreenwood1990;
        _melToLinear = &melToLinearGreenwood1990;
    
        _linearToMelMatrix = &linearToMelMatrixGreenwood1990;
        _melToLinearMatrix = &melToLinearMatrixGreenwood1990;

        break;

    }
    
    Real highMel = _linearToMel( _highFrequency );
    Real lowMel = _linearToMel( _lowFrequency );
    
    LOUDIA_DEBUG("MELBANDS: lowMel: " << lowMel << ", highMel: " << highMel);

    Real stepMel = (highMel - lowMel) / (_bandCount + 1.0);
    Real stepSpectrum = Real(_fftSize) / _sampleRate;
    
    // start Mel frequencies of filters
    MatrixXR starts(_bandCount, 1);
    for (int i=0; i<starts.rows(); i++) {
        starts(i, 0) = (Real(i) * stepMel + lowMel);
    }

    MatrixXR startsLinear;
    _melToLinearMatrix(starts, &startsLinear);
    startsLinear *= stepSpectrum;

    // stop Mel frequencies of filters
    MatrixXR stops(_bandCount, 1);
    for (int i=0; i<stops.rows(); i++) {
        stops(i, 0) = (Real(i + 2) * stepMel + lowMel);
    }

    MatrixXR stopsLinear;
    _melToLinearMatrix(stops, &stopsLinear);
    stopsLinear *= stepSpectrum;


    // center Mel frequencies of filters
    MatrixXR centers(_bandCount, 1);
    for (int i=0; i<centers.rows(); i++) {
        centers(i, 0) = (Real(i + 1) * stepMel + lowMel);
    }

    _centersLinear = startsLinear + (stopsLinear - startsLinear) / 2.0;
    //melToLinearMatrixGreenwood1990(centers, &centersLinear);
    //centersLinear *= stepSpectrum;
    
    // start bins of filters
    MatrixXI startBins = startsLinear.array().ceil().cast<Integer>();

    // stop bins of filters
    MatrixXI stopBins = stopsLinear.array().ceil().cast<Integer>();

    std::vector<MatrixXR> weights;
    
    // fill in the weights
    for (int i=0; i < startBins.rows(); i++) {
        int startBin = startBins(i, 0);
        int stopBin = stopBins(i, 0);
        
        int filterLength = stopBin - startBin;
        
        LOUDIA_DEBUG("MELBANDS: filterLength: " << filterLength);
        MatrixXR newFilter = MatrixXR::Constant(1, 1, 1.0);
        if (filterLength != 0){
        newFilter.resize(filterLength, 1);
        
        Real start = startsLinear(i, 0);
        Real center = _centersLinear(i, 0);
        Real stop = stopsLinear(i, 0);
        
        triangleWindow(&newFilter, start - startBin, stop  - startBin, center  - startBin);
        }
        
        weights.push_back(newFilter);
    }

    _bands.setStartsWeights(startBins, weights);

    LOUDIA_DEBUG("MELBANDS: Finished set up...");
    }

    void MelBands::process(const MatrixXR& spectrum, MatrixXR* bands) {  
    _bands.process(spectrum, bands);
    }

    void MelBands::triangleWindow(MatrixXR* window, Real start, Real stop, Real center, Real height) {
    int filterLength = (*window).rows();

    if (center == -1) {
        LOUDIA_DEBUG("MELBANDS: Triangle window setting the center by default");
        center = start + (stop - start) / 2.0;
    }

    if ((center < 0) || (center > filterLength) || (center == start) || (center == stop)) {
        // Throw an exception invalid filter center
    }

    LOUDIA_DEBUG("MELBANDS: Creating triangle window");
    LOUDIA_DEBUG("MELBANDS: Triangle start: " << start << ", stop: " << stop << ", center: " << center << ", height: " << height);
    
    for (int i = 0; i < filterLength; i++) {
        if (i <= center) {
        (*window)(i,0) =  height * (Real(i) - start) / (center - start);
        } else {
        (*window)(i,0) =  height * (Real(1.0) - ((Real(i) - center) / (stop - center)));
        }
    }
    
    LOUDIA_DEBUG("MELBANDS: Triangle window created: [" << (*window).transpose() << "]");
    }

    void MelBands::reset(){
    // Initial values
    _bands.reset();
    }

    void MelBands::starts(MatrixXI* result) const {
    return _bands.starts( result );
    }

    void MelBands::centers(MatrixXR* result) const {
    (*result) = _centersLinear;
    }

    std::vector<MatrixXR> MelBands::weights() const {
    return _bands.weights();
    }

    void MelBands::bandWeights(int band, MatrixXR* bandWeights) const {
    return _bands.bandWeights( band, bandWeights );
    }

    int MelBands::bands() const {
    return _bands.bands();
    }

    Real MelBands::lowFrequency() const{
    return _lowFrequency;
    }
    
    void MelBands::setLowFrequency( Real frequency, bool callSetup ){
    _lowFrequency = frequency;
    if ( callSetup ) setup();
    }

    Real MelBands::highFrequency() const{
    return _highFrequency;
    }
    
    void MelBands::setHighFrequency( Real frequency, bool callSetup ){
    _highFrequency = frequency;
    if ( callSetup ) setup();
    }

    Real MelBands::sampleRate() const{
    return _sampleRate;
    }
    
    void MelBands::setSampleRate( Real frequency, bool callSetup ){
    _sampleRate = frequency;
    if ( callSetup ) setup();
    }

    int MelBands::bandCount() const {
    return _bandCount;
    }

    void MelBands::setBandCount( int count, bool callSetup ) {
    _bandCount = count;
    if ( callSetup ) setup();
    }

    int MelBands::fftSize() const{
    return _fftSize;
    }

    void MelBands::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }

    MelBands::ScaleType MelBands::scaleType() const{
    return _scaleType;
    }

    void MelBands::setScaleType( MelBands::ScaleType type, bool callSetup ) {
    _scaleType = type;
    if ( callSetup ) setup();
    }

    Real linearToMelGreenwood1990(Real linearFreq) {
    return log10((linearFreq / 165.4) + 1.0) / 2.1;
    }

    Real melToLinearGreenwood1990(Real melFreq) {
    return 165.4 * (pow(10.0, 2.1 * melFreq) - 1.0);
    }

    void linearToMelMatrixGreenwood1990(const MatrixXR& linearFreq, MatrixXR* melFreq) {
    LOUDIA_DEBUG("MELBANDS: Scaling (Greenwood 1990) linearFreq: " << linearFreq);

    (*melFreq) = ((linearFreq / 165.4).array() + 1.0).logN(10) / 2.1;
    }

    void melToLinearMatrixGreenwood1990(const MatrixXR& melFreq, MatrixXR* linearFreq) {
    LOUDIA_DEBUG("MELBANDS: Scaling (Greenwood 1990) melFreq: " << melFreq);

    (*linearFreq) = 165.4 * ((melFreq * 2.1).array().expN(10.0) - 1.0);
    }


    Real linearToMelStevens1937(Real linearFreq) {
    return log((linearFreq / 700.0) + 1.0) * 1127.01048;
    }

    Real melToLinearStevens1937(Real melFreq) {
    return (exp(melFreq / 1127.01048) - 1.0) * 700.0;
    }

    void linearToMelMatrixStevens1937(const MatrixXR& linearFreq, MatrixXR* melFreq) {
    (*melFreq) = ((linearFreq / 700.0).array() + 1.0).log() * 1127.01048;
    }

    void melToLinearMatrixStevens1937(const MatrixXR& melFreq, MatrixXR* linearFreq) {
    (*linearFreq) = ((melFreq / 1127.01048).array().exp() - 1.0) * 700.0;
    }


    Real linearToMelFant1968(Real linearFreq) {
    return (1000.0 / log(2.0)) * log(1.0 + linearFreq / 1000.0);
    }

    Real melToLinearFant1968(Real melFreq) {
    return 1000.0 * (exp(melFreq * log(2.0) / 1000.0) - 1.0);
    }

    void linearToMelMatrixFant1968(const MatrixXR& linearFreq, MatrixXR* melFreq) {
    (*melFreq) = (1000.0 / log(2.0)) * ((linearFreq / 1000.0).array() + 1.0).log();
    }

    void melToLinearMatrixFant1968(const MatrixXR& melFreq, MatrixXR* linearFreq) {
    (*linearFreq) = 1000.0 * ((melFreq * log(2.0) / 1000.0).array().exp() - 1.0);
    }

}