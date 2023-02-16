#pragma once

namespace Loudia
{
        /**
    * @class Bands
    *
    * @brief Algorithm to calculate the sum of values in a given set of weighted bands.
    *
    * This class represents an object to perform calculatinos of band values.
    * The value of a band corresponds to the sum of the values of the input array in the band's
    * positions multiplied with the band's weights.
    *
    * In this implementation the positions of a given band are defined by the index of the first
    * array cell of the band and the size of the weights array.
    * The full configuration of the bands algorithm is defined using a single column matrix for the
    * starts of the bands and a vector of single row matrices for the weights each band.
    *
    * Note that the number of rows of the starts matrix and the size of the vector of weights must
    * be the same, and this will be the number of bands.
    *
    * @author Ricard Marxer
    *
    * @sa MelBands
    */
    class Bands {
    protected:
        // Internal parameters
        MatrixXI _starts;
        std::vector<MatrixXR> _weights;

        // Internal variables

    public:
        /**
            Constructs a Bands object with a single band covering the entire array.
        */
        Bands();

        /**
            Constructs a Bands object with the specified @a starts and @a
            weights setting.
            @param starts single column matrix of Integers that determine the
            first array cell of each band.
            @param weights vector of single row matrices of Reals that determine the
            values of the weights of each band.
        */
        Bands(MatrixXI starts, std::vector<MatrixXR> weights);

        /**
            Destroys the Bands algorithm and frees its resources.
        */
        ~Bands();

        /**
            Calculates the bands of @a frames using the specified starts and weights properties.
            @param frames matrix of Real values.
            @param bands pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a frames and as many columns as the number of bands
            (rows in the starts matrix and elements in the weights vector).
            Note that if the output matrix is not of the required size it will be resized,
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR&  frames, MatrixXR* bands);

        void setup();
        void reset();

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
            Determines the @a starts positions and @a weights of the bands.
            Note that the number of rows of @a starts and the size of @a weights should be the same and
            will determine the number of bands.
        */
        void setStartsWeights(const MatrixXI& starts, std::vector<MatrixXR> weights, bool callSetup = true);
    };

    /**
    * @class BarkBands
    *
    * @brief Algorithm to get the band values of Bark-scale frequency warpped magnitude spectrums.
    *
    * This class represents an object to calculate the Bark bands on vectors
    * of Real values.  This method is a special case of the Bands algorithm.
    *
    * The method consists in a set rectangular windows spaced evenly on a
    * Bark-frequency scale.
    *
    * The sampleRate and FFT size of the input spectrum are specified using setSampleRate() and
    * setFftSize().
    *
    * The first and last bands are specified using setLowBand() and
    * setHighBand().
    *
    *
    * @author Ricard Marxer
    *
    * @sa Bands, MelBands
    */
    class BarkBands {
    protected:
        Real _lowBand;
        Real _highBand;
        Real _sampleRate;
        int _fftSize;

        Bands _bands;
        MatrixXR _centersLinear;

    public:
        /**
            Constructs a Bark bands object with the specified @a lowBand, @a highBand,
            @a bandCount, @a sampleRate, @a fftSize and @a scaleType settings.
            @param lowBand band of the lowest Bark band,
            must be greater than zero 0 and lower than half the sampleRate.
            @param highBand band of the highest Bark band,
            must be greater than zero 0 and lower than half the sampleRate.
            @param bandCount number of Bark bands.
            @param sampleRate sampleRate frequency of the input signal.
            @param fftSize size of the FFT.
        */
        BarkBands(int lowBand = 0, int highBand = 23, Real sampleRate = 44100.0, int fftSize = 1024);

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
            Return the low band.
            The default is 0.
            @sa lowBand, highBand, setLowBand, setHighBand
        */
        Real lowBand() const;

        /**
            Specifies the low @a band of the spectral whitening.
            The given @a band must be in the range of 0 to the sampleRate / 2.
            @sa lowBand, highBand, setHighBand
        */
        void setLowBand( Real band, bool callSetup = true );

        /**
            Return the high band.
            The default is 23.
            @sa lowBand, setLowBand, setHighBand
        */
        Real highBand() const;

        /**
            Specifies the high @a band.
            The given @a band must be in the range of 0 to the sampleRate / 2.
            @sa lowBand, highBand, setLowBand
        */
        void setHighBand( Real band, bool callSetup = true );

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

    };

    Bands::Bands() : 
    _starts(1, 1)
    { 
    _weights.push_back(MatrixXR::Constant(1, 1, 0.0));

    setup();
    }


    Bands::Bands(MatrixXI starts, std::vector<MatrixXR> weights) {
    LOUDIA_DEBUG("BANDS: Constructor starts: " << starts);

    if ( starts.rows() != (int)weights.size() ) {
        // Throw an exception
    }

    if ( starts.cols() != 1 ) {
        // Throw an exception
    }


    for (int i = 0; i < (int)weights.size(); i++){
        if ( weights[i].cols() != 1 ) {
        // Throw an exception
        }
    }

    _starts = starts;
    _weights = weights;

    setup();
    }

    Bands::~Bands() {}

    void Bands::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("BANDS: Setting up...");

    reset();
    LOUDIA_DEBUG("BANDS: Finished set up...");
    }


    void Bands::process(const MatrixXR& spectrum, MatrixXR* bands){

    (*bands).resize(spectrum.rows(), _starts.rows());

    for (int j = 0; j < spectrum.rows(); j++) {
        for (int i = 0; i < _starts.rows(); i++ ) {
        (*bands)(j, i) = spectrum.block(j, _starts(i, 0), 1, _weights[i].rows()).row(0).dot(_weights[i].col(0));
        }
    }
    }

    void Bands::reset(){
    // Initial values
    }

    vector<MatrixXR> Bands::weights() const {
    return _weights;
    }

    void Bands::bandWeights(int band, MatrixXR* bandWeights) const {
    (*bandWeights) =  _weights[ band ];
    }

    void Bands::starts(MatrixXI* result) const {
    (*result) = _starts;
    }

    void Bands::setStartsWeights(const MatrixXI& starts, std::vector<MatrixXR> weights, bool callSetup ) {
    _weights = weights;
    _starts = starts;
    
    if ( callSetup ) setup();
    }

    int Bands::bands() const {
    return _starts.rows();
    }

    BarkBands::BarkBands(int lowBand, int highBand, Real sampleRate, int fftSize) 
    {
    
    LOUDIA_DEBUG("BARKBANDS: Constructor lowBand: " << _lowBand << 
            ", highBand: " << _highBand << 
            ", sampleRate: " << _sampleRate << 
            ", fftSize: " << _fftSize);

    if ( lowBand >= highBand ) {
        // Throw an exception, highBand must be higher than lowBand
    }
    
    setLowBand( lowBand, false );
    setHighBand( highBand, false );
    setSampleRate( sampleRate, false );
    setFftSize( fftSize, false );
    
    setup();
    
    LOUDIA_DEBUG("BARKBANDS: Constructed");
    
    }

    void BarkBands::setup(){
    LOUDIA_DEBUG("BARKBANDS: Setting up...");

    // In some cases the first boundary is set to 0
    MatrixXR startFreqs(25, 1);
    startFreqs << 20, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500;

    MatrixXR centerFreqs(24, 1);
    centerFreqs << 60, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500;
    
    MatrixXI startBins = ((startFreqs.block(_lowBand, 0, _highBand - _lowBand + 2, 1) / _sampleRate) * _fftSize).cast<int>();
    
    std::vector<MatrixXR> weights;
    for (int i = 0; i < startBins.rows() - 1; i++) {
        MatrixXR bandWeights = MatrixXR::Ones(startBins(i+1) - startBins(i), 1);
        weights.push_back(bandWeights);
    }

    _bands.setStartsWeights(startBins.block(0, 0, _highBand - _lowBand + 1, 1), weights);

    LOUDIA_DEBUG("BARKBANDS: Finished set up...");
    }

    void BarkBands::process(const MatrixXR& spectrum, MatrixXR* bands) {  
    _bands.process(spectrum, bands);
    }

    void BarkBands::reset(){
    // Initial values
    _bands.reset();
    }

    void BarkBands::starts(MatrixXI* result) const {
    return _bands.starts( result );
    }

    void BarkBands::centers(MatrixXR* result) const {
    (*result) = _centersLinear;
    }

    std::vector<MatrixXR> BarkBands::weights() const {
    return _bands.weights();
    }

    void BarkBands::bandWeights(int band, MatrixXR* bandWeights) const {
    return _bands.bandWeights( band, bandWeights );
    }

    int BarkBands::bands() const {
    return _bands.bands();
    }

    Real BarkBands::lowBand() const{
    return _lowBand;
    }
    
    void BarkBands::setLowBand( Real band, bool callSetup ){
    _lowBand = band;
    if ( callSetup ) setup();
    }

    Real BarkBands::highBand() const{
    return _highBand;
    }
    
    void BarkBands::setHighBand( Real band, bool callSetup ){
    _highBand = band;
    if ( callSetup ) setup();
    }

    Real BarkBands::sampleRate() const{
    return _sampleRate;
    }
    
    void BarkBands::setSampleRate( Real frequency, bool callSetup ){
    _sampleRate = frequency;
    if ( callSetup ) setup();
    }

    int BarkBands::fftSize() const{
    return _fftSize;
    }

    void BarkBands::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }
}    