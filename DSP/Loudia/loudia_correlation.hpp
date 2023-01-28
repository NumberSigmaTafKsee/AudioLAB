#pragma once

namespace Loudia
{
    /**
    * @class Autocorrelation
    *
    * @brief Algorithm to perform an autocorrelation of vectors of Real values.
    *
    * This class represents an object to perform a correlation of a vector with itself.
    *
    * The correlation can be performed using two methods:
    * -# Direct method 
    * -# FFT method
    *
    * The Direct method consists in applying the correlation formula directly
    * in the time domain.
    *
    * The FFT method consists in performing an Fast Fourier Transform of the
    * vector and multiply it by its conjugate.
    * Finally the algorithm applies an IFFT to the result of the 
    * multiplication in order to obtain the autocorrelation for all
    * the time lags.
    *
    * The Direct method performs faster than the FFT method only
    * on vectors of small sizes. The decision point for selecting one of
    * the two methods depends on the platform.
    *
    * The method performed can be specified using setUseFft().
    *
    * @author Ricard Marxer
    *
    * @sa Correlation, PitchACF
    */
    class Autocorrelation {
    protected:
        // Internal parameters
        int _inputSize;
        int _minLag;
        int _maxLag;
        bool _useFft;

        // Internal variables
        int _calcMinLag;
        int _calcMaxLag;
        FFT _fft;
        IFFT _ifft;
        MatrixXC _tempFft;
        MatrixXR _temp;

    public:
        /**
            Constructs an Autocorrelation object with the specified @a inputSize, 
            @a maxLag and @a minLag settings.
            
            @param inputSize size of the inputs arrays to be autocorrelated,
            must be > 0.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param maxLag maximum lag to be calculated
            @param minLag minimum lag to be calculated
            @param useFft determines whether or not to use the FFT method
        */
        Autocorrelation(int inputSize, int maxLag, int minLag, bool useFft);
        Autocorrelation(int inputSize, int maxLag, int minLag);
        Autocorrelation(int inputSize, int maxLag);
        Autocorrelation(int inputSize = 1024);

        /**
            Destroys the Autocorrelation algorithm and frees its resources.
        */
        ~Autocorrelation();

        void setup();
        void reset();

        /**
            Performs an autocorrelation on each of the rows of @a frames.
            Puts the resulting autocorrelations in the rows of @a autocorrelation.
            
            @param frames matrix of Real values.  The number of columns of @a 
            frames must be equal to the inputSize property.
            
            @param autocorrelation pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a frames and maxLag - minLag columns.
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& frames, MatrixXR* autocorrelation);

        /**
            Returns the size of the input arrays to be autocorrelated.
            The default is 1024.
            
            @sa setInputSize()
        */
        int inputSize() const;  

        /**
            Specifies the @a size of the input.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 the algorithm will perform faster.
            
            @sa inputSize()
        */
        void setInputSize( int size, bool callSetup = true );

        /**
            Returns the minimum lag to be calculated.
            The default is 0.
            
            @sa maxLag(), setMinLag(), setMaxLag()
        */
        int minLag() const;  
        
        /**
            Specifies the minimum @a lag of the autocorrelation.
            The given @a lag will be constratined between -inputSize + 1 and inputSize.
            
            @sa minLag(), maxLag(), setMaxLag()
        */  
        void setMinLag( int lag, bool callSetup = true );

        /**
            Returns the maximum lag to be calculated.
            The default is inputSize.
            
            @sa minLag(), setMinLag(), setMaxLag()
        */
        int maxLag() const;  
        
        /**
            Specifies the maximum @a lag of the autocorrelation.
            The given @a lag will be constratined between -inputSize + 1 and inputSize.
            
            @sa minLag(), maxLag(), setMinLag()
        */
        void setMaxLag( int lag, bool callSetup = true );

        /**
            Returns @c true if the FFT method should be used for the autocorrelation.
            The default is True for inputSize larger than 128; otherwise it is False.
            
            @sa setUseFft()
        */
        bool useFft() const;  
        
        /**
            Specifies whether the autocorrelation should be performed using the FFT method.
            
            @sa useFft()
        */
        void setUseFft( bool useFft, bool callSetup = true );
    };

    /**
    * @class Correlation
    *
    * @brief Algorithm to perform the correlation between two vectors of Real values.
    *
    * This class represents an object to perform a correlation of two vectors.
    *
    * The correlation can be performed using two methods:
    * -# Direct method 
    * -# FFT method
    *
    * The Direct method consists in applying the correlation formula directly
    * in the time domain.
    *
    * The FFT method consists in performing an Fast Fourier Transform of each
    * of the vectors and multiply the first by the conjugate of the second.
    * Finally the algorithm applies
    * an IFFT to the result of the multiplication in order to obtain the 
    * autocorrelation for all the time lags.
    *
    * The Direct method performs faster than the FFT method only
    * on vectors of small sizes. The decision point for selecting one of
    * the two methods depends on the platform.
    *
    * The method performed can be specified using setUseFft().
    *
    * @author Ricard Marxer
    *
    * @sa Autocorrelation
    */
    class Correlation {
    protected:
        // Internal parameters
        int _inputSizeA;
        int _inputSizeB;
        int _minLag;
        int _maxLag;
        bool _useFft;
        int _fftSize;

        // Internal variables
        FFT _fft;
        IFFT _ifft;

        MatrixXC _fftA;
        MatrixXC _fftB;
        MatrixXR _result;

    public:
        /**
            Constructs an Autocorrelation object with the specified @a inputSize, 
            @a maxLag and @a minLag settings.
            
            @param inputSizeA size of the first input arrays to be autocorrelated,
            must be > 0.
            The algorithm performs faster for sizes which are a power of 2.
            @param inputSizeB size of the second input arrays to be autocorrelated,
            must be > 0.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param maxLag maximum lag to be calculated
            @param minLag minimum lag to be calculated
            @param useFft determines whether or not to use the FFT method
        */
        //Correlation();
        Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag, bool useFft);
        Correlation(int inputSizeA = 1024, int inputSizeB = 1024, int maxLag = std::numeric_limits<int>::max(), int minLag = -std::numeric_limits<int>::max());
        
        /**
            Destroys the Correlation algorithm and frees its resources.
        */
        ~Correlation();

        void setup();
        void reset();

        /**
            Performs a Correlation between each of the rows of @a framesA and
            each of the rows of @b framesB respectively.
            Puts the resulting correlations in the rows of @a correlation.
            
            @param framesA matrix of Real values.  The number of columns of @a 
            framesA must be equal to the inputSizeA property.
            @param framesB matrix of Real values.  The number of columns of @a 
            framesB must be equal to the inputSizeB property.
            Note that @a framesA and @a framesB should have the same number of rows.
            @param correlation pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a framesA (and @a framesB) and maxLag - minLag columns.
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& framesA, const MatrixXR& framesB, MatrixXR* correlation);

        /**
            Returns the size of the first input arrays to be correlated.
            The default is 1024.
            
            @sa inputSizeB(), setInputSizeA(), setInputSizeB()
        */
        int inputSizeA() const;  

        /**
            Specifies the @a size of the first of the input arrays.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 the algorithm will perform faster.
            
            @sa inputSizeA(), inputSizeB(), setInputSizeB()
        */
        void setInputSizeA( int size, bool callSetup = true );

        /**
            Returns the size of the second of the input arrays to be correlated.
            The default is 1024.
            
            @sa inputSizeA(), setInputSizeA(), setInputSizeB()
        */
        int inputSizeB() const;  

        /**
            Specifies the @a size of the second of the input arrays.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 the algorithm will perform faster.
            
            @sa inputSizeA(), inputSizeB(), setInputSizeA()
        */
        void setInputSizeB( int size, bool callSetup = true );

        /**
            Returns the minimum lag to be calculated.
            The default is -max(_inputSizeA, _inputSizeB) + 1.
            
            @sa maxLag(), setMinLag(), setMaxLag()
        */
        int minLag() const;  
        
        /**
            Specifies the minimum @a lag of the ocorrelation.
            The given @a lag will be constratined between - max( inputSizeA, inputSizeB ) + 1 and min( inputSizeA, inputSizeB ).
            Note that the lag should be smaller than the maxLag.
            @sa minLag(), maxLag(), setMaxLag()
        */  
        void setMinLag( int lag, bool callSetup = true );

        /**
            Returns the maximum lag to be calculated.
            The default is inputSize.
            
            @sa minLag(), setMinLag(), setMaxLag()
        */
        int maxLag() const;  
        
        /**
            Specifies the maximum @a lag of the correlation.
            The given @a lag will be constratined between - max( inputSizeA, inputSizeB ) + 1 and min( inputSizeA, inputSizeB ).
            Note that the lag should be larger than the maxLag.
            @sa minLag(), maxLag(), setMinLag()
        */
        void setMaxLag( int lag, bool callSetup = true );

        /**
            Returns @c true if the FFT method should be used for the correlation.
            The default is True for inputSize larger than 128; otherwise it is False.
            
            @sa setUseFft()
        */
        bool useFft() const;  
        
        /**
            Specifies whether the autocorrelation should be performed using the FFT method.
            
            @sa useFft()
        */
        void setUseFft( bool useFft, bool callSetup = true );
    };

    Autocorrelation::Autocorrelation(int inputSize)
    {
    LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize);

    setInputSize( inputSize, false );
    setMinLag( 0, false );
    setMaxLag( inputSize, false );
    setUseFft( (_maxLag - _minLag) > 128, false );

    setup();
    }

    Autocorrelation::Autocorrelation(int inputSize, int maxLag)
    {
    LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
            << " maxLag: " << maxLag);

    setInputSize( inputSize, false );
    setMinLag( 0, false );
    setMaxLag( maxLag, false );
    setUseFft( (_maxLag - _minLag) > 128, false );

    setup();
    }


    Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag)
    {
    LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
            << " minLag: " << minLag
            << " maxLag: " << maxLag);

    setInputSize( inputSize, false );
    setMinLag( minLag, false );
    setMaxLag( maxLag, false );
    setUseFft( (_maxLag - _minLag) > 128, false );

    setup();
    }


    Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag, bool useFft)
    {
    LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
            << " minLag: " << minLag
            << " maxLag: " << maxLag
            << " useFft: " << useFft);
    
    setInputSize( inputSize, false );
    setMinLag( minLag, false );
    setMaxLag( maxLag, false );
    setUseFft( useFft, false );
    
    setup();
    }

    Autocorrelation::~Autocorrelation(){}

    void Autocorrelation::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("AUTOCORRELATION: Setting up...");

    _calcMinLag = min(_inputSize + 1, max(-_inputSize + 1, _minLag));  
    _calcMaxLag = min(_inputSize + 1, max(-_inputSize + 1, _maxLag));

    if ( _useFft ) {
        _fft.setFftSize( nextPowerOf2( _inputSize * 2, false ) );
        _fft.setZeroPhase( false, false );
        
        _ifft.setFftSize( nextPowerOf2( _inputSize * 2, false ) );
        _ifft.setZeroPhase( false, false );

        _fft.setup();
        _ifft.setup();
    }
    
    reset();
    
    LOUDIA_DEBUG("AUTOCORRELATION: Finished setup.");
    }

    void Autocorrelation::process(const MatrixXR& frames, MatrixXR* autocorrelation){
    const int rows = frames.rows();

    (*autocorrelation).resize(rows, _maxLag - _minLag);
    
    (*autocorrelation).setZero();
    
    if ( _useFft ) {
        _fft.process(frames, &_tempFft);
        
        _tempFft.array() *= _tempFft.conjugate().array();
        
        _ifft.process(_tempFft, &_temp);

        (*autocorrelation).block(0, _calcMinLag - _minLag, rows, _calcMaxLag - _calcMinLag) = _temp.block(0, 0, rows, _calcMaxLag - _calcMinLag);

    } else {
        correlate(frames, frames, &_temp, _calcMinLag, _calcMaxLag);
        
        (*autocorrelation).block(0, _calcMinLag - _minLag, rows, _calcMaxLag - _calcMinLag) = _temp;
        
    }
    }

    void Autocorrelation::reset(){
    // Initial values
    }

    int Autocorrelation::inputSize() const {
    return _inputSize;
    }
    
    void Autocorrelation::setInputSize( int size, bool callSetup ) {
    _inputSize = size;
    if ( callSetup ) setup();
    }

    int Autocorrelation::minLag() const {
    return _minLag;
    }
    
    void Autocorrelation::setMinLag( int lag, bool callSetup ) {
    _minLag = lag;
    if ( callSetup ) setup();
    }

    int Autocorrelation::maxLag() const {
    return _maxLag;
    }
    
    void Autocorrelation::setMaxLag( int lag, bool callSetup ) {
    _maxLag = lag;
    if ( callSetup ) setup();
    }

    bool Autocorrelation::useFft() const {
    return _useFft;
    }  

    void Autocorrelation::setUseFft( bool useFft, bool callSetup ) {
    _useFft = useFft;
    if ( callSetup ) setup();
    }

    /*
    Correlation::Correlation()
    {
    int inputSizeA = 1024;
    int inputSizeB = 1024;
    setInputSizeA( inputSizeA, false );
    setInputSizeB( inputSizeB, false );
    setMinLag( -std::numeric_limits<int>::max(), false );
    setMaxLag( std::numeric_limits<int>::max(), false );
    setUseFft( (_maxLag - _minLag) > 128, false );
    setup();
    }
    */
    Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag)
    {
    LOUDIA_DEBUG("CORRELATION: Construction inputSizeA: " << inputSizeA
                << " inputSizeB: " << inputSizeB
                << " minLag: " << minLag
                << " maxLag: " << maxLag);

    setInputSizeA( inputSizeA, false );
    setInputSizeB( inputSizeB, false );
    setMinLag( minLag, false );
    setMaxLag( maxLag, false );
    setUseFft( (_maxLag - _minLag) > 128, false );

    setup();
    }


    Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag, bool useFft)
    {
    LOUDIA_DEBUG("CORRELATION: Construction inputSizeA: " << inputSizeA
                << " inputSizeB: " << inputSizeB
                << " minLag: " << minLag
                << " maxLag: " << maxLag
                << " useFft: " << useFft);
    
    setInputSizeA( inputSizeA, false );
    setInputSizeB( inputSizeB, false );
    setMinLag( minLag, false );
    setMaxLag( maxLag, false );
    setUseFft( useFft, false );
    
    setup();
    }

    Correlation::~Correlation(){}

    void Correlation::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("CORRELATION: Setting up...");

    if ( _useFft ) {
        _fftSize = nextPowerOf2( ( (_inputSizeA + _inputSizeB) - 1) * 2 );
        
        _fft.setFftSize( _fftSize, false );
        _fft.setZeroPhase( false, false );
        _fft.setup();

        _ifft.setFftSize( _fftSize, false );
        _ifft.setZeroPhase( false, false );
        _ifft.setup();
    }
    
    reset();
    
    LOUDIA_DEBUG("CORRELATION: Finished setup.");
    }

    void Correlation::process(const MatrixXR& inputA, const MatrixXR& inputB, MatrixXR* correlation){
    const int rows = inputA.rows();

    if ( rows != inputB.rows() ) {
        // Thorw ValueError rows of A and B must be the same
    }
    
    (*correlation).resize(rows, _maxLag - _minLag);
    
    if ( _useFft ) {
        
        _fft.process(inputA, &_fftA);
        _fft.process(inputB, &_fftB);
            
        _ifft.process(_fftA.array() * _fftB.conjugate().array(), &_result);
        
        // TODO: use Eigen rowwise().shift(_fftSize - 2) when it will exist
        for(int i = _minLag;  i < _maxLag; i++ ){
        (*correlation).col(i - _minLag) = _result.col(((_fftSize-2) + (i - _minLag)) % _fftSize);
        }

    } else {
        
        correlate(inputA, inputB, correlation, _minLag, _maxLag);
        
    }
    }

    void Correlation::reset(){
    // Initial values
    }

    int Correlation::inputSizeA() const {
    return _inputSizeA;
    }
    
    void Correlation::setInputSizeA( int size, bool callSetup ) {
    _inputSizeA = size;
    if ( callSetup ) setup();
    }

    int Correlation::inputSizeB() const {
    return _inputSizeB;
    }
    
    void Correlation::setInputSizeB( int size, bool callSetup ) {
    _inputSizeB = size;
    if ( callSetup ) setup();
    }

    int Correlation::minLag() const {
    return _minLag;
    }
    
    void Correlation::setMinLag( int lag, bool callSetup ) {
    if ( lag >= _maxLag ) {
        // Thorw ValueError, "The minLag should be smaller than the maxLag."
    }

    _minLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
    _minLag = min( min(_inputSizeA, _inputSizeB), _minLag);  
    if ( callSetup ) setup();
    }

    int Correlation::maxLag() const {
    return _maxLag;
    }
    
    void Correlation::setMaxLag( int lag, bool callSetup ) {
    if ( lag <= _minLag ) {
        // Thorw ValueError, "The maxLag should be larger than the minLag."
    }
    
    _maxLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
    _maxLag = min( min(_inputSizeA, _inputSizeB), _maxLag);  
    if ( callSetup ) setup();
    }

    bool Correlation::useFft() const {
    return _useFft;
    }  

    void Correlation::setUseFft( bool useFft, bool callSetup ) {
    _useFft = useFft;
    if ( callSetup ) setup();
    }
}