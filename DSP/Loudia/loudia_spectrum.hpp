#pragma once

namespace Loudia
{
    /**
    * @class FFT
    *
    * @brief Algorithm to perform a Fast Fourier Transform of a vector of Real values.
    *
    * This class represents an object to perform Fast Fourier Transforms (FFT) on Real data.
    * The FFT is a fast implementation of a Discrete Fourier Transform (DFT).
    * The algorithm takes as input N-point vectors of Real values (N being the frame size)
    * and returns (M / 2 + 1) point vectors of Complex values (M being the FFT size).
    *
    * Note that the algorithm works fastest when M is a power of 2.
    *
    * Since the input is Real valued, the FFT will be symmetric
    * and only half of the output is needed.
    * This processing unit will only return the (M / 2 + 1)-point array 
    * corresponding to positive frequencies of the FFT.
    *
    * When M is different than N the input data is zero padded at the end.
    * Alternatively the algorithm can perform an N/2 rotation and zero pad the center
    * before the FFT to allow a zero phase transform.
    * This is done by using the setZeroPhase() method.
    *
    * @author Ricard Marxer
    *
    * @sa FFTComplex, IFFT, IFFTComplex
    */
    class FFT{
        protected:
        int _fftSize;
        bool _zeroPhase;

        int _halfSize;
        
        Real* _in;
        fftwf_complex* _out;
        fftwf_plan _fftplan;
        

        public:
        /**
            Constructs an FFT object with the specified @a fftSize and @a
            zeroPhase setting.
            
            @param fftSize size of the FFT transform must be > 0, 
            it is the target size of the transform.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param zeroPhase determines whether
            or not to perform the zero phase transform.
        */
        FFT(int fftSize = 1024, bool zeroPhase = true);
        
        /**
            Destroys the FFT algorithm and frees its resources.
        */
        ~FFT();
        
        /**
            Performs a Fast Fourier Transform on each of the rows of @a frames and
            puts the resulting FFT in the rows of @a fft.
            
            @param frames matrix of Real values.  The number of columns of @a frames must
            be smaller or equal to the fftSize property.
            
            @param fft pointer to a matrix of Complex values for the output.  The matrix should
            have the same number of rows as @a frames and (fftSize / 2) + 1 columns. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& frames, MatrixXC* fft);

        void setup();
        void reset();

        /**
            Returns the size of the FFT to be performed.  The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT to be performed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Returns the zero phase setting.  The default is True.
            
            @sa setZeroPhase()
        */
        bool zeroPhase() const;

        /**
            Specifies the @a zeroPhase setting.
            
            @sa zeroPhase()
        */
        void setZeroPhase( bool zeroPhase, bool callSetup = true );

    };    

    /**
    * @class FFTComplex
    *
    * @brief Algorithm to perform a Fast Fourier Transform of a vector of Complex values.
    *
    * This class represents an object to perform Fast Fourier Transforms (FFT) on Real data.
    * The FFT is a fast implementation of a Discrete Fourier Transform (DFT).
    * The algorithm takes as input N point vectors of Real values (N being the frame size) 
    * and returns M point vectors of Complex values (M being the FFT size).
    *
    * Note that the algorithm works fastest when M is a power of 2.
    *
    * When M is different than N the input data is zero padded at the end.
    * Alternatively the algorithm can perform an N/2 rotation and zero pad the center
    * before the FFT to allow a zero phase transform.
    * This is done by using the setZeroPhase() method.
    *
    * @author Ricard Marxer
    *
    * @sa FFTComplex, IFFT, IFFTComplex
    */
    class FFTComplex{
    protected:
        int _frameSize;
        int _fftSize;
        bool _zeroPhase;

        fftwf_complex* _in;
        
        fftwf_complex* _out;

        fftwf_plan _fftplan;
        
        template <typename FrameMatrixType>
        void process(const FrameMatrixType& frames, MatrixXC* fft);


    public:
        /**
            Constructs an FFT object with the specified @a fftSize and @a
            zeroPhase setting.
            
            @param frameSize size of the frame must be > 0, 
            it is the size of the input frames.
            
            @param fftSize size of the FFT transform must be > 0, 
            it is the target size of the transform.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param zeroPhase determines whether
            or not to perform the zero phase transform.
        */
        FFTComplex(int frameSize, int fftSize, bool zeroPhase = true);
        
        /**
            Destroys the algorithm and frees its resources.
        */
        ~FFTComplex();
        
        /**
            Performs a Fast Fourier Transform on each of the rows of @a frames and
            puts the resulting FFT in the rows of @a fft.
            
            @param frames matrix of Real values.  The number of columns of @a frames must
            be equal to the frameSize property.
            
            @param fft pointer to a matrix of Complex values for the output.  The matrix should
            have the same number of rows as @a frames and (fftSize / 2) + 1 columns. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& frames, MatrixXC* fft);
        void process(const MatrixXR& frames, MatrixXC* fft);
        
        void setup();
        void reset();

        /**
            Returns the size of the FFT to be performed.  The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT to be performed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Returns the size of the frame to be processed.
            The default is 1024.
            
            @sa setFrameSize()
        */
        int frameSize() const;

        /**
            Specifies the @a size of the frame to be processed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa frameSize()
        */
        void setFrameSize( int size, bool callSetup = true );

        /**
            Returns the zero phase setting.  The default is True.
            
            @sa setZeroPhase()
        */
        bool zeroPhase() const;

        /**
            Specifies the @a zeroPhase setting.
            
            @sa zeroPhase()
        */
        void setZeroPhase( bool zeroPhase, bool callSetup = true );
    };

    /**
    * @class IFFT
    *
    * @brief Algorithm to perform an Inverse Fast Fourier Transform of a vector of Complex 
    * values representing the positive frequencies half of a symmetric FFT.
    *
    * The IFFT is a fast implementation of an Inverse Discrete Fourier Transform (IDFT).
    * The algorithm takes as input (M / 2 + 1) point vectors of Complex 
    * values (M being the FFT size), and returns N point vectors of Real 
    * values (N being the frame size).
    *
    * The input of the IFFT is assumed to be the positive frequencies half of an M point
    * magnitude symmetric and phase antisymmetric FFT.  Therefore the result is a Real value 
    * vector.
    *
    * Note that N can be smaller than M.
    * In this case the last ( M - N ) coefficients
    * will be discarded, since it assumes that zero padding has been made
    * at the end of the frame prior to the forward FFT transform.
    *
    * Alternatively the algorithm can undo the center zeropadding and
    * the N/2 rotation if done durnig the FFT forward transform.
    * This is specified by using the setZeroPhase() method.
    *
    * @author Ricard Marxer
    *
    * @sa FFT, FFTComplex, IFFTComplex
    */
    class IFFT{
    protected:
        int _fftSize;
        bool _zeroPhase;

        int _halfSize;

        fftwf_complex* _in;
        Real* _out;

        fftwf_plan _fftplan;
    

    public:
        /**
            Constructs an IFFT object with the specified @a fftSize and @a
            zeroPhase setting.
            
            @param fftSize size of the IFFT transform must be > 0, 
            it is the target size of the transform.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param zeroPhase specifies whether
            or not the zero phase method was performed.
        */
        IFFT(int fftSize = 1024, bool zeroPhase = true);
        
        /**
            Destroys the IFFT algorithm and frees its resources.
        */
        ~IFFT();
        
        /**
            Performs a Inverse Fast Fourier Transform on each of the rows of @a fft and
            puts the resulting IFFT in the rows of @a frames.
            
            @param fft matrix of Complex values.  The number of columns of @a fft must
            be equal to the (fftSize / 2) + 1, 
            where fftSize is parameter of the constructor or specified by setFftSize().
            
            @param frame pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a fft and fftSize columns.  
            Note that if the zeroPhase setting is true, the resulting IFFT transforms
            will be rotated to compensate for Zero Phase method that may have been performed
            when the FFT had been done.
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& fft, MatrixXR* frame);
        
        void setup();
        void reset();

        /**
            Returns the size of the FFT to be performed.  The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the IFFT to be performed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Returns the zero phase setting.  The default is True.
            
            @sa setZeroPhase()
        */
        bool zeroPhase() const;

        /**
            Specifies the @a zeroPhase setting.
            
            @sa zeroPhase()
        */
        void setZeroPhase( bool zeroPhase, bool callSetup = true );
    };

    /**
    * @class IFFTComplex
    *
    * @brief Algorithm to perform an Inverse Fast Fourier Transform of a vector of Complex 
    * values representing the full FFT.
    *
    * The IFFT is a fast implementation of an Inverse Discrete Fourier Transform (IDFT).
    * The algorithm takes as input M point vectors of Complex 
    * values (M being the FFT size), and returns N point vectors of Real 
    * values (N being the frame size).
    *
    * Note that N can be smaller than M.
    * In this case the last ( M - N ) coefficients
    * will be discarded, since it assumes that zero padding has been made
    * at the end of the frame prior to the forward FFT transfor.
    *
    * Alternatively the algorithm can undo the center zeropadding and
    * the N/2 rotation if done durnig the FFT forward transform.
    * This is specified by using the setZeroPhase() method.
    *
    * @author Ricard Marxer
    *
    * @sa FFT, FFTComplex, IFFT
    */
    class IFFTComplex{
    protected:
        int _fftSize;
        int _frameSize;
        bool _zeroPhase;

        fftwf_complex* _in;
        fftwf_complex* _out;

        fftwf_plan _fftplan;
        
        template <typename FrameMatrixType>
        void process(const FrameMatrixType& ffts, MatrixXC* frames);

    public:
        IFFTComplex(int fftSize, int frameSize, bool zeroPhase = true);
        ~IFFTComplex();
        
        void process(const MatrixXC& ffts, MatrixXC* frames);
        void process(const MatrixXR& ffts, MatrixXC* frames);
        
        void setup();
        void reset();

        /**
            Returns the size of the FFT to be processed.
            The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT to be processed.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Returns the size of the target frame.
            The default is 1024.
            
            @sa setFrameSize()
        */
        int frameSize() const;

        /**
            Specifies the @a size of the target frame.
            The given @a size must be higher than 0.
            Note that if @a size is a power of 2 will perform faster.
            
            @sa frameSize()
        */
        void setFrameSize( int size, bool callSetup = true );

        /**
            Returns the zero phase setting.
            The default is true.
            
            @sa setZeroPhase()
        */
        bool zeroPhase() const;

        /**
            Specifies the @a zeroPhase setting.
            
            @sa zeroPhase()
        */
        void setZeroPhase( bool zeroPhase, bool callSetup = true );
    };

    /**
    * @class DCT
    *
    * @brief Algorithm to perform a Discrete Cosine Transform of a vector of Real values.
    *
    * This class represents an object to perform Discrete Cosine Transform (DCT) on Real data.
    * The algorithm takes as input N-point vectors of Real values 
    * and returns M-point vectors of Real values.
    *
    * 5 types of DCT are implemented:
    * -# Type I
    * -# Type II
    * -# Type III
    * -# Type IV
    * -# Octave's Implementation
    *
    * The DCT type can be selected using the 
    * setDCTType() taking as argument a DCTType.
    *
    *
    * @author Ricard Marxer
    *
    * @sa FFT
    */
    class DCT {
    public:
    /**
        @enum DCTType
        @brief Specifies the type of the DCT.
        @sa dctType
    */
    enum DCTType {
        I = 0 /**< DCT Type-I */,
        II = 1 /**< DCT Type-II */,
        III = 2 /**< DCT Type-III */,
        IV = 3 /**< DCT Type-IV */,
        OCTAVE = 4 /**< Octave's implementation */
    };

    protected:
        // Internal parameters
        int _inputSize;
        int _dctSize;
        Real _scale;

        DCTType _dctType;

        // Internal variables
        MatrixXR _dctMatrix;

        void type1Matrix(MatrixXR* dctMatrix);

        void type2Matrix(MatrixXR* dctMatrix);

        void typeOctaveMatrix(MatrixXR* dctMatrix);

    public:
        /**
            Constructs a DCT object with the given @a inputSize, @a dctSize,
            @a scale, @a dctType parameters.
        */
        DCT(int inputSize = 1024, int dctSize = 1024, bool scale = false, DCTType dctType = OCTAVE);
        
        /**
            Destroys the algorithm and frees its resources.
        */
        ~DCT();

        void reset();
        void setup();

        /**
            Performs a Discrete Cosine Transform on each of the rows of @a frames and
            puts the resulting DCT in the rows of @a dct.
            
            @param frames matrix of Real values.  The number of columns of @a frames must
            be equal to the inputSize.
            
            @param dct pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a frames and dctSize columns. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& frames, MatrixXR* dct);

        /**
            Returns the type of the DCT
            
            By default it is OCTAVE.
        */
        DCTType dctType() const;

        /**
            Specifies the type of the DCT
        */
        void setDctType( DCTType type, bool callSetup = true );

        /**
            Returns the size of the input arrays.
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
            Returns the output size of the DCT
            
            Note that the result will when performing
            the DCT at the most inputSize coefficients
            will be outputed.
            
            By default it is 1024.
        */
        int dctSize() const;

        /**
            Specifies the output size of the DCT
            
            Note that the result will when performing
            the DCT at the most inputSize coefficients
            will be outputed.
            
            By default it is 1024.
        */
        void setDctSize( int size, bool callSetup = true );
    };


    FFT::FFT(int fftSize, bool zeroPhase) :
    _in( NULL ),
    _out( NULL ),
    _fftplan( NULL )
    {
    LOUDIA_DEBUG("FFT: Constructor fftSize: " << fftSize 
            << ", zeroPhase: " << zeroPhase);

    setFftSize( fftSize, false );
    setZeroPhase( zeroPhase, false );
    
    setup();
    
    LOUDIA_DEBUG("FFT: Constructed");
    }

    FFT::~FFT(){
    LOUDIA_DEBUG("FFT: Destroying...");
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    LOUDIA_DEBUG("FFT: Destroyed out");
    }

    void FFT::setup(){
    LOUDIA_DEBUG("FFT: Setting up...");
    
    // Free the ressources if needed 
    // before setting them up
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    
    _halfSize = ( _fftSize / 2 ) + 1;
    
    // Allocate the ressources needed
    _in = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
    _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize);
    
    _fftplan = fftwf_plan_dft_r2c_1d( _fftSize, _in, _out,
                                        FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
        
    LOUDIA_DEBUG("FFT: Finished set up...");
    }

    void FFT::process(const MatrixXR& frames, MatrixXC* ffts){
    const int cols = frames.cols();
    const int rows = frames.rows();

    if(_fftSize < cols){
        // Throw exception, the FFT size must be greater or equal than the input size
    }
    
    (*ffts).resize(rows, _halfSize);

    for (int i = 0; i < rows; i++){    
        // Fill the buffer with zeros
        Eigen::Map<MatrixXR>(_in, 1, _fftSize) = MatrixXR::Zero(1, _fftSize);
        
        // Put the data in _in
        if(_zeroPhase){

        int half_plus = (int)ceil((Real)cols / 2.0);
        int half_minus = (int)floor((Real)cols / 2.0);
        
        // Put second half of the frame at the beginning 
        Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, half_plus) = frames.row(i).block(0, half_minus, 1, half_plus);
        
        // and first half of the frame at the end
        Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) = frames.row(i).block(0, 0, 1, half_minus);


        }else{
        // Put all of the frame at the beginning
        Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, cols) = frames.row(i);
        }

        // Process the data
        fftwf_execute(_fftplan);

        // Take the data from _out
        (*ffts).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _halfSize);
    }
    }

    void FFT::reset(){
    }

    int FFT::fftSize() const{
    return _fftSize;
    }

    void FFT::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }

    bool FFT::zeroPhase() const{
    return _zeroPhase;
    }

    void FFT::setZeroPhase( bool zeroPhase, bool callSetup ) {
    _zeroPhase = zeroPhase;
    if ( callSetup ) setup();
    }

    FFTComplex::FFTComplex(int frameSize, int fftSize, bool zeroPhase) :
    _in( NULL ),
    _out( NULL ),
    _fftplan( NULL )
    {
    LOUDIA_DEBUG("FFTComplex: Constructor frameSize: " << frameSize 
            << ", fftSize: " << fftSize 
            << ", zeroPhase: " << zeroPhase);

    if(_fftSize < _frameSize){
        // Throw exception, the FFTComplex size must be greater or equal than the input size
    }

    setFrameSize( frameSize, false );
    setFftSize( fftSize, false );
    setZeroPhase( zeroPhase, false );
    
    setup();
    
    LOUDIA_DEBUG("FFTComplex: Constructed");
    }

    FFTComplex::~FFTComplex(){
    LOUDIA_DEBUG("FFT: Destroying...");
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    LOUDIA_DEBUG("FFT: Destroyed out");
    }

    void FFTComplex::setup(){
    LOUDIA_DEBUG("FFTComplex: Setting up...");

    // Free the ressources if needed 
    // before setting them up
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    
    _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
    _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
        
    _fftplan = fftwf_plan_dft_1d( _fftSize, _in, _out,
                                    FFTW_FORWARD, FFTW_ESTIMATE );
    
    LOUDIA_DEBUG("FFTComplex: Finished set up...");
    }

    template<typename FrameMatrixType>
    void FFTComplex::process(const FrameMatrixType& frames, MatrixXC* ffts){
    (*ffts).resize(frames.rows(), _fftSize);

    for (int i = 0; i < frames.rows(); i++){    
        // Fill the buffer with zeros
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize) = MatrixXC::Zero(1, _fftSize);
        
        // Put the data in _in
        if(_zeroPhase){

        int half_plus = (int)ceil((Real)_frameSize / 2.0);
        int half_minus = (int)floor((Real)_frameSize / 2.0);

        // Put second half of the frame at the beginning 
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, 0, 1, half_plus) = frames.row(i).block(0, half_minus, 1, half_plus).template cast<Complex>();
        
        // and first half of the frame at the end
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) = frames.row(i).block(0, 0, 1, half_minus).template cast<Complex>();


        }else{

        // Put all of the frame at the beginning
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, 0, 1, _frameSize) = frames.row(i).template cast<Complex>();
        }
        // Process the data
        fftwf_execute(_fftplan);

        // Take the data from _out
        (*ffts).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _fftSize);
    }
    }

    void FFTComplex::process(const MatrixXR& frames, MatrixXC* ffts){
    process<MatrixXR>(frames, ffts);
    }

    void FFTComplex::process(const MatrixXC& frames, MatrixXC* ffts){
    process<MatrixXC>(frames, ffts);
    }

    void FFTComplex::reset(){
    }

    int FFTComplex::fftSize() const{
    return _fftSize;
    }

    void FFTComplex::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }

    int FFTComplex::frameSize() const{
    return _frameSize;
    }

    void FFTComplex::setFrameSize( int size, bool callSetup ) {
    _frameSize = size;
    if ( callSetup ) setup();
    }

    bool FFTComplex::zeroPhase() const{
    return _zeroPhase;
    }

    void FFTComplex::setZeroPhase( bool zeroPhase, bool callSetup ) {
    _zeroPhase = zeroPhase;
    if ( callSetup ) setup();
    }

    IFFT::IFFT(int fftSize, bool zeroPhase) :
    _in( NULL ),
    _out( NULL ),
    _fftplan( NULL )
    {
    LOUDIA_DEBUG("IFFT: Constructor fftSize: " << fftSize
            << ", zeroPhase: " << zeroPhase);

    setFftSize( fftSize, false );
    setZeroPhase( zeroPhase, false );
    
    setup();
    
    LOUDIA_DEBUG("IFFT: Constructed");
    }

    IFFT::~IFFT(){
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    }

    void IFFT::setup(){
    LOUDIA_DEBUG("IFFT: Setting up...");

    // Free the ressources if needed 
    // before setting them up
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    
    _halfSize = ( _fftSize / 2 ) + 1;

    // Allocate the ressources needed  
    _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize );
    _out = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
    
    _fftplan = fftwf_plan_dft_c2r_1d( _fftSize, _in, _out,
                                        FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
        
    LOUDIA_DEBUG("IFFT: Finished set up...");
    }

    void IFFT::process(const MatrixXC& ffts, MatrixXR* frames){
    const int rows = ffts.rows();
    
    (*frames).resize(rows, _fftSize);

    for (int i = 0; i < rows; i++){    
        // Fill the buffer with zeros
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _halfSize) = ffts.row(i);
        
        // Process the data
        fftwf_execute(_fftplan);

        // Take the data from _out
        if(_zeroPhase){

        int half_plus = (int)ceil((Real)_fftSize / 2.0);
        int half_minus = (int)floor((Real)_fftSize / 2.0);
        
        // Take second half of the frame from the beginning 
        (*frames).row(i).block(0, half_minus, 1, half_plus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, 0, 1, half_plus) / _fftSize;
        
        // and first half of the frame from the end
        (*frames).row(i).block(0, 0, 1, half_minus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) / _fftSize;

        }else{
        
        // Take all of the frame from the beginning
        (*frames).row(i) = Eigen::Map<MatrixXR>(_out, 1, _fftSize) / _fftSize;

        }
    }
    }

    void IFFT::reset(){
    }

    int IFFT::fftSize() const{
    return _fftSize;
    }

    void IFFT::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }

    bool IFFT::zeroPhase() const{
    return _zeroPhase;
    }

    void IFFT::setZeroPhase( bool zeroPhase, bool callSetup ) {
    _zeroPhase = zeroPhase;
    if ( callSetup ) setup();
    }

    IFFTComplex::IFFTComplex(int frameSize, int fftSize, bool zeroPhase) :
    _in( NULL ),
    _out( NULL ),
    _fftplan( NULL )
    {
    LOUDIA_DEBUG("IFFTComplex: Constructor frameSize: " << frameSize 
            << ", fftSize: " << fftSize 
            << ", zeroPhase: " << zeroPhase);

    if(_fftSize < _frameSize){
        // Throw exception, the IFFTComplex size must be greater or equal than the input size
    }
    
    setup();
    
    LOUDIA_DEBUG("IFFTComplex: Constructed");
    }

    IFFTComplex::~IFFTComplex(){
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    }

    void IFFTComplex::setup(){
    LOUDIA_DEBUG("IFFTComplex: Setting up...");
    
    // Free the ressources if needed 
    // before setting them up
    if ( _fftplan ) {
        LOUDIA_DEBUG("FFT: Destroying plan");
        fftwf_destroy_plan( _fftplan );
    }

    if ( _in ) {
        LOUDIA_DEBUG("FFT: Destroying in");
        fftwf_free( _in ); 
    }

    if ( _out ) {
        LOUDIA_DEBUG("FFT: Destroying out");
        fftwf_free( _out );
    }
    
    _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
    _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
        
    _fftplan = fftwf_plan_dft_1d( _fftSize, _in, _out,
                                    FFTW_BACKWARD, FFTW_ESTIMATE );
    
    LOUDIA_DEBUG("IFFTComplex: Finished set up...");
    }

    template<typename FrameMatrixType>
    void IFFTComplex::process(const FrameMatrixType& ffts, MatrixXC* frames){
    const int rows = ffts.rows();

    (*frames).resize(rows, _fftSize);

    for (int i = 0; i < rows; i++){
        
        // Fill the buffer with zeros
        Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize) = ffts.row(i).template cast<Complex>();
        
        // Process the data
        fftwf_execute(_fftplan);

        // Take the data from _out
        if(_zeroPhase){

        int half_plus = (int)ceil((Real)_frameSize / 2.0);
        int half_minus = (int)floor((Real)_frameSize / 2.0);
        
        // Take second half of the frame from the beginning 
        (*frames).row(i).block(0, half_minus, 1, half_plus) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _fftSize).block(0, 0, 1, half_plus) / _fftSize;
        
        // and first half of the frame from the end
        (*frames).row(i).block(0, 0, 1, half_minus) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) / _fftSize;

        }else{

        // Take all of the frame from the beginning
        (*frames).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _fftSize).block(0, 0, 1, _frameSize) / _fftSize;
        }
    }
    }

    void IFFTComplex::process(const MatrixXR& ffts, MatrixXC* frames){
    process<MatrixXR>(ffts, frames);
    }

    void IFFTComplex::process(const MatrixXC& ffts, MatrixXC* frames){
    process<MatrixXC>(ffts, frames);
    }

    void IFFTComplex::reset(){
    }

    int IFFTComplex::fftSize() const{
    return _fftSize;
    }

    void IFFTComplex::setFftSize( int size, bool callSetup ) {
    _fftSize = size;
    if ( callSetup ) setup();
    }

    int IFFTComplex::frameSize() const{
    return _frameSize;
    }

    void IFFTComplex::setFrameSize( int size, bool callSetup ) {
    _frameSize = size;
    if ( callSetup ) setup();
    }

    bool IFFTComplex::zeroPhase() const{
    return _zeroPhase;
    }

    void IFFTComplex::setZeroPhase( bool zeroPhase, bool callSetup ) {
    _zeroPhase = zeroPhase;
    if ( callSetup ) setup();
    }


    DCT::DCT(int inputSize, int dctSize, bool scale, DCTType dctType) :
            _scale( scale )
    {
        LOUDIA_DEBUG("DCT: Construction inputSize: " << inputSize
                    << ", dctSize: " << dctSize
                    << ", dctType: " << dctType);

        if (inputSize < dctSize) {
            // TODO: Throw an exception since dctSize is the number of coefficients to output and it cannot output more
            return;
        }

        setInputSize( inputSize, false );
        setDctSize( dctSize, false );
        setDctType( dctType, false );

        setup();
    }

    DCT::~DCT(){}

    void DCT::setup(){
        // Prepare the buffers
        LOUDIA_DEBUG("DCT: Setting up...");
        
        _dctMatrix.resize(_inputSize, _inputSize);

        switch(_dctType) {
        case I:
            type1Matrix( &_dctMatrix );
            break;

        case II:
            type2Matrix( &_dctMatrix );
            break;

        case III:
            // Throw ImplementationError not implemented yet
            break;

        case IV:
            // Throw ImplementationError not implemented yet
            break;

        case OCTAVE:
            typeOctaveMatrix( &_dctMatrix );
            break;

        }


        reset();
        LOUDIA_DEBUG("DCT: Finished setup.");
    }

    void DCT::type1Matrix(MatrixXR* dctMatrix) {
        int size = (*dctMatrix).rows();

        Real norm = 1.0;
        if ( _scale ) norm = sqrt(Real(2.0)/Real(size - 1));

        for(int i=0; i < size; i++){
            (*dctMatrix)(i, size - 1) = norm * 0.5 * pow((Real)-1, (Real)i);
            for(int j=1; j < size-1; j++){
                (*dctMatrix)(i,j) = norm * cos(Real(j * i) * M_PI / Real(size - 1));
            }
        }

        (*dctMatrix).col(0).setConstant(norm * 0.5);

    }

    void DCT::type2Matrix(MatrixXR* dctMatrix) {
        int size = (*dctMatrix).rows();

        // In MATLAB the first column is scaled differently
        Real norm0 = 1.0;
        Real norm = 1.0;
        if ( _scale ) {
            norm0 = sqrt(Real(1.0)/Real(size));
            norm = sqrt(Real(2.0)/Real(size));
        }

        for(int i=0; i < size; i++){
            for(int j=0; j < size; j++){
                if (j==0) {
                    (*dctMatrix)(i,j) = norm0 * cos(Real(j) * M_PI / Real(size) * (Real(i) + 0.5));
                } else {
                    (*dctMatrix)(i,j) = norm * cos(Real(j) * M_PI / Real(size) * (Real(i) + 0.5));
                }
            }
        }
    }

    void DCT::typeOctaveMatrix(MatrixXR* dctMatrix) {
        int size = (*dctMatrix).rows();

        Real norm = 1.0;
        if ( _scale ) norm = sqrt(2.0/Real(size));

        for(int i=0; i < size; i++){
            for(int j=1; j < size; j++){
                (*dctMatrix)(i,j) = norm * cos(Real(j) * M_PI / Real(2 * size) * (Real(2 * i - 1)));
            }
        }

        (*dctMatrix).col(0).setConstant( norm * sqrt(0.5) );
    }

    void DCT::process(const MatrixXR& input, MatrixXR* dctCoeffs){
        (*dctCoeffs).resize(input.rows(), _dctSize);

        for ( int i = 0 ; i < input.rows(); i++) {
            (*dctCoeffs).row(i) = input.row(i) * _dctMatrix.block(0, 0, input.cols(), _dctSize);
        }
    }

    void DCT::reset(){
        // Initial values
    }

    DCT::DCTType DCT::dctType() const{
        return _dctType;
    }

    void DCT::setDctType( DCTType type, bool callSetup ) {
        if (type == _dctType) {
        return;
        }
    
        _dctType = type;
        
        if ( callSetup ) setup();
    }

    int DCT::inputSize() const{
        return _inputSize;
    }

    void DCT::setInputSize( int size, bool callSetup ) {
        if (size == _inputSize) {
        return;
        }
        
        _inputSize = size;
        if ( callSetup ) setup();
    }

    int DCT::dctSize() const{
        return _dctSize;
    }

    void DCT::setDctSize( int size, bool /*callSetup*/ ) {
        _dctSize = size;
    }
}