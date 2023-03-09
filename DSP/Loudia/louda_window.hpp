#pragma once

namespace Loudia
{
        /**
    * @class Window
    *
    * @brief Algorithm to create and apply several type of windows on vectors of Real 
    * or Complex values.
    *
    * This class represents an object to apply a window on frames of Real or Complex data.
    * The algorithm takes as input N-point vectors of Real (or Complex) values 
    * and returns N-point vectors of Real (or Complex) values where the samples are weighted
    * by a weighting window.
    *
    * 5 types of windows are implemented:
    * -# Rectangular
    * -# Hann or Hanning
    * -# Hamming
    * -# Cosine
    * -# Blackmann
    * -# Blackmann Harris
    * -# Nuttall
    * -# Blackman Nuttall
    *
    * The Window type can be selected using the 
    * setWindowType() method.
    *
    * Additionally a Custom window can be specified using 
    * the setWindow() method.
    *
    * @author Ricard Marxer
    *
    * @sa FFT
    */
    class Window{
    public:
        /**
            @enum WindowType
            @brief Specifies the type of the window.
            @sa windowType
        */
        enum WindowType {
            RECTANGULAR         = 0 /**< Rectangular window */,
            HANN                = 1 /**< Hann window */,
            HANNING             = 2 /**< Alias for a Hann window */,
            HAMMING             = 3 /**< Hamming window */,
            COSINE              = 4 /**< Cosine window */,
            BLACKMAN            = 5 /**< Blackman window */,
            BLACKMANHARRIS      = 6 /**< Blackman-Harris window */,
            NUTTALL             = 7 /**< Nuttall window */,
            BLACKMANNUTTALL     = 8 /**< Blackman-Nuttall window */,
            CUSTOM              = 9 /**< Custom window. Note that this window type must be select
                                    when setting the window using setWindow()*/
        };

    protected:
        int _inputSize;
        WindowType _windowType;
        MatrixXR _window;
        
        MatrixXR hann(int length);
        MatrixXR hamming(int length);
        MatrixXR cosine(int length);

        MatrixXR blackmanType(int length, Real a0, Real a1, Real a2, Real a3);
        MatrixXR blackman(int length);
        MatrixXR nuttall(int length);
        MatrixXR blackmanHarris(int length);
        MatrixXR blackmanNuttall(int length);

        template<typename FrameMatrixType, typename WindowedMatrixType>
        void process(const FrameMatrixType& frames, WindowedMatrixType* windowedFrames);
    
    public: 
        /**
            Constructs a Window object with the given @a inputSize and @a windowType parameters
            given.
        */
        Window(int inputSize = 1024, WindowType windowType = RECTANGULAR);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~Window();

        void setup();
        void reset();

        /**
            Applies the window on each of the rows of @a frames and
            puts the result in the rows of @a windowedFrames.
            
            @param frames matrix of Real (or Complex) values.  The number of columns of @a frames must
            be equal to the inputSize property.
            
            @param windowedFrames pointer to a matrix of Real (or Complex) values for the output.
            The matrix should have the same number of rows and columns as @a frames. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */  
        void process(const MatrixXC& frames, MatrixXC* windowedFrames);
        void process(const MatrixXR& frames, MatrixXR* windowedFrames);
        void process(const MatrixXR& frames, MatrixXC* windowedFrames);

        /**
            Returns the input size of the algorithm.
            
            By default it is 1024.
        */
        int inputSize() const;
        
        /**
            Specifies the input size of the algorithm.
        */
        void setInputSize( int size, bool callSetup = true );
        
        /**
            Return the type of the window
            
            By default it is RECTANGULAR.
        */
        WindowType windowType() const;
        
        /**
            Specify the type of the window.
        */
        void setWindowType( WindowType type, bool callSetup = true );

        /**
            Return the single row matrix of Real values representing the window.
            
            The number of cols of the window will be equal to inputSize.
            
            By default it is a single row matrix with all values set to 1.0.
        */  
        const MatrixXR& window() const;

        /**
            Specify the single row matrix of Real values representing the window.
            
            The number of cols of the window must be equal to inputSize.
            
            Note that when the window is set, using setWindow(),
            the window type is automatically set to CUSTOM.
            By default it is a single row matrix with all values set to 1.0.
        */
        void setWindow( const MatrixXR& window, bool callSetup = true );
    };

    Window::Window(int inputSize, Window::WindowType windowType)
    {
    LOUDIA_DEBUG("WINDOW: Constructor inputSize: " << inputSize << 
            ", windowType: " << windowType);

    setInputSize( inputSize, false );
    setWindow( MatrixXR::Ones(1, inputSize), false );
    setWindowType( windowType, false );

    setup();

    LOUDIA_DEBUG("WINDOW: Constructed");
    }

    Window::~Window(){}

    void Window::setup(){
    LOUDIA_DEBUG("WINDOW: Setting up...");
    
    switch(_windowType){

    case RECTANGULAR:
        _window = MatrixXR::Ones(1, _inputSize);
        break;

    case HANN:
    case HANNING:
        _window = hann(_inputSize);
        break;

    case HAMMING:
        _window = hamming(_inputSize);
        break;

    case COSINE:
        _window = hamming(_inputSize);
        break;
        
    case BLACKMAN:
        _window = blackman(_inputSize);
        break;

    case NUTTALL:
        _window = nuttall(_inputSize);
        break;

    case BLACKMANHARRIS:
        _window = blackmanHarris(_inputSize);
        break;

    case BLACKMANNUTTALL:
        _window = blackmanNuttall(_inputSize);
        break;
    case CUSTOM:
        break;
        
    default:
        // Throw ValueError unknown window type
        break;
    }
    
    LOUDIA_DEBUG("WINDOW: Finished set up...");
    }

    MatrixXR Window::hann(int length){
    MatrixXR result(1, length);

    for(int i = 0; i < length; i++ ){
        result(0, i) = 0.5 * (1 - cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0)));
    }

    return result;
    }

    MatrixXR Window::hamming(int length){
    MatrixXR result(1, length);

    for(int i = 0; i < length; i++ ){
        result(0, i) = 0.53836 - 0.46164 * cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0));
    }

    return result;
    }

    MatrixXR Window::cosine(int length){
    MatrixXR result(1, length);

    for(int i = 0; i < length; i++ ){
        result(0, i) = sin((M_PI * (Real)i) / ((Real)length - 1.0));
    }

    return result;
    }

    MatrixXR Window::blackmanType(int length, Real a0, Real a1, Real a2, Real a3){
    MatrixXR result(1, length);

    Real pi_length_1 = M_PI / ((Real)length - 1.0);

    for(int i = 0; i < length; i++ ){
        result(0, i) = a0 \
                    - a1 * cos(2.0 * (Real)i * pi_length_1) \
                    + a2 * cos(4.0 * (Real)i * pi_length_1) \
                    - a3 * cos(6.0 * (Real)i * pi_length_1);
    }

    return result;
    }

    MatrixXR Window::blackman(int length){
    Real a0 = (1 - 0.16) / 2.0;
    Real a1 = 0.5;
    Real a2 = 0.16 / 2.0;
    Real a3 = 0.0;
    
    return blackmanType(length, a0, a1, a2, a3);
    }


    MatrixXR Window::nuttall(int length){
    Real a0 = 0.355768;
    Real a1 = 0.487396;
    Real a2 = 0.144232;
    Real a3 = 0.012604;
    
    return blackmanType(length, a0, a1, a2, a3);
    }


    MatrixXR Window::blackmanHarris(int length){
    Real a0 = 0.35875;
    Real a1 = 0.48829;
    Real a2 = 0.14128;
    Real a3 = 0.01168;
    
    return blackmanType(length, a0, a1, a2, a3);
    }


    MatrixXR Window::blackmanNuttall(int length){
    Real a0 = 0.3635819;
    Real a1 = 0.4891775;
    Real a2 = 0.1365995;
    Real a3 = 0.0106411;
    
    return blackmanType(length, a0, a1, a2, a3);
    }

    template<typename FrameMatrixType, typename WindowedMatrixType>
    void Window::process(const FrameMatrixType& frames, WindowedMatrixType* windowedFrames){
    (*windowedFrames).resize(frames.rows(), _inputSize);

    for (int i = 0; i < frames.rows(); i++){
        // Process and set
        (*windowedFrames).row(i) = (frames.row(i).array() * _window.array()).template cast<typename WindowedMatrixType::Scalar>();
    }

    }

    void Window::process(const MatrixXC& frames, MatrixXC* windowedFrames){
    process<MatrixXC, MatrixXC>(frames, windowedFrames);
    }

    void Window::process(const MatrixXR& frames, MatrixXC* windowedFrames){
    process<MatrixXR, MatrixXC>(frames, windowedFrames);
    }

    void Window::process(const MatrixXR& frames, MatrixXR* windowedFrames){
    process<MatrixXR, MatrixXR>(frames, windowedFrames);
    }


    void Window::reset(){
    }

    int Window::inputSize() const{
    return _inputSize;
    }

    void Window::setInputSize( int size, bool callSetup ) {
    _inputSize = size;
    if ( callSetup ) setup();
    }


    Window::WindowType Window::windowType() const{
    return _windowType;
    }

    void Window::setWindowType( WindowType type, bool callSetup ) {
    _windowType = type;
    if ( callSetup ) setup();
    }


    const MatrixXR& Window::window() const{
    return _window;
    }

    void Window::setWindow( const MatrixXR& window, bool callSetup ){
    if (window.cols() != _inputSize || window.rows() != 1) {
        // Throw exception wrong window size
    }

    setWindowType(CUSTOM, false);
    _window = window;

    if ( callSetup ) setup();
    }
}