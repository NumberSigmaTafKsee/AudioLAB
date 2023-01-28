#pragma once

namespace Loudia
{
    /**
    * @class LPC
    *
    * @brief Algorithm to calculate the Linear Predictive Coding of vectors of Real values.
    *
    * This class represents an object to perform a Linear Predictive Coding on vectors 
    * of Real values.  Which is a useful technique for estimating a parametric representation
    * of a spectrum magnitude.  The algorithm estimates a set of M coefficients of a IIR filter
    * whose frequency response approximates the vector of Reals passed as input.
    *
    * This algorithm implements the Levinson-Durbin recursion for solving the
    * following linear equation system:
    *
    * R a = r
    *
    * where R is the Toeplitz matrix made of the first M - 1 autocorrelation 
    * coefficients of the input vector and r is a vector made of M - 1 
    * autocorrelation coefficients starting from the second of the input
    * vector.
    *
    * Optionally a pre-emphasis FIR filter may be applied to the input vector
    * in order to enhance estimation of higher frequencies. The pre-emphasis filter
    * consists of a 2 coefficient filter of the form b = [1, -b1] where usually:
    *
    * 0.96 <= b1 <= 0.99
    *
    * The b1 coefficient defaults to 0, but can be specified using setPreEmphasis().
    *
    * @author Ricard Marxer
    *
    * @sa MelBands, Bands, MFCC
    */
    class LPC {
    protected:
        // Internal parameters
        int _inputSize;
        int _coefficientCount;
        Real _preEmphasis;
        
        // Internal variables
        MatrixXR _pre;
        MatrixXR _preRow;
        MatrixXR _temp;
        MatrixXR _acorr;

        Filter _preFilter;
        Autocorrelation _acorrelation;

    public:
        /**
            Constructs an LPC object with the specified @a inputSize, 
            @a coefficientCount and @a preEmphasis settings.
            
            @param inputSize size of the inputs arrays,
            must be > 0.
            The algorithm performs faster for sizes which are a power of 2.
            
            @param coefficientCount number of coefficients to be estimated
            @param preEmphasis second coefficient of the FIR pre-emphasis filter
        */
        LPC(int inputSize = 1024, int coefficientCount = 15, Real preEmphasis = 0.0);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~LPC();

        void setup();
        void reset();

        /**
            Performs an LPC on each of the rows of @a frames.
            Puts the resulting LPC coefficients in the rows of @a lpcCoefficients, the
            reflection coefficients in @a reflectionCoefficients and the error in @a error.
            
            @param frames matrix of Real values.  The number of columns of @a 
            frames must be equal to the input size specified using setInputSize().
            
            @param lpcCoefficients pointer to a matrix of Real values for the LPC coefficients.
            The matrix should have the same number of rows as @a frames and coefficientCount columns.
            @param reflectionCoefficients pointer to a matrix of 
            Real values for the reflection coefficients.
            The matrix should have the same number of rows as @a frames 
            and coefficientCount + 1 columns.
            
            @param error pointer to a matrix of 
            Real values for the LPC error gain.
            The matrix should have the same number of rows as @a frames 
            and 1 single column.
            Note that if the output matrices are not of the required sizes they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& frames, MatrixXR* lpcCoefficients,
                    MatrixXR* reflectionCoefficients, MatrixXR* error);

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
            Returns the number of coefficients to be calculated.
            The default is 15.
            
            @sa setCoefficientCount()
        */
        int coefficientCount() const;

        /**
            Specifies the @a count of coefficients to be calculated.
            The given @a count must be in the range between 0 and (input size - 1).
                
            @sa coefficientCount()
        */
        void setCoefficientCount( int count, bool callSetup = true );

        /**
            Returns the second coefficient of the FIR preemphasis filter.
            The default is 0.0.
            
            @sa setPreEmphasis()
        */
        Real preEmphasis() const;

        /**
            Specifies the second @a coefficient of the FIR preemphasis filter.
            
            @sa preEmphasis()
        */
        void setPreEmphasis( Real coefficient, bool callSetup = true );

    };

    class LPCResidual {
    protected:
        // Internal parameters
        int _frameSize;
        
        // Internal variables
        MatrixXR _result;

        Filter _filter;
    
    public:
        LPCResidual(int frameSize);

        ~LPCResidual();

        void setup();

        void process(const MatrixXR& frame, const MatrixXR& lpcCoeffs, MatrixXR* residual);

        void reset();

    };


    LPC::LPC(int inputSize, int coefficientCount, Real preEmphasis) 
    {
    LOUDIA_DEBUG("LPC: Constructor inputSize: " << inputSize 
            << ", coefficientCount: " << coefficientCount
            << ", preEmphasis: " << preEmphasis);

    if ( coefficientCount > inputSize ) {
        // Thorw ValueError, the number of coefficients must be smaller or equal than the frame size.
    }
    
    setInputSize( inputSize, false ); 
    setCoefficientCount( coefficientCount, false );
    setPreEmphasis( preEmphasis, false );
    
    setup();
    }

    LPC::~LPC() {}


    void LPC::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("LPC: Setting up...");
    
    if ( _preEmphasis != 0.0 ) {
        MatrixXR preCoeffs(2, 1);
        preCoeffs << 1, -_preEmphasis;
        _preFilter.setB( preCoeffs );
    }
    
    _acorrelation.setInputSize( _inputSize, false );
    _acorrelation.setMaxLag( _coefficientCount + 1, false );
    _acorrelation.setUseFft( true, false );
    _acorrelation.setup();
    
    reset();
    
    LOUDIA_DEBUG("LPC: Finished set up...");
    }


    void LPC::process(const MatrixXR& frame, MatrixXR* lpcCoeffs, MatrixXR* reflectionCoeffs, MatrixXR* error){
    LOUDIA_DEBUG("LPC: Processing...");
    const int rows = frame.rows();
    const int cols = frame.cols();

    if ( cols != _inputSize ) {
        // Throw ValueError, the frames passed are the wrong size
    }
    
    _pre.resize(rows, cols);
    
    if ( _preEmphasis != 0.0 ) {
        for ( int row = 0; row < rows; row++) {
        _preFilter.process( frame.transpose(), &_preRow );
        _pre.row( row ) = _preRow.transpose();
        }
    } else {
        _pre = frame;
    }
    
    LOUDIA_DEBUG("LPC: Processing autocorrelation");
    
    _acorrelation.process(_pre, &_acorr);
    
    LOUDIA_DEBUG("LPC: Processing Levinson-Durbin recursion");

    (*lpcCoeffs).resize(rows, _coefficientCount);
    (*reflectionCoeffs).resize(rows, _coefficientCount - 1);
    (*error).resize(rows, 1);
    
    // Initial values of the LPC coefficients
    (*lpcCoeffs).setZero();
    (*lpcCoeffs).col(0).setOnes();
    
    // Initial value of the Error
    (*error).col(0) = _acorr.col(0);

    (*reflectionCoeffs).setZero();

    for ( int row = 0; row < rows; row++) {  
        Real gamma;
        
        if ((_acorr.array() == 0.).all())
        continue;

        for ( int i = 1; i < _coefficientCount; i++ ) {
        gamma = _acorr(row, i);

        // Use the Eigen reverse()      
    //       if ( i >= 2) {
    //         gamma += ((*lpcCoeffs).row(row).segment(1, i-1) * _acorr.row(row).segment(1, i-1).transpose().reverse())(0,0);
    //       }
        
        // instead of manually walking it in reverse order
        for (int j = 1; j <= i-1; ++j) {
        gamma += (*lpcCoeffs)(row, j) * _acorr(row, i-j);  
        }
        
        // Get the reflection coefficient
        (*reflectionCoeffs)(row, i-1) = - gamma / (*error)(row, 0);

        // Update the error      
        (*error)(row, 0) *= (1 - (*reflectionCoeffs)(row, i-1) * (*reflectionCoeffs).conjugate()(row, i-1));
        
        // Update the LPC coefficients
        if(i >= 2){
            _temp = (*lpcCoeffs).row(row).segment(1, i-1);
            reverseCols(&_temp);
            
            (*lpcCoeffs).row(row).segment(1, i-1) += (*reflectionCoeffs)(row, i-1) * _temp.conjugate();
        }
        
        (*lpcCoeffs)(row, i) = (*reflectionCoeffs)(row, i-1);
        }
    }
    
    LOUDIA_DEBUG("LPC: Finished Processing");
    }

    void LPC::reset(){
    // Initial values

    if ( _preEmphasis != 0.0 ) {
        _preFilter.reset( );
    }

    _acorrelation.reset( );


    }

    int LPC::inputSize() const {
    return _inputSize;
    }
    
    void LPC::setInputSize( int size, bool callSetup ) {
    _inputSize = size;
    if ( callSetup ) setup();
    }


    int LPC::coefficientCount() const {
    return _coefficientCount;
    }

    void LPC::setCoefficientCount( int count, bool callSetup ) {
    _coefficientCount = count;
    if ( callSetup ) setup();
    }

    Real LPC::preEmphasis() const {
    return _preEmphasis;
    }

    void LPC::setPreEmphasis( Real coefficient, bool callSetup ) {
    _preEmphasis = coefficient;
    if ( callSetup ) setup();
    }

    LPCResidual::LPCResidual(int frameSize) : 
    _frameSize( frameSize ),
    _filter( 1 )
    {
    LOUDIA_DEBUG("LPCRESIDUAL: Constructor frameSize: " << frameSize);
    
    setup();
    }

    LPCResidual::~LPCResidual() {}


    void LPCResidual::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("LPCRESIDUAL: Setting up...");

    reset();
    
    LOUDIA_DEBUG("LPCRESIDUAL: Finished set up...");
    }


    void LPCResidual::process(const MatrixXR& frame, const MatrixXR& lpcCoeffs, MatrixXR* residual){
    LOUDIA_DEBUG("LPCRESIDUAL: Processing...");
    const int rows = frame.rows();
    const int cols = frame.cols();

    if ( cols != _frameSize ) {
        // Throw ValueError, the frames passed are the wrong size
    }

    (*residual).resize(rows, cols);

    LOUDIA_DEBUG("LPCRESIDUAL: Processing setting filter");
    _filter.setA( lpcCoeffs.transpose() );
    
    for ( int row = 0; row < rows; row++ ) {
        
        LOUDIA_DEBUG("LPCRESIDUAL: Processing filter");
        _filter.process(frame.row( row ).transpose(), &_result);
        
        (*residual).row( row ) = _result.transpose();
    }

    LOUDIA_DEBUG("LPCRESIDUAL: Finished Processing");
    }

    void LPCResidual::reset(){
    // Initial values
    }
}