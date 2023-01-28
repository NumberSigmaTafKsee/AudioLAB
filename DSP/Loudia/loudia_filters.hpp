#pragma once

namespace Loudia
{

/**
    * @class Filter
    *
    * @brief Algorithm to apply one or several IIR filters given the Real value coefficients.
    *
    * This class represents an object to apply one or several IIR filters.
    * The coefficients must be manually set by the user.  To create and use
    * some special parametrized filters such as Low Pass, High Pass, Band Pass and Band Stop
    * filters refer to IIRFilter.
    *
    * For Complex coefficient filters, refer to FilterComplex.
    *
    * This filter implementation allows single and multiple channel filtering.
    * The number of channels is specified using the setChannelCount().
    * The a and b coefficients of the filter are specified by two matrices of
    * Real values. The rows of the matrix are the time indices of the filter
    * and the columns (if more than one) are the channels.
    *
    * Three different situations are possible, depending on the number of columns
    * in the coefficients matrices and in the input matrix:
    * - if the number of columns of the coefficients matrices 
    * are equal to the number columns in the input matrix, then
    * each column of the input matrix is filtered by a column 
    * of the coefficients matrix. This is the situation when trying to
    * filter differently all the channels of a multi-channel signal.
    * - if the coefficients matrix has one single column 
    * and the input matrix has multiple columns then each 
    * column of the input matrix is filtered by the single column of the
    * coefficients matrices.  This is the situation when trying to
    * filter equally all the channels of a multi-channel signal.
    * - if the coefficients matrices have multiple columns each and the 
    * input matrix has multiple columns then the column of the input matrix is filtered 
    * by each of the columns of the coefficients matrices.  This is the situation
    * when applying a filterbank to a single channel signal.
    *
    * Note that in all cases the number of columns in a and b coefficients matrices
    * must be the same.
    * Note that the channel count determines the number of output channels in any situation
    * and is therefore must be equal to the maximum number of channels between input and 
    * coefficient matrices.
    *
    * @author Ricard Marxer
    *
    * @sa Filter
    */
    class Filter {
    protected:
        // Internal parameters
        int _channelCount; 
        int _length;

        // Internal variables
        MatrixXR _ina;
        MatrixXR _inb;

        MatrixXR _a;
        MatrixXR _b;

        MatrixXR _z;
        MatrixXR _samples;

        void setupState();
        void setupCoeffs();

    public:
        /**
            Constructs a band pass filter object with the given @a channelCount, @a b,
            and @a a coefficients given.
        */
        Filter(int channelCount = 1);
        Filter(const MatrixXR& b, const MatrixXR& a, int channelCount);

        void setup();
        void reset();

        /**
            Performs a filtering of each of the columns of @a samples.
            Puts the resulting filtered in the columns of @a filtered.
            
            @param samples matrix of Real values.  A column represents a channel and a row 
            represents a time index.
            
            @param filtered pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows and columns as @a samples.
            
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
            
            @sa IIRFilter::process
        */
        void process(const MatrixXR& samples, MatrixXR* filtered);

        /**
            Returns the a coefficients of the filter
            
            Both the b and a coefficients matrices are normalized
            by the first row of the a coefficients matrix.
            
            Note that if the first row of the a coefficients matrix
            has elements to zero, some of the filtered samples will 
            result in NaN.
            
            Note that the number of columns in a and b must be the same,
            and that it must be equal to 1 or Filter::channelCount.
            By default it is a single element matrix of value 1.0.
            @sa setA(), b(), setB()
        */
        void a( MatrixXR* a ) const;

        /**
            Specifies the a coefficients of the filter
            
            Both the b and a coefficients matrices are normalized
            by the first row of the a coefficients matrix.
            
            Note that if the first row of the a coefficients matrix
            has elements to zero, some of the filtered samples will 
            result in NaN.
            
            Note that the number of columns in a and b must be the same,
            and that it must be equal to 1 or Filter::channelCount.
            
            @sa a(), b(), setB()
        */
        void setA( const MatrixXR& a, bool callSetup = true );

        /**
            Returns the b coefficients of the filter
            
            Both the b and a coefficients matrices are normalized
            by the first row of the a coefficients matrix.
            
            Note that if the first row of the a coefficients matrix
            has elements to zero, some of the filtered samples will 
            result in NaN.
            Note that the number of columns in a and b must be the same,
            and that it must be equal to 1 or the channel count.
            
            By default it is a single element matrix of value 1.0.
            @sa setB(), a(), setA()
        */
        void b( MatrixXR* b ) const;

        /**
            Specifies the b coefficients of the filter
            
            Both the b and a coefficients matrices are normalized
            by the first row of the a coefficients matrix.
            
            Note that if the first row of the a coefficients matrix
            has elements to zero, some of the filtered samples will 
            result in NaN.
            Note that the number of columns in a and b must be the same,
            and that it must be equal to 1 or the channel count.
            
            @sa b(), a(), setA()
        */
        void setB( const MatrixXR& b, bool callSetup = true );

        /**
            Returns the number of output channles of the filter
            
            Note that the number of channels must be equal to the 
            number of columns in a and b or to the number of columns
            in the input matrix.
            
            By default it is 1.
        */
        int channelCount() const;

        /**
            Specifies the number of output channles of the filter
            
            Note that the number of channels must be equal to the 
            number of columns in a and b or to the number of columns
            in the input matrix.
        */
        void setChannelCount( int count, bool callSetup = true );


        /**
            Returns the length of the filter, which is 
            the maximum number of rows between the a and b coefficients 
            matrices.
        */
        int length() const;
    };

    /**
    * @class BandFilter
    *
    * @brief Algorithm to create and apply several types of low pass, high pass, 
    * band pass and band stop filters.
    *
    * This class represents an object to create and apply several types of band filters.
    * Additionally the coefficients, zeros, poles and gains of the created filters
    * can be retrieved.
    *
    * 4 types of bands are implemented:
    * -# Low Pass
    * -# High Pass
    * -# Band Pass
    * -# Band Stop
    *
    * The band type can be selected using the 
    * setBandType() taking as argument a BandType.
    *
    * The critical frequencies are specified using 
    * setLowFrequency() and setHighFrequency().
    * Note that for low pass and high pass filters which have one single critical frequency
    * only setLowFrequency() has an effect.
    *
    * 4 types of filters are implemented:
    * -# Chebyshev I
    * -# Chebyshev II
    * -# Bessel
    * -# Butterworth
    *
    * The filter type can be selected using the 
    * setFilterType() taking as argument a FilterType.
    *
    * The order of the filters can be specified using setOrder().
    *
    * For Chebyshev I filters the pass band ripple can be specified using
    * setPassRipple().  Note that this method has no effect if
    * a different type of filter is used.
    * 
    * For Chebyshev II filters the stop band attenuation is specified using
    * setStopAttenuation().  Note that this method has no effect if
    * a different type of filter is used.
    *
    * @author Ricard Marxer
    *
    * @sa Filter
    */
    class BandFilter {
    public:
        /**
            @enum FilterType
            @brief Specifies the type of the filter.
            @sa filterType
        */
        enum FilterType {
            CHEBYSHEVI = 0 /**< Chebyshev Type-I filter */,
            CHEBYSHEVII = 1 /**< Chebyshev Type-II filter */,
            BUTTERWORTH = 2 /**< Butterworth filter */,
            BESSEL = 3 /**< Bessel filter */
        };

        /**
            @enum BandType
            @brief Specifies the type of the band.
            
            @sa bandType
        */
        enum BandType {
            LOWPASS = 0 /**< Low pass filter */,
            HIGHPASS = 1 /**< High pass filter */,
            BANDPASS = 2 /**< Band pass filter */,
            BANDSTOP = 3 /**< Band stop filter */
        };

    protected:
        int _order;
        Real _lowFrequency;
        Real _highFrequency;
        Real _passRipple;
        Real _stopAttenuation;
        int _channelCount;
        
        Filter _filter;

        FilterType _filterType;
        BandType _bandType;
    
    public:
        /**
            Constructs a band pass filter object with the given @a order, @a lowFrequency,
            @a highFrequency, @a filterType, @a ripplePass and @a attenuationStop parameters
            given.
        */
        BandFilter(int order = 4, Real lowFrequency = 0.0, Real highFrequency = 1.0, BandType bandType = LOWPASS, FilterType filterType = CHEBYSHEVII, Real ripplePass = 0.05, Real attenuationStop = 40.0);

        void setup();
        void reset();

        /**
            Performs a filtering of each of the columns of @a samples.
            Puts the resulting filtered in the columns of @a filtered.
            
            @param samples matrix of Real values.  A column represents a channel and a row 
            represents a time index.
            
            @param filtered pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows and columns as @a samples.
            
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
            
            @sa Filter::process
        */
        void process( const MatrixXR& samples, MatrixXR* filtered );

        /**
            Return in @a a the single column matrix @a a coefficients.
            @sa b
        */
        void a( MatrixXR* a ) const;
        
        /**
            Return in @a b the single column matrix @a b coefficients.
            @sa a
        */
        void b( MatrixXR* b ) const;

        /**
            Return the order of the filter.
            The default is 4.
            @sa setOrder
        */
        int order() const;

        /**
            Specifies the @a order of the filter.
            The given @a order must be higher than 0.
            Note that orders higher than 25 are not allowed for Bessel filters.
            
            @sa order
        */
        void setOrder( int order, bool callSetup = true );

        /**
            Return the low frequency of the filter.
            The default is 0.0.
            @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
        */  
        Real lowFrequency() const;  

        /**
            Specifies the low normalized @a frequency of the filter.
            The given @a frequency must be in the range of 0 to 1.
            
            @sa lowFrequency, highFrequency, setHighFrequency
        */
        void setLowFrequency( Real frequency, bool callSetup = true );

        /**
            Return the stop frequency of the filter.
            The default is 1.0.
            @sa lowFrequency, setLowFrequency, setHighFrequency
        */  
        Real highFrequency() const;  

        /**
            Specifies the stop normalized @a frequency of the filter.
            The given @a frequency must be in the range of 0 to 1.
            
            @sa lowFrequency, highFrequency, setLowFrequency
        */
        void setHighFrequency( Real frequency, bool callSetup = true );

        /**
            Return the filter type.
            The default is CHEBYSHEVII.
            The given @a frequency must be in the range of 0 to 1.
            @sa setFilterType, order, setOrder
        */
        FilterType filterType() const;

        /**
            Specifies the filter @a type.
            
            @sa lowFrequency, highFrequency, setLowFrequency
        */
        void setFilterType( FilterType type, bool callSetup = true );

        /**
            Return the type of the band of the filter.
            
            By default it is LOWPASS.
            @sa setBandType()
        */
        BandType bandType() const;

        /**
            Specifies the type of the band of the filter.
            
            @sa bandType()
        */
        void setBandType( BandType type, bool callSetup = true );

        /**
            Returns the ripple of the pass band in dB
            
            Note that this property only has an effect if 
            the filter type used is CHEBYSHEVI.
            By default it is 0.05.
            @sa setPassRipple(), stopAttenuation(), setStopAttenuation()
        */
        Real passRipple() const;

        /**
            Specifies the ripple of the pass band in dB
            
            Note that this property only has an effect if 
            the filter type used is CHEBYSHEVI.
            
            @sa passRipple(), stopAttenuation(), setStopAttenuation()
        */
        void setPassRipple( Real rippleDB, bool callSetup = true );

        /**
            Returns the attenuation of the stop band in dB
            
            Note that this property only has an effect if 
            the filter type used is CHEBYSHEVII.
            By default it is 40.0.
            @sa passRipple(), setPassRipple(), setStopAttenuation()
        */
        Real stopAttenuation() const;

        /**
            Specifies the attenuation of the stop band in dB
            
            Note that this property only has an effect if 
            the filter type used is CHEBYSHEVII.
            
            @sa passRipple(), setPassRipple(), stopAttenuation()
        */
        void setStopAttenuation( Real attenuationDB, bool callSetup = true );
    };

    Filter::Filter(int channelCount)
    {
    LOUDIA_DEBUG("FILTER: Constructor channelCount:" << channelCount);

    setChannelCount( channelCount, false );

    setA(MatrixXR::Ones(1, _channelCount));
    setB(MatrixXR::Ones(1, _channelCount));

    setup();
    }


    Filter::Filter(const MatrixXR& b,
                const MatrixXR& a,
                int channelCount)
    {

    LOUDIA_DEBUG("FILTER: Constructor channelCount:" << channelCount);
    LOUDIA_DEBUG("FILTER: Constructor b:" << b.transpose() << ", a:" << a.transpose());

    setChannelCount( channelCount, false );

    setA( a );
    setB( b );
        
    setup();
    }

    void Filter::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("FILTER: Setting up...");

    setupState();

    _samples.resize(1, _channelCount);
    
    reset();

    LOUDIA_DEBUG("FILTER: Finished set up...");
    }


    void Filter::setupState() { 
    if( _z.rows() != _length ){
        _z.resize(_length, _channelCount);

        reset(); // if the state has changed, we reset it
    }

    }


    void Filter::process(const MatrixXR& samples, MatrixXR* output){
    // Process will be called with a matrix where columns will be channelCount
    // and rows will be the time axis

    //DEBUG("FILTER: Entering process...");
    _samples.resize(samples.rows(), _channelCount);
    (*output).resize(samples.rows(), _channelCount);
    //DEBUG("FILTER: After resize...");  
    

    // Check that it has one column or as many as channelCount
    if ((samples.cols() != 1) && (samples.cols() != _channelCount)) {
        // TODO: Throw an exception
        LOUDIA_DEBUG("FILTER: Error in shape of 'samples'. samples.cols():" << samples.cols() << ", _channelCount:" << _channelCount);
        return;
    }

    // Set the input
    if (samples.cols() == 1) {
        // If only one column has been defined repeat it for all columns
        for (int i=0; i < _samples.cols(); i++) {
        _samples.block(0, i, samples.rows(), 1) = samples.col(0);
        }
    }else{
        // Else set it directly
        _samples.block(0, 0, samples.rows(), _channelCount) = samples;
    }
    //DEBUG("FILTER: _a: " << _a);

    //DEBUG("FILTER: After setting coeffs...");    
    for ( int i = 0; i < _samples.rows(); i++ ) {
        if ( _length > 1 ) {
        //DEBUG("output.rows(): " << (*output).rows());
        //DEBUG("output.cols(): " << (*output).cols());
        //DEBUG("_z.rows(): " << _z.rows());
        //DEBUG("_z.cols(): " << _z.cols());
        //DEBUG("_b.rows(): " << _b.rows());
        //DEBUG("_b.cols(): " << _b.cols());
        //DEBUG("_a.rows(): " << _a.rows());
        //DEBUG("_a.cols(): " << _a.cols());
        (*output).row( i ) = _z.row( 0 ).array() + (_b.row( 0 ).array() * _samples.row( i ).array());
        
        //DEBUG("FILTER: After setting output..., output: " << (*output).row( i ));
        // Fill in middle delays
        for ( int j = 0; j < _length - 1; j++ ) {      
            _z.row( j ) = _z.row( j + 1 ).array() + (_samples.row( i ).array() * _b.row( j + 1 ).array()) - ((*output).row( i ).array() * _a.row( j + 1 ).array());
        }

        // Calculate the last delay
        _z.row( _length - 2 ) = (_samples.row( i ).array() * _b.row( _length - 1 ).array()) - ((*output).row( i ).array() * _a.row( _length - 1 ).array());

        } else {
        (*output).row( i ) = _samples.row( i ) * _b.row( 0 );
        }
    }
    //DEBUG("FILTER: output: " << (*output));
    //DEBUG("FILTER: After processing...");
    }


    void Filter::setA( const MatrixXR& a, bool callSetup ){
    _ina = a;

    if ( callSetup ) {
        setupCoeffs();
        setupState();
    }
    }

    void Filter::setB( const MatrixXR& b, bool callSetup ){
    _inb = b;

    if ( callSetup ) {
        setupCoeffs();
        setupState();
    }
    }

    void Filter::setupCoeffs() {
    _length = max(_inb.rows(), _ina.rows());

    // Normalize by the first value value the denominator coefficients
    // since a[0] must be 1.0
    // TODO: throw an exception when a[0] == 0
    LOUDIA_DEBUG("FILTER: Initializing 'a' coeffs");
    _a = MatrixXR::Zero(_length, _channelCount);

    // Check that it has one column or as many as channelCount
    if ((_ina.cols() != 1) && (_ina.cols() != _channelCount)) {
        // TODO: Throw an exception
        LOUDIA_DEBUG("FILTER: Error in shape of 'a' coeffs. _ina.cols():" << _ina.cols() << ", _channelCount:" << _channelCount);
        return;
    }

    // Set the coefficients
    if (_ina.cols() == 1) {
        // If only one column has been defined repeat it for all columns
        for (int i=0; i < _a.cols(); i++) {
        _a.block(0, i, _ina.rows(), 1) = _ina.col(0);
        }
    }else{
        // Else set it directly
        _a.block(0, 0, _ina.rows(), _channelCount) = _ina;
    }

    for(int i = 0; i < _a.rows(); i++){
        _a.row(i) = _a.row(i).array() / _ina.row(0).array();
    }

    LOUDIA_DEBUG("FILTER: Setting the 'a' coefficients.");
    LOUDIA_DEBUG("FILTER: 'a' coefficients: " << _a.transpose());

    _b = MatrixXR::Zero(_length, _channelCount);

    // Check that it has one column or as many as channelCount
    if ((_inb.cols() != 1) && (_inb.cols() != _channelCount)) {
        // TODO: Throw an exception
        LOUDIA_DEBUG("FILTER: Error in shape of 'b' coeffs. b.cols():" << _inb.cols() << ", _channelCount:" << _channelCount);
        return;
    }

    // Set the coefficients
    if (_inb.cols() == 1) {
        // If only one column has been defined repeat it for all columns
        for (int i=0; i < _b.cols(); i++) {
        _b.block(0, i, _inb.rows(), 1) = _inb.col(0);
        }
    }else{
        // Else set it directly
        _b.block(0, 0, _inb.rows(), _channelCount) = _inb;
    }

    for(int i = 0; i < _b.rows(); i++){
        _b.row(i) = _b.row(i).array() / _ina.row(0).array();
    }  
    
    LOUDIA_DEBUG("FILTER: Setting the 'b' coefficients.");
    LOUDIA_DEBUG("FILTER: 'b' coefficients: " << _b.transpose());  

    }

    void Filter::a(MatrixXR* a) const {
    (*a) = _a;
    }

    void Filter::b(MatrixXR* b) const {
    (*b) = _b;
    }


    void Filter::reset(){
    // Initial values
    _z = MatrixXR::Zero(_length, _channelCount);
    }

    int Filter::channelCount() const {
    return _channelCount;
    }

    void Filter::setChannelCount( int count, bool callSetup ) {
    _channelCount = count;
    
    if ( callSetup ) setup();  
    }

    int Filter::length() const {
    return _length;
    }

    BandFilter::BandFilter( int order, Real lowFrequency, Real highFrequency, BandType bandType, FilterType filterType, Real passRipple, Real stopAttenuation ) : 
    _channelCount( 1 ),
    _filter( _channelCount )
    {
    LOUDIA_DEBUG("BANDFILTER: Constructor order: " << order 
            << ", lowFrequency: " << lowFrequency
            << ", highFrequency: " << highFrequency
            << ", passRipple: " << passRipple
            << ", stopAttenuation: " << stopAttenuation );

    if ( order < 1 ) {
        // Throw an exception
    }

    setOrder( order, false );
    setLowFrequency( lowFrequency, false );
    setHighFrequency( highFrequency, false );
    setPassRipple( passRipple, false );
    setStopAttenuation( stopAttenuation, false );
    setFilterType( filterType, false ); 
    setBandType( bandType, false );
    
    setup();
    
    LOUDIA_DEBUG("BANDFILTER: Constructed");
    }

    void BandFilter::setup(){
    LOUDIA_DEBUG("BANDFILTER: Setting up...");

    _filter.setChannelCount( _channelCount, false );

    LOUDIA_DEBUG("BANDFILTER: Getting zpk");  
    // Get the lowpass z, p, k
    MatrixXC zeros, poles;
    Real gain;

    switch( _filterType ){
    case CHEBYSHEVI:
        chebyshev1(_order, _passRipple, _channelCount, &zeros, &poles, &gain);
        break;

    case CHEBYSHEVII:
        chebyshev2(_order, _stopAttenuation, _channelCount, &zeros, &poles, &gain);
        break;

    case BUTTERWORTH:
        butterworth(_order, _channelCount, &zeros, &poles, &gain);
        break;

    case BESSEL:
        bessel(_order, _channelCount, &zeros, &poles, &gain);
        break;
    }
    
    LOUDIA_DEBUG("BANDFILTER: zeros:" << zeros );
    LOUDIA_DEBUG("BANDFILTER: poles:" << poles );
    LOUDIA_DEBUG("BANDFILTER: gain:" << gain );
    
    // Convert zpk to ab coeffs
    MatrixXC a;
    MatrixXC b;
    zpkToCoeffs(zeros, poles, gain, &b, &a);

    LOUDIA_DEBUG("BANDFILTER: Calculated the coeffs");

    // Since we cannot create matrices of Nx0
    // we have created at least one Zero in 0
    if ( zeros == MatrixXC::Zero(zeros.rows(), zeros.cols()) ){
        // Now we must remove the last coefficient from b
        MatrixXC temp = b.block(0, 0, b.rows(), b.cols()-1);
        b = temp;
    }

    // Get the warped critical frequency
    Real fs = 2.0;
    Real warped = 2.0 * fs * tan( M_PI * _lowFrequency / fs );
    
    Real warpedStop = 2.0 * fs * tan( M_PI * _highFrequency / fs );
    Real warpedCenter = sqrt(warped * warpedStop);
    Real warpedBandwidth = warpedStop - warped;

    // Warpped coeffs
    MatrixXC wa;
    MatrixXC wb;

    LOUDIA_DEBUG("BANDFILTER: Create the band type filter from the analog prototype");

    switch( _bandType ){
    case LOWPASS:
        lowPassToLowPass(b, a, warped, &wb, &wa);
        break;
        
    case HIGHPASS:
        lowPassToHighPass(b, a, warped, &wb, &wa);  
        break;

    case BANDPASS:
        lowPassToBandPass(b, a, warpedCenter, warpedBandwidth, &wb, &wa);
        break;
        
    case BANDSTOP:
        lowPassToBandStop(b, a, warpedCenter, warpedBandwidth, &wb, &wa);
        break;
    }

    LOUDIA_DEBUG("BANDFILTER: Calculated the low pass to band pass");
    
    // Digital coeffs
    MatrixXR da;
    MatrixXR db;
    bilinear(wb, wa, fs, &db, &da);
    
    LOUDIA_DEBUG("BANDFILTER: setup the coeffs");

    // Set the coefficients to the filter
    _filter.setA( da.transpose() );
    _filter.setB( db.transpose() );
    
    _filter.setup();
    
    LOUDIA_DEBUG("BANDFILTER: Finished set up...");
    }

    void BandFilter::a(MatrixXR* a) const{
    _filter.a(a);
    }

    void BandFilter::b(MatrixXR* b) const{
    _filter.b(b);
    }

    void BandFilter::process(const MatrixXR& samples, MatrixXR* filtered) {
    _filter.process(samples, filtered);
    }

    void BandFilter::reset(){
    // Initial values
    _filter.reset();
    }

    int BandFilter::order() const{
    return _order;
    }

    void BandFilter::setOrder( int order, bool callSetup ){
    _order = order;
    if ( callSetup ) setup();
    }

    Real BandFilter::lowFrequency() const{
    return _lowFrequency;
    }
    
    void BandFilter::setLowFrequency( Real frequency, bool callSetup ){
    _lowFrequency = frequency;
    if ( callSetup ) setup();
    }

    Real BandFilter::highFrequency() const{
    return _highFrequency;
    }
    
    void BandFilter::setHighFrequency( Real frequency, bool callSetup ){
    _highFrequency = frequency;
    if ( callSetup ) setup();
    }

    BandFilter::FilterType BandFilter::filterType() const{
    return _filterType;
    }

    void BandFilter::setFilterType( FilterType type, bool callSetup ){
    _filterType = type;
    if ( callSetup ) setup();
    }

    BandFilter::BandType BandFilter::bandType() const{
    return _bandType;
    }

    void BandFilter::setBandType( BandType type, bool callSetup ){
    _bandType = type;
    if ( callSetup ) setup();
    }

    Real BandFilter::passRipple() const{
    return _passRipple;
    }

    void BandFilter::setPassRipple( Real rippleDB, bool callSetup ){
    _passRipple = rippleDB;
    if ( callSetup ) setup();
    }

    Real BandFilter::stopAttenuation() const{
    return _stopAttenuation;
    }

    void BandFilter::setStopAttenuation( Real attenuationDB, bool callSetup ){
    _stopAttenuation = attenuationDB;
    if ( callSetup ) setup();
    }
}
