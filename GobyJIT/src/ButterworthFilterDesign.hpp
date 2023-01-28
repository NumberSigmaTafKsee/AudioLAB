// Butterworth-Filter-Design
// https://github.com/ruohoruotsi/Butterworth-Filter-Design
// It has some interesting code but the higher order eq should be reduced

#include <vector>
#include <complex>


typedef std::complex <double> complex_double;


namespace FX::Filters::ButterworthDesign
{

    // A biquad filter expression:
    // y[n] = b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2];

    class Biquad {
        
    public:
        Biquad();
        ~Biquad();
        
        
        // Coefficients are public, as client classes and test harnesses need direct access.
        // Special accessors to better encapsulate data, would be less readable,
        // so skipping for pedagogical reasons. Another idea would be to make this a struct.
        //
        double b0, b1, b2, a1, a2;  // second order section variable
        double a3, a4, b3, b4;      // fourth order section variables
        
        
        // Coefficients for a DF2T fourth order section (Used for EQ filters)
        void DF2TFourthOrderSection(double B0, double B1, double B2, double B3, double B4,
                                    double A0, double A1, double A2, double A3, double A4);
        
        // Coefficients for a DF2T biquad section.
        void DF2TBiquad(double B0, double B1, double B2,
                        double A0, double A1, double A2);
    };



    class BiquadChain {
        
    public:
        
        BiquadChain(int count);
        
        BiquadChain();
        ~BiquadChain();
        
        void resize(int count);
        void reset();
        
        // Process the biquad filter chain on the input buffer, write to output buffer
        // buffers arrays contain count elements with a single stride separating successive elements.
        void processBiquad(const float * input, float * output, const int stride, const int count, const Biquad * coeffs);
        
        // Extension for fourth order sections
        void processFourthOrderSections(const float * input, float * output, int stride, int count, const Biquad * coeffs);
        
    private:
        
        int numFilters;
        double _xn1, _xn2;
        std::vector <double> _yn, _yn1, _yn2;
        
        // Fourth order section variables
        double xn3, xn4;
        std::vector <double> _yn3, _yn4;
        
        void allocate(int count);
    };


    class Butterworth {
        
    public:
        
        Butterworth(){
        }
        ~Butterworth(){
        }
        
        static vector <complex_double> prototypeAnalogLowPass(int filterOrder);
        
        bool loPass(double fs, double f1, double f2, int filterOrder,
                    vector <Biquad> & coeffs, double & overallGain){
            return coefficients(kLoPass, fs, f1, f2, filterOrder, coeffs, overallGain);
        }
        
        bool hiPass(double fs, double f1, double f2, int filterOrder,
                    vector <Biquad> & coeffs, double & overallGain){
            return coefficients(kHiPass, fs, f1, f2, filterOrder, coeffs, overallGain);
        }
        
        bool bandPass(double fs, double f1, double f2, int filterOrder,
                    vector <Biquad> & coeffs, double & overallGain){
            return coefficients(kBandPass, fs, f1, f2, filterOrder, coeffs, overallGain);
        }
        
        bool bandStop(double fs, double f1, double f2, int filterOrder,
                    vector <Biquad> & coeffs, double & overallGain){
            return coefficients(kBandStop, fs, f1, f2, filterOrder, coeffs, overallGain);
        }
        
        
        enum FILTER_TYPE {
            kLoPass     = 10000,
            kHiPass     = 10001,
            kBandPass   = 10002,
            kBandStop   = 10003,
            kLoShelf    = 10004,
            kHiShelf    = 10005, // high order EQ
            kParametric = 10006  // high order EQ
        };
        
        
        //******************************************************************************
        
        // generic coeffs
        bool coefficients(FILTER_TYPE filter, const double fs, const double freq1_cutoff, const double freq2_cutoff, const int filterOrder,
                        vector <Biquad> & coeffs, double & overallGain);
        
        // high order EQ
        bool coefficientsEQ(FILTER_TYPE filter, const double fs, const double f1, const double f2, const int filterOrder,
                            vector <Biquad> & coeffs, double overallGain);
        
        
    private:
        
        double blt(complex_double & sz);
        bool s2Z();
        bool zp2SOS();
        
        // Private helper methods to convert the lowpass analog prototype to the desired filter type
        void convert2lopass();
        void convert2hipass();
        void convert2bandpass();
        void convert2bandstop();
        
        // Internal state used during computation of coefficients
        double f1, f2;
        
        int numPoles, numZeros;
        std::vector <complex_double> zeros;
        std::vector <complex_double> poles;
        
        double Wc; // Omega cutoff == passband edge freq
        double bw; // Bandwidth
        
        double gain;
        double preBLTgain;
        
        int nba;
        double * ba;
    };

    #pragma mark Biquad

    Biquad::Biquad(){
    }
    Biquad::~Biquad(){
    }


    void Biquad::DF2TFourthOrderSection(double B0, double B1, double B2, double B3, double B4,
                                        double A0, double A1, double A2, double A3, double A4){
        b0 = B0  / A0;
        b1 = B1  / A0;
        b2 = B2  / A0;
        b3 = B3  / A0;
        b4 = B4  / A0;
        
        a1 = (-A1) / A0;  // The negation conforms the Direct-Form II Transposed discrete-time
        a2 = (-A2) / A0;  // filter (DF2T) coefficients to the expectations of the process function.
        a3 = (-A3) / A0;
        a4 = (-A4) / A0;
    }


    void Biquad::DF2TBiquad(double B0, double B1, double B2,
                            double A0, double A1, double A2){
        b0 = B0  / A0;
        b1 = B1  / A0;
        b2 = B2  / A0;
        a1 = (-A1) / A0;  // The negation conforms the Direct-Form II Transposed discrete-time
        a2 = (-A2) / A0;  // filter (DF2T) coefficients to the expectations of the process function.
    }



    #pragma mark - BiquadChain

    void BiquadChain::allocate(int count){
        
        numFilters = count;
        _yn.resize(numFilters);
        _yn1.resize(numFilters);
        _yn2.resize(numFilters);
        
        // Fourth order sections
        _yn3.resize(numFilters);
        _yn4.resize(numFilters);
    }


    BiquadChain::BiquadChain() : numFilters(0){
    }


    BiquadChain::BiquadChain(int count){
        allocate(count);
        reset();
    }


    void BiquadChain::resize(int count){
        allocate(count);
    }


    void BiquadChain::reset(){
        
        _xn1 = 0;
        _xn2 = 0;
        
        for(int i = 0; i < numFilters; i++){
            _yn[i] = 0;
            _yn1[i] = 0;
            _yn2[i] = 0;
            
            // Fourth order sections
            _yn3[i] = 0;
            _yn4[i] = 0;
        }
    }


    void BiquadChain::processBiquad(const float * input, float * output, const int stride, const int count, const Biquad * coeffs){
        
        double * yn = &_yn[0];
        double * yn1 = &_yn1[0];
        double * yn2 = &_yn2[0];
        
        for(int n = 0; n < count; n++){
            double xn = *input;
            
            yn[0] = coeffs[0].b0 * xn + coeffs[0].b1 * _xn1 + coeffs[0].b2 * _xn2
            + coeffs[0].a1 * yn1[0] + coeffs[0].a2 * yn2[0];
            
            for(int i = 1; i < numFilters; i++){
                yn[i] = coeffs[i].b0 * yn[i - 1] + coeffs[i].b1 * yn1[i - 1] + coeffs[i].b2 * yn2[i - 1]
                + coeffs[i].a1 * yn1[i] + coeffs[i].a2 * yn2[i];
            }
            
            // Shift delay line elements.
            for(int i = 0; i < numFilters; i++){
                yn2[i] = yn1[i];
                yn1[i] = yn[i];
            }
            _xn2 = _xn1;
            _xn1 = xn;
            
            // Store result and stride
            *output = yn[numFilters - 1];
            
            input += stride;
            output += stride;
        }
    }


    void BiquadChain::processFourthOrderSections(const float * input, float * output, const int stride, const int count, const Biquad * coeffs){
        
        double * yn = &_yn[0];
        double * yn1 = &_yn1[0];
        double * yn2 = &_yn2[0];
        double * yn3 = &_yn3[0];
        double * yn4 = &_yn4[0];
        
        for(int n = 0; n < count; n++){
            double xn = *input;
            
            yn[0] = coeffs[0].b0 * xn
            + coeffs[0].b1 * _xn1
            + coeffs[0].b2 * _xn2
            + coeffs[0].b3 * xn3
            + coeffs[0].b4 * xn4 +
            
            coeffs[0].a1 * yn1[0]
            + coeffs[0].a2 * yn2[0]
            + coeffs[0].a3 * yn3[0]
            + coeffs[0].a4 * yn4[0];
            
            for(int i = 1; i < numFilters; i++){
                yn[i] = coeffs[i].b0 * yn[i - 1]
                + coeffs[i].b1 * yn1[i - 1]
                + coeffs[i].b2 * yn2[i - 1]
                + coeffs[i].b3 * yn3[i - 1]
                + coeffs[i].b4 * yn4[i - 1] +
                
                coeffs[i].a1 * yn1[i]
                + coeffs[i].a2 * yn2[i]
                + coeffs[i].a3 * yn3[i]
                + coeffs[i].a4 * yn4[i];
            }
            
            // Shift delay line elements.
            for(int i = 0; i < numFilters; i++){
                yn4[i] = yn3[i];
                yn3[i] = yn2[i];
                yn2[i] = yn1[i];
                yn1[i] = yn[i];
            }
            
            xn4 = xn3;
            xn3 = _xn2;
            _xn2 = _xn1;
            _xn1 = xn;
            
            // Store result and stride
            *output = yn[numFilters - 1];
            
            input += stride;
            output += stride;
        }
    }

    #define LOG_OUTPUT 0    // Enable output logging

    #pragma mark - Public

    //******************************************************************************
    // Lowpass analogue prototype. Places Butterworth poles evenly around
    // the S-plane unit circle.
    //
    // Reference: MATLAB buttap(filterOrder)
    //******************************************************************************

    vector <complex_double>
    Butterworth::prototypeAnalogLowPass(int filterOrder){
        
        vector <complex_double> poles;
        
        for(uint32_t k = 0; k < (filterOrder + 1) / 2; k++){
            double theta = (double)(2 * k + 1) * M_PI / (2 * filterOrder);
            double real = -sin(theta);
            double imag = cos(theta);
            poles.push_back(complex_double(real,  imag));
            poles.push_back(complex_double(real, -imag)); // conjugate
        }
        
        return poles;
    }


    //******************************************************************************
    // Generate Butterworth coefficients, the main public method
    //
    // Reference: MATLAB butter(n, Wn, ...)
    //            http://en.wikipedia.org/wiki/Butterworth_filter
    //******************************************************************************

    bool Butterworth::coefficients(FILTER_TYPE filter, const double fs, const double freq1_cutoff,
                                const double freq2_cutoff, const int filterOrder,
                                vector <Biquad> & coeffs,
                                double & overallGain){
        
        //******************************************************************************
        // Init internal state based on filter design requirements
        
        zeros.resize(2 * filterOrder);
        poles.resize(2 * filterOrder);
        
        f1 = freq1_cutoff;
        f2 = freq2_cutoff;
        
        Wc = 0;  // Omega cutoff = passband edge freq
        bw = 0;
        
        
        //******************************************************************************
        // Prewarp
        
        f1 = 2 * tan(M_PI * f1 / fs);
        f2 = 2 * tan(M_PI * f2 / fs);
        
    #if DEBUG && LOG_OUTPUT
        cout << endl << endl;
        cout << "[Butterworth Filter Design] prewarped f1 = " << f1 << endl;
        cout << "[Butterworth Filter Design] prewarped f2 = " << f2 << endl;
    #endif
        
        
        //******************************************************************************
        // Design basic S-plane poles-only analogue LP prototype
        
        // Get zeros & poles of prototype analogue low pass.
        vector <complex_double> tempPoles = prototypeAnalogLowPass(filterOrder);
        
        // Copy tmppole into poles
        int index = 0;
        for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
            poles[index] = *itr;
            index++;
        }
        
        numPoles = (int)tempPoles.size();
        numZeros = 0;       // butterworth LP prototype has no zeros.
        gain = 1.0;         // always 1 for the butterworth prototype lowpass.
        
        
        //******************************************************************************
        // Convert prototype to target filter type (LP/HP/BP/BS) - S-plane
        
        // Re-orient BP/BS corner frequencies if necessary
        if(f1 > f2){
            double temp = f2;
            f2 = f1;
            f1 = temp;
        }
        
        // Cutoff Wc = f2
        switch(filter){
                
            case kLoPass:
                convert2lopass();
                break;
                
            case kHiPass:
                convert2hipass();
                break;
                
            case kBandPass:
                convert2bandpass();
                break;
                
            case kBandStop:
                convert2bandstop();
                break;
                
            default: {
    #if LOG_OUTPUT
                cout << "[Butterworth Filter Design] Unknown Filter Type" << endl;
    #endif
                return false;
            }
        }
        
        
        //******************************************************************************
        // SANITY CHECK: Ensure poles are in the left half of the S-plane
        for(uint32_t i = 0; i < numPoles; i++){
            if(poles[i].real() > 0){
    #if LOG_OUTPUT
                cerr << "[Butterworth Filter Design] Error: poles must be in the left half plane" << endl;
    #endif
                return false;
            }
        }
        
        
        //******************************************************************************
        // Map zeros & poles from S-plane to Z-plane
        
        nba = 0;
        ba = new double[2 * std::max(numPoles, numZeros) + 5];
        preBLTgain = gain;
        
        if(!s2Z()){
    #if LOG_OUTPUT
            cerr << "[Butterworth Filter Design] Error: s2Z failed" << endl << endl;
    #endif
            return false;
        }
        
        
        //******************************************************************************
        // Split up Z-plane poles and zeros into SOS
        
        if(!zp2SOS()){
    #if LOG_OUTPUT
            cerr << "[Butterworth Filter Design] Error: zp2SOS failed" << endl << endl;
    #endif
            return false;
        }
        
        // correct the overall gain
        if(filter == kLoPass || filter == kBandPass){ // pre-blt is okay for S-plane
            ba[0] = preBLTgain * (preBLTgain / gain); // 2nd term is how much BLT boosts,
        }
        else if(filter == kHiPass || kBandStop){ // HF gain != DC gain
            ba[0] = 1 / ba[0];
        }
        
        
        //******************************************************************************
        // Init biquad chain with coefficients from SOS
        
        overallGain = ba[0];
        int numFilters = filterOrder / 2;
        if(filter == kBandPass || filter == kBandStop){
            numFilters = filterOrder; // we have double the # of biquad sections
            
            // IOHAVOC filterOrder is never used again? figure this out FIXME
            // filterOrder *= 2;
        }
        
        coeffs.resize(numFilters);
        for(uint32_t i = 0; i < numFilters; i++){
            (coeffs)[i].DF2TBiquad(1.0,                     // b0
                                ba[4 * i + 1],           // b1
                                ba[4 * i + 2],           // b2
                                1.0,                     // a0
                                ba[4 * i + 3],           // a1
                                ba[4 * i + 4]);          // a2
        }
        
    #if LOG_OUTPUT
        
        // ba[0] contains your gain factor
        cout << endl;
        cout << "[Butterworth Filter Design] preBLTgain = " << preBLTgain << endl;
        cout << "[Butterworth Filter Design] gain = " << gain << endl;
        
        cout << "[Butterworth Filter Design] ba[0] = " << ba[0] << endl;
        cout << "[Butterworth Filter Design] coeff size = " << nba << endl << endl;
        
        for(uint32_t i = 0; i < (nba - 1) / 4; i++){
            // b0,b1,b2: a0,a1,a2:= 1.0, ba[4*i+1], ba[4*i+2], 1.0, ba[4*i+3], ba[4*i+4]
            cout << "[Butterworth Filter Design] Biquads:= 1.0 " << ba[4 * i + 1] << " "
            << ba[4 * i + 2] << " "
            << ba[4 * i + 3] << " "
            << ba[4 * i + 4] << endl;
        }
    #endif
        
        return true;
    }


    //******************************************************************************
    // High-Order Equalizer Filters: Low-Shelf, High-Shelf & Parametric Boost-Cut
    //
    // Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer
    //            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
    //            http://eceweb1.rutgers.edu/~orfanidi/hpeq/
    ////******************************************************************************

    bool Butterworth::coefficientsEQ(FILTER_TYPE filter, double fs, double f1,
                                    double f2, int filterOrder,
                                    vector <Biquad> & coeffs,
                                    double overallGain){
        
        // Convert band edges to radians/second
        double w1 = 2.0 * M_PI * (f1 / fs);
        double w2 = 2.0 * M_PI * (f2 / fs);
        
        // Compute center frequency w0 & bandwidth Dw
        // for parametric case in radians/sample
        double Dw = w2 - w1;
        double w0 = acos(sin(w1 + w2) / (sin(w1) + sin(w2)));
        if(w2 == M_PI){
            w0 = M_PI;
        }
        
        // Setup gain
        double G0 = 1.0;        // Reference gain == 0 dB
        double G = overallGain;
        double GB = 0.75 * G;   // Setup Peak-to-Bandwidth Gain Ratio. What about the 3-dB point??
        
        G = pow(10, (G  / 20.0));     // G  = 10^(G/20);
        GB = pow(10, (GB / 20.0));    // GB = 10^(GB/20);
        
        //  Do not proceed with design if G == G0
        if(G == G0){
            return true;
        }
        
        int L = filterOrder / 2;
        double c0 = cos(w0);
        
        if(w0 == 0){
            c0 = 1.0;
        }
        if(w0 == M_PI / 2){
            c0 = 0.0;
        }
        if(w0 == M_PI){
            c0 = -1.0;
        }
        
        double WB = tan(Dw / 2.0);
        double epsilon = sqrt(((G * G) - (GB * GB)) / ((GB * GB) - (G0 * G0)));
        
        double g = pow(G,  (1.0 / filterOrder));
        double g0 = pow(G0, (1.0 / filterOrder));
        
        // Butterworth
        double b = WB / pow(epsilon, (1.0 / filterOrder));
        
        // Ensure size 'L' of coeff vector is correct!
        coeffs.resize(L);
        for(uint32_t i = 1; i <= L; i++){
            double phi = (2 * i - 1.0) * M_PI / (2.0 * filterOrder);
            double si = sin(phi);
            double D = b * b + 2.0 * b * si + 1.0;
            
            if(filter == kLoShelf || filter == kHiShelf){ // Compute SOS coefficients
                
                double b0h = (g * g * b * b + 2 * g0 * g * b * si + g0 * g0) / D;
                double b1h = (filter == kHiShelf) ? -2.0 * (g * g * b * b - g0 * g0) / D : 2.0 * (g * g * b * b - g0 * g0) / D;
                double b2h = (g * g * b * b - 2 * g0 * g * b * si + g0 * g0) / D;
                double a1h = (filter == kHiShelf) ? -2.0 * (b * b - 1.0) / D : 2.0 * (b * b - 1.0) / D;
                double a2h = (b * b - 2 * b * si + 1.0) / D;
                
                // High-order HP/LP shelving filter coefficients can be expressed as 2nd-order sections (SOS)
                // i.e. biquads
                (coeffs)[i - 1].DF2TBiquad(b0h,         // b0
                                        b1h,         // b1
                                        b2h,         // b2
                                        1.0,         // a0
                                        a1h,         // a1
                                        a2h);        // a2
            }
            else if(filter == kParametric){ // Compute 4th order sections
                
                double b0 = (g * g * b * b + g0 * g0 + 2 * g * g0 * si * b) / D;
                double b1 = -4 * c0 * (g0 * g0 + g * g0 * si * b) / D;
                double b2 = 2 * ((1 + 2 * c0 * c0) * g0 * g0 - g * g * b * b) / D;
                double b3 = -4 * c0 * (g0 * g0 - g * g0 * si * b) / D;
                double b4 = (g * g * b * b + g0 * g0 - 2 * g * g0 * si * b) / D;
                double a1 = -4 * c0 * (1 + si * b) / D;
                double a2 = 2 * (1 + 2 * c0 * c0 - b * b) / D;
                double a3 = -4 * c0 * (1 - si * b) / D;
                double a4 = (b * b - 2 * si * b + 1) / D;
                
                // Parameteric EQ filter coefficients (like band pass & band stop) are twice the
                // specified filter order. This is normal and by design. Unlike bandpass & bandstop
                // though, the realization via the Bilinear Transform (BLT) renders 4th order sections.
                // So rather than split 4th order sections into 2nd order sections (biquads),
                // with fancy polynomial root factoring, we use them as is.
                // There are no stability issues for sections this size.
                (coeffs)[i - 1].DF2TFourthOrderSection(b0,          // b0
                                                    b1,          // b1
                                                    b2,          // b2
                                                    b3,          // b3
                                                    b4,          // b4
                                                    1.0,         // a0
                                                    a1,          // a1
                                                    a2,          // a2
                                                    a3,          // a3
                                                    a4);         // a4
            }
        }
        
        return true;
    }




    #pragma mark - Filter design utility methods

    //******************************************************************************
    //
    // Z = (2 + S) / (2 - S) is the S-plane to Z-plane bilinear transform
    //
    // Reference: http://en.wikipedia.org/wiki/Bilinear_transform
    //
    //******************************************************************************

    double Butterworth::blt(complex_double & sz){
        
        complex_double two(2.0, 0);
        complex_double s = sz;      // sz is an input(S-plane) & output(Z-plane) arg
        sz = (two + s) / (two - s);
        
        // return the gain
        return abs((two - s));
    }


    //******************************************************************************
    //
    // Convert poles & zeros from S-plane to Z-plane via Bilinear Tranform (BLT)
    //
    //******************************************************************************

    bool Butterworth::s2Z(){
        
        // blt zeros
        for(uint32_t i = 0; i < numZeros; i++){
            gain /= blt(zeros[i]);
        }
        
        // blt poles
        for(uint32_t i = 0; i < numPoles; i++){
            gain *= blt(poles[i]);
        }
        
        return true;
    }


    //******************************************************************************
    //
    // Convert filter poles and zeros to second-order sections
    //
    // Reference: http://www.mathworks.com/help/signal/ref/zp2sos.html
    //
    //******************************************************************************

    bool Butterworth::zp2SOS(){
        
        int filterOrder = std::max(numZeros, numPoles);
        complex_double * zerosTempVec = new complex_double[filterOrder];
        complex_double * polesTempVec = new complex_double[filterOrder];
        
        // Copy
        for(uint32_t i = 0; i < numZeros; i++){
            zerosTempVec[i] = zeros[i];
        }
        
        // Add zeros at -1, so if S-plane degenerate case where
        // numZeros = 0 will map to -1 in Z-plane.
        for(uint32_t i = numZeros; i < filterOrder; i++){
            zerosTempVec[i] = complex_double(-1, 0);
        }
        
        // Copy
        for(uint32_t i = 0; i < numPoles; i++){
            polesTempVec[i] = poles[i];
        }
        
        ba[0] = gain; // store gain
        
        int numSOS = 0;
        for(uint32_t i = 0; i + 1 < filterOrder; i += 2, numSOS++){
            ba[4 * numSOS + 1] = -(zerosTempVec[i] + zerosTempVec[i + 1]).real();
            ba[4 * numSOS + 2] =  (zerosTempVec[i] * zerosTempVec[i + 1]).real();
            ba[4 * numSOS + 3] = -(polesTempVec[i] + polesTempVec[i + 1]).real();
            ba[4 * numSOS + 4] =  (polesTempVec[i] * polesTempVec[i + 1]).real();
        }
        
        // Odd filter order thus one pair of poles/zeros remains
        if(filterOrder % 2 == 1){
            ba[4 * numSOS + 1] = -zerosTempVec[filterOrder - 1].real();
            ba[4 * numSOS + 2] = ba[4 * numSOS + 4] = 0;
            ba[4 * numSOS + 3] = -polesTempVec[filterOrder - 1].real();
            numSOS++;
        }
        
        // Set output param
        nba = 1 + 4 * numSOS;
        
        delete[] zerosTempVec;
        delete[] polesTempVec;
        return true;
    }




    #pragma mark - Analog lowpss prototype conversion methods

    //******************************************************************************
    // Convert analog lowpass prototype poles to lowpass
    //******************************************************************************

    void Butterworth::convert2lopass(){
        Wc = f2;                                   // critical frequency
        
        gain *= pow(Wc, numPoles);
        
        numZeros = 0;                              // poles only
        for(uint32_t i = 0; i < numPoles; i++){    // scale poles by the cutoff Wc
            poles[i] = Wc * poles[i];
        }
    }


    //******************************************************************************
    // Convert lowpass poles & zeros to highpass
    // with Wc = f2, use:  hp_S = Wc / lp_S;
    //******************************************************************************

    void Butterworth::convert2hipass(){
        Wc = f2;   // Critical frequency
        
        // Calculate gain
        complex_double prodz(1.0, 0.0);
        complex_double prodp(1.0, 0.0);
        
        for(uint32_t i = 0; i < numZeros; i++){
            prodz *= -zeros[i];
        }
        
        for(uint32_t i = 0; i < numPoles; i++){
            prodp *= -poles[i];
        }
        
        gain *= prodz.real() / prodp.real();
        
        // Convert LP poles to HP
        for(uint32_t i = 0; i < numPoles; i++){
            if(abs(poles[i])){
                poles[i] = complex_double(Wc) / poles[i]; //  hp_S = Wc / lp_S;
                
            }
        }
        
        // Init with zeros, no non-zero values to convert
        numZeros = numPoles;
        for(uint32_t i = 0; i < numZeros; i++){
            zeros[i] = complex_double(0.0);
        }
    }


    //******************************************************************************
    // Convert lowpass poles to bandpass
    // use:  bp_S = 0.5 * lp_S * BW +
    //				   0.5 * sqrt ( BW^2 * lp_S^2 - 4*Wc^2 )
    // where   BW = W2 - W1
    //		    Wc^2 = W2 * W1
    //******************************************************************************

    void Butterworth::convert2bandpass(){
        bw = f2 - f1;
        Wc = sqrt(f1 * f2);
        
        // Calculate bandpass gain
        gain *= pow(bw, numPoles - numZeros);
        
        // Convert LP poles to BP: these two sets of for-loops result in an ordered
        // list of poles and their complex conjugates
        vector <complex_double> tempPoles;
        for(uint32_t i = 0; i < numPoles; i++){    // First set of poles + conjugates
            if(abs(poles[i])){
                complex_double firstterm = 0.5 * poles[i] * bw;
                complex_double secondterm = 0.5 * sqrt((bw * bw) * (poles[i] * poles[i]) - 4 * Wc * Wc);
                tempPoles.push_back(firstterm + secondterm);
            }
        }
        
        for(uint32_t i = 0; i < numPoles; i++){    // Second set of poles + conjugates
            if(abs(poles[i])){
                complex_double firstterm = 0.5 * poles[i] * bw;
                complex_double secondterm = 0.5 * sqrt((bw * bw) * (poles[i] * poles[i]) - 4 * Wc * Wc);
                tempPoles.push_back(firstterm - secondterm);      // complex conjugate
            }
        }
        
        // Init zeros, no non-zero values to convert
        numZeros = numPoles;
        for(uint32_t i = 0; i < numZeros; i++){
            zeros[i] = complex_double(0.0);
        }
        
        // Copy converted poles to output array
        int index = 0;
        for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
            poles[index] = *itr;
            index++;
        }
        numPoles = (int)tempPoles.size();
    }


    //******************************************************************************
    // Convert lowpass poles to bandstop
    // use:  bs_S = 0.5 * BW / lp_S +
    //				   0.5 * sqrt ( BW^2 / lp_S^2 - 4*Wc^2 )
    // where   BW = W2 - W1
    //		 Wc^2 = W2 * W1
    //******************************************************************************

    void Butterworth::convert2bandstop(){
        bw = f2 - f1;
        Wc = sqrt(f1 * f2);
        
        // Compute gain
        complex_double prodz(1.0, 0.0);
        complex_double prodp(1.0, 0.0);
        for(uint32_t i = 0; i < numZeros; i++){
            prodz *= -zeros[i];
        }
        
        for(uint32_t i = 0; i < numPoles; i++){
            prodp *= -poles[i];
        }
        
        gain *= prodz.real() / prodp.real();
        
        // Convert LP zeros to band stop
        numZeros = numPoles;
        vector <complex_double> ztmp;
        for(uint32_t i = 0; i < numZeros; i++){
            ztmp.push_back(complex_double(0.0,  Wc));
            ztmp.push_back(complex_double(0.0, -Wc));         // complex conjugate
        }
        
        vector <complex_double> tempPoles;
        for(uint32_t i = 0; i < numPoles; i++){    // First set of poles + conjugates
            if(abs(poles[i])){
                complex_double term1 = 0.5 * bw / poles[i];
                complex_double term2 = 0.5 * sqrt((bw * bw) / (poles[i] * poles[i]) - (4 * Wc * Wc));
                tempPoles.push_back(term1 + term2);
            }
        }
        
        for(uint32_t i = 0; i < numPoles; i++){    // Second set of poles + conjugates
            if(abs(poles[i])){
                complex_double term1 = 0.5 * bw / poles[i];
                complex_double term2 = 0.5 * sqrt((bw * bw) / (poles[i] * poles[i]) - (4 * Wc * Wc));
                tempPoles.push_back(term1 - term2);       // complex conjugate
            }
        }
        
        // Copy converted zeros to output array
        int index = 0;
        for(vector <complex_double>::iterator itr = ztmp.begin(); itr != ztmp.end(); itr++){
            zeros[index] = *itr;
            index++;
        }
        numZeros = (int)ztmp.size();
        
        // Copy converted poles to output array
        index = 0;
        for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
            poles[index] = *itr;
            index++;
        }
        numPoles = (int)tempPoles.size();
    }
}