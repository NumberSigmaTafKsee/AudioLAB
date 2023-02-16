#pragma once

#include "qmdsp_math.hpp"
#include "qmdsp_fft.hpp"

namespace qmdsp
{
    struct CQConfig {
        double FS;         // samplerate
        double min;        // minimum frequency
        double max;        // maximum frequency
        int BPO;           // bins per octave
        double CQThresh;   // threshold
    };

    class ConstantQ
    {
    public:
        ConstantQ(CQConfig config);
        ~ConstantQ();

        void process(const double* FFTRe, const double* FFTIm,
                    double* CQRe, double* CQIm);

        double* process(const double* FFTData);

        void sparsekernel();

        double getQ() { return m_dQ; }
        int getK() { return m_uK; }
        int getFFTLength() { return m_FFTLength; }
        int getHop() { return m_hop; }

    private:
        void initialise(CQConfig config);
        void deInitialise();
            
        double* m_CQdata;
        double m_FS;
        double m_FMin;
        double m_FMax;
        double m_dQ;
        double m_CQThresh;
        int m_hop;
        int m_BPO;
        int m_FFTLength;
        int m_uK;

        struct SparseKernel {
            std::vector<int> is;
            std::vector<int> js;
            std::vector<double> imag;
            std::vector<double> real;
        };

        SparseKernel *m_sparseKernel;
    };

    struct ChromaConfig {
        double FS;
        double min;
        double max;
        int BPO;
        double CQThresh;
        MathUtilities::NormaliseType normalise;
    };

    class Chromagram 
    {
    public: 
        Chromagram( ChromaConfig Config );
        ~Chromagram();

        /**
        * Process a time-domain input signal of length equal to
        * getFrameSize().
        * 
        * The returned buffer contains the chromagram values indexed by
        * bin, with the number of values corresponding to the BPO field
        * in the ChromaConfig supplied at construction. It is owned by
        * the Chromagram object and is reused from one process call to
        * the next.
        */
        double *process(const double *data);
        
        /**
        * Process a frequency-domain input signal generated from a
        * time-domain signal of length equal to getFrameSize() that has
        * been windowed and "fftshifted" to place the zero index in the
        * centre of the frame. The real and imag buffers must each
        * contain the full getFrameSize() frequency bins.
        * 
        * The returned buffer contains the chromagram values indexed by
        * bin, with the number of values corresponding to the BPO field
        * in the ChromaConfig supplied at construction. It is owned by
        * the Chromagram object and is reused from one process call to
        * the next.
        */
        double *process(const double *real, const double *imag);
        
        void unityNormalise(double* src);

        // Complex arithmetic
        double kabs( double real, double imag );
            
        // Results
        int getK() { return m_uK;}
        int getFrameSize() { return m_frameSize; }
        int getHopSize()   { return m_hopSize; }
        
    private:
        int initialise( ChromaConfig Config );
        int deInitialise();

        Window<double> *m_window;
        double *m_windowbuf;
            
        double* m_chromadata;
        double m_FMin;
        double m_FMax;
        int m_BPO;
        int m_uK;

        MathUtilities::NormaliseType m_normalise;

        int m_frameSize;
        int m_hopSize;

        FFTReal* m_FFT;
        ConstantQ* m_ConstantQ;

        double* m_FFTRe;
        double* m_FFTIm;
        double* m_CQRe;
        double* m_CQIm;

        bool m_skGenerated;
    };


    //----------------------------------------------------------------------------

    ConstantQ::ConstantQ( CQConfig config ) :
        m_sparseKernel(0)
    {
        initialise(config);
    }

    ConstantQ::~ConstantQ()
    {
        deInitialise();
    }

    static double squaredModule(const double & xx, const double & yy) {
        return xx*xx + yy*yy;
    }

    void ConstantQ::sparsekernel()
    {
        SparseKernel *sk = new SparseKernel();

        double* windowRe = new double [ m_FFTLength ];
        double* windowIm = new double [ m_FFTLength ];
        double* transfWindowRe = new double [ m_FFTLength ];
        double* transfWindowIm = new double [ m_FFTLength ];
            
        // for each bin value K, calculate temporal kernel, take its fft
        // to calculate the spectral kernel then threshold it to make it
        // sparse and add it to the sparse kernels matrix
        
        double squareThreshold = m_CQThresh * m_CQThresh;

        FFT fft(m_FFTLength);
            
        for (int j = m_uK - 1; j >= 0; --j) {
            
            for (int i = 0; i < m_FFTLength; ++i) {
                windowRe[i] = 0;
                windowIm[i] = 0;
            }

            // Compute a complex sinusoid windowed with a hamming window
            // of the right length
            
            int windowLength = (int)ceil
                (m_dQ * m_FS / (m_FMin * pow(2, (double)j / (double)m_BPO)));

            int origin = m_FFTLength/2 - windowLength/2;

            for (int i = 0; i < windowLength; ++i) {
                double angle = (2.0 * M_PI * m_dQ * i) / windowLength;
                windowRe[origin + i] = cos(angle);
                windowIm[origin + i] = sin(angle);
            }

            // Shape with hamming window
            Window<double> hamming(HammingWindow, windowLength);
            hamming.cut(windowRe + origin);
            hamming.cut(windowIm + origin);

            // Scale
            for (int i = 0; i < windowLength; ++i) {
                windowRe[origin + i] /= windowLength;
            }
            for (int i = 0; i < windowLength; ++i) {
                windowIm[origin + i] /= windowLength;
            }

            // Input is expected to have been fftshifted, so do the
            // same to the input to the fft that contains the kernel
            for (int i = 0; i < m_FFTLength/2; ++i) {
                double temp = windowRe[i];
                windowRe[i] = windowRe[i + m_FFTLength/2];
                windowRe[i + m_FFTLength/2] = temp;
            }
            for (int i = 0; i < m_FFTLength/2; ++i) {
                double temp = windowIm[i];
                windowIm[i] = windowIm[i + m_FFTLength/2];
                windowIm[i + m_FFTLength/2] = temp;
            }
        
            fft.process(false, windowRe, windowIm, transfWindowRe, transfWindowIm);

            // convert to sparse form
            for (int i = 0; i < m_FFTLength; i++) {

                // perform thresholding
                double mag = squaredModule(transfWindowRe[i], transfWindowIm[i]);
                if (mag <= squareThreshold) continue;
                    
                // Insert non-zero position indexes
                sk->is.push_back(i);
                sk->js.push_back(j);

                // take conjugate, normalise and add to array for sparse kernel
                sk->real.push_back( transfWindowRe[i] / m_FFTLength);
                sk->imag.push_back(-transfWindowIm[i] / m_FFTLength);
            }
        }

        delete [] windowRe;
        delete [] windowIm;
        delete [] transfWindowRe;
        delete [] transfWindowIm;

        m_sparseKernel = sk;
    }

    void ConstantQ::initialise( CQConfig Config )
    {
        m_FS = Config.FS;             // Sample rate
        m_FMin = Config.min;          // Minimum frequency
        m_FMax = Config.max;          // Maximum frequency
        m_BPO = Config.BPO;           // Bins per octave
        m_CQThresh = Config.CQThresh; // Threshold for sparse kernel generation

        // Q value for filter bank
        m_dQ = 1/(pow(2,(1/(double)m_BPO))-1);

        // No. of constant Q bins
        m_uK = (int)ceil(m_BPO * log(m_FMax/m_FMin)/log(2.0));

        // Length of fft required for this Constant Q filter bank
        m_FFTLength = MathUtilities::nextPowerOfTwo(int(ceil(m_dQ * m_FS/m_FMin)));

        // Hop from one frame to next
        m_hop = m_FFTLength / 8;

        // allocate memory for cqdata
        m_CQdata = new double [2*m_uK];
    }

    void ConstantQ::deInitialise()
    {
        delete [] m_CQdata;
        delete m_sparseKernel;
    }

    //-----------------------------------------------------------------------------
    double* ConstantQ::process( const double* fftdata )
    {
        if (!m_sparseKernel) {
            std::cerr << "ERROR: ConstantQ::process: Sparse kernel has not been initialised" << std::endl;
            return m_CQdata;
        }

        SparseKernel *sk = m_sparseKernel;

        for (int row=0; row < 2 * m_uK; row++) {
            m_CQdata[ row ] = 0;
            m_CQdata[ row+1 ] = 0;
        }
        const int *fftbin = &(sk->is[0]);
        const int *cqbin = &(sk->js[0]);
        const double *real = &(sk->real[0]);
        const double *imag = &(sk->imag[0]);
        const int sparseCells = int(sk->real.size());
            
        for (int i = 0; i < sparseCells; i++) {
            const int row = cqbin[i];
            const int col = fftbin[i];
            if (col == 0) continue;
            const double & r1 = real[i];
            const double & i1 = imag[i];
            const double & r2 = fftdata[ (2*m_FFTLength) - 2*col - 2 ];
            const double & i2 = fftdata[ (2*m_FFTLength) - 2*col - 2 + 1 ];
            // add the multiplication
            m_CQdata[ 2*row  ] += (r1*r2 - i1*i2);
            m_CQdata[ 2*row+1] += (r1*i2 + i1*r2);
        }

        return m_CQdata;
    }

    void ConstantQ::process(const double *FFTRe, const double* FFTIm,
                            double *CQRe, double *CQIm)
    {
        if (!m_sparseKernel) {
            std::cerr << "ERROR: ConstantQ::process: Sparse kernel has not been initialised" << std::endl;
            return;
        }

        SparseKernel *sk = m_sparseKernel;

        for (int row = 0; row < m_uK; row++) {
            CQRe[ row ] = 0;
            CQIm[ row ] = 0;
        }

        const int *fftbin = &(sk->is[0]);
        const int *cqbin = &(sk->js[0]);
        const double *real = &(sk->real[0]);
        const double *imag = &(sk->imag[0]);
        const int sparseCells = int(sk->real.size());
            
        for (int i = 0; i<sparseCells; i++) {
            const int row = cqbin[i];
            const int col = fftbin[i];
            if (col == 0) continue;
            const double & r1 = real[i];
            const double & i1 = imag[i];
            const double & r2 = FFTRe[ m_FFTLength - col ];
            const double & i2 = FFTIm[ m_FFTLength - col ];
            // add the multiplication
            CQRe[ row ] += (r1*r2 - i1*i2);
            CQIm[ row ] += (r1*i2 + i1*r2);
        }
    }


    //----------------------------------------------------------------------------

    Chromagram::Chromagram( ChromaConfig Config ) :
        m_skGenerated(false)
    {
        initialise( Config );
    }

    int Chromagram::initialise( ChromaConfig Config )
    {       
        m_FMin = Config.min;                // min freq
        m_FMax = Config.max;                // max freq
        m_BPO  = Config.BPO;                // bins per octave
        m_normalise = Config.normalise;     // if frame normalisation is required

        // Extend range to a full octave
        double octaves = log(m_FMax / m_FMin) / log(2.0);
        m_FMax = m_FMin * pow(2.0, ceil(octaves));

        // Create array for chroma result
        m_chromadata = new double[ m_BPO ];

        // Create Config Structure for ConstantQ operator
        CQConfig ConstantQConfig;

        // Populate CQ config structure with parameters
        // inherited from the Chroma config
        ConstantQConfig.FS = Config.FS;
        ConstantQConfig.min = m_FMin;
        ConstantQConfig.max = m_FMax;
        ConstantQConfig.BPO = m_BPO;
        ConstantQConfig.CQThresh = Config.CQThresh;
            
        // Initialise ConstantQ operator
        m_ConstantQ = new ConstantQ( ConstantQConfig );

        // No. of constant Q bins
        m_uK = m_ConstantQ->getK();

        // Initialise working arrays
        m_frameSize = m_ConstantQ->getFFTLength();
        m_hopSize = m_ConstantQ->getHop();

        // Initialise FFT object    
        m_FFT = new FFTReal(m_frameSize);

        m_FFTRe = new double[ m_frameSize ];
        m_FFTIm = new double[ m_frameSize ];
        m_CQRe  = new double[ m_uK ];
        m_CQIm  = new double[ m_uK ];

        m_window = 0;
        m_windowbuf = 0;

        return 1;
    }

    Chromagram::~Chromagram()
    {
        deInitialise();
    }

    int Chromagram::deInitialise()
    {
        delete[] m_windowbuf;
        delete m_window;

        delete [] m_chromadata;

        delete m_FFT;

        delete m_ConstantQ;

        delete [] m_FFTRe;
        delete [] m_FFTIm;
        delete [] m_CQRe;
        delete [] m_CQIm;
        return 1;
    }

    //----------------------------------------------------------------------------------
    // returns the absolute value of complex number xx + i*yy
    double Chromagram::kabs(double xx, double yy)
    {
        double ab = sqrt(xx*xx + yy*yy);
        return(ab);
    }
    //-----------------------------------------------------------------------------------


    void Chromagram::unityNormalise(double *src)
    {
        double min, max;
        double val = 0;

        MathUtilities::getFrameMinMax( src, m_BPO, & min, &max );

        for (int i = 0; i < m_BPO; i++) {
            val = src[ i ] / max;
            src[ i ] = val;
        }
    }


    double *Chromagram::process(const double *data)
    {
        if (!m_skGenerated) {
            // Generate CQ Kernel 
            m_ConstantQ->sparsekernel();
            m_skGenerated = true;
        }

        if (!m_window) {
            m_window = new Window<double>(HammingWindow, m_frameSize);
            m_windowbuf = new double[m_frameSize];
        }

        for (int i = 0; i < m_frameSize; ++i) {
            m_windowbuf[i] = data[i];
        }
        m_window->cut(m_windowbuf);

        // The frequency-domain version expects pre-fftshifted input - so
        // we must do the same here
        for (int i = 0; i < m_frameSize/2; ++i) {
            double tmp = m_windowbuf[i];
            m_windowbuf[i] = m_windowbuf[i + m_frameSize/2];
            m_windowbuf[i + m_frameSize/2] = tmp;
        }

        m_FFT->forward(m_windowbuf, m_FFTRe, m_FFTIm);

        return process(m_FFTRe, m_FFTIm);
    }

    double *Chromagram::process(const double *real, const double *imag)
    {
        if (!m_skGenerated) {
            // Generate CQ Kernel 
            m_ConstantQ->sparsekernel();
            m_skGenerated = true;
        }

        // initialise chromadata to 0
        for (int i = 0; i < m_BPO; i++) m_chromadata[i] = 0;

        // Calculate ConstantQ frame
        m_ConstantQ->process( real, imag, m_CQRe, m_CQIm );
            
        // add each octave of cq data into Chromagram
        const int octaves = m_uK / m_BPO;
        for (int octave = 0; octave < octaves; octave++) {
            int firstBin = octave*m_BPO;
            for (int i = 0; i < m_BPO; i++) {
                m_chromadata[i] += kabs( m_CQRe[ firstBin + i ],
                                        m_CQIm[ firstBin + i ]);
            }
        }

        MathUtilities::normalise(m_chromadata, m_BPO, m_normalise);

        return m_chromadata;
    }
}