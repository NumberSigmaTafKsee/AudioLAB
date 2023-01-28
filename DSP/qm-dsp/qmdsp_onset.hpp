#pragma once

#include "qmdsp_math.hpp"
#include "qmdsp_phasevocoder.hpp"


#define DF_HFC (1)
#define DF_SPECDIFF (2)
#define DF_PHASEDEV (3)
#define DF_COMPLEXSD (4)
#define DF_BROADBAND (5)

namespace qmdsp
{

    struct DFConfig{
        int stepSize; // DF step in samples
        int frameLength; // DF analysis window - usually 2*step. Must be even!
        int DFType; // type of detection function ( see defines )
        double dbRise; // only used for broadband df (and required for it)
        bool adaptiveWhitening; // perform adaptive whitening
        double whiteningRelaxCoeff; // if < 0, a sensible default will be used
        double whiteningFloor; // if < 0, a sensible default will be used
    };

    class DetectionFunction  
    {
    public:
        double* getSpectrumMagnitude();
        DetectionFunction( DFConfig config );
        virtual ~DetectionFunction();

        /**
        * Process a single time-domain frame of audio, provided as
        * frameLength samples.
        */
        double processTimeDomain(const double* samples);

        /**
        * Process a single frequency-domain frame, provided as
        * frameLength/2+1 real and imaginary component values.
        */
        double processFrequencyDomain(const double* reals, const double* imags);

    private:
        void whiten();
        double runDF();

        double HFC(int length, double* src);
        double specDiff(int length, double* src);
        double phaseDev(int length, double *srcPhase);
        double complexSD(int length, double *srcMagnitude, double *srcPhase);
        double broadband(int length, double *srcMagnitude);
            
    private:
        void initialise( DFConfig Config );
        void deInitialise();

        int m_DFType;
        int m_dataLength;
        int m_halfLength;
        int m_stepSize;
        double m_dbRise;
        bool m_whiten;
        double m_whitenRelaxCoeff;
        double m_whitenFloor;

        double* m_magHistory;
        double* m_phaseHistory;
        double* m_phaseHistoryOld;
        double* m_magPeaks;

        double* m_windowed; // Array for windowed analysis frame
        double* m_magnitude; // Magnitude of analysis frame ( frequency domain )
        double* m_thetaAngle;// Phase of analysis frame ( frequency domain )
        double* m_unwrapped; // Unwrapped phase of analysis frame

        Window<double> *m_window;
        PhaseVocoder* m_phaseVoc;   // Phase Vocoder
    };

    struct PPWinThresh
    {
        int pre;
        int post;

        PPWinThresh(int x, int y) :
            pre(x),
            post(y)
        {
        }
    };

    struct QFitThresh
    {
        double a;
        double b;
        double c;

        QFitThresh(double x, double y, double z) :
            a(x),
            b(y),
            c(z)
        {
        }
    };

    struct PPickParams
    {
        int length; // detection function length
        double tau; // time resolution of the detection function
        int alpha; // alpha-norm parameter
        double cutoff;// low-pass filter cutoff freq
        int LPOrd; // low-pass filter order
        double* LPACoeffs; // low-pass filter denominator coefficients
        double* LPBCoeffs; // low-pass filter numerator coefficients
        PPWinThresh WinT;// window size in frames for adaptive thresholding [pre post]:
        QFitThresh QuadThresh;
        float delta; // delta threshold used as an offset when computing the smoothed detection function

        PPickParams() :
            length(0),
            tau(0),
            alpha(0),
            cutoff(0),
            LPOrd(0),
            LPACoeffs(NULL),
            LPBCoeffs(NULL),
            WinT(0,0),
            QuadThresh(0,0,0),
            delta(0)
        {
        }
    };

    class PeakPicking  
    {
    public:
        PeakPicking( PPickParams Config );
        virtual ~PeakPicking();
            
        void process( double* src, int len, std::vector<int> &onsets  );

    private:
        void initialise( PPickParams Config  );
        void deInitialise();
        int  quadEval( std::vector<double> &src, std::vector<int> &idx );
            
        DFProcConfig m_DFProcessingParams;

        int m_DFLength ;
        double Qfilta ;
        double Qfiltb;
        double Qfiltc;

        double* m_workBuffer;
            
        DFProcess*  m_DFSmoothing;
    };

    //////////////////////////////////////////////////////////////////////
    // Construction/Destruction
    //////////////////////////////////////////////////////////////////////

    DetectionFunction::DetectionFunction( DFConfig config ) :
        m_window(0)
    {
        m_magHistory = NULL;
        m_phaseHistory = NULL;
        m_phaseHistoryOld = NULL;
        m_magPeaks = NULL;

        initialise( config );
    }

    DetectionFunction::~DetectionFunction()
    {
        deInitialise();
    }


    void DetectionFunction::initialise( DFConfig Config )
    {
        m_dataLength = Config.frameLength;
        m_halfLength = m_dataLength/2 + 1;

        m_DFType = Config.DFType;
        m_stepSize = Config.stepSize;
        m_dbRise = Config.dbRise;

        m_whiten = Config.adaptiveWhitening;
        m_whitenRelaxCoeff = Config.whiteningRelaxCoeff;
        m_whitenFloor = Config.whiteningFloor;
        if (m_whitenRelaxCoeff < 0) m_whitenRelaxCoeff = 0.9997;
        if (m_whitenFloor < 0) m_whitenFloor = 0.01;

        m_magHistory = new double[ m_halfLength ];
        memset(m_magHistory,0, m_halfLength*sizeof(double));
                    
        m_phaseHistory = new double[ m_halfLength ];
        memset(m_phaseHistory,0, m_halfLength*sizeof(double));

        m_phaseHistoryOld = new double[ m_halfLength ];
        memset(m_phaseHistoryOld,0, m_halfLength*sizeof(double));

        m_magPeaks = new double[ m_halfLength ];
        memset(m_magPeaks,0, m_halfLength*sizeof(double));

        m_phaseVoc = new PhaseVocoder(m_dataLength, m_stepSize);

        m_magnitude = new double[ m_halfLength ];
        m_thetaAngle = new double[ m_halfLength ];
        m_unwrapped = new double[ m_halfLength ];

        m_window = new Window<double>(HanningWindow, m_dataLength);
        m_windowed = new double[ m_dataLength ];
    }

    void DetectionFunction::deInitialise()
    {
        delete [] m_magHistory ;
        delete [] m_phaseHistory ;
        delete [] m_phaseHistoryOld ;
        delete [] m_magPeaks ;

        delete m_phaseVoc;

        delete [] m_magnitude;
        delete [] m_thetaAngle;
        delete [] m_windowed;
        delete [] m_unwrapped;

        delete m_window;
    }

    double DetectionFunction::processTimeDomain(const double *samples)
    {
        m_window->cut(samples, m_windowed);

        m_phaseVoc->processTimeDomain(m_windowed, 
                                    m_magnitude, m_thetaAngle, m_unwrapped);

        if (m_whiten) whiten();

        return runDF();
    }

    double DetectionFunction::processFrequencyDomain(const double *reals,
                                                    const double *imags)
    {
        m_phaseVoc->processFrequencyDomain(reals, imags,
                                        m_magnitude, m_thetaAngle, m_unwrapped);

        if (m_whiten) whiten();

        return runDF();
    }

    void DetectionFunction::whiten()
    {
        for (int i = 0; i < m_halfLength; ++i) {
            double m = m_magnitude[i];
            if (m < m_magPeaks[i]) {
                m = m + (m_magPeaks[i] - m) * m_whitenRelaxCoeff;
            }
            if (m < m_whitenFloor) m = m_whitenFloor;
            m_magPeaks[i] = m;
            m_magnitude[i] /= m;
        }
    }

    double DetectionFunction::runDF()
    {
        double retVal = 0;

        switch( m_DFType )
        {
        case DF_HFC:
            retVal = HFC( m_halfLength, m_magnitude);
            break;
            
        case DF_SPECDIFF:
            retVal = specDiff( m_halfLength, m_magnitude);
            break;
            
        case DF_PHASEDEV:
            // Using the instantaneous phases here actually provides the
            // same results (for these calculations) as if we had used
            // unwrapped phases, but without the possible accumulation of
            // phase error over time
            retVal = phaseDev( m_halfLength, m_thetaAngle);
            break;
            
        case DF_COMPLEXSD:
            retVal = complexSD( m_halfLength, m_magnitude, m_thetaAngle);
            break;

        case DF_BROADBAND:
            retVal = broadband( m_halfLength, m_magnitude);
            break;
        }
            
        return retVal;
    }

    double DetectionFunction::HFC(int length, double *src)
    {
        double val = 0;
        for (int i = 0; i < length; i++) {
            val += src[ i ] * ( i + 1);
        }
        return val;
    }

    double DetectionFunction::specDiff(int length, double *src)
    {
        double val = 0.0;
        double temp = 0.0;
        double diff = 0.0;

        for (int i = 0; i < length; i++) {
            
            temp = fabs( (src[ i ] * src[ i ]) - (m_magHistory[ i ] * m_magHistory[ i ]) );
                    
            diff= sqrt(temp);

            // (See note in phaseDev below.)

            val += diff;

            m_magHistory[ i ] = src[ i ];
        }

        return val;
    }


    double DetectionFunction::phaseDev(int length, double *srcPhase)
    {
        double tmpPhase = 0;
        double tmpVal = 0;
        double val = 0;

        double dev = 0;

        for (int i = 0; i < length; i++) {
            tmpPhase = (srcPhase[ i ]- 2*m_phaseHistory[ i ]+m_phaseHistoryOld[ i ]);
            dev = MathUtilities::princarg( tmpPhase );

            // A previous version of this code only counted the value here
            // if the magnitude exceeded 0.1.  My impression is that
            // doesn't greatly improve the results for "loud" music (so
            // long as the peak picker is reasonably sophisticated), but
            // does significantly damage its ability to work with quieter
            // music, so I'm removing it and counting the result always.
            // Same goes for the spectral difference measure above.
                    
            tmpVal  = fabs(dev);
            val += tmpVal ;

            m_phaseHistoryOld[ i ] = m_phaseHistory[ i ] ;
            m_phaseHistory[ i ] = srcPhase[ i ];
        }
            
        return val;
    }


    double DetectionFunction::complexSD(int length, double *srcMagnitude, double *srcPhase)
    {
        double val = 0;
        double tmpPhase = 0;
        double tmpReal = 0;
        double tmpImag = 0;
    
        double dev = 0;
        ComplexData meas = ComplexData( 0, 0 );
        ComplexData j = ComplexData( 0, 1 );

        for (int i = 0; i < length; i++) {
            
            tmpPhase = (srcPhase[ i ]- 2*m_phaseHistory[ i ]+m_phaseHistoryOld[ i ]);
            dev= MathUtilities::princarg( tmpPhase );
                    
            meas = m_magHistory[i] - ( srcMagnitude[ i ] * exp( j * dev) );

            tmpReal = real( meas );
            tmpImag = imag( meas );

            val += sqrt( (tmpReal * tmpReal) + (tmpImag * tmpImag) );
                    
            m_phaseHistoryOld[ i ] = m_phaseHistory[ i ] ;
            m_phaseHistory[ i ] = srcPhase[ i ];
            m_magHistory[ i ] = srcMagnitude[ i ];
        }

        return val;
    }

    double DetectionFunction::broadband(int length, double *src)
    {
        double val = 0;
        for (int i = 0; i < length; ++i) {
            double sqrmag = src[i] * src[i];
            if (m_magHistory[i] > 0.0) {
                double diff = 10.0 * log10(sqrmag / m_magHistory[i]);
                if (diff > m_dbRise) val = val + 1;
            }
            m_magHistory[i] = sqrmag;
        }
        return val;
    }        

    double* DetectionFunction::getSpectrumMagnitude()
    {
        return m_magnitude;
    }
    //////////////////////////////////////////////////////////////////////
    // Construction/Destruction
    //////////////////////////////////////////////////////////////////////

    PeakPicking::PeakPicking( PPickParams Config )
    {
        m_workBuffer = NULL;
        initialise( Config );
    }

    PeakPicking::~PeakPicking()
    {
        deInitialise();
    }

    void PeakPicking::initialise( PPickParams Config )
    {
        m_DFLength = Config.length ;
        Qfilta = Config.QuadThresh.a ;
        Qfiltb = Config.QuadThresh.b ;
        Qfiltc = Config.QuadThresh.c ;
            
        m_DFProcessingParams.length = m_DFLength; 
        m_DFProcessingParams.LPOrd = Config.LPOrd; 
        m_DFProcessingParams.LPACoeffs = Config.LPACoeffs; 
        m_DFProcessingParams.LPBCoeffs = Config.LPBCoeffs; 
        m_DFProcessingParams.winPre  = Config.WinT.pre;
        m_DFProcessingParams.winPost = Config.WinT.post; 
        m_DFProcessingParams.AlphaNormParam = Config.alpha;
        m_DFProcessingParams.isMedianPositive = false;
        m_DFProcessingParams.delta = Config.delta; //add the delta threshold as an adjustable parameter

        m_DFSmoothing = new DFProcess( m_DFProcessingParams );

        m_workBuffer = new double[ m_DFLength ];
        memset( m_workBuffer, 0, sizeof(double)*m_DFLength);
    }

    void PeakPicking::deInitialise()
    {
        delete [] m_workBuffer;
        delete m_DFSmoothing;
        m_workBuffer = NULL;
    }

    void PeakPicking::process( double* src, int len, vector<int> &onsets )
    {
        if (len < 4) return;

        vector <double> m_maxima;   

        // Signal conditioning 
        m_DFSmoothing->process( src, m_workBuffer );
            
        for (int i = 0; i < len; i++) {
            m_maxima.push_back( m_workBuffer[ i ] );                
        }
            
        quadEval( m_maxima, onsets );

        for( int b = 0; b < (int)m_maxima.size(); b++) {
            src[ b ] = m_maxima[ b ];
        }
    }

    int PeakPicking::quadEval( vector<double> &src, vector<int> &idx )
    {
        int maxLength;

        vector <int> m_maxIndex;
        vector <int> m_onsetPosition;
            
        vector <double> m_maxFit;
        vector <double> m_poly;
        vector <double> m_err;

        m_poly.push_back(0);
        m_poly.push_back(0);
        m_poly.push_back(0);

        for (int t = -2; t < 3; t++) {
            m_err.push_back( (double)t );
        }

        for (int i = 2; i < int(src.size()) - 2; i++) {
            if ((src[i] > src[i-1]) && (src[i] > src[i+1]) && (src[i] > 0) ) {
                m_maxIndex.push_back(i);
            }
        }

        maxLength = int(m_maxIndex.size());

        double selMax = 0;

        for (int j = 0; j < maxLength ; j++) {
            for (int k = -2; k <= 2; ++k) {
                selMax = src[ m_maxIndex[j] + k ] ;
                m_maxFit.push_back(selMax);                 
            }

            TPolyFit::PolyFit2(m_err, m_maxFit, m_poly);

            double f = m_poly[0];
            double h = m_poly[2];

            if (h < -Qfilta || f > Qfiltc) {
                idx.push_back(m_maxIndex[j]);
            }
                    
            m_maxFit.clear();
        }

        return 1;
    }
}