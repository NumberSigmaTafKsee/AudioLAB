#pragma once

#include "qmdsp_math.hpp"
#include <stdexcept>
#include <climits>

namespace qmdsp
{
    class Filter
    {
    public:
        struct Parameters {
            std::vector<double> a;
            std::vector<double> b;
        };

        /**
        * Construct an IIR filter with numerators b and denominators
        * a. The filter will have order b.size()-1. To make an FIR
        * filter, leave the vector a in the param struct empty.
        * Otherwise, a and b must have the same number of values.
        */
        Filter(Parameters params);
        
        ~Filter();

        void reset();

        /**
        * Filter the input sequence \arg in of length \arg n samples, and
        * write the resulting \arg n samples into \arg out. There must be
        * enough room in \arg out for \arg n samples to be written.
        */
        void process(const double *const QM_R__ in,
                    double *const QM_R__ out,
                    const int n);

        int getOrder() const { return m_order; }
        
    private:
        int m_order;
        int m_sz;
        std::vector<double> m_a;
        std::vector<double> m_b;
        std::vector<double> m_bufa;
        std::vector<double> m_bufb;
        int m_offa;
        int m_offb;
        int m_offmax;
        bool m_fir;

        Filter(const Filter &); // not supplied
        Filter &operator=(const Filter &); // not supplied
    };
        

    /**
    * Zero-phase digital filter, implemented by processing the data
    * through a filter specified by the given filter parameters (see
    * Filter) and then processing it again in reverse.
    */
    class FiltFilt  
    {
    public:
        FiltFilt(Filter::Parameters);
        virtual ~FiltFilt();

        void process(const double *const QM_R__ src,
                    double *const QM_R__ dst,
                    const int length);

    private:
        Filter m_filter;
        int m_ord;
    };


    struct DFProcConfig
    {
        int length; 
        int LPOrd; 
        double *LPACoeffs; 
        double *LPBCoeffs; 
        int winPre;
        int winPost; 
        double AlphaNormParam;
        bool isMedianPositive;
        float delta; //delta threshold used as an offset when computing the smoothed detection function

        DFProcConfig() :
            length(0),
            LPOrd(0),
            LPACoeffs(NULL),
            LPBCoeffs(NULL),
            winPre(0),
            winPost(0),
            AlphaNormParam(0),
            isMedianPositive(false),
            delta(0)
        {
        }
    };

    class DFProcess  
    {
    public:
        DFProcess( DFProcConfig Config );
        virtual ~DFProcess();

        void process( double* src, double* dst );
            
    private:
        void initialise( DFProcConfig Config );
        void deInitialise();
        void removeDCNormalize( double *src, double*dst );
        void medianFilter( double* src, double* dst );

        int m_length;
        int m_FFOrd;

        int m_winPre;
        int m_winPost;

        double m_alphaNormParam;

        double* filtSrc;
        double* filtDst;

        double* m_filtScratchIn;
        double* m_filtScratchOut;

        FiltFilt* m_FiltFilt;

        bool m_isMedianPositive;
        float m_delta; //add delta threshold
    };


    class Framer  
    {
    public:
        Framer();
        virtual ~Framer();

        void setSource(double* src, int64_t length);
        void configure(int frameLength, int hop);
        
        int getMaxNoFrames();
        void getFrame(double* dst);

        void resetCounters();

    private:
        int64_t m_sampleLen;          // DataLength (samples)
        int m_framesRead;             // Read Frames Index
    
        double* m_srcBuffer;
        double* m_dataFrame;          // Analysis Frame Buffer
        double* m_strideFrame;        // Stride Frame Buffer
        int m_frameLength;            // Analysis Frame Length
        int m_stepSize;               // Analysis Frame Stride

        int m_maxFrames;

        int64_t m_srcIndex;
    };

    //////////////////////////////////////////////////////////////////////
    // Construction/Destruction
    //////////////////////////////////////////////////////////////////////

    DFProcess::DFProcess( DFProcConfig config )
    {
        filtSrc = NULL;
        filtDst = NULL;     
        m_filtScratchIn = NULL;
        m_filtScratchOut = NULL;

        m_FFOrd = 0;

        initialise( config );
    }

    DFProcess::~DFProcess()
    {
        deInitialise();
    }

    void DFProcess::initialise( DFProcConfig config )
    {
        m_length = config.length;
        m_winPre = config.winPre;
        m_winPost = config.winPost;
        m_alphaNormParam = config.AlphaNormParam;

        m_isMedianPositive = config.isMedianPositive;

        filtSrc = new double[ m_length ];
        filtDst = new double[ m_length ];

        Filter::Parameters params;
        params.a = std::vector<double>
            (config.LPACoeffs, config.LPACoeffs + config.LPOrd + 1);
        params.b = std::vector<double>
            (config.LPBCoeffs, config.LPBCoeffs + config.LPOrd + 1);
        
        m_FiltFilt = new FiltFilt(params);
            
        //add delta threshold
        m_delta = config.delta;
    }

    void DFProcess::deInitialise()
    {
        delete [] filtSrc;
        delete [] filtDst;
        delete [] m_filtScratchIn;
        delete [] m_filtScratchOut;
        delete m_FiltFilt;
    }

    void DFProcess::process(double *src, double* dst)
    {
        if (m_length == 0) return;

        removeDCNormalize( src, filtSrc );

        m_FiltFilt->process( filtSrc, filtDst, m_length );

        medianFilter( filtDst, dst );
    }


    void DFProcess::medianFilter(double *src, double *dst)
    {
        int i,k,j,l;
        int index = 0;

        double val = 0;

        double* y = new double[ m_winPost + m_winPre + 1];
        memset( y, 0, sizeof( double ) * ( m_winPost + m_winPre + 1) );

        double* scratch = new double[ m_length ];

        for( i = 0; i < m_winPre; i++) {
            
            if (index >= m_length) {
                break;
            }

            k = i + m_winPost + 1;

            for( j = 0; j < k; j++) {
                y[ j ] = src[ j ];
            }
            scratch[ index ] = MathUtilities::median( y, k );
            index++;
        }

        for(  i = 0; i + m_winPost + m_winPre < m_length; i ++) {
            
            if (index >= m_length) {
                break;
            }
                            
            l = 0;
            for(  j  = i; j < ( i + m_winPost + m_winPre + 1); j++) {
                y[ l ] = src[ j ];
                l++;
            }

            scratch[index] = MathUtilities::median( y, (m_winPost + m_winPre + 1 ));
            index++;
        }

        for( i = std::max( m_length - m_winPost, 1); i < m_length; i++) {
            
            if (index >= m_length) {
                break;
            }

            k = std::max( i - m_winPre, 1);

            l = 0;
            for( j = k; j < m_length; j++) {
                y[ l ] = src[ j ];
                l++;
            }
                    
            scratch[index] = MathUtilities::median( y, l);
            index++;
        }

        for( i = 0; i < m_length; i++ ) {
            //add a delta threshold used as an offset when computing the smoothed detection function
            //(helps to discard noise when detecting peaks) 
            val = src[ i ] - scratch[ i ] - m_delta;
                    
            if( m_isMedianPositive ) {
                if( val > 0 ) {
                    dst[ i ]  = val;
                } else {
                    dst[ i ]  = 0;
                }
            } else {
                dst[ i ]  = val;
            }
        }
            
        delete [] y;
        delete [] scratch;
    }


    void DFProcess::removeDCNormalize( double *src, double*dst )
    {
        double DFmax = 0;
        double DFMin = 0;
        double DFAlphaNorm = 0;

        MathUtilities::getFrameMinMax( src, m_length, &DFMin, &DFmax );

        MathUtilities::getAlphaNorm( src, m_length, m_alphaNormParam, &DFAlphaNorm );

        for (int i = 0; i < m_length; i++) {
            dst[ i ] = ( src[ i ] - DFMin ) / DFAlphaNorm; 
        }
    }


    FiltFilt::FiltFilt(Filter::Parameters parameters) :
        m_filter(parameters)
    {
        m_ord = m_filter.getOrder();
    }

    FiltFilt::~FiltFilt()
    {
    }

    void FiltFilt::process(const double *const QM_R__ src,
                        double *const QM_R__ dst,
                        const int length)
    {       
        int i;

        if (length == 0) return;

        int nFilt = m_ord + 1;
        int nFact = 3 * (nFilt - 1);
        int nExt = length + 2 * nFact;

        double *filtScratchIn = new double[ nExt ];
        double *filtScratchOut = new double[ nExt ];
            
        for (i = 0; i < nExt; i++) {
            filtScratchIn[ i ] = 0.0;
            filtScratchOut[ i ] = 0.0;
        }

        // Edge transients reflection
        double sample0 = 2 * src[ 0 ];
        double sampleN = 2 * src[ length - 1 ];

        int index = 0;
        for (i = nFact; i > 0; i--) {
            if (i < length) {
                filtScratchIn[index] = sample0 - src[ i ];
            }
            ++index;
        }
        index = 0;
        for (i = 0; i < nFact; i++) {
            if (i + 1 < length) {
                filtScratchIn[(nExt - nFact) + index] =
                    sampleN - src[ (length - 2) - i ];
            }
            ++index;
        }

        for (i = 0; i < length; i++) {
            filtScratchIn[ i + nFact ] = src[ i ];
        }
        
        ////////////////////////////////
        // Do 0Ph filtering
        m_filter.process(filtScratchIn, filtScratchOut, nExt);
            
        // reverse the series for FILTFILT 
        for (i = 0; i < nExt; i++) { 
            filtScratchIn[ i ] = filtScratchOut[ nExt - i - 1];
        }

        // clear filter state
        m_filter.reset();
        
        // do FILTER again 
        m_filter.process(filtScratchIn, filtScratchOut, nExt);

        // reverse the series to output
        for (i = 0; i < length; i++) {
            dst[ i ] = filtScratchOut[ nExt - nFact - i - 1 ];
        }

        delete [] filtScratchIn;
        delete [] filtScratchOut;
    }


    Filter::Filter(Parameters params)
    {
        if (params.a.empty()) {
            m_fir = true;
            if (params.b.empty()) {
                throw logic_error("Filter must have at least one pair of coefficients");
            }
        } else {
            m_fir = false;
            if (params.a.size() != params.b.size()) {
                throw logic_error("Inconsistent numbers of filter coefficients");
            }
        }
        
        m_sz = int(params.b.size());
        m_order = m_sz - 1;

        m_a = params.a;
        m_b = params.b;
        
        // We keep some empty space at the start of the buffer, and
        // encroach gradually into it as we add individual sample
        // calculations at the start. Then when we run out of space, we
        // move the buffer back to the end and begin again. This is
        // significantly faster than moving the whole buffer along in
        // 1-sample steps every time.

        m_offmax = 20;
        m_offa = m_offmax;
        m_offb = m_offmax;

        if (!m_fir) {
            m_bufa.resize(m_order + m_offmax);
        }

        m_bufb.resize(m_sz + m_offmax);
    }

    Filter::~Filter()
    {
    }

    void
    Filter::reset()
    {
        m_offb = m_offmax;
        m_offa = m_offmax;

        if (!m_fir) {
            m_bufa.assign(m_bufa.size(), 0.0);
        }

        m_bufb.assign(m_bufb.size(), 0.0);
    }

    void
    Filter::process(const double *const QM_R__ in,
                    double *const QM_R__ out,
                    const int n)
    {
        for (int s = 0; s < n; ++s) {

            if (m_offb > 0) {
                --m_offb;
            } else {
                for (int i = m_sz - 2; i >= 0; --i) {
                    m_bufb[i + m_offmax + 1] = m_bufb[i];
                }
                m_offb = m_offmax;
            }
            m_bufb[m_offb] = in[s];

            double b_sum = 0.0;
            for (int i = 0; i < m_sz; ++i) {
                b_sum += m_b[i] * m_bufb[i + m_offb];
            }

            double outval;

            if (m_fir) {

                outval = b_sum;

            } else {

                double a_sum = 0.0;
                for (int i = 0; i < m_order; ++i) {
                    a_sum += m_a[i + 1] * m_bufa[i + m_offa];
                }

                outval = b_sum - a_sum;

                if (m_offa > 0) {
                    --m_offa;
                } else {
                    for (int i = m_order - 2; i >= 0; --i) {
                        m_bufa[i + m_offmax + 1] = m_bufa[i];
                    }
                    m_offa = m_offmax;
                }
                m_bufa[m_offa] = outval;
            }
            
            out[s] = outval;
        }
    }


    Framer::Framer() :
        m_sampleLen(0),
        m_framesRead(0),
        m_srcBuffer(0),
        m_dataFrame(0),
        m_strideFrame(0),
        m_frameLength(0),
        m_stepSize(0),
        m_maxFrames(0),
        m_srcIndex(0)
    {
    }

    Framer::~Framer()
    {
        delete[] m_dataFrame;
        delete[] m_strideFrame;
    }

    void Framer::configure(int frameLength, int hop)
    {
        m_frameLength = frameLength;
        m_stepSize = hop;

        resetCounters();

        delete[] m_dataFrame;  
        m_dataFrame = new double[ m_frameLength ];

        delete [] m_strideFrame;        
        m_strideFrame = new double[ m_stepSize ];
    }

    void Framer::getFrame(double *dst)
    {
        if ((m_srcIndex + int64_t(m_frameLength)) < m_sampleLen) {

            for (int i = 0; i < m_frameLength; i++) {
                dst[i] = m_srcBuffer[m_srcIndex++]; 
            }
            m_srcIndex -= (m_frameLength - m_stepSize);

        } else { // m_srcIndex is within m_frameLength of m_sampleLen

            int rem = int(m_sampleLen - m_srcIndex);
            int zero = m_frameLength - rem;

            for (int i = 0; i < rem; i++) {
                dst[i] = m_srcBuffer[m_srcIndex++];
            }
                    
            for (int i = 0; i < zero; i++ ) {
                dst[rem + i] = 0.0;
            }

            m_srcIndex -= (rem - m_stepSize);
        }

        m_framesRead++;
    }

    void Framer::resetCounters()
    {
        m_framesRead = 0;
        m_srcIndex = 0;
    }

    int Framer::getMaxNoFrames()
    {
        return m_maxFrames;
    }

    void Framer::setSource(double *src, int64_t length)
    {
        m_srcBuffer = src;
        m_sampleLen = length;

        int64_t maxFrames = length / int64_t(m_stepSize);
        if (maxFrames * int64_t(m_stepSize) < length) {
            ++maxFrames;
        }
        if (maxFrames > INT_MAX) maxFrames = INT_MAX;
        m_maxFrames = maxFrames;
    }    
}