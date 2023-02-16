#pragma once

#include "qmdsp_math.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm>

namespace qmdsp
{
    /**
    * Decimator carries out a fast downsample by a power-of-two
    * factor. Only a limited number of factors are supported, from two to
    * whatever getHighestSupportedFactor() returns. This is much faster
    * than Resampler but has a worse signal-noise ratio.
    */
    class Decimator  
    {
    public:
        /**
        * Construct a Decimator to operate on input blocks of length
        * inLength, with decimation factor decFactor.  inLength should be
        * a multiple of decFactor.  Output blocks will be of length
        * inLength / decFactor.
        *
        * decFactor must be a power of two.  The highest supported factor
        * is obtained through getHighestSupportedFactor(); for higher
        * factors, you will need to chain more than one decimator.
        */
        Decimator(int inLength, int decFactor);
        virtual ~Decimator();

        /**
        * Process inLength samples (as supplied to constructor) from src
        * and write inLength / decFactor samples to dst.  Note that src
        * and dst may be the same or overlap (an intermediate buffer is
        * used).
        */
        void process( const double* src, double* dst );

        /**
        * Process inLength samples (as supplied to constructor) from src
        * and write inLength / decFactor samples to dst.  Note that src
        * and dst may be the same or overlap (an intermediate buffer is
        * used).
        */
        void process( const float* src, float* dst );

        int getFactor() const { return m_decFactor; }
        static int getHighestSupportedFactor() { return 8; }

        void resetFilter();

    private:
        void deInitialise();
        void initialise( int inLength, int decFactor );
        void doAntiAlias( const double* src, double* dst, int length );
        void doAntiAlias( const float* src, double* dst, int length );

        int m_inputLength;
        int m_outputLength;
        int m_decFactor;

        double Input;
        double Output ;

        double o1,o2,o3,o4,o5,o6,o7;

        double a[ 9 ];
        double b[ 9 ];
            
        double* decBuffer;
    };

    /**
    * DecimatorB carries out a fast downsample by a power-of-two
    * factor. It only knows how to decimate by a factor of 2, and will
    * use repeated decimation for higher factors. A Butterworth filter of
    * order 6 is used for the lowpass filter.
    */
    class DecimatorB
    {
    public:
        void process( const double* src, double* dst );
        void process( const float* src, float* dst );

        /**
        * Construct a DecimatorB to operate on input blocks of length
        * inLength, with decimation factor decFactor.  inLength should be
        * a multiple of decFactor.  Output blocks will be of length
        * inLength / decFactor.
        *
        * decFactor must be a power of two.
        */
        DecimatorB(int inLength, int decFactor);
        virtual ~DecimatorB();

        int getFactor() const { return m_decFactor; }

    private:
        void deInitialise();
        void initialise(int inLength, int decFactor);
        void doAntiAlias(const double* src, double* dst, int length, int filteridx);
        void doProcess();

        int m_inputLength;
        int m_outputLength;
        int m_decFactor;

        std::vector<std::vector<double> > m_o;

        double m_a[7];
        double m_b[7];
            
        double *m_aaBuffer;
        double *m_tmpBuffer;
    };


    /**
    * Resampler resamples a stream from one integer sample rate to
    * another (arbitrary) rate, using a kaiser-windowed sinc filter.  The
    * results and performance are pretty similar to libraries such as
    * libsamplerate, though this implementation does not support
    * time-varying ratios (the ratio is fixed on construction).
    *
    * See also Decimator, which is faster and rougher but supports only
    * power-of-two downsampling factors.
    */
    class Resampler
    {
    public:
        /**
        * Construct a Resampler to resample from sourceRate to
        * targetRate.
        */
        Resampler(int sourceRate, int targetRate);

        /**
        * Construct a Resampler to resample from sourceRate to
        * targetRate, using the given filter parameters.
        */
        Resampler(int sourceRate, int targetRate,
                double snr, double bandwidth);

        virtual ~Resampler();

        /**
        * Read n input samples from src and write resampled data to
        * dst. The return value is the number of samples written, which
        * will be no more than ceil((n * targetRate) / sourceRate). The
        * caller must ensure the dst buffer has enough space for the
        * samples returned.
        */
        int process(const double *src, double *dst, int n);

        /**
        * Read n input samples from src and return resampled data by
        * value.
        */
        std::vector<double> process(const double *src, int n);

        /**
        * Return the number of samples of latency at the output due by
        * the filter. (That is, the output will be delayed by this number
        * of samples relative to the input.)
        */
        int getLatency() const { return m_latency; }

        /**
        * Carry out a one-off resample of a single block of n
        * samples. The output is latency-compensated.
        */
        static std::vector<double> resample
        (int sourceRate, int targetRate, const double *data, int n);

    private:
        int m_sourceRate;
        int m_targetRate;
        int m_gcd;
        int m_filterLength;
        int m_latency;
        double m_peakToPole;
        
        struct Phase {
            int nextPhase;
            std::vector<double> filter;
            int drop;
        };

        Phase *m_phaseData;
        int m_phase;
        std::vector<double> m_buffer;
        int m_bufferOrigin;

        void initialise(double, double);
        double reconstructOne();
    };


    //////////////////////////////////////////////////////////////////////
    // Construction/Destruction
    //////////////////////////////////////////////////////////////////////

    Decimator::Decimator( int inLength, int decFactor )
    {
        m_inputLength = 0;
        m_outputLength = 0;
        m_decFactor = 1;

        initialise( inLength, decFactor );
    }

    Decimator::~Decimator()
    {
        deInitialise();
    }

    void Decimator::initialise( int inLength, int decFactor)
    {
        m_inputLength = inLength;
        m_decFactor = decFactor;
        m_outputLength = m_inputLength / m_decFactor;

        decBuffer = new double[ m_inputLength ];

        // If adding new factors here, add them to
        // getHighestSupportedFactor in the header as well

        if(m_decFactor == 8) {

            //////////////////////////////////////////////////
            b[0] = 0.060111378492136;
            b[1] = -0.257323420830598;
            b[2] = 0.420583503165928;
            b[3] = -0.222750785197418;
            b[4] = -0.222750785197418;
            b[5] = 0.420583503165928;
            b[6] = -0.257323420830598;
            b[7] = 0.060111378492136;

            a[0] = 1;
            a[1] = -5.667654878577432;
            a[2] = 14.062452278088417;
            a[3] = -19.737303840697738;
            a[4] = 16.889698874608641;
            a[5] = -8.796600612325928;
            a[6] = 2.577553446979888;
            a[7] = -0.326903916815751;
            //////////////////////////////////////////////////

        } else if( m_decFactor == 4 ) {
            
            //////////////////////////////////////////////////
            b[ 0 ] = 0.10133306904918619;
            b[ 1 ] = -0.2447523353702363;
            b[ 2 ] = 0.33622528590120965;
            b[ 3 ] = -0.13936581560633518;
            b[ 4 ] = -0.13936581560633382;
            b[ 5 ] = 0.3362252859012087;
            b[ 6 ] = -0.2447523353702358;
            b[ 7 ] = 0.10133306904918594;

            a[ 0 ] = 1;
            a[ 1 ] = -3.9035590278139427;
            a[ 2 ] = 7.5299379980621133;
            a[ 3 ] = -8.6890803793177511;
            a[ 4 ] = 6.4578667096099176;
            a[ 5 ] = -3.0242979431223631;
            a[ 6 ] = 0.83043385136748382;
            a[ 7 ] = -0.094420800837809335;
            //////////////////////////////////////////////////

        } else if( m_decFactor == 2 ) {
            
            //////////////////////////////////////////////////
            b[ 0 ] = 0.20898944260075727;
            b[ 1 ] = 0.40011234879814367;
            b[ 2 ] = 0.819741973072733;
            b[ 3 ] = 1.0087419911682323;
            b[ 4 ] = 1.0087419911682325;
            b[ 5 ] = 0.81974197307273156;
            b[ 6 ] = 0.40011234879814295;
            b[ 7 ] = 0.20898944260075661;

            a[ 0 ] = 1;
            a[ 1 ] = 0.0077331184208358217;
            a[ 2 ] = 1.9853971155964376;
            a[ 3 ] = 0.19296739275341004;
            a[ 4 ] = 1.2330748872852182;
            a[ 5 ] = 0.18705341389316466;
            a[ 6 ] = 0.23659265908013868;
            a[ 7 ] = 0.032352924250533946;

        } else {
            
            if ( m_decFactor != 1 ) {
                std::cerr << "WARNING: Decimator::initialise: unsupported decimation factor " << m_decFactor << ", no antialiasing filter will be used" << std::endl;
            }

            //////////////////////////////////////////////////
            b[ 0 ] = 1;
            b[ 1 ] = 0;
            b[ 2 ] = 0;
            b[ 3 ] = 0;
            b[ 4 ] = 0;
            b[ 5 ] = 0;
            b[ 6 ] = 0;
            b[ 7 ] = 0;

            a[ 0 ] = 1;
            a[ 1 ] = 0;
            a[ 2 ] = 0;
            a[ 3 ] = 0;
            a[ 4 ] = 0;
            a[ 5 ] = 0;
            a[ 6 ] = 0;
            a[ 7 ] = 0;
        }

        resetFilter();
    }

    void Decimator::deInitialise()
    {
        delete [] decBuffer;
    }

    void Decimator::resetFilter()
    {
        Input = Output = 0;

        o1=o2=o3=o4=o5=o6=o7=0;
    }

    void Decimator::doAntiAlias(const double *src, double *dst, int length)
    {
        for (int i = 0; i < length; i++ ) {
            
            Input = (double)src[ i ];

            Output = Input * b[ 0 ] + o1;

            o1 = Input * b[ 1 ] - Output * a[ 1 ] + o2;
            o2 = Input * b[ 2 ] - Output * a[ 2 ] + o3;
            o3 = Input * b[ 3 ] - Output * a[ 3 ] + o4;
            o4 = Input * b[ 4 ] - Output * a[ 4 ] + o5;
            o5 = Input * b[ 5 ] - Output * a[ 5 ] + o6;
            o6 = Input * b[ 6 ] - Output * a[ 6 ] + o7;
            o7 = Input * b[ 7 ] - Output * a[ 7 ] ;

            dst[ i ] = Output;
        }
    }

    void Decimator::doAntiAlias(const float *src, double *dst, int length)
    {
        for (int i = 0; i < length; i++ ) {
            
            Input = (double)src[ i ];

            Output = Input * b[ 0 ] + o1;

            o1 = Input * b[ 1 ] - Output * a[ 1 ] + o2;
            o2 = Input * b[ 2 ] - Output * a[ 2 ] + o3;
            o3 = Input * b[ 3 ] - Output * a[ 3 ] + o4;
            o4 = Input * b[ 4 ] - Output * a[ 4 ] + o5;
            o5 = Input * b[ 5 ] - Output * a[ 5 ] + o6;
            o6 = Input * b[ 6 ] - Output * a[ 6 ] + o7;
            o7 = Input * b[ 7 ] - Output * a[ 7 ] ;

            dst[ i ] = Output;
        }
    }

    void Decimator::process(const double *src, double *dst)
    {
        if (m_decFactor == 1) {
            for (int i = 0; i < m_outputLength; i++ ) {
                dst[i] = src[i];
            }
            return;
        }
            
        doAntiAlias( src, decBuffer, m_inputLength );

        int idx = 0;

        for (int i = 0; i < m_outputLength; i++ ) {
            dst[ idx++ ] = decBuffer[ m_decFactor * i ];
        }
    }

    void Decimator::process(const float *src, float *dst)
    {
        if (m_decFactor == 1) {
            for (int i = 0; i < m_outputLength; i++ ) {
                dst[i] = src[i];
            }
            return;
        }

        doAntiAlias( src, decBuffer, m_inputLength );

        int idx = 0;

        for (int i = 0; i < m_outputLength; i++ ) {
            dst[ idx++ ] = decBuffer[ m_decFactor * i ];
        }
    }


    DecimatorB::DecimatorB(int inLength, int decFactor)
    {
        m_inputLength = 0;
        m_outputLength = 0;
        m_decFactor = 1;
        m_aaBuffer = 0;
        m_tmpBuffer = 0;

        initialise(inLength, decFactor);
    }

    DecimatorB::~DecimatorB()
    {
        deInitialise();
    }

    void DecimatorB::initialise(int inLength, int decFactor)
    {
        m_inputLength = inLength;
        m_decFactor = decFactor;
        m_outputLength = m_inputLength / m_decFactor;

        if (m_decFactor < 2 || !MathUtilities::isPowerOfTwo(m_decFactor)) {
            std::cerr << "ERROR: DecimatorB::initialise: Decimation factor must be a power of 2 and at least 2 (was: " << m_decFactor << ")" << std::endl;
            m_decFactor = 0;
            return;
        }

        if (m_inputLength % m_decFactor != 0) {
            std::cerr << "ERROR: DecimatorB::initialise: inLength must be a multiple of decimation factor (was: " << m_inputLength << ", factor is " << m_decFactor << ")" << std::endl;
            m_decFactor = 0;
            return;
        }        

        m_aaBuffer = new double[m_inputLength];
        m_tmpBuffer = new double[m_inputLength];

        // Order 6 Butterworth lowpass filter
        // Calculated using e.g. MATLAB butter(6, 0.5, 'low')

        m_b[0] = 0.029588223638661;
        m_b[1] = 0.177529341831965;
        m_b[2] = 0.443823354579912;
        m_b[3] = 0.591764472773216;
        m_b[4] = 0.443823354579912;
        m_b[5] = 0.177529341831965;
        m_b[6] = 0.029588223638661;

        m_a[0] = 1.000000000000000;
        m_a[1] = 0.000000000000000;
        m_a[2] = 0.777695961855673;
        m_a[3] = 0.000000000000000;
        m_a[4] = 0.114199425062434;
        m_a[5] = 0.000000000000000;
        m_a[6] = 0.001750925956183;

        for (int factor = m_decFactor; factor > 1; factor /= 2) {
            m_o.push_back(vector<double>(6, 0.0));
        }
    }

    void DecimatorB::deInitialise()
    {
        delete [] m_aaBuffer;
        delete [] m_tmpBuffer;
    }

    void DecimatorB::doAntiAlias(const double *src, double *dst, int length,
                                int filteridx)
    {
        std::vector<double> &o = m_o[filteridx];

        for (int i = 0; i < length; i++) {

            double input = src[i];
            double output = input * m_b[0] + o[0];

            o[0] = input * m_b[1] - output * m_a[1] + o[1];
            o[1] = input * m_b[2] - output * m_a[2] + o[2];
            o[2] = input * m_b[3] - output * m_a[3] + o[3];
            o[3] = input * m_b[4] - output * m_a[4] + o[4];
            o[4] = input * m_b[5] - output * m_a[5] + o[5];
            o[5] = input * m_b[6] - output * m_a[6];

            dst[i] = output;
        }
    }

    void DecimatorB::doProcess()
    {
        int filteridx = 0;
        int factorDone = 1;

        while (factorDone < m_decFactor) {

            doAntiAlias(m_tmpBuffer, m_aaBuffer,
                        m_inputLength / factorDone,
                        filteridx);

            filteridx ++;
            factorDone *= 2;

            for (int i = 0; i < m_inputLength / factorDone; ++i) {
                m_tmpBuffer[i] = m_aaBuffer[i * 2];
            }
        }
    }

    void DecimatorB::process(const double *src, double *dst)
    {
        if (m_decFactor == 0) return;

        for (int i = 0; i < m_inputLength; ++i) {
            m_tmpBuffer[i] = src[i];
        }

        doProcess();
        
        for (int i = 0; i < m_outputLength; ++i) {
            dst[i] = m_tmpBuffer[i];
        }
    }

    void DecimatorB::process(const float *src, float *dst)
    {
        if (m_decFactor == 0) return;

        for (int i = 0; i < m_inputLength; ++i) {
            m_tmpBuffer[i] = src[i];
        }

        doProcess();
        
        for (int i = 0; i < m_outputLength; ++i) {
            dst[i] = m_tmpBuffer[i];
        }
    }    


    Resampler::Resampler(int sourceRate, int targetRate) :
        m_sourceRate(sourceRate),
        m_targetRate(targetRate)
    {
    #ifdef DEBUG_RESAMPLER
        std::cerr << "Resampler::Resampler(" <<  sourceRate << "," << targetRate << ")" << std::endl;
    #endif
        initialise(100, 0.02);
    }

    Resampler::Resampler(int sourceRate, int targetRate, 
                        double snr, double bandwidth) :
        m_sourceRate(sourceRate),
        m_targetRate(targetRate)
    {
        initialise(snr, bandwidth);
    }

    Resampler::~Resampler()
    {
        delete[] m_phaseData;
    }

    void
    Resampler::initialise(double snr, double bandwidth)
    {
        int higher = std::max(m_sourceRate, m_targetRate);
        int lower = std::min(m_sourceRate, m_targetRate);

        m_gcd = MathUtilities::gcd(lower, higher);
        m_peakToPole = higher / m_gcd;

        if (m_targetRate < m_sourceRate) {
            // antialiasing filter, should be slightly below nyquist
            m_peakToPole = m_peakToPole / (1.0 - bandwidth/2.0);
        }

        KaiserWindow::Parameters params =
            KaiserWindow::parametersForBandwidth(snr, bandwidth, higher / m_gcd);

        params.length =
            (params.length % 2 == 0 ? params.length + 1 : params.length);
        
        params.length =
            (params.length > 200001 ? 200001 : params.length);

        m_filterLength = params.length;

        std::vector<double> filter;

        KaiserWindow kw(params);
        SincWindow sw(m_filterLength, m_peakToPole * 2);

        filter = std::vector<double>(m_filterLength, 0.0);
        for (int i = 0; i < m_filterLength; ++i) filter[i] = 1.0;
        sw.cut(filter.data());
        kw.cut(filter.data());
        
        int inputSpacing = m_targetRate / m_gcd;
        int outputSpacing = m_sourceRate / m_gcd;

    #ifdef DEBUG_RESAMPLER
        std::cerr << "resample " << m_sourceRate << " -> " << m_targetRate
            << ": inputSpacing " << inputSpacing << ", outputSpacing "
            << outputSpacing << ": filter length " << m_filterLength
            << std::endl;
    #endif

        // Now we have a filter of (odd) length flen in which the lower
        // sample rate corresponds to every n'th point and the higher rate
        // to every m'th where n and m are higher and lower rates divided
        // by their gcd respectively. So if x coordinates are on the same
        // scale as our filter resolution, then source sample i is at i *
        // (targetRate / gcd) and target sample j is at j * (sourceRate /
        // gcd).

        // To reconstruct a single target sample, we want a buffer (real
        // or virtual) of flen values formed of source samples spaced at
        // intervals of (targetRate / gcd), in our example case 3.  This
        // is initially formed with the first sample at the filter peak.
        //
        // 0  0  0  0  a  0  0  b  0
        //
        // and of course we have our filter
        //
        // f1 f2 f3 f4 f5 f6 f7 f8 f9
        //
        // We take the sum of products of non-zero values from this buffer
        // with corresponding values in the filter
        //
        // a * f5 + b * f8
        //
        // Then we drop (sourceRate / gcd) values, in our example case 4,
        // from the start of the buffer and fill until it has flen values
        // again
        //
        // a  0  0  b  0  0  c  0  0
        //
        // repeat to reconstruct the next target sample
        //
        // a * f1 + b * f4 + c * f7
        //
        // and so on.
        //
        // Above I said the buffer could be "real or virtual" -- ours is
        // virtual. We don't actually store all the zero spacing values,
        // except for padding at the start; normally we store only the
        // values that actually came from the source stream, along with a
        // phase value that tells us how many virtual zeroes there are at
        // the start of the virtual buffer.  So the two examples above are
        //
        // 0 a b  [ with phase 1 ]
        // a b c  [ with phase 0 ]
        //
        // Having thus broken down the buffer so that only the elements we
        // need to multiply are present, we can also unzip the filter into
        // every-nth-element subsets at each phase, allowing us to do the
        // filter multiplication as a simply vector multiply. That is, rather
        // than store
        //
        // f1 f2 f3 f4 f5 f6 f7 f8 f9
        // 
        // we store separately
        //
        // f1 f4 f7
        // f2 f5 f8
        // f3 f6 f9
        //
        // Each time we complete a multiply-and-sum, we need to work out
        // how many (real) samples to drop from the start of our buffer,
        // and how many to add at the end of it for the next multiply.  We
        // know we want to drop enough real samples to move along by one
        // computed output sample, which is our outputSpacing number of
        // virtual buffer samples. Depending on the relationship between
        // input and output spacings, this may mean dropping several real
        // samples, one real sample, or none at all (and simply moving to
        // a different "phase").

        m_phaseData = new Phase[inputSpacing];

        for (int phase = 0; phase < inputSpacing; ++phase) {

            Phase p;

            p.nextPhase = phase - outputSpacing;
            while (p.nextPhase < 0) p.nextPhase += inputSpacing;
            p.nextPhase %= inputSpacing;
            
            p.drop = int(std::ceil(std::max(0.0, double(outputSpacing - phase))
                            / inputSpacing));

            int filtZipLength = int(std::ceil(double(m_filterLength - phase)
                                        / inputSpacing));

            for (int i = 0; i < filtZipLength; ++i) {
                p.filter.push_back(filter[i * inputSpacing + phase]);
            }

            m_phaseData[phase] = p;
        }

    #ifdef DEBUG_RESAMPLER
        int cp = 0;
        int totDrop = 0;
        for (int i = 0; i < inputSpacing; ++i) {
            std::cerr << "phase = " << cp << ", drop = " << m_phaseData[cp].drop
                << ", filter length = " << m_phaseData[cp].filter.size()
                << ", next phase = " << m_phaseData[cp].nextPhase << std::endl;
            totDrop += m_phaseData[cp].drop;
            cp = m_phaseData[cp].nextPhase;
        }
        std::cerr << "total drop = " << totDrop << std::endl;
    #endif

        // The May implementation of this uses a pull model -- we ask the
        // resampler for a certain number of output samples, and it asks
        // its source stream for as many as it needs to calculate
        // those. This means (among other things) that the source stream
        // can be asked for enough samples up-front to fill the buffer
        // before the first output sample is generated.
        // 
        // In this implementation we're using a push model in which a
        // certain number of source samples is provided and we're asked
        // for as many output samples as that makes available. But we
        // can't return any samples from the beginning until half the
        // filter length has been provided as input. This means we must
        // either return a very variable number of samples (none at all
        // until the filter fills, then half the filter length at once) or
        // else have a lengthy declared latency on the output. We do the
        // latter. (What do other implementations do?)
        //
        // We want to make sure the first "real" sample will eventually be
        // aligned with the centre sample in the filter (it's tidier, and
        // easier to do diagnostic calculations that way). So we need to
        // pick the initial phase and buffer fill accordingly.
        // 
        // Example: if the inputSpacing is 2, outputSpacing is 3, and
        // filter length is 7,
        // 
        //    x     x     x     x     a     b     c ... input samples
        // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 ... 
        //          i        j        k        l    ... output samples
        // [--------|--------] <- filter with centre mark
        //
        // Let h be the index of the centre mark, here 3 (generally
        // int(filterLength/2) for odd-length filters).
        //
        // The smallest n such that h + n * outputSpacing > filterLength
        // is 2 (that is, ceil((filterLength - h) / outputSpacing)), and
        // (h + 2 * outputSpacing) % inputSpacing == 1, so the initial
        // phase is 1.
        //
        // To achieve our n, we need to pre-fill the "virtual" buffer with
        // 4 zero samples: the x's above. This is int((h + n *
        // outputSpacing) / inputSpacing). It's the phase that makes this
        // buffer get dealt with in such a way as to give us an effective
        // index for sample a of 9 rather than 8 or 10 or whatever.
        //
        // This gives us output latency of 2 (== n), i.e. output samples i
        // and j will appear before the one in which input sample a is at
        // the centre of the filter.

        int h = int(m_filterLength / 2);
        int n = std::ceil(double(m_filterLength - h) / outputSpacing);
        
        m_phase = (h + n * outputSpacing) % inputSpacing;

        int fill = (h + n * outputSpacing) / inputSpacing;
        
        m_latency = n;

        m_buffer = std::std::vector<double>(fill, 0);
        m_bufferOrigin = 0;

    #ifdef DEBUG_RESAMPLER
        std::cerr << "initial phase " << m_phase << " (as " << (m_filterLength/2) << " % " << inputSpacing << ")"
                << ", latency " << m_latency << std::endl;
    #endif
    }

    double
    Resampler::reconstructOne()
    {
        Phase &pd = m_phaseData[m_phase];
        double v = 0.0;
        int n = pd.filter.size();

        if (n + m_bufferOrigin > (int)m_buffer.size()) {
            cerr << "ERROR: n + m_bufferOrigin > m_buffer.size() [" << n << " + "
                << m_bufferOrigin << " > " << m_buffer.size() << "]" << endl;
            throw std::logic_error("n + m_bufferOrigin > m_buffer.size()");
        }

        const double *const QM_R__ buf(m_buffer.data() + m_bufferOrigin);
        const double *const QM_R__ filt(pd.filter.data());

        for (int i = 0; i < n; ++i) {
            // NB gcc can only vectorize this with -ffast-math
            v += buf[i] * filt[i];
        }

        m_bufferOrigin += pd.drop;
        m_phase = pd.nextPhase;
        return v;
    }

    int
    Resampler::process(const double *src, double *dst, int n)
    {
        m_buffer.insert(m_buffer.end(), src, src + n);

        int maxout = int(ceil(double(n) * m_targetRate / m_sourceRate));
        int outidx = 0;

    #ifdef DEBUG_RESAMPLER
        cerr << "process: buf siz " << m_buffer.size() << " filt siz for phase " << m_phase << " " << m_phaseData[m_phase].filter.size() << endl;
    #endif

        double scaleFactor = (double(m_targetRate) / m_gcd) / m_peakToPole;

        while (outidx < maxout &&
            m_buffer.size() >= m_phaseData[m_phase].filter.size() + m_bufferOrigin) {
            dst[outidx] = scaleFactor * reconstructOne();
            outidx++;
        }

        if (m_bufferOrigin > (int)m_buffer.size()) {
            cerr << "ERROR: m_bufferOrigin > m_buffer.size() [" 
                << m_bufferOrigin << " > " << m_buffer.size() << "]" << endl;
            throw std::logic_error("m_bufferOrigin > m_buffer.size()");
        }

        m_buffer = std::vector<double>(m_buffer.begin() + m_bufferOrigin, m_buffer.end());
        m_bufferOrigin = 0;
        
        return outidx;
    }
        
    std::vector<double>
    Resampler::process(const double *src, int n)
    {
        int maxout = int(std::ceil(double(n) * m_targetRate / m_sourceRate));
        std::vector<double> out(maxout, 0.0);
        int got = process(src, out.data(), n);
        assert(got <= maxout);
        if (got < maxout) out.resize(got);
        return out;
    }

    std::vector<double>
    Resampler::resample(int sourceRate, int targetRate, const double *data, int n)
    {
        Resampler r(sourceRate, targetRate);

        int latency = r.getLatency();

        // latency is the output latency. We need to provide enough
        // padding input samples at the end of input to guarantee at
        // *least* the latency's worth of output samples. that is,

        int inputPad = int(std::ceil((double(latency) * sourceRate) / targetRate));

        // that means we are providing this much input in total:

        int n1 = n + inputPad;

        // and obtaining this much output in total:

        int m1 = int(std::ceil((double(n1) * targetRate) / sourceRate));

        // in order to return this much output to the user:

        int m = int(std::ceil((double(n) * targetRate) / sourceRate));
        
    #ifdef DEBUG_RESAMPLER
        std::cerr << "n = " << n << ", sourceRate = " << sourceRate << ", targetRate = " << targetRate << ", m = " << m << ", latency = " << latency << ", inputPad = " << inputPad << ", m1 = " << m1 << ", n1 = " << n1 << ", n1 - n = " << n1 - n << std::endl;
    #endif

        vector<double> pad(n1 - n, 0.0);
        vector<double> out(m1 + 1, 0.0);

        int gotData = r.process(data, out.data(), n);
        int gotPad = r.process(pad.data(), out.data() + gotData, pad.size());
        int got = gotData + gotPad;
        
    #ifdef DEBUG_RESAMPLER
        std::cerr << "resample: " << n << " in, " << pad.size() << " padding, " << got << " out (" << gotData << " data, " << gotPad << " padding, latency = " << latency << ")" << std::endl;
    #endif
    #ifdef DEBUG_RESAMPLER_VERBOSE
        int printN = 50;
        std::cerr << "first " << printN << " in:" << std::endl;
        for (int i = 0; i < printN && i < n; ++i) {
            if (i % 5 == 0) cerr << endl << i << "... ";
            std::cerr << data[i] << " ";
        }
        std::cerr << std::endl;
    #endif

        int toReturn = got - latency;
        if (toReturn > m) toReturn = m;

        std::vector<double> sliced(out.begin() + latency, 
                            out.begin() + latency + toReturn);

    #ifdef DEBUG_RESAMPLER_VERBOSE
        std::cerr << "first " << printN << " out (after latency compensation), length " << sliced.size() << ":";
        for (int i = 0; i < printN && i < sliced.size(); ++i) {
            if (i % 5 == 0) std::cerr << endl << i << "... ";
            std::cerr << sliced[i] << " ";
        }
        std::cerr << std::endl;
    #endif

        return sliced;
    }

}