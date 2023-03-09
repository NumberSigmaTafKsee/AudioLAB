#pragma once

#include "qmdsp_math.hpp"
#include "qmdsp_filter.hpp"
#include "qmdsp_fft.hpp"
#include "qmdsp_resampler.hpp"
#include <vector>
#include <cstddef>

namespace qmdsp
{
    /**
    * This class takes an input audio signal and a sequence of beat
    * locations (calculated e.g. by TempoTrackV2) and estimates which of
    * the beat locations are downbeats (first beat of the bar).
    * 
    * The input audio signal is expected to have been downsampled to a
    * very low sampling rate (e.g. 2700Hz).  A utility function for
    * downsampling and buffering incoming block-by-block audio is
    * provided.
    */
    class DownBeat
    {
    public:
        /**
        * Construct a downbeat locator that will operate on audio at the
        * downsampled by the given decimation factor from the given
        * original sample rate, plus beats extracted from the same audio
        * at the given original sample rate with the given frame
        * increment.
        *
        * decimationFactor must be a power of two no greater than 64, and
        * dfIncrement must be a multiple of decimationFactor.
        */
        DownBeat(float originalSampleRate,
                size_t decimationFactor,
                size_t dfIncrement);
        ~DownBeat();

        void setBeatsPerBar(int bpb);

        /**
        * Estimate which beats are down-beats.
        * 
        * audio contains the input audio stream after downsampling, and
        * audioLength contains the number of samples in this downsampled
        * stream.
        *
        * beats contains a series of beat positions expressed in
        * multiples of the df increment at the audio's original sample
        * rate, as described to the constructor.
        *
        * The returned downbeat array contains a series of indices to the
        * beats array.
        */
        void findDownBeats(const float *audio, // downsampled
                        size_t audioLength, // after downsampling
                        const std::vector<double> &beats,
                        std::vector<int> &downbeats);

        /**
        * Return the beat spectral difference function.  This is
        * calculated during findDownBeats, so this function can only be
        * meaningfully called after that has completed.  The returned
        * vector contains one value for each of the beat times passed in
        * to findDownBeats, less one.  Each value contains the spectral
        * difference between region prior to the beat's nominal position
        * and the region following it.
        */
        void getBeatSD(std::vector<double> &beatsd) const;
        
        /**
        * For your downsampling convenience: call this function
        * repeatedly with input audio blocks containing dfIncrement
        * samples at the original sample rate, to decimate them to the
        * downsampled rate and buffer them within the DownBeat class.
        *     
        * Call getBufferedAudio() to retrieve the results after all
        * blocks have been processed.
        */
        void pushAudioBlock(const float *audio);
        
        /**
        * Retrieve the accumulated audio produced by pushAudioBlock calls.
        */
        const float *getBufferedAudio(size_t &length) const;

        /**
        * Clear any buffered downsampled audio data.
        */
        void resetAudioBuffer();

    private:
        typedef std::vector<int> i_vec_t;
        typedef std::vector<std::vector<int> > i_mat_t;
        typedef std::vector<double> d_vec_t;
        typedef std::vector<std::vector<double> > d_mat_t;

        void makeDecimators();
        double measureSpecDiff(d_vec_t oldspec, d_vec_t newspec);

        int m_bpb;
        float m_rate;
        size_t m_factor;
        size_t m_increment;
        Decimator *m_decimator1;
        Decimator *m_decimator2;
        float *m_buffer;
        float *m_decbuf;
        size_t m_bufsiz;
        size_t m_buffill;
        size_t m_beatframesize;
        double *m_beatframe;
        FFTReal *m_fft;
        double *m_fftRealOut;
        double *m_fftImagOut;
        d_vec_t m_beatsd;
    };

    struct WinThresh
    {
        int pre;
        int post;
    };

    struct TTParams
    {
        int winLength; //Analysis window length
        int lagLength; //Lag & Stride size
        int alpha; //alpha-norm parameter
        int LPOrd; // low-pass Filter order
        double* LPACoeffs; //low pass Filter den coefficients
        double* LPBCoeffs; //low pass Filter num coefficients
        WinThresh WinT;//window size in frames for adaptive thresholding [pre post]:
    };


    class TempoTrack  
    {
    public:
        TempoTrack( TTParams Params );
        virtual ~TempoTrack();

        std::vector<int> process( std::vector <double> DF,
                                std::vector <double> *tempoReturn = 0);

            
    private:
        void initialise( TTParams Params );
        void deInitialise();

        int beatPredict( int FSP, double alignment, double period, int step);
        int phaseMM( double* DF, double* weighting, int winLength, double period );
        void createPhaseExtractor( double* Filter, int winLength,  double period,  int fsp, int lastBeat );
        int findMeter( double* ACF,  int len, double period );
        void constDetect( double* periodP, int currentIdx, int* flag );
        void stepDetect( double* periodP, double* periodG, int currentIdx, int* flag );
        void createCombFilter( double* Filter, int winLength, int TSig, double beatLag );
        double tempoMM( double* ACF, double* weight, int sig );
            
        int m_dataLength;
        int m_winLength;
        int m_lagLength;

        double m_rayparam;
        double m_sigma;
        double m_DFWVNnorm;

        std::vector<int>  m_beats; // Vector of detected beats

        double m_lockedTempo;

        double* m_tempoScratch;
        double* m_smoothRCF; // Smoothed Output of Comb Filterbank (m_tempoScratch)
            
        // Processing Buffers 
        double* m_rawDFFrame; // Original Detection Function Analysis Frame
        double* m_smoothDFFrame; // Smoothed Detection Function Analysis Frame
        double* m_frameACF; // AutoCorrelation of Smoothed Detection Function 

        //Low Pass Coefficients for DF Smoothing
        double* m_ACoeffs;
        double* m_BCoeffs;
            
        // Objetcs/operators declaration
        Framer m_DFFramer;
        DFProcess* m_DFConditioning;
        Correlation m_correlator;
        // Config structure for DFProcess
        DFProcConfig m_DFPParams;

        // also want to smooth m_tempoScratch 
        DFProcess* m_RCFConditioning;
        // Config structure for RCFProcess
        DFProcConfig m_RCFPParams;
    };    



    //!!! Question: how far is this actually sample rate dependent?  I
    // think it does produce plausible results for e.g. 48000 as well as
    // 44100, but surely the fixed window sizes and comb filtering will
    // make it prefer double or half time when run at e.g. 96000?

    class TempoTrackV2
    {
    public:
        /**
        * Construct a tempo tracker that will operate on beat detection
        * function data calculated from audio at the given sample rate
        * with the given frame increment.
        *
        * Currently the sample rate and increment are used only for the
        * conversion from beat frame location to bpm in the tempo array.
        */
        TempoTrackV2(float sampleRate, int dfIncrement);
        ~TempoTrackV2();

        // Returned beat periods are given in df increment units; inputtempo and tempi in bpm
        void calculateBeatPeriod(const std::vector<double> &df,
                                std::vector<double> &beatPeriod,
                                std::vector<double> &tempi) {
            calculateBeatPeriod(df, beatPeriod, tempi, 120.0, false);
        }

        // Returned beat periods are given in df increment units; inputtempo and tempi in bpm
        // MEPD 28/11/12 Expose inputtempo and constraintempo parameters
        // Note, if inputtempo = 120 and constraintempo = false, then functionality is as it was before
        void calculateBeatPeriod(const std::vector<double> &df,
                                std::vector<double> &beatPeriod,
                                std::vector<double> &tempi,
                                double inputtempo, bool constraintempo);

        // Returned beat positions are given in df increment units
        void calculateBeats(const std::vector<double> &df,
                            const std::vector<double> &beatPeriod,
                            std::vector<double> &beats) {
            calculateBeats(df, beatPeriod, beats, 0.9, 4.0);
        }

        // Returned beat positions are given in df increment units
        // MEPD 28/11/12 Expose alpha and tightness parameters
        // Note, if alpha = 0.9 and tightness = 4, then functionality is as it was before
        void calculateBeats(const std::vector<double> &df,
                            const std::vector<double> &beatPeriod,
                            std::vector<double> &beats,
                            double alpha, double tightness);

    private:
        typedef std::vector<int> i_vec_t;
        typedef std::vector<std::vector<int> > i_mat_t;
        typedef std::vector<double> d_vec_t;
        typedef std::vector<std::vector<double> > d_mat_t;

        float m_rate;
        int m_increment;

        void adapt_thresh(d_vec_t &df);
        double mean_array(const d_vec_t &dfin, int start, int end);
        void filter_df(d_vec_t &df);
        void get_rcf(const d_vec_t &dfframe, const d_vec_t &wv, d_vec_t &rcf);
        void viterbi_decode(const d_mat_t &rcfmat, const d_vec_t &wv,
                            d_vec_t &bp, d_vec_t &tempi);
        double get_max_val(const d_vec_t &df);
        int get_max_ind(const d_vec_t &df);
        void normalise_vec(d_vec_t &df);
    };    


    DownBeat::DownBeat(float originalSampleRate,
                    size_t decimationFactor,
                    size_t dfIncrement) :
        m_bpb(0),
        m_rate(originalSampleRate),
        m_factor(decimationFactor),
        m_increment(dfIncrement),
        m_decimator1(0),
        m_decimator2(0),
        m_buffer(0),
        m_decbuf(0),
        m_bufsiz(0),
        m_buffill(0),
        m_beatframesize(0),
        m_beatframe(0)
    {
        // beat frame size is next power of two up from 1.3 seconds at the
        // downsampled rate (happens to produce 4096 for 44100 or 48000 at
        // 16x decimation, which is our expected normal situation)
        m_beatframesize = MathUtilities::nextPowerOfTwo
            (int((m_rate / decimationFactor) * 1.3));
        if (m_beatframesize < 2) {
            m_beatframesize = 2;
        }
        m_beatframe = new double[m_beatframesize];
        m_fftRealOut = new double[m_beatframesize];
        m_fftImagOut = new double[m_beatframesize];
        m_fft = new FFTReal(m_beatframesize);
    }

    DownBeat::~DownBeat()
    {
        delete m_decimator1;
        delete m_decimator2;
        if (m_buffer) free(m_buffer);
        delete[] m_decbuf;
        delete[] m_beatframe;
        delete[] m_fftRealOut;
        delete[] m_fftImagOut;
        delete m_fft;
    }

    void
    DownBeat::setBeatsPerBar(int bpb)
    {
        m_bpb = bpb;
    }

    void
    DownBeat::makeDecimators()
    {
        if (m_factor < 2) return;
        size_t highest = Decimator::getHighestSupportedFactor();
        if (m_factor <= highest) {
            m_decimator1 = new Decimator(m_increment, m_factor);
            return;
        }
        m_decimator1 = new Decimator(m_increment, highest);
        m_decimator2 = new Decimator(m_increment / highest, m_factor / highest);
        m_decbuf = new float[m_increment / highest];
    }

    void
    DownBeat::pushAudioBlock(const float *audio)
    {
        if (m_buffill + (m_increment / m_factor) > m_bufsiz) {
            if (m_bufsiz == 0) m_bufsiz = m_increment * 16;
            else m_bufsiz = m_bufsiz * 2;
            if (!m_buffer) {
                m_buffer = (float *)malloc(m_bufsiz * sizeof(float));
            } else {
                m_buffer = (float *)realloc(m_buffer, m_bufsiz * sizeof(float));
            }
        }
        if (!m_decimator1 && m_factor > 1) {
            makeDecimators();
        }
        if (m_decimator2) {
            m_decimator1->process(audio, m_decbuf);
            m_decimator2->process(m_decbuf, m_buffer + m_buffill);
        } else if (m_decimator1) {
            m_decimator1->process(audio, m_buffer + m_buffill);
        } else {
            // just copy across (m_factor is presumably 1)
            for (size_t i = 0; i < m_increment; ++i) {
                (m_buffer + m_buffill)[i] = audio[i];
            }
        }
        m_buffill += m_increment / m_factor;
    }
        
    const float *
    DownBeat::getBufferedAudio(size_t &length) const
    {
        length = m_buffill;
        return m_buffer;
    }

    void
    DownBeat::resetAudioBuffer()
    {
        if (m_buffer) {
            free(m_buffer);
        }
        m_buffer = 0;
        m_buffill = 0;
        m_bufsiz = 0;
    }

    void
    DownBeat::findDownBeats(const float *audio,
                            size_t audioLength,
                            const d_vec_t &beats,
                            i_vec_t &downbeats)
    {
        // FIND DOWNBEATS BY PARTITIONING THE INPUT AUDIO FILE INTO BEAT SEGMENTS
        // WHERE THE AUDIO FRAMES ARE DOWNSAMPLED  BY A FACTOR OF 16 (fs ~= 2700Hz)
        // THEN TAKING THE JENSEN-SHANNON DIVERGENCE BETWEEN BEAT SYNCHRONOUS SPECTRAL FRAMES

        // IMPLEMENTATION (MOSTLY) FOLLOWS:
        //  DAVIES AND PLUMBLEY "A SPECTRAL DIFFERENCE APPROACH TO EXTRACTING DOWNBEATS IN MUSICAL AUDIO"
        //  EUSIPCO 2006, FLORENCE, ITALY

        d_vec_t newspec(m_beatframesize / 2); // magnitude spectrum of current beat
        d_vec_t oldspec(m_beatframesize / 2); // magnitude spectrum of previous beat

        m_beatsd.clear();

        if (audioLength == 0) return;

        for (size_t i = 0; i + 1 < beats.size(); ++i) {

            // Copy the extents of the current beat from downsampled array
            // into beat frame buffer

            size_t beatstart = (beats[i] * m_increment) / m_factor;
            size_t beatend = (beats[i+1] * m_increment) / m_factor;
            if (beatend >= audioLength) beatend = audioLength - 1;
            if (beatend < beatstart) beatend = beatstart;
            size_t beatlen = beatend - beatstart;

            // Also apply a Hanning window to the beat frame buffer, sized
            // to the beat extents rather than the frame size.  (Because
            // the size varies, it's easier to do this by hand than use
            // our Window abstraction.)

            for (size_t j = 0; j < beatlen && j < m_beatframesize; ++j) {
                double mul = 0.5 * (1.0 - cos(TWO_PI * (double(j) / double(beatlen))));
                m_beatframe[j] = audio[beatstart + j] * mul;
            }

            for (size_t j = beatlen; j < m_beatframesize; ++j) {
                m_beatframe[j] = 0.0;
            }

            // Now FFT beat frame
            
            m_fft->forward(m_beatframe, m_fftRealOut, m_fftImagOut);
            
            // Calculate magnitudes

            for (size_t j = 0; j < m_beatframesize/2; ++j) {
                newspec[j] = sqrt(m_fftRealOut[j] * m_fftRealOut[j] +
                                m_fftImagOut[j] * m_fftImagOut[j]);
            }

            // Preserve peaks by applying adaptive threshold

            MathUtilities::adaptiveThreshold(newspec);

            // Calculate JS divergence between new and old spectral frames

            if (i > 0) { // otherwise we have no previous frame
                m_beatsd.push_back(measureSpecDiff(oldspec, newspec));
            }

            // Copy newspec across to old

            for (size_t j = 0; j < m_beatframesize/2; ++j) {
                oldspec[j] = newspec[j];
            }
        }

        // We now have all spectral difference measures in specdiff

        int timesig = m_bpb;
        if (timesig == 0) timesig = 4;

        d_vec_t dbcand(timesig); // downbeat candidates

        for (int beat = 0; beat < timesig; ++beat) {
            dbcand[beat] = 0;
        }

    // look for beat transition which leads to greatest spectral change
    for (int beat = 0; beat < timesig; ++beat) {
        int count = 0;
        for (int example = beat-1; example < (int)m_beatsd.size(); example += timesig) {
            if (example < 0) continue;
            dbcand[beat] += (m_beatsd[example]) / timesig;
            ++count;
        }
        if (count > 0) {
            dbcand[beat] /= count;
        }
    }

        // first downbeat is beat at index of maximum value of dbcand
        int dbind = MathUtilities::getMax(dbcand);

        // remaining downbeats are at timesig intervals from the first
        for (int i = dbind; i < (int)beats.size(); i += timesig) {
            downbeats.push_back(i);
        }
    }

    double
    DownBeat::measureSpecDiff(d_vec_t oldspec, d_vec_t newspec)
    {
        // JENSEN-SHANNON DIVERGENCE BETWEEN SPECTRAL FRAMES

        int SPECSIZE = 512;   // ONLY LOOK AT FIRST 512 SAMPLES OF SPECTRUM. 
        if (SPECSIZE > int(oldspec.size())/4) {
            SPECSIZE = int(oldspec.size())/4;
        }
        double SD = 0.;
        double sd1 = 0.;

        double sumnew = 0.;
        double sumold = 0.;
    
        for (int i = 0;i < SPECSIZE;i++) {
            
            newspec[i] +=EPS;
            oldspec[i] +=EPS;
            
            sumnew+=newspec[i];
            sumold+=oldspec[i];
        } 
        
        for (int i = 0;i < SPECSIZE;i++) {
            
            newspec[i] /= (sumnew);
            oldspec[i] /= (sumold);
            
            // IF ANY SPECTRAL VALUES ARE 0 (SHOULDN'T BE ANY!) SET THEM TO 1
            if (newspec[i] == 0) {
                newspec[i] = 1.;
            }
            
            if (oldspec[i] == 0) {
                oldspec[i] = 1.;
            }
            
            // JENSEN-SHANNON CALCULATION
            sd1 = 0.5*oldspec[i] + 0.5*newspec[i];  
            SD = SD + (-sd1*log(sd1)) +
                (0.5*(oldspec[i]*log(oldspec[i]))) +
                (0.5*(newspec[i]*log(newspec[i])));
        }
        
        return SD;
    }

    void
    DownBeat::getBeatSD(vector<double> &beatsd) const
    {
        for (int i = 0; i < (int)m_beatsd.size(); ++i) {
            beatsd.push_back(m_beatsd[i]);
        }
    }    


    //////////////////////////////////////////////////////////////////////
    // Construction/Destruction
    //////////////////////////////////////////////////////////////////////

    TempoTrack::TempoTrack( TTParams Params )
    {
        m_tempoScratch = NULL;
        m_rawDFFrame = NULL;
        m_smoothDFFrame = NULL;
        m_frameACF = NULL;
        m_smoothRCF = NULL;

        m_dataLength = 0;
        m_winLength = 0;
        m_lagLength = 0;

        m_rayparam = 0;
        m_sigma = 0;
        m_DFWVNnorm = 0;

        initialise( Params );
    }

    TempoTrack::~TempoTrack()
    {
        deInitialise();
    }

    void TempoTrack::initialise( TTParams Params )
    {       
        m_winLength = Params.winLength;
        m_lagLength = Params.lagLength;

        m_rayparam   = 43.0;
        m_sigma = sqrt(3.9017);
        m_DFWVNnorm = exp( ( log( 2.0 ) / m_rayparam ) * ( m_winLength + 2 ) );

        m_rawDFFrame = new double[ m_winLength ];
        m_smoothDFFrame = new double[ m_winLength ];
        m_frameACF = new double[ m_winLength ];
        m_tempoScratch = new double[ m_lagLength ];
        m_smoothRCF = new double[ m_lagLength ];

        m_DFFramer.configure( m_winLength, m_lagLength );
            
        m_DFPParams.length = m_winLength;
        m_DFPParams.AlphaNormParam = Params.alpha;
        m_DFPParams.LPOrd = Params.LPOrd;
        m_DFPParams.LPACoeffs = Params.LPACoeffs;
        m_DFPParams.LPBCoeffs = Params.LPBCoeffs;
        m_DFPParams.winPre = Params.WinT.pre;
        m_DFPParams.winPost = Params.WinT.post;
        m_DFPParams.isMedianPositive = true;
            
        m_DFConditioning = new DFProcess( m_DFPParams );

        // these are parameters for smoothing m_tempoScratch
        m_RCFPParams.length = m_lagLength;
        m_RCFPParams.AlphaNormParam = Params.alpha;
        m_RCFPParams.LPOrd = Params.LPOrd;
        m_RCFPParams.LPACoeffs = Params.LPACoeffs;
        m_RCFPParams.LPBCoeffs = Params.LPBCoeffs;
        m_RCFPParams.winPre = Params.WinT.pre;
        m_RCFPParams.winPost = Params.WinT.post;
        m_RCFPParams.isMedianPositive = true;

        m_RCFConditioning = new DFProcess( m_RCFPParams );
    }

    void TempoTrack::deInitialise()
    {       
        delete [] m_rawDFFrame;
        delete [] m_smoothDFFrame;
        delete [] m_smoothRCF;  
        delete [] m_frameACF;
        delete [] m_tempoScratch;
        delete m_DFConditioning;
        delete m_RCFConditioning;
    }

    void TempoTrack::createCombFilter(double* Filter, int winLength, int /* TSig */, double beatLag)
    {
        int i;

        if( beatLag == 0 ) {
            for( i = 0; i < winLength; i++ ) {    
                Filter[ i ] =
                    ( ( i + 1 ) / pow( m_rayparam, 2.0) ) *
                    exp( ( -pow(( i + 1 ),2.0 ) /
                        ( 2.0 * pow( m_rayparam, 2.0))));
            }
        } else {   
            m_sigma = beatLag/4;
            for( i = 0; i < winLength; i++ ) {
                double dlag = (double)(i+1) - beatLag;
                Filter[ i ] =  exp(-0.5 * pow(( dlag / m_sigma), 2.0) ) /
                    (sqrt(TWO_PI) * m_sigma);
            }
        }
    }

    double TempoTrack::tempoMM(double* ACF, double* weight, int tsig)
    {
        double period = 0;
        double maxValRCF = 0.0;
        int maxIndexRCF = 0;

        double* pdPeaks;

        int maxIndexTemp;
        double maxValTemp;
        int count; 
            
        int numelem,i,j;
        int a, b;

        for( i = 0; i < m_lagLength; i++ ) {
            m_tempoScratch[ i ] = 0.0;
        }

        if( tsig == 0 ) {
            //if time sig is unknown, use metrically unbiased version of Filterbank
            numelem = 4;
        } else {
            numelem = tsig;
        }

    #ifdef DEBUG_TEMPO_TRACK
        std::cerr << "tempoMM: m_winLength = " << m_winLength << ", m_lagLength = " << m_lagLength << ", numelem = " << numelem << std::endl;
    #endif

        for(i=1;i<m_lagLength-1;i++) {
            //first and last output values are left intentionally as zero
            for (a=1;a<=numelem;a++) {
                for(b=(1-a);b<a;b++) {
                    if( tsig == 0 ) {                                       
                        m_tempoScratch[i] += ACF[a*(i+1)+b-1] * (1.0 / (2.0 * (double)a-1)) * weight[i];
                    } else {
                        m_tempoScratch[i] += ACF[a*(i+1)+b-1] * 1 * weight[i];
                    }
                }
            }
        }


        //////////////////////////////////////////////////
        // MODIFIED BEAT PERIOD EXTRACTION //////////////
        /////////////////////////////////////////////////

        // find smoothed version of RCF ( as applied to Detection Function)
        m_RCFConditioning->process( m_tempoScratch, m_smoothRCF);

        if (tsig != 0) { // i.e. in context dependent state

    //     NOW FIND MAX INDEX OF ACFOUT
            for( i = 0; i < m_lagLength; i++) {
                if( m_tempoScratch[ i ] > maxValRCF) {
                    maxValRCF = m_tempoScratch[ i ];
                    maxIndexRCF = i;
                }
            }

        } else { // using rayleigh weighting

            vector <vector<double> > rcfMat;
            
            double sumRcf = 0.;
            
            double maxVal = 0.;
            // now find the two values which minimise rcfMat
            double minVal = 0.;
            int p_i = 1; // periodicity for row i;
            int p_j = 1; //periodicity for column j;
            
            for ( i=0; i<m_lagLength; i++) {
                m_tempoScratch[i] =m_smoothRCF[i];
            }       

            // normalise m_tempoScratch so that it sums to zero.
            for ( i=0; i<m_lagLength; i++) {
                sumRcf += m_tempoScratch[i];
            }       
            
            for( i=0; i<m_lagLength; i++) {
                m_tempoScratch[i] /= sumRcf;
            }       
            
            // create a matrix to store m_tempoScratchValues modified by log2 ratio
            for ( i=0; i<m_lagLength; i++) {
                rcfMat.push_back  ( vector<double>() ); // adds a new row...
            }
            
            for (i=0; i<m_lagLength; i++) {
                for (j=0; j<m_lagLength; j++) {
                    rcfMat[i].push_back (0.);
                }
            }
            
            // the 'i' and 'j' indices deliberately start from '1' and not '0'
            for ( i=1; i<m_lagLength; i++) {
                for (j=1; j<m_lagLength; j++) {
                    double log2PeriodRatio = log( static_cast<double>(i)/
                                                static_cast<double>(j) ) /
                        log(2.0);
                    rcfMat[i][j] = ( abs(1.0-abs(log2PeriodRatio)) );
                    rcfMat[i][j] += ( 0.01*( 1./(m_tempoScratch[i]+m_tempoScratch[j]) ) );
                }
            }
                    
            // set diagonal equal to maximum value in rcfMat 
            // we don't want to pick one strong middle peak - we need a combination of two peaks.
            
            for ( i=1; i<m_lagLength; i++) {
                for (j=1; j<m_lagLength; j++) {
                    if (rcfMat[i][j] > maxVal) {       
                        maxVal = rcfMat[i][j];
                    }
                }
            }
            
            for ( i=1; i<m_lagLength; i++) {
                rcfMat[i][i] = maxVal;
            }
            
            // now find the row and column number which minimise rcfMat
            minVal = maxVal;
                    
            for ( i=1; i<m_lagLength; i++) {
                for ( j=1; j<m_lagLength; j++) {
                    if (rcfMat[i][j] < minVal) {       
                        minVal = rcfMat[i][j];
                        p_i = i;
                        p_j = j;
                    }
                }
            }
            
            
            // initially choose p_j (arbitrary) - saves on an else statement
            int beatPeriod = p_j;
            if (m_tempoScratch[p_i] > m_tempoScratch[p_j]) {
                beatPeriod = p_i;
            }
                    
            // now write the output
            maxIndexRCF = static_cast<int>(beatPeriod);
        }


        double locked = 5168.f / maxIndexRCF;
        if (locked >= 30 && locked <= 180) {
            m_lockedTempo = locked;
        }

    #ifdef DEBUG_TEMPO_TRACK
        std::cerr << "tempoMM: locked tempo = " << m_lockedTempo << std::endl;
    #endif

        if( tsig == 0 ) {
            tsig = 4;
        }

    #ifdef DEBUG_TEMPO_TRACK
        std::cerr << "tempoMM: maxIndexRCF = " << maxIndexRCF << std::endl;
    #endif
            
        if( tsig == 4 ) {
            
    #ifdef DEBUG_TEMPO_TRACK
            std::cerr << "tsig == 4" << std::endl;
    #endif

            pdPeaks = new double[ 4 ];
            for( i = 0; i < 4; i++ ){ pdPeaks[ i ] = 0.0;}

            pdPeaks[ 0 ] = ( double )maxIndexRCF + 1;

            maxIndexTemp = 0;
            maxValTemp = 0.0;
            count = 0;

            for( i = (2 * maxIndexRCF + 1) - 1; i < (2 * maxIndexRCF + 1) + 2; i++ ) {
                if( ACF[ i ] > maxValTemp ) {
                    maxValTemp = ACF[ i ];
                    maxIndexTemp = count;
                }
                count++;
            }
            pdPeaks[ 1 ] = (double)( maxIndexTemp + 1 + ( (2 * maxIndexRCF + 1 ) - 2 ) + 1 )/2;

            maxIndexTemp = 0;
            maxValTemp = 0.0;
            count = 0;

            for( i = (3 * maxIndexRCF + 2 ) - 2; i < (3 * maxIndexRCF + 2 ) + 3; i++ ) {
                if( ACF[ i ] > maxValTemp ) {
                    maxValTemp = ACF[ i ];
                    maxIndexTemp = count;
                }
                count++;
            }
            pdPeaks[ 2 ] = (double)( maxIndexTemp + 1 + ( (3 * maxIndexRCF + 2) - 4 ) + 1 )/3;

            maxIndexTemp = 0;
            maxValTemp = 0.0;
            count = 0;

            for( i = ( 4 * maxIndexRCF + 3) - 3; i < ( 4 * maxIndexRCF + 3) + 4; i++ ) {
                if( ACF[ i ] > maxValTemp ) {
                    maxValTemp = ACF[ i ];
                    maxIndexTemp = count;
                }
                count++;
            }

            pdPeaks[ 3 ] = (double)( maxIndexTemp + 1 + ( (4 * maxIndexRCF + 3) - 9 ) + 1 )/4 ;


            period = MathUtilities::mean( pdPeaks, 4 );

        } else {
            
    #ifdef DEBUG_TEMPO_TRACK
            std::cerr << "tsig != 4" << std::endl;
    #endif

            pdPeaks = new double[ 3 ];
            for( i = 0; i < 3; i++ ) {
                pdPeaks[ i ] = 0.0;
            }

            pdPeaks[ 0 ] = ( double )maxIndexRCF + 1;

            maxIndexTemp = 0;
            maxValTemp = 0.0;
            count = 0;

            for( i = (2 * maxIndexRCF + 1) - 1; i < (2 * maxIndexRCF + 1) + 2; i++ ) {
                if( ACF[ i ] > maxValTemp ) {
                    maxValTemp = ACF[ i ];
                    maxIndexTemp = count;
                }
                count++;
            }
            pdPeaks[ 1 ] = (double)( maxIndexTemp + 1 + ( (2 * maxIndexRCF + 1 ) - 2 ) + 1 )/2;

            maxIndexTemp = 0;
            maxValTemp = 0.0;
            count = 0;

            for( i = (3 * maxIndexRCF + 2 ) - 2; i < (3 * maxIndexRCF + 2 ) + 3; i++ ) {
                if( ACF[ i ] > maxValTemp ) {
                    maxValTemp = ACF[ i ];
                    maxIndexTemp = count;
                }
                count++;
            }
            pdPeaks[ 2 ] = (double)( maxIndexTemp + 1 + ( (3 * maxIndexRCF + 2) - 4 ) + 1 )/3;


            period = MathUtilities::mean( pdPeaks, 3 );
        }

        delete [] pdPeaks;

        return period;
    }

    void TempoTrack::stepDetect( double* periodP, double* periodG, int currentIdx, int* flag )
    {
        double stepthresh = 1 * 3.9017;

        if( *flag ) {
            if(abs(periodG[ currentIdx ] - periodP[ currentIdx ]) > stepthresh) {
                // do nuffin'
            }
        } else {
            if(fabs(periodG[ currentIdx ]-periodP[ currentIdx ]) > stepthresh) {
                *flag = 3;
            }
        }
    }

    void TempoTrack::constDetect( double* periodP, int currentIdx, int* flag )
    {
        double constthresh = 2 * 3.9017;

        if( fabs( 2 * periodP[ currentIdx ] - periodP[ currentIdx - 1] - periodP[ currentIdx - 2] ) < constthresh) {
            *flag = 1;
        } else {
            *flag = 0;
        }
    }

    int TempoTrack::findMeter(double *ACF, int len, double period)
    {
        int i;
        int p = (int)MathUtilities::round( period );
        int tsig;

        double Energy_3 = 0.0;
        double Energy_4 = 0.0;

        double temp3A = 0.0;
        double temp3B = 0.0;
        double temp4A = 0.0;
        double temp4B = 0.0;

        double* dbf = new double[ len ]; int t = 0;
        for( int u = 0; u < len; u++ ){ dbf[ u ] = 0.0; }

        if( (double)len < 6 * p + 2 ) {
            
            for( i = ( 3 * p - 2 ); i < ( 3 * p + 2 ) + 1; i++ ) {
                temp3A += ACF[ i ];
                dbf[ t++ ] = ACF[ i ];
            }
            
            for( i = ( 4 * p - 2 ); i < ( 4 * p + 2 ) + 1; i++ ) {
                temp4A += ACF[ i ];
            }

            Energy_3 = temp3A;
            Energy_4 = temp4A;

        } else {
            
            for( i = ( 3 * p - 2 ); i < ( 3 * p + 2 ) + 1; i++ ) {
                temp3A += ACF[ i ];
            }
            
            for( i = ( 4 * p - 2 ); i < ( 4 * p + 2 ) + 1; i++ ) {
                temp4A += ACF[ i ];
            }

            for( i = ( 6 * p - 2 ); i < ( 6 * p + 2 ) + 1; i++ ) {
                temp3B += ACF[ i ];
            }
            
            for( i = ( 2 * p - 2 ); i < ( 2 * p + 2 ) + 1; i++ ) {
                temp4B += ACF[ i ];
            }

            Energy_3 = temp3A + temp3B;
            Energy_4 = temp4A + temp4B;
        }

        if (Energy_3 > Energy_4) {
            tsig = 3;
        } else {
            tsig = 4;
        }

        return tsig;
    }

    void TempoTrack::createPhaseExtractor(double *Filter, int /* winLength */, double period, int fsp, int lastBeat)
    {       
        int p = (int)MathUtilities::round( period );
        int predictedOffset = 0;

    #ifdef DEBUG_TEMPO_TRACK
        std::cerr << "TempoTrack::createPhaseExtractor: period = " << period << ", p = " << p << std::endl;
    #endif

        if (p > 10000) {
            std::cerr << "TempoTrack::createPhaseExtractor: WARNING! Highly implausible period value " << p << "!" << std::endl;
            period = 5168 / 120;
        }

        double* phaseScratch = new double[ p*2 + 2 ];
        for (int i = 0; i < p*2 + 2; ++i) phaseScratch[i] = 0.0;

            
        if ( lastBeat != 0 ) {
            
            lastBeat = (int)MathUtilities::round((double)lastBeat );///(double)winLength);

            predictedOffset = lastBeat + p - fsp;

            if (predictedOffset < 0) {
                lastBeat = 0;
            }
        }

        if ( lastBeat != 0 ) {
            
            int mu = p;
            double sigma = (double)p/8;
            double PhaseMin = 0.0;
            double PhaseMax = 0.0;
            int scratchLength = p*2;
            double temp = 0.0;

            for(  int i = 0; i < scratchLength; i++ ) {
                phaseScratch[ i ] = exp( -0.5 * pow( ( i - mu ) / sigma, 2 ) ) / ( sqrt(TWO_PI) *sigma );
            }

            MathUtilities::getFrameMinMax( phaseScratch, scratchLength, &PhaseMin, &PhaseMax );
                            
            for(int i = 0; i < scratchLength; i ++) {
                temp = phaseScratch[ i ];
                phaseScratch[ i ] = (temp - PhaseMin)/PhaseMax;
            }

    #ifdef DEBUG_TEMPO_TRACK
            std::cerr << "predictedOffset = " << predictedOffset << std::endl;
    #endif

            int index = 0;
            for (int i = p - ( predictedOffset - 1); i < p + ( p - predictedOffset) + 1; i++) {
    #ifdef DEBUG_TEMPO_TRACK
                std::cerr << "assigning to filter index " << index << " (size = " << p*2 << ")" << " value " << phaseScratch[i] << " from scratch index " << i << std::endl;
    #endif
                Filter[ index++ ] = phaseScratch[ i ];
            }
        } else {
            for( int i = 0; i < p; i ++) {
                Filter[ i ] = 1;
            }
        }
            
        delete [] phaseScratch;
    }

    int TempoTrack::phaseMM(double *DF, double *weighting, int winLength, double period)
    {
        int alignment = 0;
        int p = (int)MathUtilities::round( period );

        double temp = 0.0;

        double* y = new double[ winLength ];
        double* align = new double[ p ];

        for( int i = 0; i < winLength; i++ ) {   
            y[ i ] = (double)( -i + winLength  )/(double)winLength;
            y[ i ] = pow(y [i ],2.0); // raise to power 2.
        }

        for( int o = 0; o < p; o++ ) { 
            temp = 0.0;
            for (int i = 1 + (o - 1); i < winLength; i += (p + 1)) {
                temp = temp + DF[ i ] * y[ i ]; 
            }
            align[ o ] = temp * weighting[ o ];       
        }


        double valTemp = 0.0;
        for(int i = 0; i < p; i++) {
            if( align[ i ] > valTemp ) {
                valTemp = align[ i ];
                alignment = i;
            }
        }

        delete [] y;
        delete [] align;

        return alignment;
    }

    int TempoTrack::beatPredict(int FSP0, double alignment, double period, int step )
    {
        int beat = 0;

        int p = (int)MathUtilities::round( period );
        int align = (int)MathUtilities::round( alignment );
        int FSP = (int)MathUtilities::round( FSP0 );

        int FEP = FSP + ( step );

        beat = FSP + align;

        m_beats.push_back( beat );

        while( beat + p < FEP ) {
            beat += p;
            m_beats.push_back( beat );
        }

        return beat;
    }



    vector<int> TempoTrack::process( vector <double> DF,
                                    vector <double> *tempoReturn )
    {
        m_dataLength = DF.size();
            
        m_lockedTempo = 0.0;

        double period = 0.0;
        int stepFlag = 0;
        int constFlag = 0;
        int FSP = 0;
        int tsig = 0;
        int lastBeat = 0;

        vector <double> causalDF;

        causalDF = DF;

        //Prepare Causal Extension DFData
    //    int DFCLength = m_dataLength + m_winLength;
            
        for( int j = 0; j < m_winLength; j++ ) {
            causalDF.push_back( 0 );
        }
            
            
        double* RW = new double[ m_lagLength ];
        for (int clear = 0; clear < m_lagLength; clear++){ RW[ clear ] = 0.0;}

        double* GW = new double[ m_lagLength ];
        for (int clear = 0; clear < m_lagLength; clear++){ GW[ clear ] = 0.0;}

        double* PW = new double[ m_lagLength ];
        for(int clear = 0; clear < m_lagLength; clear++){ PW[ clear ] = 0.0;}

        m_DFFramer.setSource( &causalDF[0], m_dataLength );

        int TTFrames = m_DFFramer.getMaxNoFrames();

    #ifdef DEBUG_TEMPO_TRACK
        std::cerr << "TTFrames = " << TTFrames << std::endl;
    #endif
            
        double* periodP = new double[ TTFrames ];
        for(int clear = 0; clear < TTFrames; clear++){ periodP[ clear ] = 0.0;}
            
        double* periodG = new double[ TTFrames ];
        for(int clear = 0; clear < TTFrames; clear++){ periodG[ clear ] = 0.0;}
            
        double* alignment = new double[ TTFrames ];
        for(int clear = 0; clear < TTFrames; clear++){ alignment[ clear ] = 0.0;}

        m_beats.clear();

        createCombFilter( RW, m_lagLength, 0, 0 );

        int TTLoopIndex = 0;

        for( int i = 0; i < TTFrames; i++ ) {
            
            m_DFFramer.getFrame( m_rawDFFrame );

            m_DFConditioning->process( m_rawDFFrame, m_smoothDFFrame );

            m_correlator.doAutoUnBiased( m_smoothDFFrame, m_frameACF, m_winLength );
                    
            periodP[ TTLoopIndex ] = tempoMM( m_frameACF, RW, 0 );

            if( GW[ 0 ] != 0 ) {
                periodG[ TTLoopIndex ] = tempoMM( m_frameACF, GW, tsig );
            } else {
                periodG[ TTLoopIndex ] = 0.0;
            }

            stepDetect( periodP, periodG, TTLoopIndex, &stepFlag );

            if( stepFlag == 1) {
                constDetect( periodP, TTLoopIndex, &constFlag );
                stepFlag = 0;
            } else {
                stepFlag -= 1;
            }

            if( stepFlag < 0 ) {
                stepFlag = 0;
            }

            if( constFlag != 0) {
                
                tsig = findMeter( m_frameACF, m_winLength, periodP[ TTLoopIndex ] );
            
                createCombFilter( GW, m_lagLength, tsig, periodP[ TTLoopIndex ] );
                            
                periodG[ TTLoopIndex ] = tempoMM( m_frameACF, GW, tsig ); 

                period = periodG[ TTLoopIndex ];

    #ifdef DEBUG_TEMPO_TRACK
                std::cerr << "TempoTrack::process: constFlag == " << constFlag << ", TTLoopIndex = " << TTLoopIndex << ", period from periodG = " << period << std::endl;
    #endif

                createPhaseExtractor( PW, m_winLength, period, FSP, 0 ); 

                constFlag = 0;

            } else {
                
                if( GW[ 0 ] != 0 ) {
                    period = periodG[ TTLoopIndex ];

    #ifdef DEBUG_TEMPO_TRACK
                    std::cerr << "TempoTrack::process: GW[0] == " << GW[0] << ", TTLoopIndex = " << TTLoopIndex << ", period from periodG = " << period << std::endl;
    #endif

                    if (period > 10000) {
                        std::cerr << "TempoTrack::process: WARNING!  Highly implausible period value " << period << "!" << std::endl;
                        std::cerr << "periodG contains (of " << TTFrames << " frames): " << std::endl;
                        for (int i = 0; i < TTLoopIndex + 3 && i < TTFrames; ++i) {
                            std::cerr << i << " -> " << periodG[i] << std::endl;
                        }
                        std::cerr << "periodP contains (of " << TTFrames << " frames): " << std::endl;
                        for (int i = 0; i < TTLoopIndex + 3 && i < TTFrames; ++i) {
                            std::cerr << i << " -> " << periodP[i] << std::endl;
                        }
                        period = 5168 / 120;
                    }

                    createPhaseExtractor( PW, m_winLength, period, FSP, lastBeat ); 

                }
                else
                {
                    period = periodP[ TTLoopIndex ];

    #ifdef DEBUG_TEMPO_TRACK
                    std::cerr << "TempoTrack::process: GW[0] == " << GW[0] << ", TTLoopIndex = " << TTLoopIndex << ", period from periodP = " << period << std::endl;
    #endif

                    createPhaseExtractor( PW, m_winLength, period, FSP, 0 ); 
                }
            }

            alignment[ TTLoopIndex ] = phaseMM( m_rawDFFrame, PW, m_winLength, period ); 

            lastBeat = beatPredict(FSP, alignment[ TTLoopIndex ], period, m_lagLength );

            FSP += (m_lagLength);

            if (tempoReturn) tempoReturn->push_back(m_lockedTempo);

            TTLoopIndex++;
        }


        delete [] periodP;
        delete [] periodG;
        delete [] alignment;

        delete [] RW;
        delete [] GW;
        delete [] PW;

        return m_beats;
    }


    #define   EPS 0.0000008 // just some arbitrary small number

    TempoTrackV2::TempoTrackV2(float rate, int increment) :
        m_rate(rate), m_increment(increment) {
    }

    TempoTrackV2::~TempoTrackV2() { }

    void
    TempoTrackV2::filter_df(d_vec_t &df)
    {
        int df_len = int(df.size());
        
        d_vec_t a(3);
        d_vec_t b(3);
        d_vec_t lp_df(df_len);

        //equivalent in matlab to [b,a] = butter(2,0.4);
        a[0] = 1.0000;
        a[1] = -0.3695;
        a[2] = 0.1958;
        b[0] = 0.2066;
        b[1] = 0.4131;
        b[2] = 0.2066;

        double inp1 = 0.;
        double inp2 = 0.;
        double out1 = 0.;
        double out2 = 0.;


        // forwards filtering
        for (int i = 0; i < df_len; i++) {
            lp_df[i] =  b[0]*df[i] + b[1]*inp1 + b[2]*inp2 - a[1]*out1 - a[2]*out2;
            inp2 = inp1;
            inp1 = df[i];
            out2 = out1;
            out1 = lp_df[i];
        }

        // copy forwards filtering to df...
        // but, time-reversed, ready for backwards filtering
        for (int i = 0; i < df_len; i++) {
            df[i] = lp_df[df_len - i - 1];
        }

        for (int i = 0; i < df_len; i++) {
            lp_df[i] = 0.;
        }

        inp1 = 0.; inp2 = 0.;
        out1 = 0.; out2 = 0.;

        // backwards filetering on time-reversed df
        for (int i = 0; i < df_len; i++) {
            lp_df[i] =  b[0]*df[i] + b[1]*inp1 + b[2]*inp2 - a[1]*out1 - a[2]*out2;
            inp2 = inp1;
            inp1 = df[i];
            out2 = out1;
            out1 = lp_df[i];
        }

        // write the re-reversed (i.e. forward) version back to df
        for (int i = 0; i < df_len; i++) {
            df[i] = lp_df[df_len - i - 1];
        }
    }


    // MEPD 28/11/12
    // This function now allows for a user to specify an inputtempo (in BPM)
    // and a flag "constraintempo" which replaces the general rayleigh weighting for periodicities
    // with a gaussian which is centered around the input tempo
    // Note, if inputtempo = 120 and constraintempo = false, then functionality is
    // as it was before
    void
    TempoTrackV2::calculateBeatPeriod(const vector<double> &df,
                                    vector<double> &beat_period,
                                    vector<double> &tempi,
                                    double inputtempo, bool constraintempo)
    {
        // to follow matlab.. split into 512 sample frames with a 128 hop size
        // calculate the acf,
        // then the rcf.. and then stick the rcfs as columns of a matrix
        // then call viterbi decoding with weight vector and transition matrix
        // and get best path

        int wv_len = 128;

        // MEPD 28/11/12
        // the default value of inputtempo in the beat tracking plugin is 120
        // so if the user specifies a different inputtempo, the rayparam will be updated
        // accordingly.
        // note: 60*44100/512 is a magic number
        // this might (will?) break if a user specifies a different frame rate for the onset detection function
        double rayparam = (60*44100/512)/inputtempo;

        // make rayleigh weighting curve
        d_vec_t wv(wv_len);

        // check whether or not to use rayleigh weighting (if constraintempo is false)
        // or use gaussian weighting it (constraintempo is true)
        if (constraintempo) {
            for (int i = 0; i < wv_len; i++) {
                // MEPD 28/11/12
                // do a gaussian weighting instead of rayleigh
                wv[i] = exp( (-1.*pow((double(i)-rayparam),2.)) / (2.*pow(rayparam/4.,2.)) );
            }
        } else {
            for (int i = 0; i < wv_len; i++) {
                // MEPD 28/11/12
                // standard rayleigh weighting over periodicities
                wv[i] = (double(i) / pow(rayparam,2.)) * exp((-1.*pow(-double(i),2.)) / (2.*pow(rayparam,2.)));
            }
        }

        // beat tracking frame size (roughly 6 seconds) and hop (1.5 seconds)
        int winlen = 512;
        int step = 128;

        // matrix to store output of comb filter bank, increment column of matrix at each frame
        d_mat_t rcfmat;
        int col_counter = -1;
        int df_len = int(df.size());

        // main loop for beat period calculation
        for (int i = 0; i+winlen < df_len; i+=step) {
            
            // get dfframe
            d_vec_t dfframe(winlen);
            for (int k=0; k < winlen; k++) {
                dfframe[k] = df[i+k];
            }
            // get rcf vector for current frame
            d_vec_t rcf(wv_len);
            get_rcf(dfframe,wv,rcf);

            rcfmat.push_back( d_vec_t() ); // adds a new column
            col_counter++;
            for (int j = 0; j < wv_len; j++) {
                rcfmat[col_counter].push_back( rcf[j] );
            }
        }

        // now call viterbi decoding function
        viterbi_decode(rcfmat,wv,beat_period,tempi);
    }


    void
    TempoTrackV2::get_rcf(const d_vec_t &dfframe_in, const d_vec_t &wv, d_vec_t &rcf)
    {
        // calculate autocorrelation function
        // then rcf
        // just hard code for now... don't really need separate functions to do this

        // make acf

        d_vec_t dfframe(dfframe_in);

        MathUtilities::adaptiveThreshold(dfframe);

        int dfframe_len = int(dfframe.size());
        int rcf_len = int(rcf.size());
        
        d_vec_t acf(dfframe_len);

        for (int lag = 0; lag < dfframe_len; lag++) {
            double sum = 0.;
            double tmp = 0.;

            for (int n = 0; n < (dfframe_len - lag); n++) {
                tmp = dfframe[n] * dfframe[n + lag];
                sum += tmp;
            }
            acf[lag] = double(sum/ (dfframe_len - lag));
        }

        // now apply comb filtering
        int numelem = 4;

        for (int i = 2; i < rcf_len; i++) { // max beat period
            for (int a = 1; a <= numelem; a++) { // number of comb elements
                for (int b = 1-a; b <= a-1; b++) { // general state using normalisation of comb elements
                    rcf[i-1] += ( acf[(a*i+b)-1]*wv[i-1] ) / (2.*a-1.);     // calculate value for comb filter row
                }
            }
        }

        // apply adaptive threshold to rcf
        MathUtilities::adaptiveThreshold(rcf);

        double rcfsum =0.;
        for (int i = 0; i < rcf_len; i++) {
            rcf[i] += EPS ;
            rcfsum += rcf[i];
        }

        // normalise rcf to sum to unity
        for (int i = 0; i < rcf_len; i++) {
            rcf[i] /= (rcfsum + EPS);
        }
    }

    void
    TempoTrackV2::viterbi_decode(const d_mat_t &rcfmat, const d_vec_t &wv, d_vec_t &beat_period, d_vec_t &tempi)
    {
        // following Kevin Murphy's Viterbi decoding to get best path of
        // beat periods through rfcmat
        
        int wv_len = int(wv.size());
        
        // make transition matrix
        d_mat_t tmat;
        for (int i = 0; i < wv_len; i++) {
            tmat.push_back ( d_vec_t() ); // adds a new column
            for (int j = 0; j < wv_len; j++) {
                tmat[i].push_back(0.); // fill with zeros initially
            }
        }

        // variance of Gaussians in transition matrix
        // formed of Gaussians on diagonal - implies slow tempo change
        double sigma = 8.;
        // don't want really short beat periods, or really long ones
        for (int i = 20; i  < wv_len - 20; i++) {
            for (int j = 20; j < wv_len - 20; j++) {
                double mu = double(i);
                tmat[i][j] = exp( (-1.*pow((j-mu),2.)) / (2.*pow(sigma,2.)) );
            }
        }

        // parameters for Viterbi decoding... this part is taken from
        // Murphy's matlab

        d_mat_t delta;
        i_mat_t psi;
        for (int i = 0; i < int(rcfmat.size()); i++) {
            delta.push_back(d_vec_t());
            psi.push_back(i_vec_t());
            for (int j = 0; j < int(rcfmat[i].size()); j++) {
                delta[i].push_back(0.); // fill with zeros initially
                psi[i].push_back(0); // fill with zeros initially
            }
        }

        int T = int(delta.size());

        if (T < 2) return; // can't do anything at all meaningful

        int Q = int(delta[0].size());

        // initialize first column of delta
        for (int j = 0; j < Q; j++) {
            delta[0][j] = wv[j] * rcfmat[0][j];
            psi[0][j] = 0;
        }

        double deltasum = 0.;
        for (int i = 0; i < Q; i++) {
            deltasum += delta[0][i];
        }
        for (int i = 0; i < Q; i++) {
            delta[0][i] /= (deltasum + EPS);
        }

        for (int t=1; t < T; t++)
        {
            d_vec_t tmp_vec(Q);

            for (int j = 0; j < Q; j++) {
                for (int i = 0; i < Q; i++) {
                    tmp_vec[i] = delta[t-1][i] * tmat[j][i];
                }

                delta[t][j] = get_max_val(tmp_vec);

                psi[t][j] = get_max_ind(tmp_vec);

                delta[t][j] *= rcfmat[t][j];
            }

            // normalise current delta column
            double deltasum = 0.;
            for (int i = 0; i < Q; i++) {
                deltasum += delta[t][i];
            }
            for (int i = 0; i < Q; i++) {
                delta[t][i] /= (deltasum + EPS);
            }
        }

        i_vec_t bestpath(T);
        d_vec_t tmp_vec(Q);
        for (int i = 0; i < Q; i++) {
            tmp_vec[i] = delta[T-1][i];
        }

        // find starting point - best beat period for "last" frame
        bestpath[T-1] = get_max_ind(tmp_vec);

        // backtrace through index of maximum values in psi
        for (int t=T-2; t>0 ;t--) {
            bestpath[t] = psi[t+1][bestpath[t+1]];
        }

        // weird but necessary hack -- couldn't get above loop to terminate at t >= 0
        bestpath[0] = psi[1][bestpath[1]];

        int lastind = 0;
        for (int i = 0; i < T; i++) {
            int step = 128;
            for (int j = 0; j < step; j++) {
                lastind = i*step+j;
                beat_period[lastind] = bestpath[i];
            }
    //        std::cerr << "bestpath[" << i << "] = " << bestpath[i] << " (used for beat_periods " << i*step << " to " << i*step+step-1 << ")" << std::endl;
        }

        // fill in the last values...
        for (int i = lastind; i < int(beat_period.size()); i++) {
            beat_period[i] = beat_period[lastind];
        }

        for (int i = 0; i < int(beat_period.size()); i++) {
            tempi.push_back((60. * m_rate / m_increment)/beat_period[i]);
        }
    }

    double
    TempoTrackV2::get_max_val(const d_vec_t &df)
    {
        double maxval = 0.;
        int df_len = int(df.size());
        
        for (int i = 0; i < df_len; i++) {
            if (maxval < df[i]) {
                maxval = df[i];
            }
        }

        return maxval;
    }

    int
    TempoTrackV2::get_max_ind(const d_vec_t &df)
    {
        double maxval = 0.;
        int ind = 0;
        int df_len = int(df.size());
        
        for (int i = 0; i < df_len; i++) {
            if (maxval < df[i]) {
                maxval = df[i];
                ind = i;
            }
        }

        return ind;
    }

    void
    TempoTrackV2::normalise_vec(d_vec_t &df)
    {
        double sum = 0.;
        int df_len = int(df.size());
        
        for (int i = 0; i < df_len; i++) {
            sum += df[i];
        }

        for (int i = 0; i < df_len; i++) {
            df[i]/= (sum + EPS);
        }
    }

    // MEPD 28/11/12
    // this function has been updated to allow the "alpha" and "tightness" parameters
    // of the dynamic program to be set by the user
    // the default value of alpha = 0.9 and tightness = 4
    void
    TempoTrackV2::calculateBeats(const vector<double> &df,
                                const vector<double> &beat_period,
                                vector<double> &beats, double alpha, double tightness)
    {
        if (df.empty() || beat_period.empty()) return;

        int df_len = int(df.size());

        d_vec_t cumscore(df_len); // store cumulative score
        i_vec_t backlink(df_len); // backlink (stores best beat locations at each time instant)
        d_vec_t localscore(df_len); // localscore, for now this is the same as the detection function

        for (int i = 0; i < df_len; i++) {
            localscore[i] = df[i];
            backlink[i] = -1;
        }

        //double tightness = 4.;
        //double alpha = 0.9;
        // MEPD 28/11/12
        // debug statements that can be removed.
    //    std::cerr << "alpha" << alpha << std::endl;
    //    std::cerr << "tightness" << tightness << std::endl;

        // main loop
        for (int i = 0; i < df_len; i++) {
            
            int prange_min = -2*beat_period[i];
            int prange_max = round(-0.5*beat_period[i]);

            // transition range
            int txwt_len = prange_max - prange_min + 1;
            d_vec_t txwt (txwt_len);
            d_vec_t scorecands (txwt_len);

            for (int j = 0; j < txwt_len; j++) {
                
                double mu = double(beat_period[i]);
                txwt[j] = exp( -0.5*pow(tightness * log((round(2*mu)-j)/mu),2));

                // IF IN THE ALLOWED RANGE, THEN LOOK AT CUMSCORE[I+PRANGE_MIN+J
                // ELSE LEAVE AT DEFAULT VALUE FROM INITIALISATION:  D_VEC_T SCORECANDS (TXWT.SIZE());

                int cscore_ind = i + prange_min + j;
                if (cscore_ind >= 0) {
                    scorecands[j] = txwt[j] * cumscore[cscore_ind];
                }
            }

            // find max value and index of maximum value
            double vv = get_max_val(scorecands);
            int xx = get_max_ind(scorecands);

            cumscore[i] = alpha*vv + (1.-alpha)*localscore[i];
            backlink[i] = i+prange_min+xx;

    //        std::cerr << "backlink[" << i << "] <= " << backlink[i] << std::endl;
        }

        // STARTING POINT, I.E. LAST BEAT.. PICK A STRONG POINT IN cumscore VECTOR
        d_vec_t tmp_vec;
        for (int i = df_len - beat_period[beat_period.size()-1] ; i < df_len; i++) {
            tmp_vec.push_back(cumscore[i]);
        }

        int startpoint = get_max_ind(tmp_vec) +
            df_len - beat_period[beat_period.size()-1] ;

        // can happen if no results obtained earlier (e.g. input too short)
        if (startpoint >= int(backlink.size())) {
            startpoint = int(backlink.size()) - 1;
        }

        // USE BACKLINK TO GET EACH NEW BEAT (TOWARDS THE BEGINNING OF THE FILE)
        //  BACKTRACKING FROM THE END TO THE BEGINNING.. MAKING SURE NOT TO GO BEFORE SAMPLE 0
        i_vec_t ibeats;
        ibeats.push_back(startpoint);
    //    std::cerr << "startpoint = " << startpoint << std::endl;
        while (backlink[ibeats.back()] > 0) {
    //        std::cerr << "backlink[" << ibeats.back() << "] = " << backlink[ibeats.back()] << std::endl;
            int b = ibeats.back();
            if (backlink[b] == b) break; // shouldn't happen... haha
            ibeats.push_back(backlink[b]);
        }

        // REVERSE SEQUENCE OF IBEATS AND STORE AS BEATS
        for (int i = 0; i < int(ibeats.size()); i++) {
            beats.push_back(double(ibeats[ibeats.size() - i - 1]));
        }
    }    
}