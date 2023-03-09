#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <complex>
#include <cassert>

#define TWO_PI          (2. * M_PI)

#define EPS             2.2204e-016

/* aliases to math.h functions */
#define EXP                             exp
#define COS                             cos
#define SIN                             sin
#define ABS                             fabs
#define POW                             powf
#define SQRT                    sqrtf
#define LOG10                   log10f
#define LOG                             logf
#define FLOOR                   floorf
#define TRUNC                   truncf

typedef std::complex<double> ComplexData;

/* aliases to complex.h functions */
/** sample = EXPC(complex) */
#define EXPC                    cexpf
/** complex = CEXPC(complex) */
#define CEXPC                   cexp
/** sample = ARGC(complex) */
#define ARGC                    cargf
/** sample = ABSC(complex) norm */
#define ABSC                    cabsf
/** sample = REAL(complex) */
#define REAL                    crealf
/** sample = IMAG(complex) */
#define IMAG                    cimagf


#define ISNAN(x) (sizeof(x) == sizeof(double) ? ISNANd(x) : ISNANf(x))
static inline int ISNANf(float x) { return x != x; }
static inline int ISNANd(double x) { return x != x; }
          
#define ISINF(x) (sizeof(x) == sizeof(double) ? ISINFd(x) : ISINFf(x))
static inline int ISINFf(float x) { return !ISNANf(x) && ISNANf(x - x); }
static inline int ISINFd(double x) { return !ISNANd(x) && ISNANd(x - x); }


namespace qmdsp
{
    /**
    * Static helper functions for simple mathematical calculations.
    */
    class MathUtilities  
    {
    public: 
        /**
        * Round x to the nearest integer.
        */
        static double round( double x );

        /**
        * Return through min and max pointers the highest and lowest
        * values in the given array of the given length.
        */
        static void getFrameMinMax( const double* data, int len,
                                    double* min, double* max );

        /**
        * Return the mean of the given array of the given length.
        */
        static double mean( const double* src, int len );

        /**
        * Return the mean of the subset of the given vector identified by
        * start and count.
        */
        static double mean( const std::vector<double> &data,
                            int start, int count );
        
        /**
        * Return the sum of the values in the given array of the given
        * length.
        */
        static double sum( const double* src, int len );

        /**
        * Return the median of the values in the given array of the given
        * length. If the array is even in length, the returned value will
        * be half-way between the two values adjacent to median.
        */
        static double median( const double* src, int len );

        /**
        * The principle argument function. Map the phase angle ang into
        * the range [-pi,pi).
        */
        static double princarg( double ang );

        /**
        * Floating-point division modulus: return x % y.
        */
        static double mod( double x, double y);

        /**
        * The alpha norm is the alpha'th root of the mean alpha'th power
        * magnitude. For example if alpha = 2 this corresponds to the RMS
        * of the input data, and when alpha = 1 this is the mean
        * magnitude.
        */
        static void getAlphaNorm(const double *data, int len, int alpha, double* ANorm);

        /**
        * The alpha norm is the alpha'th root of the mean alpha'th power
        * magnitude. For example if alpha = 2 this corresponds to the RMS
        * of the input data, and when alpha = 1 this is the mean
        * magnitude.
        */
        static double getAlphaNorm(const std::vector <double> &data, int alpha );

        enum NormaliseType {
            NormaliseNone,
            NormaliseUnitSum,
            NormaliseUnitMax
        };

        static void normalise(double *data, int length,
                            NormaliseType n = NormaliseUnitMax);

        static void normalise(std::vector<double> &data,
                            NormaliseType n = NormaliseUnitMax);

        /**
        * Calculate the L^p norm of a vector. Equivalent to MATLAB's
        * norm(data, p).
        */
        static double getLpNorm(const std::vector<double> &data,
                                int p);

        /**
        * Normalise a vector by dividing through by its L^p norm. If the
        * norm is below the given threshold, the unit vector for that
        * norm is returned. p may be 0, in which case no normalisation
        * happens and the data is returned unchanged.
        */
        static std::vector<double> normaliseLp(const std::vector<double> &data,
                                            int p,
                                            double threshold = 1e-6);
        
        /**
        * Threshold the input/output vector data against a moving-mean
        * average filter.
        */
        static void adaptiveThreshold(std::vector<double> &data);

        static void circShift( double* data, int length, int shift);

        static int getMax( double* data, int length, double* max = 0 );
        static int getMax( const std::vector<double> &data, double* max = 0 );
        static int compareInt(const void * a, const void * b);

        /** 
        * Return true if x is 2^n for some integer n >= 0.
        */
        static bool isPowerOfTwo(int x);

        /**
        * Return the next higher integer power of two from x, e.g. 1300
        * -> 2048, 2048 -> 2048.
        */
        static int nextPowerOfTwo(int x);

        /**
        * Return the next lower integer power of two from x, e.g. 1300 ->
        * 1024, 2048 -> 2048.
        */
        static int previousPowerOfTwo(int x);

        /**
        * Return the nearest integer power of two to x, e.g. 1300 -> 1024,
        * 12 -> 16 (not 8; if two are equidistant, the higher is returned).
        */
        static int nearestPowerOfTwo(int x);

        /**
        * Return x!
        */
        static double factorial(int x); // returns double in case it is large

        /**
        * Return the greatest common divisor of natural numbers a and b.
        */
        static int gcd(int a, int b);
    };

    class Correlation  
    {
    public:
        Correlation();
        virtual ~Correlation();

        void doAutoUnBiased( double* src, double* dst, int length );
    };
    class CosineDistance
    {
    public:
        CosineDistance() { }
        ~CosineDistance() { }

        double distance(const std::vector<double> &v1,
                        const std::vector<double> &v2);

    protected:
        double dist, dDenTot, dDen1, dDen2, dSum1;
    };

    /**
    * Helper methods for calculating Kullback-Leibler divergences.
    */
    class KLDivergence
    {
    public:
        KLDivergence() { }
        ~KLDivergence() { }

        /**
        * Calculate a symmetrised Kullback-Leibler divergence of Gaussian
        * models based on mean and variance vectors.  All input vectors
        * must be of equal size.
        */
        double distanceGaussian(const std::vector<double> &means1,
                                const std::vector<double> &variances1,
                                const std::vector<double> &means2,
                                const std::vector<double> &variances2);

        /**
        * Calculate a Kullback-Leibler divergence of two probability
        * distributions.  Input vectors must be of equal size.  If
        * symmetrised is true, the result will be the symmetrised
        * distance (equal to KL(d1, d2) + KL(d2, d1)).
        */
        double distanceDistribution(const std::vector<double> &d1,
                                    const std::vector<double> &d2,
                                    bool symmetrised);
    };

    template <typename T>
    class MedianFilter
    {
    public:
        MedianFilter(int size, float percentile = 50.f) :
            m_size(size),
            m_frame(new T[size]),
            m_sorted(new T[size]),
            m_sortend(m_sorted + size - 1) {
            setPercentile(percentile);
            reset();
        }

        ~MedianFilter() { 
            delete[] m_frame;
            delete[] m_sorted;
        }

        void setPercentile(float p) {
            m_index = int((m_size * p) / 100.f);
            if (m_index >= m_size) m_index = m_size-1;
            if (m_index < 0) m_index = 0;
        }

        void push(T value) {
            if (value != value) {
                std::cerr << "WARNING: MedianFilter::push: attempt to push NaN, pushing zero instead" << std::endl;
                // we do need to push something, to maintain the filter length
                value = T();
            }
            drop(m_frame[0]);
            const int sz1 = m_size-1;
            for (int i = 0; i < sz1; ++i) m_frame[i] = m_frame[i+1];
            m_frame[m_size-1] = value;
            put(value);
        }

        T get() const {
            return m_sorted[m_index];
        }

        int getSize() const {
            return m_size; 
        }

        T getAt(float percentile) {
            int ix = int((m_size * percentile) / 100.f);
            if (ix >= m_size) ix = m_size-1;
            if (ix < 0) ix = 0;
            return m_sorted[ix];
        }

        void reset() {
            for (int i = 0; i < m_size; ++i) m_frame[i] = 0;
            for (int i = 0; i < m_size; ++i) m_sorted[i] = 0;
        }

        static std::vector<T> filter(int size, const std::vector<T> &in) {
            std::vector<T> out;
            MedianFilter<T> f(size);
            for (int i = 0; i < int(in.size()); ++i) {
                f.push(in[i]);
                T median = f.get();
                if (i >= size/2) out.push_back(median);
            }
            while (out.size() < in.size()) {
                f.push(T());
                out.push_back(f.get());
            }
            return out;
        }

    private:
        const int m_size;
        T *const m_frame;
        T *const m_sorted;
        T *const m_sortend;
        int m_index;

        void put(T value) {
            // precondition: m_sorted contains m_size-1 values, packed at start
            // postcondition: m_sorted contains m_size values, one of which is value
            T *point = std::lower_bound(m_sorted, m_sortend, value);
            const int n = m_sortend - point;
            for (int i = n; i > 0; --i) point[i] = point[i-1];
            *point = value;
        }

        void drop(T value) {
            // precondition: m_sorted contains m_size values, one of which is value
            // postcondition: m_sorted contains m_size-1 values, packed at start
            T *point = std::lower_bound(m_sorted, m_sortend + 1, value);
            if (*point != value) {
                std::cerr << "WARNING: MedianFilter::drop: *point is " << *point
                        << ", expected " << value << std::endl;
            }
            const int n = m_sortend - point;
            for (int i = 0; i < n; ++i) point[i] = point[i+1];
            *m_sortend = T(0);
        }

        MedianFilter(const MedianFilter &); // not provided
        MedianFilter &operator=(const MedianFilter &); // not provided
    };

    class TPolyFit
    {
        typedef std::vector<std::vector<double> > Matrix;
    public:

        static double PolyFit2 (const std::vector<double> &x,  // does the work
                                const std::vector<double> &y,
                                std::vector<double> &coef);

                    
    private:
        TPolyFit &operator = (const TPolyFit &);   // disable assignment
        TPolyFit();                                // and instantiation
        TPolyFit(const TPolyFit&);                 // and copying
    
        static void Square (const Matrix &x,              // Matrix multiplication routine
                            const std::vector<double> &y,
                            Matrix &a,                    // A = transpose X times X
                            std::vector<double> &g,         // G = Y times X
                            const int nrow, const int ncol);
        // Forms square coefficient matrix

        static bool GaussJordan (Matrix &b,                  // square matrix of coefficients
                                const std::vector<double> &y, // constant std::vector
                                std::vector<double> &coef);   // solution std::vector
        // returns false if matrix singular

        static bool GaussJordan2(Matrix &b,
                                const std::vector<double> &y,
                                Matrix &w,
                                std::vector<std::vector<int> > &index);
    };

    // some utility functions

    struct NSUtility
    {
        static void swap(double &a, double &b) {
            double t = a; a = b; b = t;
        }
        
        // fills a std::vector with zeros.
        static void zeroise(std::vector<double> &array, int n) { 
            array.clear();
            for(int j = 0; j < n; ++j) array.push_back(0);
        }
        
        // fills a std::vector with zeros.
        static void zeroise(std::vector<int> &array, int n) {
            array.clear();
            for(int j = 0; j < n; ++j) array.push_back(0);
        }
        
        // fills a (m by n) matrix with zeros.
        static void zeroise(std::vector<std::vector<double> > &matrix, int m, int n) {
            std::vector<double> zero;
            zeroise(zero, n);
            matrix.clear();
            for(int j = 0; j < m; ++j) matrix.push_back(zero);
        }
        
        // fills a (m by n) matrix with zeros.
        static void zeroise(std::vector<std::vector<int> > &matrix, int m, int n) {
            std::vector<int> zero;
            zeroise(zero, n);
            matrix.clear();
            for(int j = 0; j < m; ++j) matrix.push_back(zero);
        }
        
        static double sqr(const double &x) {return x * x;}
    };


    // main PolyFit routine

    double TPolyFit::PolyFit2 (const std::vector<double> &x,
                            const std::vector<double> &y,
                            std::vector<double> &coefs)
    // nterms = coefs.size()
    // npoints = x.size()
    {
        int i, j;
        double xi, yi, yc, srs, sum_y, sum_y2;
        Matrix xmatr;        // Data matrix
        Matrix a;
        std::vector<double> g;      // Constant std::vector
        const int npoints(x.size());
        const int nterms(coefs.size());
        double correl_coef;
        NSUtility::zeroise(g, nterms);
        NSUtility::zeroise(a, nterms, nterms);
        NSUtility::zeroise(xmatr, npoints, nterms);
        if (nterms < 1) {
            std::cerr << "ERROR: PolyFit called with less than one term" << std::endl;
            return 0;
        }
        if(npoints < 2) {
            std::cerr << "ERROR: PolyFit called with less than two points" << std::endl;
            return 0;
        }
        if(npoints != (int)y.size()) {
            std::cerr << "ERROR: PolyFit called with x and y of unequal size" << std::endl;
            return 0;
        }
        for(i = 0; i < npoints; ++i) {
            //      { setup x matrix }
            xi = x[i];
            xmatr[i][0] = 1.0;         //     { first column }
            for(j = 1; j < nterms; ++j)
                xmatr[i][j] = xmatr [i][j - 1] * xi;
        }
        Square (xmatr, y, a, g, npoints, nterms);
        if(!GaussJordan (a, g, coefs)) {
            return -1;
        }
        sum_y = 0.0;
        sum_y2 = 0.0;
        srs = 0.0;
        for(i = 0; i < npoints; ++i) {
            yi = y[i];
            yc = 0.0;
            for(j = 0; j < nterms; ++j) {
                yc += coefs [j] * xmatr [i][j];
            }
            srs += NSUtility::sqr (yc - yi);
            sum_y += yi;
            sum_y2 += yi * yi;
        }

        // If all Y values are the same, avoid dividing by zero
        correl_coef = sum_y2 - NSUtility::sqr (sum_y) / npoints;
        // Either return 0 or the correct value of correlation coefficient
        if (correl_coef != 0) {
            correl_coef = srs / correl_coef;
        }
        if (correl_coef >= 1) {
            correl_coef = 0.0;
        } else {
            correl_coef = sqrt (1.0 - correl_coef);
        }
        return correl_coef;
    }


    //------------------------------------------------------------------------

    // Matrix multiplication routine
    // A = transpose X times X
    // G = Y times X

    // Form square coefficient matrix

    void TPolyFit::Square (const Matrix &x,
                        const std::vector<double> &y,
                        Matrix &a,
                        std::vector<double> &g,
                        const int nrow,
                        const int ncol)
    {
        int i, k, l;
        for(k = 0; k < ncol; ++k) {
            for(l = 0; l < k + 1; ++l) {
                a [k][l] = 0.0;
                for(i = 0; i < nrow; ++i) {
                    a[k][l] += x[i][l] * x [i][k];
                    if(k != l) {
                        a[l][k] = a[k][l];
                    }
                }
            }
            g[k] = 0.0;
            for(i = 0; i < nrow; ++i) {
                g[k] += y[i] * x[i][k];
            }
        }
    }
    //---------------------------------------------------------------------------------


    bool TPolyFit::GaussJordan (Matrix &b,
                                const std::vector<double> &y,
                                std::vector<double> &coef)
    //b square matrix of coefficients
    //y constant std::vector
    //coef solution std::vector
    //ncol order of matrix got from b.size()


    {
    /*
    { Gauss Jordan matrix inversion and solution }
    { B (n, n) coefficient matrix becomes inverse }
    { Y (n) original constant std::vector }
    { W (n, m) constant std::vector(s) become solution std::vector }
    { DETERM is the determinant }
    { ERROR = 1 if singular }
    { INDEX (n, 3) }
    { NV is number of constant vectors }
    */

        int ncol(b.size());
        int irow, icol;
        std::vector<std::vector<int> >index;
        Matrix w;

        NSUtility::zeroise(w, ncol, ncol);
        NSUtility::zeroise(index, ncol, 3);

        if (!GaussJordan2(b, y, w, index)) {
            return false;
        }

        // Interchange columns
        int m;
        for (int i = 0; i <  ncol; ++i) {
            m = ncol - i - 1;
            if(index [m][0] != index [m][1]) {
                irow = index [m][0];
                icol = index [m][1];
                for(int k = 0; k < ncol; ++k) {
                    NSUtility::swap (b[k][irow], b[k][icol]);
                }
            }
        }

        for(int k = 0; k < ncol; ++k) {
            if(index [k][2] != 0) {
                std::cerr << "ERROR: Error in PolyFit::GaussJordan: matrix is singular" << std::endl;
                return false;
            }
        }

        for( int i = 0; i < ncol; ++i) {
            coef[i] = w[i][0];
        }
    
        return true;
    }   // end;     { procedure GaussJordan }
    //----------------------------------------------------------------------------------------------


    bool TPolyFit::GaussJordan2(Matrix &b,
                                const std::vector<double> &y,
                                Matrix &w,
                                std::vector<std::vector<int> > &index)
    {
        //GaussJordan2;         // first half of GaussJordan
        // actual start of gaussj
    
        double big, t;
        double pivot;
        double determ;
        int irow = 0, icol = 0;
        int ncol(b.size());
        int nv = 1;                  // single constant std::vector

        for(int i = 0; i < ncol; ++i) {
            w[i][0] = y[i];      // copy constant std::vector
            index[i][2] = -1;
        }

        determ = 1.0;

        for (int i = 0; i < ncol; ++i) {
            // Search for largest element
            big = 0.0;

            for (int j = 0; j < ncol; ++j) {
                if (index[j][2] != 0) {
                    for (int k = 0; k < ncol; ++k) {
                        if (index[k][2] > 0) {
                            std::cerr << "ERROR: Error in PolyFit::GaussJordan2: matrix is singular" << std::endl;
                            return false;
                        }

                        if (index[k][2] < 0 && fabs(b[j][k]) > big) {
                            irow = j;
                            icol = k;
                            big = fabs(b[j][k]);
                        }
                    } //   { k-loop }
                }
            }  // { j-loop }
            
            index [icol][2] = index [icol][2] + 1;
            index [i][0] = irow;
            index [i][1] = icol;

            // Interchange rows to put pivot on diagonal
            // GJ3
            if (irow != icol) {
                determ = -determ;
                for (int m = 0; m < ncol; ++m) {
                    NSUtility::swap (b [irow][m], b[icol][m]);
                }
                if (nv > 0) {
                    for (int m = 0; m < nv; ++m) {
                        NSUtility::swap (w[irow][m], w[icol][m]);
                    }
                }
            } // end GJ3

            // divide pivot row by pivot column
            pivot = b[icol][icol];
            determ *= pivot;
            b[icol][icol] = 1.0;

            for (int m = 0; m < ncol; ++m) {
                b[icol][m] /= pivot;
            }
            if (nv > 0) {
                for (int m = 0; m < nv; ++m) {
                    w[icol][m] /= pivot;
                }
            }

            // Reduce nonpivot rows
            for (int n = 0; n < ncol; ++n) {
                if (n != icol) {
                    t = b[n][icol];
                    b[n][icol] = 0.0;
                    for (int m = 0; m < ncol; ++m) {
                        b[n][m] -= b[icol][m] * t;
                    }
                    if (nv > 0) {
                        for (int m = 0; m < nv; ++m) {
                            w[n][m] -= w[icol][m] * t;
                        }
                    }
                }
            }
        } // { i-loop }
        
        return true;
    }
    /**
    * Kaiser window: A windower whose bandwidth and sidelobe height
    * (signal-noise ratio) can be specified. These parameters are traded
    * off against the window length.
    */
    class KaiserWindow
    {
    public:
        struct Parameters {
            int length;
            double beta;
        };

        /**
        * Construct a Kaiser windower with the given length and beta
        * parameter.
        */
        KaiserWindow(Parameters p) : m_length(p.length), m_beta(p.beta) { init(); }

        /**
        * Construct a Kaiser windower with the given attenuation in dB
        * and transition width in samples.
        */
        static KaiserWindow byTransitionWidth(double attenuation,
                                            double transition) {
            return KaiserWindow
                (parametersForTransitionWidth(attenuation, transition));
        }

        /**
        * Construct a Kaiser windower with the given attenuation in dB
        * and transition bandwidth in Hz for the given samplerate.
        */
        static KaiserWindow byBandwidth(double attenuation,
                                        double bandwidth,
                                        double samplerate) {
            return KaiserWindow
                (parametersForBandwidth(attenuation, bandwidth, samplerate));
        }

        /**
        * Obtain the parameters necessary for a Kaiser window of the
        * given attenuation in dB and transition width in samples.
        */
        static Parameters parametersForTransitionWidth(double attenuation,
                                                    double transition);

        /**
        * Obtain the parameters necessary for a Kaiser window of the
        * given attenuation in dB and transition bandwidth in Hz for the
        * given samplerate.
        */
        static Parameters parametersForBandwidth(double attenuation,
                                                double bandwidth,
                                                double samplerate) {
            return parametersForTransitionWidth
                (attenuation, (bandwidth * 2 * M_PI) / samplerate);
        } 

        int getLength() const {
            return m_length;
        }

        const double *getWindow() const { 
            return m_window.data();
        }

        void cut(double *src) const { 
            cut(src, src); 
        }

        void cut(const double *src, double *dst) const {
            for (int i = 0; i < m_length; ++i) {
                dst[i] = src[i] * m_window[i];
            }
        }

    private:
        int m_length;
        double m_beta;
        std::vector<double> m_window;

        void init();
    };

    /**
    * A window containing values of the sinc function, i.e. sin(x)/x with
    * sinc(0) == 1, with x == 0 at the centre.
    */
    class SincWindow
    {
    public:
        /**
        * Construct a windower of the given length, containing the values
        * of sinc(x) with x=0 in the middle, i.e. at sample (length-1)/2
        * for odd or (length/2)+1 for even length, such that the distance
        * from -pi to pi (the nearest zero crossings either side of the
        * peak) is p samples.
        */
        SincWindow(int length, double p) : m_length(length), m_p(p) { init(); }

        int getLength() const {
            return m_length;
        }

        const double *getWindow() const { 
            return m_window.data();
        }

        void cut(double *src) const { 
            cut(src, src); 
        }

        void cut(const double *src, double *dst) const {
            for (int i = 0; i < m_length; ++i) {
                dst[i] = src[i] * m_window[i];
            }
        }

    private:
        int m_length;
        double m_p;
        std::vector<double> m_window;

        void init();
    };

    enum WindowType {
        RectangularWindow,
        BartlettWindow,
        HammingWindow,
        HanningWindow,
        BlackmanWindow,
        BlackmanHarrisWindow,

        FirstWindow = RectangularWindow,
        LastWindow = BlackmanHarrisWindow
    };

    /**
    * Various shaped windows for sample frame conditioning, including
    * cosine windows (Hann etc) and triangular and rectangular windows.
    */
    template <typename T>
    class Window
    {
    public:
        /**
        * Construct a windower of the given type and size. 
        *
        * Note that the cosine windows are periodic by design, rather
        * than symmetrical. (A window of size N is equivalent to a
        * symmetrical window of size N+1 with the final element missing.)
        */
        Window(WindowType type, int size) : m_type(type), m_size(size) { encache(); }
        Window(const Window &w) : m_type(w.m_type), m_size(w.m_size) { encache(); }
        Window &operator=(const Window &w) {
            if (&w == this) return *this;
            m_type = w.m_type;
            m_size = w.m_size;
            encache();
            return *this;
        }
        virtual ~Window() { delete[] m_cache; }
        
        void cut(T *src) const { cut(src, src); }
        void cut(const T *src, T *dst) const {
            for (int i = 0; i < m_size; ++i) {
                dst[i] = src[i] * m_cache[i];
            }
        }

        WindowType getType() const { return m_type; }
        int getSize() const { return m_size; }

        std::vector<T> getWindowData() const {
            std::vector<T> d;
            for (int i = 0; i < m_size; ++i) {
                d.push_back(m_cache[i]);
            }
            return d;
        }

    protected:
        WindowType m_type;
        int m_size;
        T *m_cache;
        
        void encache();
    };

    template <typename T>
    void Window<T>::encache()
    {
        int n = m_size;
        T *mult = new T[n];
        int i;
        for (i = 0; i < n; ++i) mult[i] = 1.0;

        switch (m_type) {
                    
        case RectangularWindow:
            for (i = 0; i < n; ++i) {
                mult[i] = mult[i] * 0.5;
            }
            break;
                
        case BartlettWindow:
            if (n == 2) {
                mult[0] = mult[1] = 0; // "matlab compatible"
            } else if (n == 3) {
                mult[0] = 0;
                mult[1] = mult[2] = 2./3.;
            } else if (n > 3) {
                for (i = 0; i < n/2; ++i) {
                    mult[i] = mult[i] * (i / T(n/2));
                    mult[i + n - n/2] = mult[i + n - n/2] * (1.0 - (i / T(n/2)));
                }
            }
            break;
                
        case HammingWindow:
            if (n > 1) {
                for (i = 0; i < n; ++i) {
                    mult[i] = mult[i] * (0.54 - 0.46 * cos(2 * M_PI * i / n));
                }
            }
            break;
                
        case HanningWindow:
            if (n > 1) {
                for (i = 0; i < n; ++i) {
                    mult[i] = mult[i] * (0.50 - 0.50 * cos(2 * M_PI * i / n));
                }
            }
            break;
                
        case BlackmanWindow:
            if (n > 1) {
                for (i = 0; i < n; ++i) {
                    mult[i] = mult[i] * (0.42 - 0.50 * cos(2 * M_PI * i / n)
                                        + 0.08 * cos(4 * M_PI * i / n));
                }
            }
            break;
                
        case BlackmanHarrisWindow:
            if (n > 1) {
                for (i = 0; i < n; ++i) {
                    mult[i] = mult[i] * (0.35875
                                        - 0.48829 * cos(2 * M_PI * i / n)
                                        + 0.14128 * cos(4 * M_PI * i / n)
                                        - 0.01168 * cos(6 * M_PI * i / n));
                }
            }
            break;
        }
            
        m_cache = mult;
    }
    /**
    * Convert between musical pitch (i.e. MIDI pitch number) and
    * fundamental frequency.
    */
    class Pitch
    {
    public:
        static float getFrequencyForPitch(int midiPitch,
                                        float centsOffset = 0,
                                        float concertA = 440.0);

        static int getPitchForFrequency(float frequency,
                                        float *centsOffsetReturn = 0,
                                        float concertA = 440.0);
    };


    #define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )

    /**  Variance-covariance matrix: creation  *****************************/

    /* Create m * m covariance matrix from given n * m data matrix. */
    void covcol(double** data, int n, int m, double** symmat)
    {
        double *mean;
        int i, j, j1, j2;

    /* Allocate storage for mean vector */

        mean = (double*) malloc(m*sizeof(double));

    /* Determine mean of column vectors of input data matrix */

        for (j = 0; j < m; j++)
        {
            mean[j] = 0.0;
            for (i = 0; i < n; i++)
            {
                mean[j] += data[i][j];
            }
            mean[j] /= (double)n;
        }

    /*
    printf("\nMeans of column vectors:\n");
    for (j = 0; j < m; j++)  {
    printf("%12.1f",mean[j]);  }   printf("\n");
    */

    /* Center the column vectors. */

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                data[i][j] -= mean[j];
            }
        }

    /* Calculate the m * m covariance matrix. */
        for (j1 = 0; j1 < m; j1++)
        {
            for (j2 = j1; j2 < m; j2++)
            {
                symmat[j1][j2] = 0.0;
                for (i = 0; i < n; i++)
                {
                    symmat[j1][j2] += data[i][j1] * data[i][j2];
                }
                symmat[j2][j1] = symmat[j1][j2];
            }
        }

        free(mean);

        return;

    }

    /**  Error handler  **************************************************/

    void erhand(char* err_msg)
    {
        fprintf(stderr,"Run-time error:\n");
        fprintf(stderr,"%s\n", err_msg);
        fprintf(stderr,"Exiting to system.\n");
        exit(1);
    }


    /**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */

    /* Householder reduction of matrix a to tridiagonal form.
    Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
    Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
    Springer-Verlag, 1976, pp. 489-494.
    W H Press et al., Numerical Recipes in C, Cambridge U P,
    1988, pp. 373-374.  */
    void tred2(double** a, int n, double* d, double* e)
    {
        int l, k, j, i;
        double scale, hh, h, g, f;
            
        for (i = n-1; i >= 1; i--)
        {
            l = i - 1;
            h = scale = 0.0;
            if (l > 0)
            {
                for (k = 0; k <= l; k++)
                    scale += fabs(a[i][k]);
                if (scale == 0.0)
                    e[i] = a[i][l];
                else
                {
                    for (k = 0; k <= l; k++)
                    {
                        a[i][k] /= scale;
                        h += a[i][k] * a[i][k];
                    }
                    f = a[i][l];
                    g = f>0 ? -sqrt(h) : sqrt(h);
                    e[i] = scale * g;
                    h -= f * g;
                    a[i][l] = f - g;
                    f = 0.0;
                    for (j = 0; j <= l; j++)
                    {
                        a[j][i] = a[i][j]/h;
                        g = 0.0;
                        for (k = 0; k <= j; k++)
                            g += a[j][k] * a[i][k];
                        for (k = j+1; k <= l; k++)
                            g += a[k][j] * a[i][k];
                        e[j] = g / h;
                        f += e[j] * a[i][j];
                    }
                    hh = f / (h + h);
                    for (j = 0; j <= l; j++)
                    {
                        f = a[i][j];
                        e[j] = g = e[j] - hh * f;
                        for (k = 0; k <= j; k++)
                            a[j][k] -= (f * e[k] + g * a[i][k]);
                    }
                }
            }
            else
                e[i] = a[i][l];
            d[i] = h;
        }
        d[0] = 0.0;
        e[0] = 0.0;
        for (i = 0; i < n; i++)
        {
            l = i - 1;
            if (d[i])
            {
                for (j = 0; j <= l; j++)
                {
                    g = 0.0;
                    for (k = 0; k <= l; k++)
                        g += a[i][k] * a[k][j];
                    for (k = 0; k <= l; k++)
                        a[k][j] -= g * a[k][i];
                }
            }
            d[i] = a[i][i];
            a[i][i] = 1.0;
            for (j = 0; j <= l; j++)
                a[j][i] = a[i][j] = 0.0;
        }
    }

    /**  Tridiagonal QL algorithm -- Implicit  **********************/

    void tqli(double* d, double* e, int n, double** z)
    {
        int m, l, iter, i, k;
        double s, r, p, g, f, dd, c, b;
            
        for (i = 1; i < n; i++)
            e[i-1] = e[i];
        e[n-1] = 0.0;
        for (l = 0; l < n; l++)
        {
            iter = 0;
            do
            {
                for (m = l; m < n-1; m++)
                {
                    dd = fabs(d[m]) + fabs(d[m+1]);
                    if (fabs(e[m]) + dd == dd) break;
                }
                if (m != l)
                {
                    if (iter++ == 30) erhand("No convergence in TLQI.");
                    g = (d[l+1] - d[l]) / (2.0 * e[l]);
                    r = sqrt((g * g) + 1.0);
                    g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                    s = c = 1.0;
                    p = 0.0;
                    for (i = m-1; i >= l; i--)
                    {
                        f = s * e[i];
                        b = c * e[i];
                        if (fabs(f) >= fabs(g))
                        {
                            c = g / f;
                            r = sqrt((c * c) + 1.0);
                            e[i+1] = f * r;
                            c *= (s = 1.0/r);
                        }
                        else
                        {
                            s = f / g;
                            r = sqrt((s * s) + 1.0);
                            e[i+1] = g * r;
                            s *= (c = 1.0/r);
                        }
                        g = d[i+1] - p;
                        r = (d[i] - g) * s + 2.0 * c * b;
                        p = s * r;
                        d[i+1] = g + p;
                        g = c * r - b;
                        for (k = 0; k < n; k++)
                        {
                            f = z[k][i+1];
                            z[k][i+1] = s * z[k][i] + c * f;
                            z[k][i] = c * z[k][i] - s * f;
                        }
                    }
                    d[l] = d[l] - p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            }  while (m != l);
        }
    }

    /* In place projection onto basis vectors */
    void pca_project(double** data, int n, int m, int ncomponents)
    {
        int  i, j, k, k2;
        double  **symmat, /* **symmat2, */ *evals, *interm;
            
        //TODO: assert ncomponents < m
            
        symmat = (double**) malloc(m*sizeof(double*));
        for (i = 0; i < m; i++)
            symmat[i] = (double*) malloc(m*sizeof(double));
                    
        covcol(data, n, m, symmat);
            
        /*********************************************************************
                    Eigen-reduction
        **********************************************************************/
            
        /* Allocate storage for dummy and new vectors. */
        evals = (double*) malloc(m*sizeof(double));     /* Storage alloc. for vector of eigenvalues */
        interm = (double*) malloc(m*sizeof(double));    /* Storage alloc. for 'intermediate' vector */
        //MALLOC_ARRAY(symmat2,m,m,double);    
        //for (i = 0; i < m; i++) {
        //      for (j = 0; j < m; j++) {
        //              symmat2[i][j] = symmat[i][j]; /* Needed below for col. projections */
        //      }
        //}
        tred2(symmat, m, evals, interm);  /* Triangular decomposition */
        tqli(evals, interm, m, symmat);   /* Reduction of sym. trid. matrix */
    /* evals now contains the eigenvalues,
    columns of symmat now contain the associated eigenvectors. */   

    /*
    printf("\nEigenvalues:\n");
    for (j = m-1; j >= 0; j--) {
    printf("%18.5f\n", evals[j]); }
    printf("\n(Eigenvalues should be strictly positive; limited\n");
    printf("precision machine arithmetic may affect this.\n");
    printf("Eigenvalues are often expressed as cumulative\n");
    printf("percentages, representing the 'percentage variance\n");
    printf("explained' by the associated axis or principal component.)\n");
            
    printf("\nEigenvectors:\n");
    printf("(First three; their definition in terms of original vbes.)\n");
    for (j = 0; j < m; j++) {
    for (i = 1; i <= 3; i++)  {
    printf("%12.4f", symmat[j][m-i]);  }
    printf("\n");  }
    */

    /* Form projections of row-points on prin. components. */
    /* Store in 'data', overwriting original data. */
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                interm[j] = data[i][j]; }   /* data[i][j] will be overwritten */
            for (k = 0; k < ncomponents; k++) {
                data[i][k] = 0.0;
                for (k2 = 0; k2 < m; k2++) {
                    data[i][k] += interm[k2] * symmat[k2][m-k-1]; }
            }
        }

    /*      
            printf("\nProjections of row-points on first 3 prin. comps.:\n");
            for (i = 0; i < n; i++) {
            for (j = 0; j < 3; j++)  {
            printf("%12.4f", data[i][j]);  }
            printf("\n");  }
    */

    /* Form projections of col.-points on first three prin. components. */
    /* Store in 'symmat2', overwriting what was stored in this. */
    //for (j = 0; j < m; j++) {
    //       for (k = 0; k < m; k++) {
    //               interm[k] = symmat2[j][k]; }  /*symmat2[j][k] will be overwritten*/
    //  for (i = 0; i < 3; i++) {
    //      symmat2[j][i] = 0.0;
    //              for (k2 = 0; k2 < m; k2++) {
    //                      symmat2[j][i] += interm[k2] * symmat[k2][m-i-1]; }
    //              if (evals[m-i-1] > 0.0005)   /* Guard against zero eigenvalue */
    //                      symmat2[j][i] /= sqrt(evals[m-i-1]);   /* Rescale */
    //              else
    //                      symmat2[j][i] = 0.0;    /* Standard kludge */
    //    }
    // }

    /*
    printf("\nProjections of column-points on first 3 prin. comps.:\n");
    for (j = 0; j < m; j++) {
    for (k = 0; k < 3; k++)  {
    printf("%12.4f", symmat2[j][k]);  }
    printf("\n");  }
    */


        for (i = 0; i < m; i++)
            free(symmat[i]);
        free(symmat);
    //FREE_ARRAY(symmat2,m);
        free(evals);
        free(interm);

    }



    Correlation::Correlation()
    {

    }

    Correlation::~Correlation()
    {

    }

    void Correlation::doAutoUnBiased(double *src, double *dst, int length)
    {
        double tmp = 0.0;
        double outVal = 0.0;

        int i, j;

        for (i = 0; i < length; i++) {
            for (j = i; j < length; j++) {
                tmp += src[ j-i ] * src[ j ]; 
            }

            outVal = tmp / ( length - i );

            if (outVal <= 0) {
                dst[ i ] = EPS;
            } else {
                dst[ i ] = outVal;
            }
            
            tmp = 0.0;
        }
    }


    double CosineDistance::distance(const std::vector<double> &v1,
                                const std::vector<double> &v2)
    {
        dist = 1.0; dDenTot = 0; dDen1 = 0; dDen2 = 0; dSum1 =0;
        double small = 1e-20;

        //check if v1, v2 same size
        if (v1.size() != v2.size())
        {
            std::cerr << "CosineDistance::distance: ERROR: vectors not the same size\n";
            return 1.0;
        }
        else
        {
            for(int i=0; i<int(v1.size()); i++)
            {
                dSum1 += v1[i]*v2[i];
                dDen1 += v1[i]*v1[i];
                dDen2 += v2[i]*v2[i];
            }
            dDenTot = sqrt(fabs(dDen1*dDen2)) + small;
            dist = 1-((dSum1)/dDenTot);
            return dist;
        }
    }

    double KLDivergence::distanceGaussian(const std::vector<double> &m1,
                                      const std::vector<double> &v1,
                                      const std::vector<double> &m2,
                                      const std::vector<double> &v2)
    {
        int sz = m1.size();

        double d = -2.0 * sz;
        double small = 1e-20;

        for (int k = 0; k < sz; ++k) {

            double kv1 = v1[k] + small;
            double kv2 = v2[k] + small;
            double km = (m1[k] - m2[k]) + small;

            d += kv1 / kv2 + kv2 / kv1;
            d += km * (1.0 / kv1 + 1.0 / kv2) * km;
        }

        d /= 2.0;

        return d;
    }

    double KLDivergence::distanceDistribution(const std::vector<double> &d1,
                                            const std::vector<double> &d2,
                                            bool symmetrised)
    {
        int sz = d1.size();

        double d = 0;
        double small = 1e-20;
        
        for (int i = 0; i < sz; ++i) {
            d += d1[i] * log10((d1[i] + small) / (d2[i] + small));
        }

        if (symmetrised) {
            d += distanceDistribution(d2, d1, false);
        }

        return d;
    }


    double MathUtilities::mod(double x, double y)
    {
        double a = floor( x / y );

        double b = x - ( y * a );
        return b;
    }

    double MathUtilities::princarg(double ang)
    {
        double ValOut;

        ValOut = mod( ang + M_PI, -2 * M_PI ) + M_PI;

        return ValOut;
    }

    void MathUtilities::getAlphaNorm(const double *data, int len, int alpha, double* ANorm)
    {
        int i;
        double temp = 0.0;
        double a=0.0;
            
        for( i = 0; i < len; i++) {
            temp = data[ i ];
            a  += ::pow( fabs(temp), double(alpha) );
        }
        a /= ( double )len;
        a = ::pow( a, ( 1.0 / (double) alpha ) );

        *ANorm = a;
    }

    double MathUtilities::getAlphaNorm( const vector <double> &data, int alpha )
    {
        int i;
        int len = data.size();
        double temp = 0.0;
        double a=0.0;
            
        for( i = 0; i < len; i++) {
            temp = data[ i ];
            a  += ::pow( fabs(temp), double(alpha) );
        }
        a /= ( double )len;
        a = ::pow( a, ( 1.0 / (double) alpha ) );

        return a;
    }

    double MathUtilities::round(double x)
    {
        if (x < 0) {
            return -floor(-x + 0.5);
        } else {
            return floor(x + 0.5);
        }
    }

    double MathUtilities::median(const double *src, int len)
    {
        if (len == 0) return 0;
        
        vector<double> scratch;
        for (int i = 0; i < len; ++i) scratch.push_back(src[i]);
        sort(scratch.begin(), scratch.end());

        int middle = len/2;
        if (len % 2 == 0) {
            return (scratch[middle] + scratch[middle - 1]) / 2;
        } else {
            return scratch[middle];
        }
    }

    double MathUtilities::sum(const double *src, int len)
    {
        int i ;
        double retVal =0.0;

        for(  i = 0; i < len; i++) {
            retVal += src[ i ];
        }

        return retVal;
    }

    double MathUtilities::mean(const double *src, int len)
    {
        double retVal =0.0;

        if (len == 0) return 0;

        double s = sum( src, len );
            
        retVal =  s  / (double)len;

        return retVal;
    }

    double MathUtilities::mean(const vector<double> &src,
                            int start,
                            int count)
    {
        double sum = 0.;
            
        if (count == 0) return 0;
        
        for (int i = 0; i < (int)count; ++i) {
            sum += src[start + i];
        }

        return sum / count;
    }

    void MathUtilities::getFrameMinMax(const double *data, int len, double *min, double *max)
    {
        int i;
        double temp = 0.0;

        if (len == 0) {
            *min = *max = 0;
            return;
        }
            
        *min = data[0];
        *max = data[0];

        for( i = 0; i < len; i++) {
            temp = data[ i ];

            if( temp < *min ) {
                *min =  temp ;
            }
            if( temp > *max ) {
                *max =  temp ;
            }
        }
    }

    int MathUtilities::getMax( double* pData, int Length, double* pMax )
    {
        int index = 0;
        int i;
        double temp = 0.0;
            
        double max = pData[0];

        for( i = 0; i < Length; i++) {
            temp = pData[ i ];

            if( temp > max ) {
                max =  temp ;
                index = i;
            }
        }

        if (pMax) *pMax = max;


        return index;
    }

    int MathUtilities::getMax( const std::vector<double> & data, double* pMax )
    {
        int index = 0;
        int i;
        double temp = 0.0;
            
        double max = data[0];

        for( i = 0; i < int(data.size()); i++) {

            temp = data[ i ];

            if( temp > max ) {
                max =  temp ;
                index = i;
            }
        }

        if (pMax) *pMax = max;


        return index;
    }

    void MathUtilities::circShift( double* pData, int length, int shift)
    {
        shift = shift % length;
        double temp;
        int i,n;

        for( i = 0; i < shift; i++) {
            
            temp=*(pData + length - 1);

            for( n = length-2; n >= 0; n--) {
                *(pData+n+1)=*(pData+n);
            }

            *pData = temp;
        }
    }

    int MathUtilities::compareInt (const void * a, const void * b)
    {
        return ( *(int*)a - *(int*)b );
    }

    void MathUtilities::normalise(double *data, int length, NormaliseType type)
    {
        switch (type) {

        case NormaliseNone: return;

        case NormaliseUnitSum:
        {
            double sum = 0.0;
            for (int i = 0; i < length; ++i) {
                sum += data[i];
            }
            if (sum != 0.0) {
                for (int i = 0; i < length; ++i) {
                    data[i] /= sum;
                }
            }
        }
        break;

        case NormaliseUnitMax:
        {
            double max = 0.0;
            for (int i = 0; i < length; ++i) {
                if (fabs(data[i]) > max) {
                    max = fabs(data[i]);
                }
            }
            if (max != 0.0) {
                for (int i = 0; i < length; ++i) {
                    data[i] /= max;
                }
            }
        }
        break;

        }
    }

    void MathUtilities::normalise(std::vector<double> &data, NormaliseType type)
    {
        switch (type) {

        case NormaliseNone: return;

        case NormaliseUnitSum:
        {
            double sum = 0.0;
            for (int i = 0; i < (int)data.size(); ++i) sum += data[i];
            if (sum != 0.0) {
                for (int i = 0; i < (int)data.size(); ++i) data[i] /= sum;
            }
        }
        break;

        case NormaliseUnitMax:
        {
            double max = 0.0;
            for (int i = 0; i < (int)data.size(); ++i) {
                if (fabs(data[i]) > max) max = fabs(data[i]);
            }
            if (max != 0.0) {
                for (int i = 0; i < (int)data.size(); ++i) data[i] /= max;
            }
        }
        break;

        }
    }

    double MathUtilities::getLpNorm(const std::vector<double> &data, int p)
    {
        double tot = 0.0;
        for (int i = 0; i < int(data.size()); ++i) {
            tot += abs(pow(data[i], p));
        }
        return pow(tot, 1.0 / p);
    }

    vector<double> MathUtilities::normaliseLp(const vector<double> &data,
                                            int p,
                                            double threshold)
    {
        int n = int(data.size());
        if (n == 0 || p == 0) return data;
        double norm = getLpNorm(data, p);
        if (norm < threshold) {
            return vector<double>(n, 1.0 / pow(n, 1.0 / p)); // unit vector
        }
        vector<double> out(n);
        for (int i = 0; i < n; ++i) {
            out[i] = data[i] / norm;
        }
        return out;
    }
        
    void MathUtilities::adaptiveThreshold(std::vector<double> &data)
    {
        int sz = int(data.size());
        if (sz == 0) return;

        vector<double> smoothed(sz);
            
        int p_pre = 8;
        int p_post = 7;

        for (int i = 0; i < sz; ++i) {

            int first = max(0,      i - p_pre);
            int last  = min(sz - 1, i + p_post);

            smoothed[i] = mean(data, first, last - first + 1);
        }

        for (int i = 0; i < sz; i++) {
            data[i] -= smoothed[i];
            if (data[i] < 0.0) data[i] = 0.0;
        }
    }

    bool
    MathUtilities::isPowerOfTwo(int x)
    {
        if (x < 1) return false;
        if (x & (x-1)) return false;
        return true;
    }

    int
    MathUtilities::nextPowerOfTwo(int x)
    {
        if (isPowerOfTwo(x)) return x;
        if (x < 1) return 1;
        int n = 1;
        while (x) { x >>= 1; n <<= 1; }
        return n;
    }

    int
    MathUtilities::previousPowerOfTwo(int x)
    {
        if (isPowerOfTwo(x)) return x;
        if (x < 1) return 1;
        int n = 1;
        x >>= 1;
        while (x) { x >>= 1; n <<= 1; }
        return n;
    }

    int
    MathUtilities::nearestPowerOfTwo(int x)
    {
        if (isPowerOfTwo(x)) return x;
        int n0 = previousPowerOfTwo(x);
        int n1 = nextPowerOfTwo(x);
        if (x - n0 < n1 - x) return n0;
        else return n1;
    }

    double
    MathUtilities::factorial(int x)
    {
        if (x < 0) return 0;
        double f = 1;
        for (int i = 1; i <= x; ++i) {
            f = f * i;
        }
        return f;
    }

    int
    MathUtilities::gcd(int a, int b)
    {
        int c = a % b;
        if (c == 0) {
            return b;
        } else {
            return gcd(b, c);
        }
    }


    KaiserWindow::Parameters
    KaiserWindow::parametersForTransitionWidth(double attenuation,
                                            double transition)
    {
        Parameters p;
        p.length = 1 + (attenuation > 21.0 ?
                        ceil((attenuation - 7.95) / (2.285 * transition)) :
                        ceil(5.79 / transition));
        p.beta = (attenuation > 50.0 ? 
                0.1102 * (attenuation - 8.7) :
                attenuation > 21.0 ? 
                0.5842 * pow(attenuation - 21.0, 0.4) + 0.07886 * (attenuation - 21.0) :
                0);
        return p;
    }

    static double besselTerm(double x, int i)
    {
        if (i == 0) {
            return 1;
        } else {
            double f = MathUtilities::factorial(i);
            return pow(x/2, i*2) / (f*f);
        }
    }

    static double bessel0(double x)
    {
        double b = 0.0;
        for (int i = 0; i < 20; ++i) {
            b += besselTerm(x, i);
        }
        return b;
    }

    void
    KaiserWindow::init()
    {
        double denominator = bessel0(m_beta);
        bool even = (m_length % 2 == 0);
        for (int i = 0; i < (even ? m_length/2 : (m_length+1)/2); ++i) {
            double k = double(2*i) / double(m_length-1) - 1.0;
            m_window.push_back(bessel0(m_beta * sqrt(1.0 - k*k)) / denominator);
        }
        for (int i = 0; i < (even ? m_length/2 : (m_length-1)/2); ++i) {
            m_window.push_back(m_window[int(m_length/2) - i - 1]);
        }
    }

    float
    Pitch::getFrequencyForPitch(int midiPitch,
                                float centsOffset,
                                float concertA)
    {
        float p = float(midiPitch) + (centsOffset / 100);
        return concertA * powf(2.0, (p - 69.0) / 12.0);
    }

    int
    Pitch::getPitchForFrequency(float frequency,
                                float *centsOffsetReturn,
                                float concertA)
    {
        float p = 12.0 * (log(frequency / (concertA / 2.0)) / log(2.0)) + 57.0;

        int midiPitch = int(p + 0.00001);
        float centsOffset = (p - midiPitch) * 100.0;

        if (centsOffset >= 50.0) {
            midiPitch = midiPitch + 1;
            centsOffset = -(100.0 - centsOffset);
        }
        
        if (centsOffsetReturn) *centsOffsetReturn = centsOffset;
        return midiPitch;
    }
    
    void
    SincWindow::init()
    {
        if (m_length < 1) {
            return;
        } else if (m_length < 2) {
            m_window.push_back(1);
            return;
        } else {

            int n0 = (m_length % 2 == 0 ? m_length/2 : (m_length - 1)/2);
            int n1 = (m_length % 2 == 0 ? m_length/2 : (m_length + 1)/2);
            double m = 2 * M_PI / m_p;

            for (int i = 0; i < n0; ++i) {
                double x = ((m_length / 2) - i) * m;
                m_window.push_back(sin(x) / x);
            }

            m_window.push_back(1.0);

            for (int i = 1; i < n1; ++i) {
                double x = i * m;
                m_window.push_back(sin(x) / x);
            }
        }
    }    
}    