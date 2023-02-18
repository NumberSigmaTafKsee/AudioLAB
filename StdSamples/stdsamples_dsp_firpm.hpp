

#include <utility>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <climits>
#include <algorithm>
#include <functional>
#include <cmath>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#if HAVE_MPFR
    #include "core_mpreal.hpp"
#endif


namespace pm {
    namespace pmmath {
		template<typename T> T sin(T);
		template<typename T> T cos(T);
		template<typename T> T tan(T);

		template<typename T> T asin(T);
		template<typename T> T acos(T);
		template<typename T> T atan(T);

		template<typename T> T log(T);
		template<typename T> T exp(T);
		template<typename T> T sqrt(T);
		template<typename T> T pow(T, T);

		template<typename T> T fabs(T);
		template<typename T> T fmax(T, T);
		template<typename T> T fmin(T, T);

		template<typename T> bool signbit(T);
		template<typename T> bool isfinite(T);
		template<typename T> bool isnan(T);

		template<typename T> long round(T);

		template<typename T> T const_pi(void);

		/* Specializations: double precision */
		template<> inline double sin<double>(double x) { return std::sin(x); };
		template<> inline double cos<double>(double x) { return std::cos(x); };
		template<> inline double tan<double>(double x) { return std::tan(x); };

		template<> inline double asin<double>(double x) { return std::asin(x); };
		template<> inline double acos<double>(double x) { return std::acos(x); };
		template<> inline double atan<double>(double x) { return std::atan(x); };

		template<> inline double log<double>(double x) { return std::log(x); };
		template<> inline double exp<double>(double x) { return std::exp(x); };
		template<> inline double sqrt<double>(double x) { return std::sqrt(x); };
		template<> inline double pow<double>(double x, double y) { return std::pow(x, y); };

		template<> inline double fabs<double>(double x) { return std::fabs(x); };
		template<> inline double fmax<double>(double x, double y) { return std::fmax(x, y); };
		template<> inline double fmin<double>(double x, double y) { return std::fmin(x, y); };

		template<> inline bool signbit<double>(double x) { return std::signbit(x); };
		template<> inline bool isfinite<double>(double x) { return std::isfinite(x); };
		template<> inline bool isnan<double>(double x) { return std::isnan(x); };

		template<> inline long round<double>(double x) { return std::round(x); };

		template<> inline double const_pi<double>(void) { return M_PI; };

		/* Specializations: long double precision */
		template<> inline long double sin<long double>(long double x) { return sinl(x); };
		template<> inline long double cos<long double>(long double x) { return cosl(x); };
		template<> inline long double tan<long double>(long double x) { return tanl(x); };

		template<> inline long double asin<long double>(long double x) { return asinl(x); };
		template<> inline long double acos<long double>(long double x) { return acosl(x); };
		template<> inline long double atan<long double>(long double x) { return atanl(x); };

		template<> inline long double log<long double>(long double x) { return logl(x); };
		template<> inline long double exp<long double>(long double x) { return expl(x); };
		template<> inline long double sqrt<long double>(long double x) { return sqrtl(x); };
		template<> inline long double pow<long double>(long double x, long double y) { return powl(x, y); };

		template<> inline long double fabs<long double>(long double x) { return fabsl(x); };
		template<> inline long double fmax<long double>(long double x, long double y) { return fmaxl(x, y); };
		template<> inline long double fmin<long double>(long double x, long double y) { return fminl(x, y); };

		template<> inline bool signbit<long double>(long double x) { return std::signbit(x); };
		template<> inline bool isfinite<long double>(long double x) { return std::isfinite(x); };
		template<> inline bool isnan<long double>(long double x) { return std::isnan(x); };

		template<> inline long round<long double>(long double x) { return std::round(x); };

		template<> inline long double const_pi<long double>(void) { return M_PI; };

		/* Specialization: multiple precision mpreal */
#ifdef HAVE_MPFR
		template<> inline mpfr::mpreal sin<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::sin(x); };
		template<> inline mpfr::mpreal cos<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::cos(x); };
		template<> inline mpfr::mpreal tan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::tan(x); };

		template<> inline mpfr::mpreal asin<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::asin(x); };
		template<> inline mpfr::mpreal acos<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::acos(x); };
		template<> inline mpfr::mpreal atan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::atan(x); };

		template<> inline mpfr::mpreal log<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::log(x); };
		template<> inline mpfr::mpreal exp<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::exp(x); };
		template<> inline mpfr::mpreal sqrt<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::sqrt(x); };
		template<> inline mpfr::mpreal pow<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::pow(x, y); };

		template<> inline mpfr::mpreal fabs<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::abs(x); };
		template<> inline mpfr::mpreal fmax<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::max(x, y); };
		template<> inline mpfr::mpreal fmin<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::min(x, y); };

		template<> inline bool signbit<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::signbit(x); };
		template<> inline bool isfinite<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::isfinite(x); };
		template<> inline bool isnan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::isnan(x); };

		template<> inline long round<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::round(x).toLong(); };

		template<> inline mpfr::mpreal const_pi<mpfr::mpreal>(void) { return mpfr::const_pi(); };
#endif

	} // namespace pmmath
    /**
     * The interval type being considered (in \f$\left[0, \pi\right]\f$
     * or in \f$\left[-1, 1\right]\f$)
     */
    enum space_t {
        FREQ,           /**< not done the change of variable
                        (i.e., we are in \f$\left[0, \pi\right]\f$) */
        CHEBY           /**< done the change of variable
                        (i.e., we are in \f$\left[-1, 1\right]\f$) */
    };

    /**
     * @brief A data type encapsulating information relevant to a
     * frequency band
     *
     * Contains important information concerning a frequency band
     * used during the execution of the exchange algorithm.
     */
    template<typename T>
    struct band_t {
        space_t space;          /**< the space in which we are working */
        std::function<T(space_t, T)> amplitude;
                                /**< the ideal amplitude for this band */
        T start;           /**< the left bound of the band */
        T stop;            /**< the right bound of the band */
        std::function<T(space_t, T)> weight;
                                /**< weight function value on the band */
        std::size_t xs;         /**< number of interpolation points taken in the band */
        std::vector<T> part;    /**< partition points (if any) inside the band */
    };

    /**
     * Gives the direction in which the change of variable is performed
     */
    enum convdir_t {
        FROMFREQ,               /**< apply the change of variable \f$y=\cos(x)\f$*/
        TOFREQ                  /**< apply the change of variable \f$y=\arccos(x)\f$*/
    };

    /*! Performs the change of variable on the set of bands of interest
    * @param[out] out output frequency bands
    * @param[in]  in input frequency bands
    * @param[in]  direction the direction in which the change of
    * variable is performed
    */
    template<typename T>
    void bandconv(std::vector<band_t<T>>& out, std::vector<band_t<T>>& in,
            convdir_t direction);


    /**
    * Gives the type of Chebyshev polynomial expansion to compute
    */
    enum chebkind_t {
            FIRST,          /**< Chebyshev expansion of the first kind*/
            SECOND          /**< Chebyshev expansion of the second kind*/
    };

    /*! Computes the cosines of the elements of a vector
    * @param[in] in the vector to process
    * @param[out] out the vector containing the cosines of the elements from the
    * vector in
    */
    template<typename T>
    void cos(std::vector<T>& out,
            std::vector<T> const& in);

    /*! Does a change of variable from the interval \f$\left[-1, 1\right]\f$ to the
    * interval \f$\left[a, b\right]\f$ on the elements of a given input vector
    * @param[in] in the vector of elements in \f$\left[-1, 1\right]\f$
    * @param[out] out the vector of elements after the change of variable
    * @param[in] a left bound of the desired interval
    * @param[in] b right bound of the desired interval
    */
    template<typename T>
    void chgvar(std::vector<T>& out,
            std::vector<T> const& in,
            T& a, T& b);

    /*! Function that generates equidistant nodes inside the
    * \f$\left[0,\pi\right]\f$ interval, meaning values of the
    * form \f$\frac{i\pi}{n-1}\f$, where \f$0\leq i\leq n-1\f$
    * @param[out] v the vector that will contain the equi-distributed points
    * @param[in] n the number of points which will be computed
    */
    template<typename T>
    void equipts(std::vector<T>& v, std::size_t n);


    /*! This function computes the values of the coefficients of the CI when
    * Chebyshev nodes of the second kind are used
    * @param[out] c vector used to hold the values of the computed coefficients
    * @param[in] fv vector that holds the value of the current function to
    * approximate at the Chebyshev nodes of the second kind scaled to the
    * appropriate interval (in our case it will be \f$\left[0,\pi\right]\f$)
    */
    template<typename T>
    void chebcoeffs(std::vector<T>& c,
                    std::vector<T>& fv);

    /*! Function that generates the coefficients of the derivative of a given CI
    *  @param[out] dc the vector of coefficients of the derivative of the CI
    *  @param[in] c the vector of coefficients of the CI whose derivative we
    *  want to compute
    *  @param[in] kind what kind of coefficient do we want to compute (for a
    *  Chebyshev expansion of the first or second kind)
    */
    template<typename T>
    void diffcoeffs(std::vector<T>& dc,
                    std::vector<T>& c,
                    chebkind_t kind = SECOND);

    /*! Chebyshev proxy rootfinding method for a given CI
    * @param[out] r the vector of computed roots of the CI
    * @param[in] c the Chebyshev coeffients of the polynomial whose roots
    * we want to find
    * @param[in] dom the real domain where we are looking for the roots
    * @param[in] kind the type of Chebyshev expansion (first or second)
    * @param[in] balance flag signaling if we should use balancing (in
    * the vein of [Parlett&Reinsch1969] "Balancing a Matrix for
    * Calculation of Eigenvalues and Eigenvectors") for the resulting
    * Chebyshev companion matrix
    */
    template<typename T>
    void roots(std::vector<T>& r, std::vector<T>& c,
            std::pair<T, T> const& dom,
            chebkind_t kind = SECOND,
            bool balance = true);

    /*! Procedure which computes the weights used in
    * the evaluation of the barycentric interpolation
    * formulas (see [Berrut&Trefethen2004] and [Pachon&Trefethen2009]
    * for the implementation ideas)
    * @param[out] w the computed weights
    * @param[in] x the interpolation points
    */
    template<typename T>
    void baryweights(std::vector<T>& w,
            std::vector<T>& x);

    /*! Determines the current reference error according to the
    * barycentric formula (internally it also computes the barycentric weights)
    * @param[out] delta the value of the current reference error
    * @param[in] x the current reference set (i.e., interpolation points)
    * @param[in] bands information relating to the ideal frequency response of
    * the filter
    */
    template<typename T>
    void compdelta(T& delta, std::vector<T>& x,
            std::vector<band_t<T>>& bands);

    /*! Determines the current reference error according to the
    * barycentric formula
    * @param[out] delta the value of the current reference error
    * @param[in] w the barycentric weights associated with the current reference
    * set
    * @param[in] x the current reference set (i.e. interpolation points)
    * @param[in] bands information relating to the ideal frequency response of
    * the filter
    */

    template<typename T>
    void compdelta(T& delta, std::vector<T>& w,
            std::vector<T>& x, std::vector<band_t<T>>& bands);

    /*! Computes the filter response at the current reference set
    * @param[out] C the vector of frequency responses at the reference set
    * @param[in] delta the current reference error
    * @param[in] x the current reference vector
    * @param[in] bands frequency band information for the ideal filter
    */
    template<typename T>
    void compc(std::vector<T>& C, T& delta,
            std::vector<T>& x, std::vector<band_t<T>>& bands);

    /*! Computes the frequency response of the current filter
    * @param[out] Pc the frequency response amplitude value at the current node
    * @param[in] xVal the current frequency node where we do our computation
    * (the point is given in the \f$\left[-1,1\right]\f$ interval,
    * and not the initial \f$\left[0,\pi\right]\f$)
    * @param[in] x the current reference set
    * @param[in] C the frequency responses at the current reference set
    * @param[in] w the current barycentric weights
    */
    template<typename T>
    void approx(T& Pc, T const& xVal,
            std::vector<T>& x, std::vector<T>& C,
            std::vector<T>& w);

    /*! Computes the approximation error at a given node using the current set of
    * reference points
    * @param[out] error the requested error value
    * @param[in] xVal the current frequency node where we do our computation
    * @param[in] delta the current reference error
    * @param[in] x the current reference set
    * @param[in] C the frequency response values at the x nodes
    * @param[in] w the barycentric weights
    * @param[in] bands frequency band information for the ideal filter
    */
    template<typename T>
    void comperror(T& error, T const& xVal,
            T& delta, std::vector<T>& x,
            std::vector<T>& C, std::vector<T>& w,
            std::vector<band_t<T>>& bands);

    /*! The ideal frequency response and weight information at the given frequency
    * node (it can be in the \f$\left[-1,1\right]\f$ interval,
    * and not the initial \f$\left[0,\pi\right]\f$, the difference is made with
    * information from the bands parameter)
    * @param[out] D ideal frequency response
    * @param[out] W weight value for the current point
    * @param[in] xVal the current frequency node where we do our computation
    * @param[in] bands frequency band information for the ideal filter
    */
    template<typename T>
    void idealvals(T& D, T& W,
            T const& xVal, std::vector<band_t<T>>& bands);

/** @enum filter_t marker to distinguish
     * between the two categories of filters (digital
     * differentiators and Hilbert transformers) that
     * can be constructed using type III and IV FIR
     * filters. */
    enum class filter_t {
        FIR_DIFFERENTIATOR,     /**< marker for constructing digital differentiators */
        FIR_HILBERT             /**< marker for constructing Hilbert transformers */
    };

    /** @enum init_t flag representing the
     * initialization strategies that can be used
     * at the lowest level of the scaling approach 
     */
    enum class init_t {
        UNIFORM,                /**< uniform initialization marker */
        SCALING,                /**< reference scaling-based initialization */
        AFP                     /**< AFP algorithm-based initialization */
    };

    /** @enum status_t code to distinguish the
     * various states in which the Parks-McClellan
     * algorithm execution finished in. */
    enum class status_t {
        STATUS_SUCCESS,                     /**< successful execution */
        STATUS_FREQUENCY_INVALID_INTERVAL,  /**< invalid frequency inputs */
        STATUS_AMPLITUDE_VECTOR_MISMATCH,   /**< amplitude/frequency vector sizes mismatch */
        STATUS_AMPLITUDE_DISCONTINUITY,     /**< discontinuous amplitude values detected */
        STATUS_WEIGHT_NEGATIVE,             /**< negative weight value detected */
        STATUS_WEIGHT_VECTOR_MISMATCH,      /**< weight/frequency vector sizes mismatch */
        STATUS_WEIGHT_DISCONTINUITY,        /**< discontinuous weight value detected */
        STATUS_SCALING_INVALID,             /**< failure in performing valid reference scaling */
        STATUS_AFP_INVALID,                 /**< numerical failure in performing AFP initialization */ 
        STATUS_COEFFICIENT_SET_INVALID,     /**< invalid final coefficient set */
        STATUS_EXCHANGE_FAILURE,            /**< runtime error in producing a valid reference set */
        STATUS_CONVERGENCE_WARNING,         /**< successful execution, but with convengence warnings */
        STATUS_UNKNOWN_FAILURE              /**< unknown runtime failure */
    };

    /**
     * @brief The type of the object returned by the Parks-McClellan algorithm.
     *
     * Utility object which contains useful information about the filter computed by the
     * Parks-McClellan algorithm
     */
    template<typename T>
    struct pmoutput_t
    {
        std::vector<T> h;           /**< the final filter coefficients*/
        std::vector<T> x;           /**< the reference set used to generate the final
                                    filter (values are in \f$[-1,1]\f$ and NOT 
                                    \f$[0,\pi]\f$)*/
        std::size_t iter;           /**< number of iterations that were necessary to
                                    achieve convergence*/
        T delta;                    /**< the final reference error */
        T q;                        /**< convergence parameter value */
        status_t status;            /**< status code for the output object */
    };

    /*! An implementation of the uniform initialization approach for
    * starting the Parks-McClellan algorithm
    * @param[out] omega the initial set of references to be computed
    * @param[in] B the frequency bands of interest (i.e., stopbands and
    * passbands for example)
    * @param[in] n the size of the reference set
    */

    template<typename T>
    void uniform(std::vector<T>& omega,
            std::vector<band_t<T>>& B, std::size_t n);

    /*! An implementation of the reference scaling approach mentioned
    * in section 4 of the article.
    * @param[out] status saves diagnostic information in case the routine does not
    * complete successfully
    * @param[out] nx the reference set obtained from scaling the initial set x
    * @param[out] ncbands contains information about the bands of interest 
    * corresponding to the new reference (i.e., how many reference points are 
    * inside each band). The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY 
    * band space)
    * @param[out] nfbands contains information about the bands of interest corresponding 
    * to the new reference (i.e., how many reference points are inside each band). 
    * The bands are given inside \f$[0,\pi]\f$ (i.e., the FREQ band space)
    * @param[in] nxs the size of the new reference set
    * @param[in] x the input reference on which will be used to perform the scaling
    * @param[in] cbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY 
    * band space)
    * @param[in] fbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[0,\pi]\f$ (i.e., the FREQ 
    * band space)
    */
    template<typename T>
    void refscaling(status_t& status,
            std::vector<T>& nx, std::vector<band_t<T>>& ncbands,
            std::vector<band_t<T>>& nfbands, std::size_t nxs,
            std::vector<T>& x, std::vector<band_t<T>>& cbands,
            std::vector<band_t<T>>& fbands);

    /*! An internal routine which implements the exchange algorithm for designing 
    * FIR filters
    * @param[in] x the initial reference set
    * @param[in] cbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY band space)
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. Since there is no explicit filter type (I to IV) given as input, 
    * the h vector of the output will correspond to the coefficients of the frequency 
    * response \f$H_d(\omega)=\sum_{k=0}^{n}h_k\cos(\omega k)\f$
    *
    * The following piece of code shows how to use this function to construct a lowpass 
    * filter. It also shows how to write frequency band specifications and how to use 
    * the uniform initialization routine. It is equivalent to the sample code from 
    * firpm. <b>Important:</b> in order to determine the final coefficients of the 
    * transfer function of the filter, some preprocessing specific to the type of the 
    * filter (I to IV) has to be done on the elements of <tt>output.h</tt> (see for 
    * example the source code of the <tt>firpm</tt> functions on how this is done for 
    * each type of filter).
    * @see firpm
    * @code
    * // frequency band specification
    * std::vector<band_t<double>> fbands(2);
    * double pi = M_PI;
    *
    * fbands[0].start = 0;
    * fbands[0].stop = pi * 0.4;
    * fbands[0].weight = [] (space_t, double) -> double {return 1.0; };
    * fbands[0].space = space_t::FREQ;
    * fbands[0].amplitude = [](space_t, double) -> double { return 1.0; };
    *
    * fbands[1].start = pi * 0.5;
    * fbands[1].stop = pi;
    * fbands[1].weight = [] (space_t, double) -> double {return 10.0; };
    * fbands[1].space = BandSpace::FREQ;
    * fbands[1].amplitude = [](space_t, double) -> double { return 10.0; };
    * std::size_t degree = 100;  // filter degree
    *
    * // reference initialization code
    * std::vector<band_t<double>> cbands;
    * std::vector<double> omega;
    * std::vector<double> x;
    * uniform(omega, fbands, degree + 2u);
    * // apply the change of variable y = cos(x) so that we are working inside [-1, 1]
    * cos(x, omega);
    * bandconv(cbands, fbands, convdir_t::FROMFREQ);
    * // apply the exchange algorithm
    * pmoutput_t<double> output = exchange(x, cbands);
    * @endcode
    */

    template<typename T>
    pmoutput_t<T> exchange(std::vector<T>& x,
            std::vector<band_t<T>>& cbands,
            double eps = 0.01,
            std::size_t nmax = 4u, 
            unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters.
    * This routine is the most general and can be set to use any of the three
    * proposed initialization types. If the default parameters are used, then
    * it will use uniform initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type I, while for odd 
    * n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number 
    * of significant digits of the minimax error that are accurate at the end of 
    * the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] strategy initialization strategy. Can be UNIFORM, SCALING or AFP
    * @param[in] depth in case the SCALING initialization strategy is used, 
    * specifies the number of scaling levels to use (by default, the value is set
    * to 1, meaning a filter of length approximatively n/2 is used to construct
    * the initial reference for the requested n coefficient filter)
    * @param[in] rstrategy in case SCALING is used, specifies how to initialize
    * the smallest length filter used to perform reference scaling (UNIFORM by
    * default)
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types I and II,
    * the values are symmetrical to the middle coefficient(s))
    *
    * An example of how to use this function is given below. It designs a degree 
    * \f$100\f$ type I lowpass filter, with passband \f$[0, 0.4\pi]\f$ and stopband 
    * \f$[0.5\pi, \pi]\f$. It has unit weight inside the passband and weight 10
    * inside the stopband. More examples, including code on how to use the customized
    * reference scaling and AFP versions <tt>firpmRS, firpmAFP</tt>, are provided inside 
    * the test files.
    * @code
    * pmoutput_t<double> output = firpm<double>(200, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    * @endcode
    */

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                init_t strategy = init_t::UNIFORM,
                std::size_t depth = 0u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters. This routine uses 
    * reference scaling by default and is just a wrapper over the <tt>firpm</tt>.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final transfer function. 
    * For even n, the filter will be type I, while for odd n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the lowest level 
    * (uniform or AFP-based) if strategy is reference scaling 
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the transfer 
    * function of the final filter (in this case, for types I and II, the values are 
    * symmetrical to the middle coefficient(s))
    */

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                std::size_t depth = 1u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters. This routine 
    * uses AFP-based initialization and is just a wrapper over <tt>firpm</tt>.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final transfer function. 
    * For even n, the filter will be type I, while for odd n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the transfer 
    * function of the final filter (in this case, for types I and II, the values are 
    * symmetrical to the middle coefficient(s))
    */

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters. 
    * This routine is the most general and can be set to use any of the three
    * proposed initialization types. If the default parameters are used, then
    * it will use uniform initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for odd 
    * n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number 
    * of significant digits of the minimax error that are accurate at the end of 
    * the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] strategy the initialization strategy
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the lowest level 
    * (uniform or AFP-based) if strategy is reference scaling 
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the 
    * transfer function of the final filter (in this case, for types III and IV, 
    * the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
            std::vector<T>const& f,
            std::vector<T>const& a,
            std::vector<T>const& w,
            filter_t type,
            double eps = 0.01,
            std::size_t nmax = 4,
            init_t strategy = init_t::UNIFORM,
            std::size_t depth = 0u,
            init_t rstrategy = init_t::UNIFORM,
            unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters. 
    * This routine uses reference scaling.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for 
    * odd n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e., quantizes the 
    * number of significant digits of the minimax error that are accurate 
    * at the end of the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the 
    * lowest level (uniform or AFP-based)
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types III and 
    * IV, the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps = 0.01,
                std::size_t nmax = 4u,
                std::size_t depth = 1u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters.
    * This routine uses AFP-based initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for 
    * odd n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e quantizes the number 
    * of significant digits of the minimax error that are accurate at the end 
    * of the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types III and 
    * IV, the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps = 0.01,
                std::size_t nmax = 4u,
                unsigned long prec = 165ul);


    template<typename T>
    void bandconv(std::vector<band_t<T>>& out, std::vector<band_t<T>>& in,
            convdir_t direction)
    {
        out.resize(in.size());
        std::size_t n = in.size() - 1u;
        for (std::size_t i{0u}; i < in.size(); ++i)
        {
            out[i].weight    = in[n - i].weight;
            out[i].amplitude = in[n - i].amplitude;
            out[i].xs        = in[n - i].xs;
            out[i].part      = in[n - i].part;
            if (direction == convdir_t::FROMFREQ)
            {
                out[i].start = pmmath::cos(in[n - i].stop);
                out[i].stop  = pmmath::cos(in[n - i].start);
                out[i].space = space_t::CHEBY;
                for(std::size_t j{0}; j < out[i].part.size(); ++j)
                    out[i].part[j] = pmmath::cos(out[i].part[j]);
                std::sort(begin(out[i].part), end(out[i].part));
            } else {
                out[i].start = pmmath::acos(in[n - i].stop);
                out[i].stop  = pmmath::acos(in[n - i].start);
                out[i].space = space_t::FREQ;
                for(std::size_t j{0}; j < out[i].part.size(); ++j)
                    out[i].part[j] = pmmath::acos(out[i].part[j]);
                std::sort(begin(out[i].part), end(out[i].part));
            }
        }
    }

    /* Template instantiation */
    template void bandconv<double>(
        std::vector<band_t<double>>& out,
        std::vector<band_t<double>>& in,
            convdir_t direction);

    template void bandconv<long double>(
        std::vector<band_t<long double>>& out,
        std::vector<band_t<long double>>& in,
            convdir_t direction);

#ifdef HAVE_MPFR
    template void bandconv<mpfr::mpreal>(
        std::vector<band_t<mpfr::mpreal>>& out,
        std::vector<band_t<mpfr::mpreal>>& in,
            convdir_t direction);
#endif

template<typename T>
    void baryweights(std::vector<T>& w,
            std::vector<T>& x)
    {
        if(x.size() > 500u)
        {
            for(std::size_t i{0u}; i < x.size(); ++i)
            {
                T one = 1;
                T denom = 0.0;
                T xi = x[i];
                for(std::size_t j{0u}; j < x.size(); ++j)
                {
                    if (j != i) {
                        denom += pmmath::log(((xi - x[j] > 0) ? (xi - x[j]) : (x[j] - xi)));
                        one *= ((xi - x[j] > 0) ? 1 : -1);
                    }
                }
                w[i] = one / pmmath::exp(denom + pmmath::log(2.0)* (x.size() - 1));
            }
        }
        else
        {
            std::size_t step = (x.size() - 2) / 15 + 1;
            T one = 1u;
            for(std::size_t i{0u}; i < x.size(); ++i)
            {
                T denom = 1.0;
                T xi = x[i];
                for(std::size_t j{0u}; j < step; ++j)
                {
                    for(std::size_t k{j}; k < x.size(); k += step)
                        if (k != i)
                            denom *= ((xi - x[k]) * 2);
                }
                w[i] = one / denom;
            }
        }
    }


    template<typename T>
    void idealvals(T& D, T& W,
            T const& x, std::vector<band_t<T>>& bands)
    {
        for (auto &it : bands) {
            if (x >= it.start && x <= it.stop) {
                D = it.amplitude(it.space, x);
                W = it.weight(it.space, x);
                return;
            }
        }
    }

    template<typename T>
    void compdelta(T& delta, std::vector<T>& x,
            std::vector<band_t<T>>& bands)
    {
        std::vector<T> w(x.size());
        baryweights(w, x);

        T num, denom, D, W, buffer;
        num = denom = D = W = 0;
        for (std::size_t i{0u}; i < w.size(); ++i)
        {
            idealvals(D, W, x[i], bands);
            buffer = w[i];
            num += buffer * D;
            buffer = w[i] / W;
            if (i % 2 == 0)
                buffer = -buffer;
            denom += buffer;
        }

        delta = num / denom;
    }

    template<typename T>
    void compdelta(T& delta, std::vector<T>& w,
            std::vector<T>& x, std::vector<band_t<T>>& bands)
    {
        T num, denom, D, W, buffer;
        num = denom = D = W = 0;
        for (std::size_t i{0u}; i < w.size(); ++i)
        {
            idealvals(D, W, x[i], bands);
            buffer = w[i];
            num += buffer * D;
            buffer = w[i] / W;
            if (i % 2 == 0)
                buffer = -buffer;
            denom += buffer;
        }
        delta = num / denom;
    }


    template<typename T>
    void compc(std::vector<T>& C, T& delta,
            std::vector<T>& omega, std::vector<band_t<T>>& bands)
    {
        T D, W;
        D = W = 0;
        for (std::size_t i{0u}; i < omega.size(); ++i)
        {
            idealvals(D, W, omega[i], bands);
            if (i % 2 != 0)
                W = -W;
            C[i] = D + (delta / W);
        }
    }

    template<typename T>
    void approx(T& Pc, T const& omega,
            std::vector<T>& x, std::vector<T>& C,
            std::vector<T>& w)
    {
        T num, denom;
        T buff;
        num = denom = 0;

        Pc = omega;
        std::size_t r = x.size();
        for (std::size_t i{0u}; i < r; ++i)
        {
            if (Pc == x[i]) {
                Pc = C[i];
                return;
            }
            buff = w[i] / (Pc - x[i]);
            num += buff * C[i];
            denom += buff;
        }
        Pc = num / denom;
    }

    template<typename T>
    void comperror(T& error, T const& xVal,
            T& delta, std::vector<T>& x,
            std::vector<T>& C, std::vector<T>& w,
            std::vector<band_t<T>>& bands)
    {
        for (std::size_t i{0u}; i < x.size(); ++i)
        {
            if (xVal == x[i]) {
                if (i % 2 == 0)
                    error = delta;
                else
                    error = -delta;
                return;
            }
        }

        T D, W;
        D = W = 0;
        idealvals(D, W, xVal, bands);
        approx(error, xVal, x, C, w);
        error -= D;
        error *= W;
    }

    /* Template instantiations */

    /* double precision */

    template void baryweights<double>(std::vector<double>& w,
            std::vector<double>& x);

    template void compdelta<double>(double& delta,
            std::vector<double>& x, std::vector<band_t<double>>& bands);

    template void compdelta<double>(double& delta,
            std::vector<double>& w, std::vector<double>& x,
            std::vector<band_t<double>>& bands);

    template void compc<double>(std::vector<double>& C, double& delta,
            std::vector<double>& x, std::vector<band_t<double>>& bands);

    template void approx<double>(double& Pc, double const& xVal,
            std::vector<double>& x, std::vector<double>& C,
            std::vector<double>& w);

    template void comperror<double>(double& error, double const& xVal,
            double& delta, std::vector<double>& x,
            std::vector<double>& C, std::vector<double>& w,
            std::vector<band_t<double>>& bands);

    /* long double precision */

    template void baryweights<long double>(std::vector<long double>& w,
            std::vector<long double>& x);

    template void compdelta<long double>(long double& delta,
            std::vector<long double>& x, std::vector<band_t<long double>>& bands);

    template void compdelta<long double>(long double& delta,
            std::vector<long double>& w, std::vector<long double>& x,
            std::vector<band_t<long double>>& bands);

    template void compc<long double>(std::vector<long double>& C, long double& delta,
            std::vector<long double> &x, std::vector<band_t<long double>>& bands);

    template void approx<long double>(long double& Pc, long double const& xVal,
            std::vector<long double>& x, std::vector<long double>& C,
            std::vector<long double>& w);

    template void comperror<long double>(long double& error, long double const& xVal,
            long double& delta, std::vector<long double>& x,
            std::vector<long double>& C, std::vector<long double>& w,
            std::vector<band_t<long double>>& bands);

#ifdef HAVE_MPFR
    // separate implementation for the MPFR version; it is much faster and
    // the higher precision should usually compensate for any eventual
    // ill-conditioning
    template<> void baryweights<mpfr::mpreal>(std::vector<mpfr::mpreal>& w,
            std::vector<mpfr::mpreal>& x)
    {
        std::size_t step = (x.size() - 2u) / 15 + 1;
        mpfr::mpreal one = 1u;
        for(std::size_t i{0u}; i < x.size(); ++i)
        {
            mpfr::mpreal denom = 1.0;
            mpfr::mpreal xi = x[i];
            for(std::size_t j{0u}; j < step; ++j)
            {
                for(std::size_t k{j}; k < x.size(); k += step)
                    if (k != i)
                        denom *= ((xi - x[k]) << 1);
            }
            w[i] = one / denom;
        }
    }
    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& w, std::vector<mpfr::mpreal>& x,
            std::vector<band_t<mpfr::mpreal>>& bands);

    template void compc<mpfr::mpreal>(std::vector<mpfr::mpreal>& C, mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void approx<mpfr::mpreal>(mpfr::mpreal& Pc, mpfr::mpreal const& xVal,
            std::vector<mpfr::mpreal>& x, std::vector<mpfr::mpreal>& C,
            std::vector<mpfr::mpreal>& w);

    template void comperror<mpfr::mpreal>(mpfr::mpreal& error, mpfr::mpreal const& xVal,
            mpfr::mpreal& delta, std::vector<mpfr::mpreal>& x,
            std::vector<mpfr::mpreal>& C, std::vector<mpfr::mpreal>& w,
            std::vector<band_t<mpfr::mpreal>>& bands);
#endif

template<typename T>
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename T>
    using VectorXcd = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;

    template<typename T>
    void balance(MatrixXd<T>& A)
    {
        std::size_t n = A.rows();

        T rNorm;      // row norm
        T cNorm;      // column norm
        bool converged = false;

        T g, f, s;
        while(!converged)
        {
            converged = true;
            for(std::size_t i{0u}; i < n; ++i)
            {
                rNorm = cNorm = 0.0;
                for(std::size_t j{0u}; j < n; ++j)
                {
                    if(i == j)
                        continue;
                    cNorm += pmmath::fabs(A(j, i));
                    rNorm += pmmath::fabs(A(i, j));
                }
                if((cNorm == 0.0) || (rNorm == 0))
                    continue;

                g = rNorm / 2.0;
                f = 1.0;
                s = cNorm + rNorm;


                while(pmmath::isfinite(cNorm) && cNorm < g)
                {
                    f *= 2.0;
                    cNorm *= 4.0;
                }

                g = rNorm * 2.0;

                while(pmmath::isfinite(cNorm) && cNorm > g)
                {
                    f /= 2.0;
                    cNorm /= 4.0;
                }

                if((rNorm + cNorm) < s * f * 0.95)
                {
                    converged = false;
                    g = T(1.0) / f;
                    // multiply by D^{-1} on the left
                    A.row(i) *= g;
                    // multiply by D on the right
                    A.col(i) *= f;
                }

            }
        }
    }

    template<typename T>
    MatrixXd<T> colleague(std::vector<T> const& c,
                    chebkind_t kind,
                    bool bal)
    {
        std::vector<T> a{c};
        std::size_t n = c.size() - 1u;
        MatrixXd<T> C(n, n);

        for(std::size_t i{0u}; i < n; ++i)
            for(std::size_t j{0u}; j < n; ++j)
                C(i, j) = 0;

        T denom = -1;
        denom /= a[n];
        denom /= 2;
        for(std::size_t i{0u}; i < c.size() - 1u; ++i)
            a[i] *= denom;
        a[n - 2u] += 0.5;

        for (std::size_t i{0u}; i < n - 1; ++i)
            C(i, i + 1) = C(i + 1, i) = 0.5;
        switch(kind) {
            case FIRST: C(n - 2, n - 1) = 1.0;  break;
            default:    C(n - 2, n - 1) = 0.5;  break;
        }
        for(std::size_t i{0u}; i < n; ++i)
            C(i, 0) = a[n - i - 1u];

        if(bal)
            balance(C);    

        return C;
    }

    template<typename T>
    void cos(std::vector<T>& out,
            std::vector<T> const& in)
    {
        out.resize(in.size());
        for(std::size_t i{0u}; i < in.size(); ++i)
            out[i] = pmmath::cos(in[i]);
    }

    template<typename T>
    void chgvar(std::vector<T>& out,
            std::vector<T> const& in,
            T& a, T& b)
    {
        out.resize(in.size());
        for(std::size_t i{0u}; i < in.size(); ++i)
            out[i] = (b + a) / 2 + in[i] * (b - a) / 2;
    }

    template<typename T>
    void clenshaw(T& result, const std::vector<T>& p,
            T const& x, T const& a, T const& b)
    {
        T bn1, bn2, bn;
        T buffer;

        bn1 = 0;
        bn2 = 0;

        // compute the value of (2*x - b - a)/(b - a) 
        // in the temporary variable buffer
        buffer = (x * 2 - b - a) / (b - a);

        int n = (int)p.size() - 1;
        for(int k{n}; k >= 0; --k) {
            bn = buffer * 2;
            bn = bn * bn1 - bn2 + p[k];
            // update values
            bn2 = bn1;
            bn1 = bn;
        }

        // set the value for the result
        // (i.e., the CI value at x)
        result = bn1 - buffer * bn2;
    }

    template<typename T>
    void clenshaw(T& result,
                std::vector<T> const& p,
                T const& x,
                chebkind_t kind)
    {
        T bn1, bn2, bn;

        int n = (int)p.size() - 1;
        bn2 = 0;
        bn1 = p[n];
        for(int k{n - 1}; k >= 1; --k) {
            bn = x * 2;
            bn = bn * bn1 - bn2 + p[k];
            // update values
            bn2 = bn1;
            bn1 = bn;
        }

        if(kind == FIRST)
            result = x * bn1 - bn2 + p[0];
        else
            result = (x * 2) * bn1 - bn2 + p[0];
    }

    template<typename T>
    void equipts(std::vector<T>& v, std::size_t n)
    {
        v.resize(n);
        // store the points in the vector v as 
        // v[i] = i * pi / (n-1)
        for(std::size_t i{0u}; i < n; ++i) {
            v[i] = pmmath::const_pi<T>() * i;
            v[i] /= (n-1);
        }
    }

    // this function computes the values of the coefficients of 
    // the CI when Chebyshev nodes of the second kind are used
    template<typename T>
    void chebcoeffs(std::vector<T>& c,
            std::vector<T>& fv)
    {
        std::size_t n = fv.size();
        std::vector<T> v(n);
        equipts(v, n);

        T buffer;

        // halve the first and last coefficients
        T oldValue1 = fv[0];
        T oldValue2 = fv[n-1u];
        fv[0] /= 2;
        fv[n-1u] /= 2;

        for(std::size_t i{0u}; i < n; ++i) {
            // compute the actual value at the Chebyshev
            // node cos(i * pi / n)
            buffer = pmmath::cos(v[i]);
            clenshaw(c[i], fv, buffer, FIRST);

            if(i == 0u || i == n-1u) {
                c[i] /= (n-1u);
            } else {
                c[i] *= 2;
                c[i] /= (n-1u);
            }
        }
        fv[0] = oldValue1;
        fv[n-1u] = oldValue2;

    }

    // function that generates the coefficients of the 
    // derivative of a given CI
    template<typename T>
    void diffcoeffs(std::vector<T>& dc,
                    std::vector<T>& c,
                    chebkind_t kind)
    {
        dc.resize(c.size()-1);
        switch(kind) {
            case FIRST: {
                int n = c.size() - 2;
                dc[n] = c[n + 1] * (2 * (n + 1));
                dc[n - 1] = c[n] * (2 * n);
                for(int i{n - 2}; i >= 0; --i) {
                    dc[i] = 2 * (i + 1);
                    dc[i] = dc[i] * c[i + 1] + dc[i + 2];
                }
                dc[0] /= 2;
            };
            break;
            default: {
                int n = c.size() - 1;
                for(int i{n}; i > 0; --i)
                    dc[i - 1] = c[i] * i;
            }
            break;
        }
    }

    template<typename T>
    void roots(std::vector<T>& r, std::vector<T>& c,
            std::pair<T, T> const& dom,
            chebkind_t kind,
            bool balance)
    {
        r.clear();
        for(auto &it : c)
            if(!pmmath::isfinite(it))
                return;

        // preprocess c and remove all the high order
        // zero coefficients
        std::size_t idx{c.size() - 1u};
        while(idx > 0u && c[idx] == 0) { --idx; }

        if(idx > 1u) {
            c.resize(idx + 1u);
            MatrixXd<T> C = colleague(c, kind, balance);
            Eigen::EigenSolver<MatrixXd<T>> es(C);
            VectorXcd<T> eigs = es.eigenvalues();

            T threshold = 1e-20;
            for(Eigen::Index i{0}; i < eigs.size(); ++i) {
                if(pmmath::fabs(eigs(i).imag()) < threshold)
                    if(dom.first < eigs(i).real() && 
                    dom.second > eigs(i).real()) {
                        r.push_back(eigs(i).real());
                    }
            }

            std::sort(begin(r), end(r));
        } else if(idx == 1) {
            r.push_back(-c[0u] / c[1u]);
        }
    }

    /* Explicit instantiation */

    /* double precision */

    template void cos<double>(std::vector<double>& out,
            std::vector<double> const& in);

    template void chgvar<double>(std::vector<double>& out,
            std::vector<double> const& in,
            double& a, double& b);

    template void equipts<double>(std::vector<double>& v, std::size_t n);


    template void chebcoeffs<double>(std::vector<double>& c,
                    std::vector<double>& fv);

    template void diffcoeffs<double>(std::vector<double>& dc,
                    std::vector<double>& c,
                    chebkind_t kind);

    template void roots<double>(std::vector<double>& r, std::vector<double>& c,
            std::pair<double, double> const& dom,
            chebkind_t kind, bool balance);

    /* long double precision */

    template void cos<long double>(std::vector<long double>& out,
            std::vector<long double> const& in);

    template void chgvar<long double>(std::vector<long double>& out,
            std::vector<long double> const& in,
            long double& a, long double& b);

    template void equipts<long double>(std::vector<long double>& v, std::size_t n);


    template void chebcoeffs<long double>(std::vector<long double>& c,
                    std::vector<long double>& fv);

    template void diffcoeffs<long double>(std::vector<long double>& dc,
                    std::vector<long double>& c,
                    chebkind_t kind);

    template void roots<long double>(std::vector<long double>& r, std::vector<long double>& c,
            std::pair<long double, long double> const& dom,
            chebkind_t kind, bool balance);

#ifdef HAVE_MPFR
    template void cos<mpfr::mpreal>(std::vector<mpfr::mpreal>& out,
            std::vector<mpfr::mpreal> const& in);

    template void chgvar<mpfr::mpreal>(std::vector<mpfr::mpreal>& out,
            std::vector<mpfr::mpreal> const& in,
            mpfr::mpreal& a, mpfr::mpreal& b);

    template void equipts<mpfr::mpreal>(std::vector<mpfr::mpreal>& v, std::size_t n);


    template void chebcoeffs<mpfr::mpreal>(std::vector<mpfr::mpreal>& c,
                    std::vector<mpfr::mpreal>& fv);

    template void diffcoeffs<mpfr::mpreal>(std::vector<mpfr::mpreal>& dc,
                    std::vector<mpfr::mpreal>& c,
                    chebkind_t kind);

    template void roots<mpfr::mpreal>(std::vector<mpfr::mpreal>& r, std::vector<mpfr::mpreal>& c,
            std::pair<mpfr::mpreal, mpfr::mpreal> const& dom,
            chebkind_t kind, bool balance);
#endif

template<typename T>
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename T>
    using VectorXd = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    // global variable used to detect cycling
    bool cycle;

    template<typename T>
    void chebvand(MatrixXd<T>& A, std::size_t degree,
            std::vector<T>& meshPoints,
            std::function<T(T)>& weightFunction)
    {

        A.resize(degree + 1u, meshPoints.size());
        for(std::size_t i{0u}; i < meshPoints.size(); ++i)
        {
            T pointWeight = weightFunction(meshPoints[i]);
            A(0u, i) = 1;
            A(1u, i) = meshPoints[i];
            for(std::size_t j{2u}; j <= degree; ++j)
                A(j, i) = meshPoints[i] * A(j - 1u, i) * 2 - A(j - 2u, i);
            for(std::size_t j{0u}; j <= degree; ++j)
                A(j, i) *= pointWeight;
        }
    }

    // approximate Fekete points
    template<typename T>
    void afp(std::vector<T>& points, MatrixXd<T>& A,
            std::vector<T>& mesh)
    {
        VectorXd<T> b = VectorXd<T>::Ones(A.rows());
        b(0) = 2;
        VectorXd<T> y = A.colPivHouseholderQr().solve(b);
        points.clear();

        for(Eigen::Index i{0}; i < y.rows(); ++i)
            if(y(i) != 0.0)
                points.push_back(mesh[i]);
        std::sort(points.begin(), points.end(),
                [](const T& lhs,
                const T& rhs) {
                    return lhs < rhs;
                });
    }

    template<typename T>
    void countBand(std::vector<band_t<T>>& cb, std::vector<T>& x)
    {
        for(auto& it : cb)
            it.xs = 0u;
        std::size_t bandIt = 0u;
        for(std::size_t i{0u}; i < x.size(); ++i)
        {
            while(bandIt < cb.size() && cb[bandIt].stop < x[i])
                bandIt++;
            ++cb[bandIt].xs;
        }
    }

    template<typename T>
    void wam(std::vector<T>& wam, std::vector<band_t<T>>& cb,
            std::size_t deg)
    {
        std::vector<T> cp;
        equipts(cp, deg + 2u);
        cos(cp, cp);
        std::sort(begin(cp), end(cp));
        for(std::size_t i{0u}; i < cb.size(); ++i)
        {
            if(cb[i].start != cb[i].stop)
            {
                std::vector<T> bufferNodes;
                chgvar(bufferNodes, cp, cb[i].start, cb[i].stop);
                bufferNodes[0] = cb[i].start;
                bufferNodes[bufferNodes.size()-1u] = cb[i].stop;
                for(auto& it : bufferNodes)
                    wam.push_back(it);
            }
            else
                wam.push_back(cb[i].start);
        }
    }

    template<typename T>
    void uniform(std::vector<T>& omega,
            std::vector<band_t<T>>& B, std::size_t n)
    {
        T avgDist = 0;
        omega.resize(n);

        std::vector<T> bandwidths(B.size());
        std::vector<std::size_t> nonPointBands;
        for(std::size_t i{0u}; i < B.size(); ++i) {
            bandwidths[i] = B[i].stop - B[i].start;
            if(bandwidths[i] > 0.0)
            {
                nonPointBands.push_back(i);
                avgDist += bandwidths[i];
            }
            B[i].xs = 1u;
        }

        avgDist /= (omega.size() - B.size());
        std::size_t npSize = nonPointBands.size();

        B[nonPointBands[npSize - 1u]].xs = omega.size() - (B.size() - npSize);
        T buffer;
        buffer = bandwidths[nonPointBands[0]] / avgDist;

        if (npSize > 1) {
            B[nonPointBands[0]].xs = pmmath::round(buffer) + 1;
            B[nonPointBands[npSize - 1u]].xs -= B[nonPointBands[0]].xs;
        }

        for(std::size_t i{1u}; i < npSize - 1u; ++i) {
            buffer = bandwidths[nonPointBands[i]] / avgDist;
            B[nonPointBands[i]].xs = pmmath::round(buffer) + 1;
            B[nonPointBands[npSize - 1u]].xs -= B[nonPointBands[i]].xs;
        }


        std::size_t startIndex = 0ul;
        for(std::size_t i{0ul}; i < B.size(); ++i) {
            if(B[i].xs > 1u)
                buffer = bandwidths[i] / (B[i].xs - 1);
            omega[startIndex] = B[i].start;
            omega[startIndex + B[i].xs - 1] = B[i].stop;
            for(std::size_t j{1ul}; j < B[i].xs - 1; ++j)
                omega[startIndex + j] = omega[startIndex + j - 1] + buffer;
            startIndex += B[i].xs;
        }
    }

    template<typename T>
    void refscaling(status_t& status,
            std::vector<T>& newX, std::vector<band_t<T>>& newChebyBands,
            std::vector<band_t<T>>& newFreqBands, std::size_t newXSize,
            std::vector<T>& x, std::vector<band_t<T>>& chebyBands,
            std::vector<band_t<T>>& freqBands)
    {
        std::vector<std::size_t> newDistribution(chebyBands.size());
        for(std::size_t i{0u}; i < chebyBands.size(); ++i)
            newDistribution[i] = 0u;
        std::size_t multipointBands = 0u;
        std::size_t offset = 0u;
        int twoInt = 0;
        for(std::size_t i{0u}; i < chebyBands.size(); ++i)
        {
            newX.push_back(x[offset]);
            if(chebyBands[i].xs > 2u)
            {
                ++multipointBands;
                for(std::size_t j{1u}; j < chebyBands[i].xs - 2u; ++j)
                {
                    newX.push_back((x[offset + j] + x[offset + j + 1]) / 2);
                    newX.push_back(x[offset + j]);
                }
                newX.push_back(x[offset + chebyBands[i].xs - 2u]);
                newX.push_back(x[offset + chebyBands[i].xs - 1u]);
                twoInt += 2;
            }
            else if(chebyBands[i].xs == 2u)
            {
                ++multipointBands;
                ++twoInt;
                newX.push_back(x[offset + 1u]);
            }
            offset += chebyBands[i].xs;
        }
        int threeInt = newXSize - newX.size() - twoInt;
        offset = 0u;
        for(std::size_t i{0u}; i < chebyBands.size(); ++i)
        {
                if(chebyBands[i].xs > 1u)
                {
                    if(threeInt > 0)
                    {
                        newX.push_back(x[offset] + (x[offset + 1] - x[offset]) / 3);
                        T secondValue = x[offset] + (x[offset + 1] - x[offset]) / 3
                            + (x[offset + 1] - x[offset]) / 3;
                        newX.push_back(secondValue);
                        threeInt--;
                        twoInt--;
                    }
                    else if (twoInt > 0)
                    {
                        newX.push_back((x[offset] + x[offset + 1]) / 2);
                        twoInt--;
                    }
                }
            offset += chebyBands[i].xs;
        }
        offset = 0;
        for(std::size_t i{0u}; i < chebyBands.size(); ++i)
        {
                if(chebyBands[i].xs > 2u)
                {
                    if(threeInt > 0)
                    {
                        newX.push_back(x[offset + chebyBands[i].xs - 2u] +
                                (x[offset + chebyBands[i].xs - 1u] -
                                    x[offset + chebyBands[i].xs - 2u]) / 3);
                        T secondValue = x[offset + chebyBands[i].xs - 2u] +
                            (x[offset + chebyBands[i].xs - 1u] -
                                x[offset + chebyBands[i].xs - 2u]) / 3 +
                            (x[offset + chebyBands[i].xs - 1u] -
                                x[offset + chebyBands[i].xs - 2u]) / 3;
                        newX.push_back(secondValue);
                        threeInt--;
                        twoInt--;
                    }
                    else if (twoInt > 0)
                    {
                        newX.push_back((x[offset + chebyBands[i].xs - 2u] +
                                    x[offset + chebyBands[i].xs - 1u]) / 2);
                        twoInt--;
                    }
                }
            offset += chebyBands[i].xs;
        }
        if(newXSize > newX.size()) {
            status = status_t::STATUS_SCALING_INVALID;
            throw std::runtime_error("ERROR: Failed to do reference scaling");
        }

        newX.resize(newXSize);
        std::sort(newX.begin(), newX.end());
        std::size_t total = 0u;
        for(std::size_t i{0u}; i < newX.size(); ++i)
        {
                for(std::size_t j{0u}; j < chebyBands.size(); ++j)
                    if(newX[i] >= chebyBands[j].start && newX[i] <= chebyBands[j].stop)
                    {
                        newDistribution[j]++;
                        ++total;
                    }
        }
        if(total != newXSize) {
            status = status_t::STATUS_SCALING_INVALID;
            throw std::runtime_error("ERROR: Failed to find reference scaling distribution");
        }

        for (std::size_t i{0u}; i < chebyBands.size(); ++i)
        {
            newFreqBands[freqBands.size() - 1u - i].xs = newDistribution[i];
            newChebyBands[i].xs = newDistribution[i];
        }
    }

    template<typename T>
    void split(std::vector<std::pair<T, T>>& subIntervals,
            std::vector<band_t<T>>& chebyBands,
            std::vector<T> &x) {
        std::vector<T> splitpts = x;
        for(std::size_t i{0u}; i < chebyBands.size(); ++i) {
            splitpts.push_back(chebyBands[i].start);
            splitpts.insert(splitpts.end(), chebyBands[i].part.begin(), chebyBands[i].part.end());
            splitpts.push_back(chebyBands[i].stop);
        }

        std::set<T> s(splitpts.begin(), splitpts.end());
        splitpts.assign(s.begin(), s.end());
        std::sort(splitpts.begin(), splitpts.end());

        std::size_t bidx{0u};
        for(std::size_t i{0u}; i < splitpts.size()-1u; ++i) {
            if(splitpts[i+1u] > chebyBands[bidx].stop) {
                ++bidx;
            } else {
                subIntervals.push_back(std::make_pair(splitpts[i], splitpts[i+1u]));
            }
        }
    }

    template<typename T>
    void extrema(status_t& status, T& convergenceOrder,
            T& delta, std::vector<T>& eigenExtrema,
            std::vector<T>& x, std::vector<band_t<T>>& chebyBands,
            std::size_t Nmax, unsigned long prec)
    {
    #ifdef HAVE_MPFR
        mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
        mpfr::mpreal::set_default_prec(prec);
    #endif

        // 1.   Split the initial [-1, 1] interval in subintervals
        //      in order that we can use a reasonable size matrix
        //      eigenvalue solver on the subintervals
        std::vector<std::pair<T, T>> subIntervals;
        std::pair<T, T> dom{std::make_pair(-1.0, 1.0)};

        split(subIntervals, chebyBands, x);

        // 2.   Compute the barycentric variables (i.e., weights)
        //      needed for the current iteration
        std::vector<T> w(x.size());
        baryweights(w, x);

        compdelta(delta, w, x, chebyBands);

        std::vector<T> C(x.size());
        compc(C, delta, x, chebyBands);

        // 3.   Use an eigenvalue solver on each subinterval to find the
        //      local extrema that are located inside the frequency bands
        std::vector<T> chebyNodes(Nmax + 1u);
        equipts(chebyNodes, Nmax + 1u);
        cos(chebyNodes, chebyNodes);

        std::vector<std::pair<T, T>> potentialExtrema;
        std::vector<T> pEx;
        T extremaErrorValueLeft;
        T extremaErrorValueRight;
        T extremaErrorValue;
        comperror(extremaErrorValue, chebyBands[0].start,
                delta, x, C, w, chebyBands);
        potentialExtrema.push_back(std::make_pair(
                chebyBands[0].start, extremaErrorValue));

        for (std::size_t i{0u}; i < chebyBands.size() - 1u; ++i)
        {
            comperror(extremaErrorValueLeft, chebyBands[i].stop,
                    delta, x, C, w, chebyBands);
            comperror(extremaErrorValueRight, chebyBands[i + 1].start,
                    delta, x, C, w, chebyBands);
            bool sgnLeft = pmmath::signbit(extremaErrorValueLeft);
            bool sgnRight = pmmath::signbit(extremaErrorValueRight);
            if (sgnLeft != sgnRight) {
                potentialExtrema.push_back(std::make_pair(
                        chebyBands[i].stop, extremaErrorValueLeft));
                potentialExtrema.push_back(std::make_pair(
                        chebyBands[i + 1u].start, extremaErrorValueRight));
            } else {
                T abs1 = pmmath::fabs(extremaErrorValueLeft);
                T abs2 = pmmath::fabs(extremaErrorValueRight);
                if(abs1 > abs2)
                    potentialExtrema.push_back(std::make_pair(
                            chebyBands[i].stop, extremaErrorValueLeft));
                else
                    potentialExtrema.push_back(std::make_pair(
                            chebyBands[i + 1u].start, extremaErrorValueRight));
            }
        }
        comperror(extremaErrorValue,
                chebyBands[chebyBands.size() - 1u].stop,
                delta, x, C, w, chebyBands);
        potentialExtrema.push_back(std::make_pair(
                chebyBands[chebyBands.size() - 1u].stop,
                extremaErrorValue));

        std::vector<std::vector<T>> pExs(subIntervals.size());
        #pragma omp parallel for
        for (std::size_t i = 0u; i < subIntervals.size(); ++i)
        {
            #ifdef HAVE_MPFR
                mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
                mpfr::mpreal::set_default_prec(prec);
            #endif
            // find the Chebyshev nodes scaled to the current subinterval
            std::vector<T> siCN(Nmax + 1u);
            chgvar(siCN, chebyNodes, subIntervals[i].first,
                    subIntervals[i].second);

            // compute the Chebyshev interpolation function values on the
            // current subinterval
            std::vector<T> fx(Nmax + 1u);
            for (std::size_t j{0u}; j < fx.size(); ++j)
                comperror(fx[j], siCN[j], delta, x, C, w,
                        chebyBands);

            // compute the values of the CI coefficients and those of its
            // derivative
            std::vector<T> c(Nmax + 1u);
            chebcoeffs(c, fx);
            std::vector<T> dc(Nmax);
            diffcoeffs(dc, c);

            // solve the corresponding eigenvalue problem and determine the
            // local extrema situated in the current subinterval
            std::vector<T> eigenRoots;
            roots(eigenRoots, dc, dom);
            if(!eigenRoots.empty()) {
                chgvar(eigenRoots, eigenRoots,
                        subIntervals[i].first, subIntervals[i].second);
                for (std::size_t j{0u}; j < eigenRoots.size(); ++j)
                    pExs[i].push_back(eigenRoots[j]);
            }
            pExs[i].push_back(subIntervals[i].first);
            pExs[i].push_back(subIntervals[i].second);
            #ifdef HAVE_MPFR
                mpfr::mpreal::set_default_prec(prevPrec);
            #endif
        }

        for(std::size_t i{0u}; i < pExs.size(); ++i)
            for(std::size_t j{0u}; j < pExs[i].size(); ++j)
                pEx.push_back(pExs[i][j]);

        std::size_t startingOffset = potentialExtrema.size();
        potentialExtrema.resize(potentialExtrema.size() + pEx.size());
        #pragma omp parallel for
        for(std::size_t i = 0u; i < pEx.size(); ++i)
        {
            #ifdef HAVE_MPFR
                mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
                mpfr::mpreal::set_default_prec(prec);
            #endif
            T valBuffer;
            comperror(valBuffer, pEx[i],
                    delta, x, C, w, chebyBands);
            potentialExtrema[startingOffset + i] = std::make_pair(pEx[i], valBuffer);
            #ifdef HAVE_MPFR
                mpfr::mpreal::set_default_prec(prevPrec);
            #endif
        }

        // sort list of potential extrema in increasing order
        std::sort(potentialExtrema.begin(), potentialExtrema.end(),
                [](const std::pair<T, T>& lhs,
                const std::pair<T, T>& rhs) {
                    return lhs.first < rhs.first;
                });

        eigenExtrema.clear();
        std::size_t extremaIt{0u};
        std::vector<std::pair<T, T>> alternatingExtrema;
        T minError = INT_MAX;
        T maxError = INT_MIN;
        T absError;

        while (extremaIt < potentialExtrema.size())
        {
            std::pair<T, T> maxErrorPoint;
            maxErrorPoint = potentialExtrema[extremaIt];
            while(extremaIt < potentialExtrema.size() - 1u &&
                (pmmath::signbit(maxErrorPoint.second) ==
                pmmath::signbit(potentialExtrema[extremaIt + 1u].second)))
            {
                ++extremaIt;
                if (pmmath::fabs(maxErrorPoint.second) < 
                    pmmath::fabs(potentialExtrema[extremaIt].second))
                    maxErrorPoint = potentialExtrema[extremaIt];
            }
            if (pmmath::isfinite(maxErrorPoint.second)) {
                alternatingExtrema.push_back(maxErrorPoint);
            }
            ++extremaIt;
        }
        std::vector<std::pair<T, T>> bufferExtrema;

        if(alternatingExtrema.size() < x.size())
        {
            status = status_t::STATUS_EXCHANGE_FAILURE;
            std::stringstream message;
            message << "ERROR: The exchange algorithm did not converge\n"
                << "TRIGGER: Not enough alternating extrema\n"
                << "POSSIBLE CAUSE: nmax too small\n";
            convergenceOrder = 2.0;
            throw std::runtime_error(message.str());
        }
        else if (alternatingExtrema.size() > x.size())
        {
            std::size_t remSuperfluous = alternatingExtrema.size() - x.size();
            if (remSuperfluous % 2 != 0)
            {
                if(remSuperfluous == 1u)
                {
                    std::vector<T> x1, x2;
                    x1.push_back(alternatingExtrema[0u].first);
                    for(std::size_t i{1u}; i < alternatingExtrema.size() - 1u; ++i)
                    {
                        x1.push_back(alternatingExtrema[i].first);
                        x2.push_back(alternatingExtrema[i].first);
                    }
                    x2.push_back(alternatingExtrema[alternatingExtrema.size() - 1u].first);
                    T delta1, delta2;
                    compdelta(delta1, x1, chebyBands);
                    compdelta(delta2, x2, chebyBands);
                    delta1 = pmmath::fabs(delta1);
                    delta2 = pmmath::fabs(delta2);
                    std::size_t sIndex{1u};
                    if(delta1 > delta2)
                        sIndex = 0u;
                    for(std::size_t i{sIndex}; i < alternatingExtrema.size() + sIndex - 1u; ++i)
                        bufferExtrema.push_back(alternatingExtrema[i]);
                    alternatingExtrema = bufferExtrema;
                    bufferExtrema.clear();
                }
                else
                {
                    T abs1 = pmmath::fabs(alternatingExtrema[0u].second);
                    T abs2 = pmmath::fabs(alternatingExtrema[alternatingExtrema.size() - 1u].second);
                    std::size_t sIndex = 0u;
                    if (abs1 < abs2)
                        sIndex = 1u;
                    for(std::size_t i{sIndex}; i < alternatingExtrema.size() + sIndex - 1u; ++i)
                        bufferExtrema.push_back(alternatingExtrema[i]);
                    alternatingExtrema = bufferExtrema;
                    bufferExtrema.clear();
                }
            }


            while (alternatingExtrema.size() > x.size())
            {
                std::size_t toRemoveIndex{0u};
                T minValToRemove;
                // change removal rule for extra extrema in case cycling is detected
                if(!cycle) {
                    minValToRemove = pmmath::fmin(pmmath::fabs(alternatingExtrema[0u].second),
                                                pmmath::fabs(alternatingExtrema[1u].second));
                } else {
                    minValToRemove = pmmath::fmax(pmmath::fabs(alternatingExtrema[0u].second),
                                                pmmath::fabs(alternatingExtrema[1u].second));
                }
                T removeBuffer;
                for (std::size_t i{1u}; i < alternatingExtrema.size() - 1u; ++i)
                {
                    if(!cycle) {
                        removeBuffer = pmmath::fmin(pmmath::fabs(alternatingExtrema[i].second),
                                    pmmath::fabs(alternatingExtrema[i + 1u].second));
                    } else {
                        removeBuffer = pmmath::fmax(pmmath::fabs(alternatingExtrema[i].second),
                                    pmmath::fabs(alternatingExtrema[i + 1u].second));
                    }
                    if (removeBuffer < minValToRemove)
                    {
                        minValToRemove = removeBuffer;
                        toRemoveIndex  = i;
                    }
                }
                for (std::size_t i{0u}; i < toRemoveIndex; ++i)
                    bufferExtrema.push_back(alternatingExtrema[i]);
                for (std::size_t i{toRemoveIndex + 2u}; i < alternatingExtrema.size(); ++i)
                    bufferExtrema.push_back(alternatingExtrema[i]);
                alternatingExtrema = bufferExtrema;
                bufferExtrema.clear();
            }
            if(cycle)
                cycle = false;
        }
        if (alternatingExtrema.size() < x.size()) {
            status = status_t::STATUS_EXCHANGE_FAILURE;
            std::stringstream message;
            message << "ERROR: The exchange algorithm did not converge\n"
                << "TRIGGER: Not enough alternating extrema\n"
                << "POSSIBLE CAUSE: nmax too small\n";
            convergenceOrder = 2.0;
            throw std::runtime_error(message.str());
        }

        for (auto& it : alternatingExtrema)
        {
            eigenExtrema.push_back(it.first);
            absError = pmmath::fabs(it.second);
            minError = pmmath::fmin(minError, absError);
            maxError = pmmath::fmax(maxError, absError);
        }

        convergenceOrder = (maxError - minError) / maxError;

        // update the extrema count in each frequency band
        std::size_t bIndex{0u};
        for(std::size_t i{0u}; i < chebyBands.size(); ++i)
        {
            chebyBands[i].xs = 0u;
        }
        for(auto &it : eigenExtrema)
        {
            if(chebyBands[bIndex].start <= it && it <= chebyBands[bIndex].stop)
            {
                ++chebyBands[bIndex].xs;
            }
            else
            {
                ++bIndex;
                ++chebyBands[bIndex].xs;
            }
        }
    #ifdef HAVE_MPFR
        mpfr::mpreal::set_default_prec(prevPrec);
    #endif
    }


    // REMARK: remember that this routine assumes that the information
    // pertaining to the reference x and the frequency bands (i.e., the
    // number of reference values inside each band) is given at the
    // beginning of the execution
    template<typename T>
    pmoutput_t<T> exchange(std::vector<T>& x,
            std::vector<band_t<T>>& chebyBands, double eps,
            std::size_t Nmax, unsigned long prec)
    {
        pmoutput_t<T> output;
        output.status = status_t::STATUS_UNKNOWN_FAILURE;

        std::size_t degree = x.size() - 2u;
        std::sort(x.begin(), x.end(),
                [](const T& lhs,
                const T& rhs) {
                    return lhs < rhs;
                });
        std::vector<T> startX{x};
        cycle = false;
        T qp1 = 0;
        T qp2 = 0;

        output.q = 1;
        output.iter = 0u;
        do {
            ++output.iter;
            extrema(output.status, output.q, output.delta,
                    output.x, startX, chebyBands, Nmax, prec);
            startX = output.x;
            if(output.iter == 1u)
                qp2 = output.q;
            else if(output.iter == 2u)
                qp1 = output.q;
            else {
                // potential cycling, tweak reference update strategy (if
                // reference set candidate is large)
                if(pmmath::fabs(qp2-output.q)/pmmath::fabs(qp2) < eps*1e-5)
                    cycle = true;
                qp2 = qp1;
                qp1 = output.q;
            }
            if(output.q > 1.0)
                break;
        } while (output.q > eps && output.iter <= 100u);
        output.status = status_t::STATUS_SUCCESS;

        if(pmmath::isnan(output.delta) || pmmath::isnan(output.q)) {
            output.status = status_t::STATUS_CONVERGENCE_WARNING;
            std::cerr << "WARNING: The exchange algorithm did not converge.\n"
                << "TRIGGER: numerical instability\n"
                << "POSSIBLE CAUSE: poor starting reference and/or "
                << "a too small value for nmax.\n";
        }

        if(output.iter >= 101u && output.q > eps) {
            output.status = status_t::STATUS_CONVERGENCE_WARNING;
            std::cerr << "WARNING: The exchange algorithm did not converge.\n"
                << "TRIGGER: exceeded iteration threshold of 100\n"
                << "POSSIBLE CAUSE: poor starting reference and/or "
                << "a too small value for nmax.\n";
        }

        output.h.resize(degree + 1u);
        std::vector<T> finalC(output.x.size());
        std::vector<T> finalAlpha(output.x.size());
        baryweights(finalAlpha, output.x);
        T finalDelta = output.delta;
        output.delta = pmmath::fabs(output.delta);
        compc(finalC, finalDelta, output.x, chebyBands);
        std::vector<T> finalChebyNodes(degree + 1);
        equipts(finalChebyNodes, degree + 1);
        cos(finalChebyNodes, finalChebyNodes);
        std::vector<T> fv(degree + 1);

        for (std::size_t i{0u}; i < fv.size(); ++i) {
            approx(fv[i], finalChebyNodes[i], output.x,
                    finalC, finalAlpha);
            if (!pmmath::isfinite(fv[i])) {
                output.status = status_t::STATUS_COEFFICIENT_SET_INVALID;
                std::stringstream message;
                message << "ERROR: Invalid frequency response generated.\n"
                    << "TRIGGER: infinite/NaN values in the final frequency response.\n"
                    << "POSSIBLE CAUSE: too small numerical precision and/or a too "
                    << "small value for nmax.";
                throw std::runtime_error(message.str());
            }
        }

        chebcoeffs(output.h, fv);

        return output;
    }

    template<typename T>
    void parseSpecification(status_t& status,
                std::vector<T> const& f,
                std::vector<T> const& a,
                std::vector<T> const& w)
    {
        if(f.size() != a.size()) {
            status = status_t::STATUS_AMPLITUDE_VECTOR_MISMATCH;
            throw std::domain_error("ERROR: Frequency and amplitude vector sizes do not match");
        }

        if(f.size() % 2 != 0) {
            status = status_t::STATUS_FREQUENCY_INVALID_INTERVAL;
            throw std::domain_error("ERROR: Frequency band edges must come in pairs");
        }

        if(f.size() != w.size() * 2u) {
            status = status_t::STATUS_WEIGHT_VECTOR_MISMATCH;
            std::stringstream message;
            message << "ERROR: Weight vector size does not match the"
                << " the number of frequency bands in the specification";
            throw std::domain_error(message.str());
        }

        for(std::size_t i{0u}; i < f.size() - 1u; ++i) {
            if(f[i] == f[i + 1u] && (a[i] != a[i + 1u])) {
                status = status_t::STATUS_AMPLITUDE_DISCONTINUITY;
                throw std::domain_error("ERROR: Adjacent bands with discontinuities are not allowed");
            }
        
            if(f[i] > f[i + 1u]) {
                status = status_t::STATUS_FREQUENCY_INVALID_INTERVAL;
                throw std::domain_error("ERROR: Frequency vector entries must be nondecreasing");
            }
        }
        for(std::size_t i{0u}; i < w.size(); ++i) {
            if(w[i] <= 0.0) {
                status = status_t::STATUS_WEIGHT_NEGATIVE;
                throw std::domain_error("ERROR: Band weights must be positive");
            }
        }

        if(f[0u] < 0.0 || f[f.size() - 1u] > 1.0) {
            status = status_t::STATUS_FREQUENCY_INVALID_INTERVAL;
            throw std::domain_error("ERROR: Normalized frequency band edges must be between 0 and 1");
        }

        bool pointBands{true};
        for(std::size_t i{0u}; i < f.size() && pointBands; i += 2u)
            if(f[i] != f[i+1u])
                pointBands = false;

        if(pointBands) {
            status = status_t::STATUS_FREQUENCY_INVALID_INTERVAL;
            throw std::domain_error("ERROR: All frequency band intervals are points");
        }
    }

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps,
                std::size_t nmax,
                init_t strategy,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec)
    {
        pmoutput_t<T> output;
        output.status = status_t::STATUS_UNKNOWN_FAILURE;

        try {
            parseSpecification(output.status, f, a, w);
            std::vector<T> h;
            std::vector<band_t<T>> fbands;
            std::vector<band_t<T>> cbands;
            std::vector<std::vector<std::size_t>> bIdx;

            std::size_t nbands{0u};
            bool newBand{true};
            for(std::size_t i{0u}; i < w.size(); ++i) {
                if(newBand) {
                    bIdx.push_back({i});
                    newBand = false;
                }
                if(i < w.size()-1u && f[2u*i+1u] == f[2u*i+2u]) {
                    if(w[i] != w[i+1u]) {
                        output.status = status_t::STATUS_WEIGHT_DISCONTINUITY;
                        throw std::domain_error("ERROR: Incompatible weights for partioned band");
                    }
            
                    bIdx[nbands].push_back(i+1u);
                } else {
                    ++nbands;
                    newBand = true;
                }
            }
            fbands.resize(bIdx.size());

            if(n % 2 != 0) {
                if(f[f.size()-1u] == 1 && a[a.size()-1u] != 0) {
                    std::cerr << "WARNING: gain at Nyquist frequency different from 0.\n"
                        << "Increasing the number of taps by one and passing to a "
                        << "type I filter" << std::endl;
                    ++n;
                }
            }
            std::size_t deg = n / 2;
            if(n % 2 == 0) {            // type I filter
                for(std::size_t i{0u}; i < fbands.size(); ++i) {
                    fbands[i].part.push_back(T(pmmath::const_pi<T>() * f[2u*bIdx[i][0u]]));
                    for(std::size_t j{0u}; j < bIdx[i].size(); ++j) {
                        fbands[i].part.push_back(T(pmmath::const_pi<T>() * f[2u*bIdx[i][j]+1u]));
                    }
                }
                for(std::size_t i{0u}; i < fbands.size(); ++i) {
                    fbands[i].start = fbands[i].part[0u];
                    fbands[i].stop  = fbands[i].part[fbands[i].part.size()-1u];
                    fbands[i].space = space_t::FREQ;

                    fbands[i].amplitude = [i, &a, &bIdx, &fbands](space_t space, T x) -> T {
                        if(space == space_t::CHEBY)
                            x = pmmath::acos(x);
                        for(std::size_t j{0u}; j < bIdx[i].size(); ++j) {
                            if(fbands[i].part[j]*(1.0-1e-12) <= x && 
                            x <= fbands[i].part[j+1u]*(1.0+1e-12)) {
                                if(a[2u*bIdx[i][j]] != a[2u*bIdx[i][j]+1u]) {
                                    return ((x-fbands[i].part[j]) * a[2u*bIdx[i][j]+1u] -
                                            (x-fbands[i].part[j+1u]) * a[2u*bIdx[i][j]]) /
                                            (fbands[i].part[j+1u] - fbands[i].part[j]);
                                } else {
                                    return a[2u*bIdx[i][j]];
                                }
                            }
                        }
                        // this should never happen
                        return a[2u*bIdx[i][0u]];
                    };
                    fbands[i].weight = [i, &w, &bIdx](space_t, T x) -> T {
                        return w[bIdx[i][0u]];
                    };
                }
            } else {                    // type II filter
                for(std::size_t i{0u}; i < fbands.size(); ++i) {
                    fbands[i].part.push_back(T(pmmath::const_pi<T>() * f[2u*bIdx[i][0u]]));
                    for(std::size_t j{0u}; j < bIdx[i].size(); ++j) {
                        if(f[2u*bIdx[i][j]+1u] == 1.0) {
                            if(f[2u*bIdx[i][j]] < 0.9999)
                                fbands[i].part.push_back(T(pmmath::const_pi<T>() * T(0.9999)));
                            else
                                fbands[i].part.push_back(T(pmmath::const_pi<T>() * 
                                                        ((f[2u*bIdx[i][j]]+1u) / 2)));
                        } else {
                            fbands[i].part.push_back(T(pmmath::const_pi<T>() * f[2u*bIdx[i][j]+1u]));
                        }
                    }
                }
                for(std::size_t i{0u}; i < fbands.size(); ++i) {
                    fbands[i].start = fbands[i].part[0u];
                    fbands[i].stop  = fbands[i].part[fbands[i].part.size()-1u];
                    fbands[i].space = space_t::FREQ;

                    fbands[i].amplitude = [i, &a, &bIdx, &fbands](space_t space, T x) -> T {
                        T nx = x;
                        if(space == space_t::CHEBY)
                            nx = pmmath::acos(x);
                        for(std::size_t j{0u}; j < bIdx[i].size(); ++j) {
                            if(fbands[i].part[j]*(1.0-1e-12) <= nx && 
                            nx <= fbands[i].part[j+1u]*(1.0+1e-12)) {
                                if(a[2u*bIdx[i][j]] != a[2u*bIdx[i][j]+1u]) {
                                    return ((nx-fbands[i].part[j]) * a[2u*bIdx[i][j]+1u] -
                                            (nx-fbands[i].part[j+1u]) * a[2u*bIdx[i][j]]) /
                                            (fbands[i].part[j+1u] - fbands[i].part[j]) / 
                                            pmmath::cos(nx/2);
                                }
                            }
                        }
                        if(space == space_t::FREQ)
                            return a[2u*bIdx[i][0u]] / pmmath::cos(x/2);
                        else
                            return a[2u*bIdx[i][0u]] / pmmath::sqrt((x+1)/2);
                    };
                    fbands[i].weight = [i, &w, &bIdx](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return pmmath::cos(x/2) * w[bIdx[i][0u]];
                        else
                            return pmmath::sqrt((x+1)/2) * w[bIdx[i][0u]];
                    };
                }
            }

            std::vector<T> x;
            bandconv(cbands, fbands, convdir_t::FROMFREQ);
            std::function<T(T)> wf = [&cbands](T x) -> T {
                for(std::size_t i{0u}; i < cbands.size(); ++i)
                    if(cbands[i].start <= x && x <= cbands[i].stop)
                        return cbands[i].weight(space_t::CHEBY, x);
                // this should never execute
                return 1.0;
            };

            if(fbands.size() > (deg + 2u) / 4)
                strategy = init_t::AFP;

            switch(strategy) {
                case init_t::UNIFORM:
                {
                    if (fbands.size() <= (deg + 2u) / 4) {
                        std::vector<T> omega;
                        uniform(omega, fbands, deg + 2u);
                        cos(x, omega);
                        bandconv(cbands, fbands, convdir_t::FROMFREQ);
                    } else {
                        // use AFP strategy for very small degrees (wrt nb of bands)
                        std::vector<T> mesh;
                        wam(mesh, cbands, deg);
                        MatrixXd<T> A;
                        chebvand(A, deg+1u, mesh, wf);
                        afp(x, A, mesh);
                        if(x.size() != deg + 2u) {
                            output.status = status_t::STATUS_AFP_INVALID;
                            std::stringstream message;
                            message << "ERROR: AFP strategy failed to produce a valid starting "
                                << "reference\n"
                                << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                            throw std::runtime_error(message.str());
                        }
                        countBand(cbands, x);
                    }
                    output = exchange(x, cbands, eps, nmax, prec);
                } break;
                case init_t::SCALING:
                {
                    std::vector<std::size_t> sdegs(depth+1u);
                    sdegs[depth] = deg;
                    for(int i{(int)depth-1}; i >= 0; --i)
                        sdegs[i] = sdegs[i+1]/2;

                    if(rstrategy == init_t::UNIFORM) {
                        if (fbands.size() <= (sdegs[0] + 2u) / 4) {
                            std::vector<T> omega;
                            uniform(omega, fbands, sdegs[0]+2u);
                            cos(x, omega);
                            bandconv(cbands, fbands, convdir_t::FROMFREQ);
                        } else {
                            // use AFP strategy for very small degrees (wrt nb of bands)
                            std::vector<T> mesh;
                            wam(mesh, cbands, sdegs[0]);
                            MatrixXd<T> A;
                            chebvand(A, sdegs[0]+1u, mesh, wf);
                            afp(x, A, mesh);
                            if(x.size() != sdegs[0] + 2u) {
                                output.status = status_t::STATUS_AFP_INVALID;
                                std::stringstream message;
                                message << "ERROR: AFP strategy failed to produce a valid "
                                    << "starting reference\n"
                                    << "POSSIBLE CAUSE: badly conditioned Chebyshev "
                                    << "Vandermonde matrix";
                                throw std::runtime_error(message.str());
                            }
                            countBand(cbands, x);
                        }
                        output = exchange(x, cbands, eps, nmax, prec);
                    } else { // AFP-based strategy
                        std::vector<T> mesh;
                        wam(mesh, cbands, sdegs[0]);
                        MatrixXd<T> A;
                        chebvand(A, sdegs[0]+1u, mesh, wf);
                        afp(x, A, mesh);
                        if(x.size() != sdegs[0] + 2u) {
                            output.status = status_t::STATUS_AFP_INVALID;
                            std::stringstream message;
                            message << "ERROR: AFP strategy failed to produce a valid "
                                << "starting reference\n"
                                << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                            throw std::runtime_error(message.str());
                        }
                        countBand(cbands, x);
                        output = exchange(x, cbands, eps, nmax, prec);
                    }
                    for(std::size_t i{1u}; i <= depth && output.q <= 0.5; ++i) {
                        x.clear();
                        refscaling(output.status, x, cbands, fbands, sdegs[i]+2u,
                                        output.x, cbands, fbands);
                        output = exchange(x, cbands, eps, nmax, prec);
                    }
                } break;
                default: { // AFP-based initialization
                    std::vector<T> mesh;
                    wam(mesh, cbands, deg);
                    MatrixXd<T> A;
                    chebvand(A, deg+1u, mesh, wf);
                    afp(x, A, mesh);
                    if(x.size() != deg + 2u) {
                        output.status = status_t::STATUS_AFP_INVALID;
                        std::stringstream message;
                        message << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                            << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                        throw std::runtime_error(message.str());
                    }
                    countBand(cbands, x);
                    output = exchange(x, cbands, eps, nmax, prec);
                }
            }

            h.resize(n+1u);
            if (output.h.size() != deg + 1u) {
                output.status = status_t::STATUS_COEFFICIENT_SET_INVALID;
                throw std::runtime_error("ERROR: final filter coefficient set is incomplete");
            }

            if(n % 2 == 0) {
                h[deg] = output.h[0];
                for(std::size_t i{0u}; i < deg; ++i)
                    h[i] = h[n-i] = output.h[deg-i] / 2u;
            } else {
                h[0] = h[n] = output.h[deg] / 4u;
                h[deg] = h[deg+1u] = (output.h[0] * 2 + output.h[1]) / 4u;
                for(std::size_t i{2u}; i < deg + 1u; ++i)
                    h[deg+1u-i] = h[deg+i] = (output.h[i-1u] + output.h[i]) / 4u;
            }
            output.h = h;
        }
        catch (std::domain_error &err) {
            std::cerr << "Invalid specification detected:" << std::endl;
            std::cerr << err.what() << std::endl;
            output.q = 2.0;
        }
        catch (std::runtime_error &err) {
            std::cerr << "Runtime error detected:" << std::endl;
            std::cerr << err.what() << std::endl;
            output.q = 2.0;
        }
        catch (...) {
            std::cerr << "Unknown exception" << std::endl;
            output.status = status_t::STATUS_UNKNOWN_FAILURE;
        }

        return output;
    }

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps,
                std::size_t nmax,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec)
    {
        if( n < 2u*f.size()) {
            std::cerr << "WARNING: too small filter length to use reference scaling." << std::endl
                << "Switching to a uniform initialization strategy." << std::endl;
            return firpm<T>(n, f, a, w, eps, nmax, init_t::UNIFORM, depth, rstrategy, prec);
        } else {
            return firpm<T>(n, f, a, w, eps, nmax, init_t::SCALING, depth, rstrategy, prec);
        }
    }

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps, std::size_t nmax,
                unsigned long prec)
    {
        return firpm<T>(n, f, a, w, eps, nmax, init_t::AFP, 0u, init_t::AFP, prec);
    }

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps,
                std::size_t nmax,
                init_t strategy,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec)
    {
        pmoutput_t<T> output;
        output.status = status_t::STATUS_UNKNOWN_FAILURE;

        try {
            parseSpecification(output.status, f, a, w);
            std::vector<T> h;
            std::vector<band_t<T>> fbands(w.size());
            std::vector<band_t<T>> cbands;
            std::vector<T> fn{f};
            std::size_t deg = n / 2u;
            T sFactor = a[1] / (f[1] * pmmath::const_pi<T>());

            if(n % 2 == 0) { // TYPE III
                if(f[0u] == 0.0) {
                    if(fn[1u] < 1e-5)
                        fn[0u] = fn[1u] / 2;
                    else
                        fn[0u] = 1e-5;
                }
                if(f[f.size() - 1u] == 1.0) {
                    if(f[f.size() - 2u] > 0.9999)
                        fn[f.size() - 1u] = (f[f.size() - 2u] + 1.0) / 2;
                    else
                        fn[f.size() - 1u] = 0.9999;
                }
                --deg;
                fbands[0u].start = pmmath::const_pi<T>() * fn[0u];
                fbands[0u].stop  = pmmath::const_pi<T>() * fn[1u];
                fbands[0u].space = space_t::FREQ;
                if(type == filter_t::FIR_DIFFERENTIATOR) {
                    fbands[0u].weight = [&w](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return (pmmath::sin(x) / x) * w[0u];
                        else
                            return (pmmath::sqrt(T(1.0) - x * x) / pmmath::acos(x)) * w[0u];
                    };
                    fbands[0u].amplitude = [sFactor](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return (x / pmmath::sin(x)) * sFactor;
                        else
                            return (pmmath::acos(x) / pmmath::sqrt(T(1.0) - x * x)) * sFactor;
                    };
                } else { // FIR_HILBERT
                    fbands[0u].weight = [&w](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return pmmath::sin(x) * w[0u];
                        else
                            return pmmath::sqrt(T(1.0) - x * x) * w[0u];
                    };
                    fbands[0u].amplitude = [&a, &fbands](space_t space, T x) -> T {
                        if(space == space_t::CHEBY)
                            x = pmmath::acos(x);
                        if(a[0u] != a[1u])
                            return (((x - fbands[0].start) * a[1u] -
                                    (x - fbands[0].stop) * a[0u]) /
                                    (fbands[0u].stop - fbands[0u].start)) / pmmath::sin(x);
                        return a[0u] / pmmath::sin(x);
                    };
                }
                for(std::size_t i{1u}; i < fbands.size(); ++i) {
                    fbands[i].start = pmmath::const_pi<T>() * fn[2u * i];
                    fbands[i].stop  = pmmath::const_pi<T>() * fn[2u * i + 1u];
                    fbands[i].space = space_t::FREQ;
                    if(type == filter_t::FIR_DIFFERENTIATOR) {
                        fbands[i].weight = [&w, i](space_t space, T x) -> T {
                            if(space == space_t::FREQ)
                                return pmmath::sin(x) * w[i];
                            else
                                return pmmath::sqrt(T(1.0) - x * x) * w[i];
                        };
                        fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                            if(a[2u * i] != a[2u * i + 1u]) {
                                if(space == space_t::CHEBY) {
                                    x = pmmath::acos(x);
                                    return ((x - fbands[i].start) * a[2u * i + 1u] -
                                            (x - fbands[i].stop) * a[2u * i]) /
                                            (fbands[i].stop - fbands[i].start);
                                }
                            }
                            return a[2u * i];
                        };
                    } else { // FIR_HILBERT
                        fbands[i].weight = [&w, i](space_t space, T x) -> T {
                            if(space == space_t::FREQ)
                                return pmmath::sin(x) * w[i];
                            else
                                return pmmath::sqrt(T(1.0) - x * x) * w[i];
                        };
                        fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                            if(space == space_t::CHEBY)
                                x = pmmath::acos(x);
                            if(a[2u * i] != a[2u * i + 1u])
                                return (((x - fbands[i].start) * a[2u * i + 1u] -
                                        (x - fbands[i].stop) * a[2u * i]) /
                                        (fbands[i].stop - fbands[i].start)) / pmmath::sin(x);
                            return a[2u * i] / pmmath::sin(x);
                        };
                    }
                }
            } else { // TYPE IV
                if(f[0u] == 0.0) {
                    if(fn[1u] < 1e-5)
                        fn[0u] = fn[1u] / 2;
                    else
                        fn[0u] = 1e-5;
                }

                fbands[0u].start = pmmath::const_pi<T>() * fn[0u];
                fbands[0u].stop  = pmmath::const_pi<T>() * fn[1u];
                fbands[0u].space = space_t::FREQ;
                if(type == filter_t::FIR_DIFFERENTIATOR) {
                    fbands[0u].weight = [&w](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return (pmmath::sin(x / 2) / x) * w[0u];
                        else
                            return (pmmath::sin(pmmath::acos(x) / 2) / pmmath::acos(x)) * w[0u];
                    };
                    fbands[0u].amplitude = [sFactor](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return (x / pmmath::sin(x / 2)) * sFactor;
                        else
                            return (pmmath::acos(x) / pmmath::sin(pmmath::acos(x) / 2)) * sFactor;
                    };
                } else { // FIR_HILBERT
                    fbands[0u].weight = [&w](space_t space, T x) -> T {
                        if(space == space_t::FREQ)
                            return pmmath::sin(x / 2) * w[0u];
                        else
                            return pmmath::sin(pmmath::acos(x) / 2) * w[0u];
                    };
                    fbands[0u].amplitude = [&fbands, &a](space_t space, T x) -> T {
                        if(space == space_t::CHEBY)
                            x = pmmath::acos(x);
                        if(a[0u] != a[1u])
                            return (((x - fbands[0].start) * a[1u] -
                                    (x - fbands[0].stop) * a[0u]) /
                                    (fbands[0].stop - fbands[0].start)) / pmmath::sin(x / 2);
                        return a[0u] / pmmath::sin(x / 2);
                    };
                }
                for(std::size_t i{1u}; i < fbands.size(); ++i) {
                    fbands[i].start = pmmath::const_pi<T>() * fn[2u * i];
                    fbands[i].stop  = pmmath::const_pi<T>() * fn[2u * i + 1u];
                    fbands[i].space = space_t::FREQ;
                    if(type == filter_t::FIR_DIFFERENTIATOR) {
                        fbands[i].weight = [&w, i](space_t space, T x) -> T {
                            if(space == space_t::FREQ)
                                return pmmath::sin(x / 2) * w[i];
                            else
                                return (pmmath::sin(pmmath::acos(x) / 2)) * w[i];
                        };
                        fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                            if(a[2u * i] != a[2u * i + 1u]) {
                                if(space == space_t::CHEBY)
                                    x = pmmath::acos(x);
                                return ((x - fbands[i].start) * a[2u * i + 1u] -
                                        (x - fbands[i].stop) * a[2u * i]) /
                                        (fbands[i].stop - fbands[i].start);
                            }
                            return a[2u * i];
                        };
                    } else { // FIR_HILBERT
                        fbands[i].weight = [&w, i](space_t space, T x) -> T {
                            if(space == space_t::FREQ)
                                return pmmath::sin(x / 2) * w[i];
                            else
                                return pmmath::sin(pmmath::acos(x) / 2) * w[i];
                        };
                        fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                            if(space == space_t::CHEBY)
                                x = pmmath::acos(x);
                            if(a[2u * i] != a[2u * i + 1u])
                                return (((x - fbands[i].start) * a[2u * i + 1u] -
                                        (x - fbands[i].stop) * a[2u * i]) /
                                        (fbands[i].stop - fbands[i].start)) / pmmath::sin(x / 2);
                            return a[2u * i] / pmmath::sin(x / 2);
                        };
                    }
                }
            }

            std::vector<T> x;
            bandconv(cbands, fbands, convdir_t::FROMFREQ);
            std::function<T(T)> wf = [&cbands](T x) -> T {
                for(std::size_t i{0u}; i < cbands.size(); ++i)
                    if(cbands[i].start <= x && x <= cbands[i].stop)
                        return cbands[i].weight(space_t::CHEBY, x);
                // this should never execute
                return 1.0;
            };

            if(fbands.size() > (deg + 2u) / 4)
                strategy = init_t::AFP;

            switch(strategy) {
                case init_t::UNIFORM:
                {
                    if (fbands.size() <= (deg + 2u) / 4) {
                        std::vector<T> omega;
                        uniform(omega, fbands, deg + 2u);
                        cos(x, omega);
                        bandconv(cbands, fbands, convdir_t::FROMFREQ);
                    } else {
                        // use AFP strategy for very small degrees (wrt nb of bands)
                        std::vector<T> mesh;
                        wam(mesh, cbands, deg);
                        MatrixXd<T> A;
                        chebvand(A, deg+1u, mesh, wf);
                        afp(x, A, mesh);
                        if(x.size() != deg + 2u) {
                            output.status = status_t::STATUS_AFP_INVALID;
                            std::stringstream message;
                            message << "ERROR: AFP strategy failed to produce a valid starting "
                                << "reference\n"
                                << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                            throw std::runtime_error(message.str());
                        }
                        countBand(cbands, x);
                    }
                    output = exchange(x, cbands, eps, nmax, prec);
                } break;
                case init_t::SCALING:
                {
                    std::vector<std::size_t> sdegs(depth+1u);
                    sdegs[depth] = deg;
                    for(int i{(int)depth-1}; i >= 0; --i)
                        sdegs[i] = sdegs[i+1]/2;

                    if(rstrategy == init_t::UNIFORM) {
                        if (fbands.size() <= (sdegs[0] + 2u) / 4) {
                            std::vector<T> omega;
                            uniform(omega, fbands, sdegs[0]+2u);
                            cos(x, omega);
                            bandconv(cbands, fbands, convdir_t::FROMFREQ);
                        } else {
                            // use AFP strategy for very small degrees (wrt nb of bands)
                            std::vector<T> mesh;
                            wam(mesh, cbands, sdegs[0]);
                            MatrixXd<T> A;
                            chebvand(A, sdegs[0]+1u, mesh, wf);
                            afp(x, A, mesh);
                            if(x.size() != sdegs[0] + 2u) {
                                output.status = status_t::STATUS_AFP_INVALID;
                                std::stringstream message;
                                message << "ERROR: AFP strategy failed to produce a valid starting "        << "reference\n"
                                    << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                                throw std::runtime_error(message.str());
                            }
                            countBand(cbands, x);
                        }
                        output = exchange(x, cbands, eps, nmax);
                    } else { // AFP-based strategy
                        std::vector<T> mesh;
                        wam(mesh, cbands, sdegs[0]);
                        MatrixXd<T> A;
                        chebvand(A, sdegs[0]+1u, mesh, wf);
                        afp(x, A, mesh);
                        if(x.size() != sdegs[0] + 2u) {
                            output.status = status_t::STATUS_AFP_INVALID;
                            std::stringstream message;
                            message << "ERROR: AFP strategy failed to produce a valid starting "            << "reference\n"
                                << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                            throw std::runtime_error(message.str());
                        }
                        countBand(cbands, x);
                        output = exchange(x, cbands, eps, nmax);
                    }
                    for(std::size_t i{1u}; i <= depth && output.q <= 0.5; ++i) {
                        x.clear();
                        refscaling(output.status, x, cbands, fbands, sdegs[i]+2u,
                                        output.x, cbands, fbands);
                        output = exchange(x, cbands, eps, nmax, prec);
                    }
                } break;
                default: { // AFP-based initialization
                    std::vector<T> mesh;
                    wam(mesh, cbands, deg);
                    MatrixXd<T> A;
                    chebvand(A, deg+1u, mesh, wf);
                    afp(x, A, mesh);
                    if(x.size() != deg + 2u) {
                        output.status = status_t::STATUS_AFP_INVALID;
                        std::stringstream message;
                        message << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                            << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix";
                        throw std::runtime_error(message.str());
                    }
                    countBand(cbands, x);
                    output = exchange(x, cbands, eps, nmax, prec);
                }
            }

            h.resize(n + 1u);
            if (output.h.size() != deg + 1u) {
                output.status = status_t::STATUS_COEFFICIENT_SET_INVALID;
                throw std::runtime_error("ERROR: final filter coefficient set is incomplete");
            }

            if(n % 2 == 0)
            {
                h[deg + 1u] = 0;
                h[deg] = (output.h[0u] * 2.0 - output.h[2]) / 4u;
                h[deg + 2u] = -h[deg];
                h[1u] = output.h[deg - 1u] / 4;
                h[2u * deg + 1u] = -h[1u];
                h[0u] =  output.h[deg] / 4;
                h[2u * (deg + 1u)] = -h[0u];
                for(std::size_t i{2u}; i < deg; ++i)
                {
                    h[deg + 1u - i] = (output.h[i - 1u] - output.h[i + 1u]) / 4;
                    h[deg + 1u + i] = -h[deg + 1u - i];
                }
            } else {
                ++deg;
                h[deg - 1u] = (output.h[0u] * 2.0 - output.h[1u]) / 4;
                h[deg] = -h[deg - 1u];
                h[0u] = output.h[deg - 1u] / 4;
                h[2u * deg - 1u] = -h[0u];
                for(std::size_t i{2u}; i < deg; ++i)
                {
                    h[deg - i] = (output.h[i - 1u] - output.h[i]) / 4;
                    h[deg + i - 1u] = -h[deg - i];
                }
            }
            output.h = h;
        }
        catch (const std::domain_error& err) {
            std::cerr << "Invalid specification detected:" << std::endl;
            std::cerr << err.what() << std::endl;
            output.q = 2.0;
        }
        catch (const std::runtime_error& err) {
            std::cerr << "Runtime error detected:" << std::endl;
            std::cerr << err.what() << std::endl;
            output.q = 2.0;
        }
        catch (...) {
            std::cerr << "Unknown exception" << std::endl;
            output.status = status_t::STATUS_UNKNOWN_FAILURE;
        }

        return output;
    }

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps,
                std::size_t nmax, std::size_t depth,
                init_t rstrategy,
                unsigned long prec)
    {
        if( n < 2u*f.size()) {
            std::cerr << "WARNING: too small filter length to use reference scaling." << std::endl
                << "Switching to a uniform initialization strategy" << std::endl;
            return firpm<T>(n, f, a, w, type, eps, nmax, init_t::UNIFORM, depth, rstrategy, prec);
        } else {
            return firpm<T>(n, f, a, w, type, eps, nmax, init_t::SCALING, depth, rstrategy, prec);
        }

    }

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps,
                std::size_t nmax,
                unsigned long prec)
    {
        return firpm<T>(n, f, a, w, type, eps,
                    nmax, init_t::AFP, 0u, init_t::AFP, prec);
    }

    
    
} // namespace pm