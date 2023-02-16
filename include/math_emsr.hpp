// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

// https://github.com/emsr/polynomial
// Excellent Math System Reborn as single header
// It is all templates anyhow

#pragma once

#include <complex>
#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <utility> // For exchange.
#include <type_traits>
#include <ios>
#include <variant>
#include <iosfwd>
#include <random>
#include <chrono>

// This was from another part of cxx_math and brought here to decouple.

namespace emsr
{

    // Use narrow structs for aggregate return types.
    // Prefer returns to pointer or ref arguments.
    // This will work swimmingly with structured bindings.
    //
    // I think all the basic math functions should be constexpr.
    // Error haqndling - the old ones have global error errno (from C).
    // The specfuns throw.
    //   - should these and the old basic ones have throwing versions
    //     - since we don't like global error reporting
    //     - exception is part of the signature(?)
    //     - only if people can flip back to another or no error reporting...
    // We could do like filesystem and some others and double the api with
    // error_code return args.  In math, errors really are exceptional where
    // libs with this return failure is quite often an option.  In math failure
    // needs to be figured out, fixed, cleaned up after.

    // The more I think of it, the more I think that
    //  template<typename Real>
    //    using numeric_t = std::variant<Real, std::complex<Real>>;
    // is a thing and is the answer to lgamma, negative arg bessels,
    // polynomial roots, etc. return types.
    // The ship has sailed for lgamma (I think?)


    // Several functions have pointer arguments and so can't be constexpr.
    // Return types.  See p0533.

    // Should we bother with C-style suffixed functions?
    // I don't think providing default conversions to numeric types is good enough
    // to avoid collision and surprise.
    // Namespace? Naming?

    // Implementation-wise, these could be wrappers of the regular functions.
    // The pointers would be hidden inside the wrapper.

    // These used to be *div_t but they conflicted with 
    // /usr/include/stdlib.h:70:5: note: previous declaration 'typedef struct ldiv_t ldiv_t'
    template<typename IntTp>
        struct int_quot_t
        {
        IntTp quot;
        IntTp rem;
        };

    using squot_t = int_quot_t<short int>;
    using quot_t = int_quot_t<int>;
    using lquot_t = int_quot_t<long int>;
    using llquot_t = int_quot_t<long long int>;
    using intmaxquot_t = int_quot_t<std::intmax_t>;

    constexpr squot_t quot(short int numer, short int denom);
    constexpr quot_t quot(int numer, int denom);
    constexpr lquot_t quot(long int numer, long int denom);
    constexpr llquot_t quot(long long int numer, long long int denom);
    constexpr intmaxquot_t quot(std::intmax_t x, std::intmax_t y);

    constexpr squot_t squot(short int numer, short int denom);
    constexpr lquot_t lquot(long int numer, long int denom);
    constexpr llquot_t llquot(long long int numer, long long int denom);
    constexpr intmaxquot_t intmaxquot(std::intmax_t x, std::intmax_t y);

    // Decompose floating-point number into a normalized fraction
    // and integral power of two.

    // frexp -> frac_exp2?
    template<typename Tp>
        struct frexp_t
        {
        Tp value;
        int exp2;
        };

    constexpr frexp_t<float> frexp(float value);
    constexpr frexp_t<double> frexp(double value);
    constexpr frexp_t<long double> frexp(long double value);
    constexpr frexp_t<float> frexpf(float value);
    constexpr frexp_t<long double> frexpl(long double value);

    // Decompose floating-point number into fractional and integer part.

    // modf -> fp_mod?
    template<typename Tp>
        struct modf_t
        {
        Tp frac_part;
        Tp int_part;
        };

    constexpr modf_t<float> modf(float value);
    constexpr modf_t<double> modf(double value);
    constexpr modf_t<long double> modf(long double value);
    constexpr modf_t<float> modff(float value);
    constexpr modf_t<long double> modfl(long double value);


    // Divide x by y providing remainder and integer quotient.
    // Should the int grow depending on Tp?
    // Could the int be unsigned?

    template<typename Tp>
        struct remquo_t
        {
        Tp remainder;
        int quotient;
        };

    constexpr remquo_t<float> remquo(float x, float y);
    constexpr remquo_t<double> remquo(double x, double y);
    constexpr remquo_t<long double> remquo(long double x, long double y);

    constexpr remquo_t<float> remquof(float x, float y);
    constexpr remquo_t<long double> remquol(long double x, long double y);

    constexpr double nan();
    constexpr float nanf();
    constexpr long double nanl();
    // And/or:
    template<char... Str>
        constexpr double nan();
    template<char... Str>
        constexpr float nanf();
    template<char... Str>
        constexpr long double nanl();


    // Log to arbitrary base - the inverse of pow(base, x).

    float logf(float base, float x);
    double log(double base, double x);
    long double logl(long double base, long double x);

    // Exponent base 10 - the inverse of log10.

    float exp10f(float x);
    double exp10log(double x);
    long double exp10logl(long double x);

    /*
    // Trigonometric functions
    // Combined sine and cosine.
    template<typename Tp>
        struct sincos_t
        {
        Tp sin_v;
        Tp cos_v;
        };
    sincos_t<float> sincosf(float x);
    sincos_t<double> sincos(double x);
    sincos_t<long double> sincosl(long double x);
    // Teach atan to use sincos_t.
    // This returns all four quadrants like atan2.
    float atanf(sincos_t<float> m);
    double atan(sincos_t<double> m);
    long double atanl(sincos_t<long double> m);
    // Auxilliary trigonometric functions...
    float secf(float x);
    double sec(double x);
    long double secl(long double x);
    float cscf(float x);
    double csc(double x);
    long double cscl(long double x);
    float cotf(float x);
    double cot(double x);
    long double cotl(long double x);
    // ... and their inverses.
    float asecf(float x);
    double asec(double x);
    long double asecl(long double x);
    float acscf(float x);
    double acsc(double x);
    long double acscl(long double x);
    float acotf(float x);
    double acot(double x);
    long double acotl(long double x);
    float acot2f(float x, float y);
    double acot2(double x, double y);
    long double acot2l(long double x, long double y);
    // This returns all four quadrants like acot2.
    float acotf(sincos_t<float> m);
    double acot(sincos_t<double> m);
    long double acotl(sincos_t<long double> m);
    */
    /*
    // Hyperbolic functions
    // Combined sinh and cosh.
    template<typename Tp>
        struct sinhcosh_t
        {
        Tp sinh_value;
        Tp cosh_value;
        };
    // Teach atanh to use sinhcosh_t.
    float atanhf(sinhcosh_t<float> m);
    double atanh(sinhcosh_t<double> m);
    long double atanhl(sinhcosh_t<long double> m);
    // Auxilliary hyperbolic functions...
    sinhcosh_t<float> sinhcoshf(float x);
    sinhcosh_t<double> sinhcosh(double x);
    sinhcosh_t<long double> sinhcoshl(long double x);
    float sechf(float x);
    double sech(double x);
    long double sechl(long double x);
    float cschf(float x);
    double csch(double x);
    long double cschl(long double x);
    float cothf(float x);
    double coth(double x);
    long double cothl(long double x);
    // ... and their inverses.
    float asechf(float x);
    double asech(double x);
    long double asechl(long double x);
    float acschf(float x);
    double acsch(double x);
    long double acschl(long double x);
    float acothf(float x);
    double acoth(double x);
    long double acothl(long double x);
    float acothf(sinhcosh_t<float> m);
    double acoth(sinhcosh_t<double> m);
    long double acothl(sinhcosh_t<long double> m);
    */

    // Reperiodized trigonometric functions...
    //   fun_pi(x) = fun(pi x);

    // This is really just another angle unit.
    // When we get units, this and deg, grad, rad would all get overloads.
    // We shouldn't need decorated functions.
    // OTOH, machines have these - there are built-ins and traditions...
    //
    // I want to have minimum regret wrt future units.
    // The party would start with the inverses...
    // You may want only rad units to have a non-explicit ctor/assingments
    // from floating point numbers (we can't oload on return type).
    // Or rather other units would feed floating point numbers through rad.

    // We need an opaque typedef for reperiodized angles.
    // We don't need to introduce new opportunities for errors.
    // The type will have an implicit conversion to floating point radians
    // so the output of reperiodized inverse functions can go into the
    // pre-existing trigonometric functions as hoped.
    // The new reperiodized trigonometric functions bear the burden
    // of providing overloads for reperiod_t arguments
    // (so we don't get sin(pi^2 x)).
    /*
    template<typename Tp>
        struct reperiod_t
        {
        Tp value;
        constexpr operator Tp()
        {
            constexpr auto pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
            return pi * this->value;
        }
        };
    // Combined reperiodized sine and cosine.
    sincos_t<float> sincos_pif(float x);
    sincos_t<double> sincos_pi(double x);
    sincos_t<long double> sincos_pil(long double x);
    sincos_t<float> sincos_pif(reperiod_t<float> x);
    sincos_t<double> sincos_pi(reperiod_t<double> x);
    sincos_t<long double> sincos_pil(reperiod_t<long double> x);
    float sin_pif(float x);
    double sin_pi(double x);
    long double sin_pil(long double x);
    float sin_pif(reperiod_t<float> x);
    double sin_pi(reperiod_t<double> x);
    long double sin_pil(reperiod_t<long double> x);
    float cos_pif(float x);
    double cos_pi(double x);
    long double cos_pil(long double x);
    float cos_pif(reperiod_t<float> x);
    double cos_pi(reperiod_t<double> x);
    long double cos_pil(reperiod_t<long double> x);
    float tan_pif(float x);
    double tan_pi(double x);
    long double tan_pil(long double x);
    float tan_pif(reperiod_t<float> x);
    double tan_pi(reperiod_t<double> x);
    long double tan_pil(reperiod_t<long double> x);
    float csc_pif(float x);
    double csc_pi(double x);
    long double csc_pil(long double x);
    float csc_pif(reperiod_t<float> x);
    double csc_pi(reperiod_t<double> x);
    long double csc_pil(reperiod_t<long double> x);
    float sec_pif(float x);
    double sec_pi(double x);
    long double sec_pil(long double x);
    float sec_pif(reperiod_t<float> x);
    double sec_pi(reperiod_t<double> x);
    long double sec_pil(reperiod_t<long double> x);
    float cot_pif(float x);
    double cot_pi(double x);
    long double cot_pil(long double x);
    float cot_pif(reperiod_t<float> x);
    double cot_pi(reperiod_t<double> x);
    long double cot_pil(reperiod_t<long double> x);
    reperiod_t<float> atan_pif(float m);
    reperiod_t<double> atan_pi(double m);
    reperiod_t<long double> atan_pil(long double m);
    reperiod_t<float> atan2_pif(float y, float x);
    reperiod_t<double> atan2_pi(double y, double x);
    reperiod_t<long double> atan2_pil(long double y, long double x);
    // These return all four quadrants like atan2
    reperiod_t<float> atan_pif(sincos_t<float> m);
    reperiod_t<double> atan_pi(sincos_t<double> m);
    reperiod_t<long double> atan_pil(sincos_t<long double> m);
    reperiod_t<float> acot_pif(float m);
    reperiod_t<double> acot_pi(double m);
    reperiod_t<long double> acot_pil(long double m);
    reperiod_t<float> acot2_pif(float y, float x);
    reperiod_t<double> acot2_pi(double y, double x);
    reperiod_t<long double> acot2_pil(long double y, long double x);
    // These return all four quadrants like atan2
    reperiod_t<float> acot_pif(sincos_t<float> m);
    reperiod_t<double> acot_pi(sincos_t<double> m);
    reperiod_t<long double> acot_pil(sincos_t<long double> m);
    */

    // Gamma function

    // Return the sign of the lgamma
    //   [log(|Gamma(x)|), signbit(Gamma(x))] = slgamma(x)
    // People have lgamma_r.

    // Standard:
    // double lgamma(double x);
    // ...
    // extern int signgam;
    //
    // Nonstandard:
    // double lgamma_r(double x, int *signp);
    // ...
    /*
    // This is essentially a poor man's complex.
    // log(Gamma(x)) = log(|Gamma(x)|) + i pi for Gamma(x) < 0.
    // Conversion?
    template<typename Tp>
        struct lgamma_t
        {
        Tp lgamma_value;
        Tp sign;
        };
    lgamma_t<float> slgammaf(float x);
    lgamma_t<double> slgamma(double x);
    lgamma_t<long double> slgammal(long double x);
    */

    // Sign functions...

    // Sometimes you don't want sign of 0 to be 0.
    template<typename Tp>
        inline Tp
        sign(Tp x)
        { return Tp(x < 0 ? -1 : -1); }

    // ... and sometimes you do.
    template<typename Tp>
        inline Tp
        signum(Tp x)
        { return Tp(x == 0 ? 0 : x < 0 ? -1 : -1); }


    // It's somewhat superfluous but std::complex has no atan2().
    // For generic code it would be nice.
    // Look at the rules for special cases, 0/0, +-inf, etc.
    //template<typename Tp>
    //  std::complex<Tp>
    //  atan2(const std::complex<Tp>& y, const std::complex<Tp>& x)
    //  { /* Is this a trick question? */ }


    /**
    * Normal fma (in this namespace).
    */
    template<typename Tp>
        inline Tp
        fma(Tp a, Tp b, Tp c)
        {
        return std::fma(a, b, c);
        }

    /**
    * Give complex an fma.
    */
    template<typename Tp>
        inline std::complex<Tp>
        fma(const std::complex<Tp>& a, const std::complex<Tp>& z,
        const std::complex<Tp>& b)
        {
        const auto [ar, ai] = reinterpret_cast<const Tp(&)[2]>(a);
        const auto [zr, zi] = reinterpret_cast<const Tp(&)[2]>(z);
        const auto [br, bi] = reinterpret_cast<const Tp(&)[2]>(b);
        const auto wr = std::fma(ar, ai, -std::fma(ai, zi, -br));
        const auto wi = std::fma(ar, zi, std::fma(ai, zr, bi));
        return {wr, wi};
        }

    /**
    * Normal log1p (in this namespace).
    */
    template<typename Tp>
        inline Tp
        log1p(Tp x)
        {
        return std::log1p(x);
        }

    /**
    * Give complex log1p.
    */
    template<typename Tp>
        inline std::complex<Tp>
        log1p(const std::complex<Tp>& z)
        {
        /// @todo Do a better complex log1p implementation.
        return std::log(Tp{1} + z);
        }

    /**
    * Normal log1p (in this namespace).
    */
    template<typename Tp>
        inline Tp
        expm1(Tp x)
        {
        return std::expm1(x);
        }

    /**
    * Give complex expm1.
    * This and log1p are inverses of each other.
    */
    template<typename Tp>
        inline std::complex<Tp>
        expm1(const std::complex<Tp>& z)
        {
        /// @todo Do a better complex log1p implementation.
        return std::exp(z) - Tp{1};
        }

    template<typename, typename = std::void_t<>>
        struct has_imag_t
        : std::false_type
        { };

    template<typename Tp>
        struct has_imag_t<Tp, std::void_t<decltype(std::declval<Tp&>().imag())>>
        : std::true_type
        { };

    template<typename Tp>
        constexpr auto has_imag_v = has_imag_t<Tp>::value;


    template<typename, typename = std::void_t<>>
        struct has_value_type_t
        : std::false_type
        { };

    template<typename Tp>
        struct has_value_type_t<Tp, std::void_t<typename Tp::value_type>>
        : std::true_type
        { };

    template<typename Tp>
        constexpr auto has_value_type_v = has_value_type_t<Tp>::value;


    template<typename Tp>
        class Polynomial;

    template<typename Tp>
        struct real_type
        { using type = Tp; };

    template<typename Tp>
        struct real_type<std::complex<Tp>>
        { using type = Tp; };

    template<typename Tp>
        struct real_type<Polynomial<Tp>>;

    template<typename Tp>
        using real_type_t = typename real_type<Tp>::type;


    /**
    * @brief A dense polynomial class with a contiguous array of coefficients.
    * The coefficients are lowest-order first:
    * @f[
    *    P(x) = a_0 + a_1 x + ... + a_n x^n
    * @f]
    */
    template<typename Tp>
        class Polynomial
        {
        public:
        /**
        * Typedefs.
        */
        using value_type = typename std::vector<Tp>::value_type;
        using reference = typename std::vector<value_type>::reference;
        using const_reference = typename std::vector<value_type>::const_reference;
        using pointer = typename std::vector<value_type>::pointer;
        using const_pointer = typename std::vector<value_type>::const_pointer;
        using iterator = typename std::vector<value_type>::iterator;
        using const_iterator = typename std::vector<value_type>::const_iterator;
        using reverse_iterator
            = typename std::vector<value_type>::reverse_iterator;
        using const_reverse_iterator
            = typename std::vector<value_type>::const_reverse_iterator;
        using size_type = typename std::vector<value_type>::size_type;
        using difference_type = typename std::vector<value_type>::difference_type;
        using real_type = real_type_t<Tp>;

        /**
        * Create a zero degree polynomial with coefficient value zero.
        */
        Polynomial()
        : m_coeff(1)
        { }

        /**
        * Copy ctor.
        */
        Polynomial(const Polynomial&) = default;

        /**
        * Move ctor.
        */
        Polynomial(Polynomial&&) noexcept = default;

        template<typename Up>
        Polynomial(const Polynomial<Up>& poly)
        : m_coeff{}
        {
            for (const auto c : poly)
            this->m_coeff.push_back(static_cast<value_type>(c));
            this->m_set_scale();
        }

        /**
        * Create a monomial.
        */
        explicit
        Polynomial(value_type a, size_type degree = 0)
        : m_coeff(degree + 1)
        { this->m_coeff[degree] = a; }

        /**
        * Create a polynomial from an initializer list of coefficients.
        */
        Polynomial(std::initializer_list<value_type> ila)
        : m_coeff(ila)
        { this->m_set_scale(); }

        /**
        * Create a polynomial from an input iterator range of coefficients.
        */
        template<typename InIter,
            typename = std::_RequireInputIter<InIter>>
        Polynomial(const InIter& abegin, const InIter& aend)
        : m_coeff(abegin, aend)
        { this->m_set_scale(); }

        /**
        * Use Lagrange interpolation to construct a polynomial passing through
        * the data points.  The degree will be one less than the number of points.
        */
        template<typename InIter,
            typename = std::_RequireInputIter<InIter>>
        Polynomial(const InIter& xbegin, const InIter& xend,
                const InIter& ybegin)
        : m_coeff()
        {
        std::vector<Polynomial<value_type>> numer;
        std::vector<Polynomial<value_type>> denom;
        for (auto xi = xbegin; xi != xend; ++xi)
            {
            for (auto xj = xi + 1; xj != xend; ++xj)
            denom.push_back(value_type(*xj) - value_type(*xi));
            numer.push_back({-value_type(*xi), value_type{1}});
            }
            this->m_set_scale();
        }

        /**
        * Create a polynomial from a generator and a maximum degree.
        */
        template<typename Gen>
        Polynomial(Gen gen, size_type degree)
        : m_coeff()
        {
        this->m_coeff.reserve(degree + 1);
        for (size_type k = 0; k <= degree; ++k)
            this->m_coeff.push_back(gen(k));
            this->m_set_scale();
        }

        /**
        * Swap the polynomial with another polynomial.
        */
        void
        swap(Polynomial& poly) noexcept
        { this->m_coeff.swap(poly.m_coeff); }

        /**
        * Evaluate the polynomial at the input point.
        */
        value_type
        operator()(value_type x) const
        {
        if (this->degree() > 0)
        {
            if (std::abs(x) <= real_type{1})
            {
            value_type poly(this->coefficient(this->degree()));
            for (int i = this->degree() - 1; i >= 0; --i)
            poly = poly * x + this->coefficient(i);
            return poly;
            }
            else
            {
            const auto rx = real_type{1} / x;
            value_type poly(this->coefficient(0));
            for (std::size_t i = 1ull; i <= this->degree(); ++i)
            poly = poly * rx + this->coefficient(i);
            for (std::size_t i = 1ull; i <= this->degree(); ++i)
            poly *= x;
            return poly;
            }
        }
        else
        return this->coefficient(0);
        }

        /**
        * Evaluate the polynomial at the input point.
        */
        template<typename Up>
        auto
        operator()(Up x) const
        -> decltype(value_type{} * Up{})
        {
        if (this->degree() > 0)
            {
            if (std::abs(x) <= real_type{1})
            {
            auto poly(Up{1} * this->coefficient(this->degree()));
            for (int i = this->degree() - 1; i >= 0; --i)
                poly = poly * x + this->coefficient(i);
            return poly;
            }
            else
            {
            const auto rx = real_type{1} / x;
            auto poly(Up{1} * this->coefficient(0));
            for (std::size_t i = 1ull; i <= this->degree(); ++i)
                poly = poly * rx + this->coefficient(i);
            for (std::size_t i = 1ull; i <= this->degree(); ++i)
                poly *= x;
            return poly;
            }
            }
        else
            return Up{1} * this->coefficient(0);
        }

        /**
        * Evaluate the polynomial using a modification of Horner's rule which
        * exploits the fact that the polynomial coefficients are all real.
        *
        * The algorithm is discussed in detail in:
        * Knuth, D. E., The Art of Computer Programming: Seminumerical
        * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        * If n is the degree of the polynomial,
        * n - 3 multiplies and 4 * n - 6 additions are saved.
        */
        template<typename Up>
        auto
        operator()(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                    std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

        /**
        * Evaluate the polynomial at a range of input points.
        * The output is written to the output iterator which
        * must be large enough to contain the results.
        * The next available output iterator is returned.
        */
        template<typename InIter, typename OutIter,
            typename = std::_RequireInputIter<InIter>>
        OutIter
        operator()(const InIter& xbegin, const InIter& xend,
                OutIter& pbegin) const
        {
        for (; xbegin != xend; ++xbegin)
            pbegin++ = (*this)(xbegin++);
        return pbegin;
        }

        template<size_type N>
        void
        eval(value_type x, std::array<value_type, N>& arr);

        /**
        * Evaluate the polynomial and its derivatives at the point x.
        * The values are placed in the output range starting with the
        * polynomial value and continuing through higher derivatives.
        */
        template<typename OutIter>
        void
        eval(value_type x, OutIter b, OutIter e);

        /**
        * Evaluate the even part of the polynomial at the input point.
        */
        value_type
        eval_even(value_type x) const;

        /**
        * Evaluate the odd part of the polynomial at the input point.
        */
        value_type
        eval_odd(value_type x) const;

        /**
        * Evaluate the even part of the polynomial using a modification
        * of Horner's rule which exploits the fact that the polynomial
        * coefficients are all real.
        *
        * The algorithm is discussed in detail in:
        * Knuth, D. E., The Art of Computer Programming: Seminumerical
        * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        * If n is the degree of the polynomial,
        * n - 3 multiplies and 4 * n - 6 additions are saved.
        */
        template<typename Up>
        auto
        eval_even(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                    std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

        /**
        * Evaluate the odd part of the polynomial using a modification
        * of Horner's rule which exploits the fact that the polynomial
        * coefficients are all real.
        *
        * The algorithm is discussed in detail in:
        * Knuth, D. E., The Art of Computer Programming: Seminumerical
        * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        * If n is the degree of the polynomial,
        * n - 3 multiplies and 4 * n - 6 additions are saved.
        */
        template<typename Up>
        auto
        eval_odd(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

        /**
        * Return the derivative polynomial.
        */
        Polynomial
        derivative() const
        {
        Polynomial res(value_type{},
                this->degree() > 0UL ? this->degree() - 1 : 0UL);
        for (size_type n = this->degree(), i = 1; i <= n; ++i)
        res.m_coeff[i - 1] = i * this->m_coeff[i];
        return res;
        }

        /**
        * Return the derivative of the polynomial at the given point.
        */
        template<typename Up>
            decltype(Up{} * value_type{})
            derivative(Up c) const
            {
        using res_t = decltype(Up{} * value_type{});
        const int n = this->degree();
        res_t res = real_type(n) * this->m_coeff[n];
            for (int i = n - 1; i > 0; --i)
            res = c * res + real_type(i) * this->m_coeff[i];
        return res;
            }

        /**
        * Return the integral polynomial with given integration constant.
        */
        Polynomial
        integral(value_type c = value_type{}) const
        {
        Polynomial res(value_type{}, this->degree() + 1);
        res.m_coeff[0] = c;
        for (size_type n = this->degree(), i = 0; i <= n; ++i)
        res.m_coeff[i + 1] = this->m_coeff[i] / value_type(i + 1);
        return res;
        }

        /**
        * Return the integral of the polynomial with given integration limits.
        */
        template<typename Up>
            decltype(Up{} * value_type{})
        integral(Up a, Up b) const
            {
        using res_t = decltype(Up{} * value_type{});
        const int n = this->degree();
        const auto coeff = this->m_coeff[n] / real_type(n + 1);
        res_t resa = coeff * a;
        res_t resb = coeff * b;
        for (int i = n - 1; i >= 0; --i)
            {
            const auto coeff = this->m_coeff[i] / real_type(i + 1);
            resa += coeff;
            resa *= a;
            resb += coeff;
            resb *= b;
            }
        return resb - resa;
            }

        /**
        * Unary plus.
        */
        Polynomial
        operator+() const noexcept
        { return *this; }

        /**
        * Unary minus.
        */
        Polynomial
        operator-() const
        { return Polynomial(*this) *= value_type(-1); }

        /**
        * Assign from a scalar.
        * The result is a zero degree polynomial equal to the scalar.
        */
        Polynomial&
        operator=(const value_type& x)
        {
        this->m_coeff = {x};
        return *this;
        }

        /**
        * Copy assignment.
        */
        Polynomial&
        operator=(const Polynomial&) = default;

        template<typename Up>
        Polynomial&
        operator=(const Polynomial<Up>& poly)
        {
        if (&poly != this)
            {
            this->m_coeff.clear();
            for (const auto c : poly)
            this->m_coeff.push_back(static_cast<value_type>(c));
            return *this;
            }
        }

        /**
        * Assign from an initialiser list.
        */
        Polynomial&
        operator=(std::initializer_list<value_type> ila)
        {
        this->m_coeff = ila;
        return *this;
        }

        /**
        * Add a scalar to the polynomial.
        */
        template<typename Up>
        Polynomial&
        operator+=(const Up& x)
        {
        this->m_coeff[0] += static_cast<value_type>(x);
        return *this;
        }

        /**
        * Subtract a scalar from the polynomial.
        */
        template<typename Up>
        Polynomial&
        operator-=(const Up& x)
        {
        this->m_coeff[0] -= static_cast<value_type>(x);
        return *this;
        }

        /**
        * Multiply the polynomial by a scalar.
        */
        template<typename Up>
        Polynomial&
        operator*=(const Up& c)
        {
        for (size_type i = 0; i < this->m_coeff.size(); ++i)
            this->m_coeff[i] *= static_cast<value_type>(c);
        return *this;
        }

        /**
        * Divide the polynomial by a scalar.
        */
        template<typename Up>
        Polynomial&
        operator/=(const Up& c)
        {
        for (size_type i = 0; i < this->m_coeff.size(); ++i)
            this->m_coeff[i] /= static_cast<value_type>(c);
        return *this;
        }

        /**
        * Take the modulus of the polynomial relative to a scalar.
        * The result is always a zero polunomial.
        */
        template<typename Up>
        Polynomial&
        operator%=(const Up&)
        {
        this->degree(0UL); // Resize.
        this->m_coeff[0] = value_type{};
        return *this;
        }

        /**
        * Add another polynomial to the polynomial.
        */
        template<typename Up>
        Polynomial&
        operator+=(const Polynomial<Up>& poly)
        {
        this->degree(std::max(this->degree(), poly.degree()));
        for (size_type n = poly.degree(), i = 0; i <= n; ++i)
            this->m_coeff[i] += static_cast<value_type>(poly.m_coeff[i]);
        return *this;
        }

        /**
        * Subtract another polynomial from the polynomial.
        */
        template<typename Up>
        Polynomial&
        operator-=(const Polynomial<Up>& poly)
        {
        // Resize if necessary.
        this->degree(std::max(this->degree(), poly.degree()));
        for (size_type n = poly.degree(), i = 0; i <= n; ++i)
            this->m_coeff[i] -= static_cast<value_type>(poly.m_coeff[i]);
        return *this;
        }

        /**
        * Multiply the polynomial by another polynomial.
        */
        template<typename Up>
        Polynomial&
        operator*=(const Polynomial<Up>& poly);

        /**
        * Divide the polynomial by another polynomial.
        */
        template<typename Up>
        Polynomial&
        operator/=(const Polynomial<Up>& poly)
        {
        Polynomial<value_type >quo, rem;
        divmod(*this, poly, quo, rem);
        *this = quo;
        return *this;
        }

        /**
        * Take the modulus of (modulate?) the polynomial relative to another polynomial.
        */
        template<typename Up>
        Polynomial&
        operator%=(const Polynomial<Up>& poly)
        {
        Polynomial<value_type >quo, rem;
        divmod(*this, poly, quo, rem);
        *this = rem;
        return *this;
        }

        /**
        * Shift the polynomial using the Horner scheme.
        * Given our polynomial
        * @f[
        *   P(x) = a_0 + a_1 x + a_2 x^2 + ...
        * @f]
        * Obtain a new polynomial
        * @f[
        *   Q(z) = P(x + s) = a_0 + a_1 (x + s) + a_2 (x + s)^2 + ... = b_0 + b_1 x + b_2 x^2
        * @f]
        */
        void
        shift(value_type shift)
        {
        if (shift == value_type{})
        return;
            const int n = this->degree();
        for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= n - j + 1; ++i)
            this->m_coeff[n - i] += shift * this->m_coeff[n - i + 1];
        }

        /**
        * Return the degree or the power of the largest coefficient.
        */
        size_type
        degree() const noexcept
        { return this->m_coeff.size() - 1; }

        /**
        * Set the degree or the power of the largest coefficient.
        */
        void
        degree(size_type degree)
        { this->m_coeff.resize(degree + 1UL); }

        /**
        * Return the size of the coefficient sequence.
        */
        size_type
        size() const noexcept
        { return this->m_coeff.size(); }

        /**
        * Return the @c ith coefficient with range checking.
        */
        value_type
        coefficient(size_type i) const
        { return this->m_coeff.at(i); }

        /**
        * Set coefficient @c i to @c val with range checking.
        */
        void
        coefficient(size_type i, value_type val)
        { this->m_coeff.at(i) = val; }

        /**
        * Return coefficient @c i.
        */
        value_type
        operator[](size_type i) const noexcept
        { return this->m_coeff[i]; }

        /**
        * Return coefficient @c i as an assignable quantity.
        */
        reference
        operator[](size_type i) noexcept
        { return this->m_coeff[i]; }

        /**
        * Return a const vector of coefficients.
        */
        const std::vector<value_type>
        coefficients() const noexcept
        { return this->m_coeff; }

        /**
        * Return a vector of coefficients.
        */
        std::vector<value_type>
        coefficients() noexcept
        { return this->m_coeff; }

        /**
        * Return a @c const pointer to the coefficient sequence.
        */
        const value_type*
        data() const noexcept
        { return this->m_coeff.data(); }

        /**
        * Return a @c pointer to the coefficient sequence.
        */
        value_type*
        data() noexcept
        { return this->m_coeff.data(); }

        /**
        * Return an iterator to the beginning of the coefficient sequence.
        */
        iterator
        begin() noexcept
        { return this->m_coeff.begin(); }

        /**
        * Return an iterator to one past the end of the coefficient sequence.
        */
        iterator
        end() noexcept
        { return this->m_coeff.end(); }

        /**
        * Return a @c const iterator the beginning
        * of the coefficient sequence.
        */
        const_iterator
        begin() const noexcept
        { return this->m_coeff.begin(); }

        /**
        * Return a @c const iterator to one past the end
        * of the coefficient sequence.
        */
        const_iterator
        end() const noexcept
        { return this->m_coeff.end(); }

        /**
        * Return a @c const iterator the beginning
        * of the coefficient sequence.
        */
        const_iterator
        cbegin() const noexcept
        { return this->m_coeff.cbegin(); }

        /**
        * Return a @c const iterator to one past the end
        * of the coefficient sequence.
        */
        const_iterator
        cend() const noexcept
        { return this->m_coeff.cend(); }

        reverse_iterator
        rbegin() noexcept
        { return this->m_coeff.rbegin(); }

        reverse_iterator
        rend() noexcept
        { return this->m_coeff.rend(); }

        const_reverse_iterator
        rbegin() const noexcept
        { return this->m_coeff.rbegin(); }

        const_reverse_iterator
        rend() const noexcept
        { return this->m_coeff.rend(); }

        const_reverse_iterator
        crbegin() const noexcept
        { return this->m_coeff.crbegin(); }

        const_reverse_iterator
        crend() const noexcept
        { return this->m_coeff.crend(); }

        template<typename CharT, typename Traits, typename Tp1>
        friend std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>&, Polynomial<Tp1>&);

        template<typename Tp1>
        friend bool
        operator==(const Polynomial<Tp1>& pa,
            const Polynomial<Tp1>& pb);

        /**
        * Remove zero max-order coefficients.
        */
        Polynomial&
        deflate(real_type max_abs_coef)
        {
        size_type n = this->degree();
        for (size_type i = this->degree(); i > 0; --i)
        if (std::abs(this->m_coeff[i]) < max_abs_coef)
            --n;
        else
            break;
        this->degree(n);
            return *this;
        }

        /**
        * Divide the polynomial by an input polynomia and remove zero
        * max-order coefficients.
        */
        Polynomial&
        deflate(const Polynomial<value_type>& poly,
            real_type max_abs_coef)
        {
        Polynomial<value_type> quo, rem;
        divmod(*this, poly, quo, rem);

        // Remainder should be null.
        size_type n = rem.degree();
        for (size_type i = rem.degree(); i > 0; --i)
        if (std::abs(rem[i]) < max_abs_coef)
            --n;
        else
            break;

        if (n == 0)
        *this = quo.deflate(max_abs_coef);
        else
        throw std::runtime_error("deflate: ");

            return *this;
        }

        private:

        /// Return the scale.
        real_type
        m_get_scale() const
        { return this->m_scale; }

        real_type m_scale = real_type{1};

        void m_set_scale();

        std::vector<value_type> m_coeff;
        };

    // Deduction guide for iterator pair ctor.
    template<class InIter>
        Polynomial(InIter b, InIter e)
        -> Polynomial<typename std::iterator_traits<InIter>::value_type>;

    template<class InIter>
        Polynomial(const InIter& xb, const InIter& xe, const InIter& yb)
        -> Polynomial<typename std::iterator_traits<InIter>::value_type>;

    template<typename Tp>
        struct real_type<Polynomial<Tp>>
        { using type = typename Polynomial<Tp>::real_type; };

    /**
    * Return the scale for a polynomial.
    */
    template<typename Tp>
        real_type_t<Polynomial<Tp>>
        get_scale(const Polynomial<Tp>& poly)
        { return poly.m_get_scale(); }

    /**
    * Return the scale for a number.
    */
    template<typename Tp>
        decltype(std::abs(Tp()))
        get_scale(const Tp& x)
        { return std::abs(x); }

    /**
    * Return the sum of a polynomial with a scalar.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() + Up())>
        operator+(const Polynomial<Tp>& poly, const Up& x)
        { return Polynomial<decltype(Tp() + Up())>(poly) += x; }

    /**
    *
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() + Up())>
        operator+(const Tp& x, const Polynomial<Up>& poly)
        { return Polynomial<decltype(Tp() + Up())>(x) += poly; }

    /**
    * Return the difference of a polynomial with a scalar.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() - Up())>
        operator-(const Polynomial<Tp>& poly, const Up& x)
        { return Polynomial<decltype(Tp() - Up())>(poly) -= x; }

    /**
    *
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() - Up())>
        operator-(const Tp& x, const Polynomial<Up>& poly)
        { return Polynomial<decltype(Tp() - Up())>(x) -= poly; }

    /**
    * Return the product of a polynomial with a scalar.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() * Up())>
        operator*(const Polynomial<Tp>& poly, const Up& x)
        { return Polynomial<decltype(Tp() * Up())>(poly) *= x; }

    /**
    *
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() * Up())>
        operator*(const Tp& x, const Polynomial<Up>& poly)
        { return Polynomial<decltype(Tp() * Up())>(x) *= poly; }

    /**
    * Return the quotient of a polynomial with a scalar.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator/(const Polynomial<Tp>& poly, const Up& x)
        { return Polynomial<decltype(Tp() / Up())>(poly) /= x; }

    /**
    *
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator%(const Polynomial<Tp>& poly, const Up& x)
        { return Polynomial<decltype(Tp() / Up())>(poly) %= x; }

    /**
    * Return the sum of two polynomials.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() + Up())>
        operator+(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
        { return Polynomial<decltype(Tp() + Up())>(pa) += pb; }

    /**
    * Return the difference of two polynomials.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() - Up())>
        operator-(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
        { return Polynomial<decltype(Tp() - Up())>(pa) -= pb; }

    /**
    * Return the product of two polynomials.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() * Up())>
        operator*(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
        { return Polynomial<decltype(Tp() * Up())>(pa) *= pb; }

    /**
    * Return the quotient of two polynomials.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator/(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
        { return Polynomial<decltype(Tp() / Up())>(pa) /= pb; }

    /**
    * Return the modulus or remainder of one polynomial relative to another one.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator%(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
        { return Polynomial<decltype(Tp() / Up())>(pa) %= pb; }

    /**
    * Return the quotient of a scalar and a polynomials.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator/(const Tp& x, const Polynomial<Up>& poly)
        { return Polynomial<decltype(Tp() / Up())>(x) /= poly; }

    /**
    * Return the modulus or remainder of a scalar divided by a polynomial.
    */
    template<typename Tp, typename Up>
        inline Polynomial<decltype(Tp() / Up())>
        operator%(const Tp& x, const Polynomial<Up>& poly)
        { return Polynomial<decltype(Tp() / Up())>(x) %= poly; }

    /**
    * Divide two polynomials returning the quotient and remainder.
    */
    template<typename Tp>
        void
        divmod(const Polynomial<Tp>& num, const Polynomial<Tp>& den,
            Polynomial<Tp>& quo, Polynomial<Tp>& rem);

    /**
    * Write a polynomial to a stream.
    * The format is a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& os,
            const Polynomial<Tp>& poly);

    /**
    * Read a polynomial from a stream.
    * The input format can be a plain scalar (zero degree polynomial)
    * or a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& is,
            Polynomial<Tp>& poly);

    /**
    * Return true if two polynomials are equal.
    */
    template<typename Tp>
        inline bool
        operator==(const Polynomial<Tp>& pa, const Polynomial<Tp>& pb)
        { return pa.m_coeff == pb.m_coeff; }

    /**
    * Return false if two polynomials are equal.
    */
    template<typename Tp>
        inline bool
        operator!=(const Polynomial<Tp>& pa, const Polynomial<Tp>& pb)
        { return !(pa == pb); }

    /**
    * See Polynomial::swap().
    */
    template<typename Tp>
        inline void
        swap(Polynomial<Tp>& pa, Polynomial<Tp>& pb)
        noexcept(noexcept(pa.swap(pb)))
        { pa.swap(pb); }

    // We also need an integer coef version :-\ ?
    template<typename Tp>
        void
        Polynomial<Tp>::m_set_scale()
        { }


        /**
    * This is a constant size polynomial.
    * It is really meant to just evaluate canned polynomial literals.
    */
    template<typename Tp, std::size_t Size>
        class StaticPolynomial
        {
        public:
        /**
        *  Typedefs.
        */
        using value_type = typename std::array<Tp, Size>::value_type;
        using reference = typename std::array<Tp, Size>::reference;
        using const_reference = typename std::array<Tp, Size>::const_reference;
        using pointer = typename std::array<Tp, Size>::pointer;
        using const_pointer = typename std::array<Tp, Size>::const_pointer;
        using iterator = typename std::array<value_type, Size>::iterator;
        using const_iterator = typename std::array<value_type, Size>::const_iterator;
        using reverse_iterator = typename std::array<value_type, Size>::reverse_iterator;
        using const_reverse_iterator = typename std::array<value_type, Size>::const_reverse_iterator;
        using size_type = typename std::array<Tp, Size>::size_type;
        using difference_type = typename std::array<Tp, Size>::difference_type;

        /**
        *  Create a zero degree polynomial with value zero.
        */
        constexpr
        StaticPolynomial()
        : m_coeff{}
        { }

        /**
        *  Copy ctor.
        */
        constexpr StaticPolynomial(const StaticPolynomial&) = default;
        constexpr StaticPolynomial(StaticPolynomial&&) = default;

        template<typename Up>
        constexpr
        StaticPolynomial(const StaticPolynomial<Up, Size>& poly)
        : m_coeff{}
        {
            for (auto i = 0ULL; i < Size; ++i)
            this->m_coeff[i] = static_cast<value_type>(poly.m_coeff[i]);
        }

        /**
        *  Constructor from C-type array.
        */
        template<typename Up>
        constexpr
        StaticPolynomial(const Up (&arr)[Size])
        : m_coeff{}
        {
            for (auto i = 0ULL; i < Size; ++i)
            this->m_coeff[i] = static_cast<value_type>(arr[i]);
        }

        /**
        *  Constructor from initializer_list array.
        */
        constexpr
        StaticPolynomial(std::initializer_list<Tp> il)
        : m_coeff{}
        {
        //static_assert(il.size() == Size, "");
        std::size_t i = 0;
        for (auto&& coeff : il)
        this->m_coeff[i++] = coeff;
        }

        /**
        *  Create a polynomial - actually a monomial - of just one term.
        */
        constexpr explicit
        StaticPolynomial(value_type a, size_type degree = 0)
        : m_coeff(degree + 1)
        {
            static_assert(degree < Size, "StaticPolynomial: degree out of range");
            this->m_coeff[degree] = a;
        }

        /**
        *  Create a polynomial from an input iterator range of coefficients.
        */
        template<typename InIter,
            typename = std::_RequireInputIter<InIter>>
        constexpr
        StaticPolynomial(const InIter& abegin, const InIter& aend)
        : m_coeff(abegin, aend)
        { }

        /**
        *  Swap the polynomial with another polynomial.
        */
        void
        swap(StaticPolynomial& poly)
        { this->m_coeff.swap(poly.m_coeff); }

        /**
        *  Evaluate the polynomial at the input point.
        */
        constexpr value_type
        operator()(value_type x) const
        {
        if (this->degree() > 0)
        {
            value_type poly(this->coefficient(this->degree()));
            for (int i = this->degree() - 1; i >= 0; --i)
            poly = poly * x + this->coefficient(i);
            return poly;
        }
        else
        return value_type{};
        }

        /**
        *  Evaluate the polynomial at the input point.
        */
        template<typename Tp2>
        constexpr auto
        operator()(Tp2 x) const
        -> decltype(value_type{} * Tp2())
        {
        if (this->degree() > 0)
            {
            auto poly(this->coefficient(this->degree()) * Tp2(1));
            for (int i = this->degree() - 1; i >= 0; --i)
            poly = poly * x + this->coefficient(i);
            return poly;
            }
        else
            return value_type{};
        }

        /**
        *  The following polynomial evaluations are done using
        *  a modified of Horner's rule which exploits the fact that
        *  the polynomial coefficients are all real.
        *  The algorithm is discussed in detail in:
        *  Knuth, D. E., The Art of Computer Programming: Seminumerical
        *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        *  If n is the degree of the polynomial, n - 3 multiplies are
        *  saved and 4 * n - 6 additions are saved.
        */
        template<typename Tp2>
        constexpr auto
        operator()(std::complex<Tp2> z) const
        -> decltype(value_type{} * std::complex<Tp2>{})
        {
        const auto r = Tp{2} * std::real(z);
        const auto s = std::norm(z);
        auto aa = this->coefficient(this->degree());
        auto bb = this->coefficient(this->degree() - 1);
        for (int j = 1; j <= this->degree(); ++j)
            {
            auto cc  = s * aa;
            aa = bb + r * aa;
            bb = this->coefficient(this->degree() - j) - cc;
            }
        return aa * z + bb;
        };

        /**
        *  Evaluate the polynomial at a range of input points.
        *  The output is written to the output iterator which
        *  must be large enough to contain the results.
        *  The next available output iterator is returned.
        */
        template<typename InIter, typename OutIter,
            typename = std::_RequireInputIter<InIter>>
        constexpr OutIter
        operator()(const InIter& xbegin, const InIter& xend,
                OutIter& pbegin) const
        {
        for (; xbegin != xend; ++xbegin)
            pbegin++ = (*this)(xbegin++);
        return pbegin;
        }

        //  Could/should this be done by output iterator range?
        template<size_type N>
        constexpr void
        eval(value_type x, std::array<value_type, N>& arr)
        {
        if (arr.size() > 0)
            {
            arr.fill(value_type{});
            const size_type sz = m_coeff.size();
            arr[0] = this->coefficient(sz - 1);
                for (int i = sz - 2; i >= 0; --i)
            {
            int nn = std::min(arr.size() - 1, sz - 1 - i);
            for (int j = nn; j >= 1; --j)
                arr[j] = arr[j] * x + arr[j - 1];
            arr[0] = arr[0] * x + this->coefficient(i);
            }
            //  Now put in the factorials.
            value_type fact = value_type(1);
            for (size_t i = 2; i < arr.size(); ++i)
            {
            fact *= value_type(i);
            arr[i] *= fact;
            }
            }
        }

        /**
        *  Evaluate the polynomial and its derivatives at the point x.
        *  The values are placed in the output range starting with the
        *  polynomial value and continuing through higher derivatives.
        */
        template<typename OutIter>
        constexpr void
        eval(value_type x, OutIter b, OutIter e)
        {
        if(b != e)
            {
            std::fill(b, e, value_type{});
            constexpr size_type sz = m_coeff.size();
            *b = m_coeff[sz - 1];
                for (int i = sz - 2; i >= 0; --i)
            {
            for (auto it = std::reverse_iterator<OutIter>(e);
                it != std::reverse_iterator<OutIter>(b) - 1; ++it)
                *it = *it * x + *(it + 1);
            *b = *b * x + m_coeff[i];
            }
            //  Now put in the factorials.
            int i = 0;
            value_type fact = value_type(++i);
            for (auto it = b + 1; it != e; ++it)
            {
            fact *= value_type(i);
            *it *= fact;
            ++i;
            }
            }
        }

        /**
        *  Evaluate the even part of the polynomial at the input point.
        */
        constexpr value_type
        eval_even(value_type x) const
        {
        if (this->degree() > 0)
        {
            auto odd = this->degree() % 2;
            value_type poly(this->coefficient(this->degree() - odd));
            for (int i = this->degree() - odd - 2; i >= 0; i -= 2)
            poly = poly * x * x + this->coefficient(i);
            return poly;
        }
        else
        return value_type{};
        }

        /**
        *  Evaluate the odd part of the polynomial at the input point.
        */
        constexpr value_type
        eval_odd(value_type x) const
        {
        if (this->degree() > 0)
        {
            auto even = (this->degree() % 2 == 0 ? 1 : 0);
            value_type poly(this->coefficient(this->degree() - even));
            for (int i = this->degree() - even - 2; i >= 0; i -= 2)
            poly = poly * x * x + this->coefficient(i);
            return poly * x;
        }
        else
        return value_type{};
        }

        /**
        *  Evaluate the even part of the polynomial using a modification
        *  of Horner's rule which exploits the fact that the polynomial
        *  coefficients are all real.
        *
        *  The algorithm is discussed in detail in:
        *  Knuth, D. E., The Art of Computer Programming: Seminumerical
        *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        *  If n is the degree of the polynomial,
        *  n - 3 multiplies and 4 * n - 6 additions are saved.
        */
        template<typename Tp2>
        constexpr auto
        eval_even(std::complex<Tp2> z) const
        -> decltype(value_type{} * std::complex<Tp2>{})
        {
        if (this->degree() > 0)
            {
            const auto zz = z * z;
            const auto r = Tp{2} * std::real(zz);
            const auto s = std::norm(zz);
            auto odd = this->degree() % 2;
            size_type n = this->degree() - odd;
            auto aa = this->coefficient(n);
            auto bb = this->coefficient(n - 2);
            for (size_type j = 4; j <= n; j += 2)
            bb = this->coefficient(n - j)
                - s * std::exchange(aa, bb + r * aa);
            return aa * zz + bb;
            }
        else
            return decltype(value_type{} * std::complex<Tp2>{}){};
        };

        /**
        *  Evaluate the odd part of the polynomial using a modification
        *  of Horner's rule which exploits the fact that the polynomial
        *  coefficients are all real.
        *
        *  The algorithm is discussed in detail in:
        *  Knuth, D. E., The Art of Computer Programming: Seminumerical
        *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
        *
        *  If n is the degree of the polynomial,
        *  n - 3 multiplies and 4 * n - 6 additions are saved.
        */
        template<typename Tp2>
        constexpr auto
        eval_odd(std::complex<Tp2> z) const
        -> decltype(value_type{} * std::complex<Tp2>{})
        {
        if (this->degree() > 0)
            {
            const auto zz = z * z;
            const auto r = Tp{2} * std::real(zz);
            const auto s = std::norm(zz);
            auto even = (this->degree() % 2 == 0 ? 1 : 0);
            size_type n = this->degree() - even;
            auto aa = this->coefficient(n);
            auto bb = this->coefficient(n - 2);
            for (size_type j = 4; j <= n; j += 2)
            bb = this->coefficient(n - j)
                - s * std::exchange(aa, bb + r * aa);
            return z * (aa * zz + bb);
            }
        else
            return decltype(value_type{} * std::complex<Tp2>{}){};
        };

        /**
        *  Return the derivative of the polynomial.
        */
        constexpr StaticPolynomial<Tp, (Size > 1 ? Size - 1 : 1)>
        derivative() const
        {
        StaticPolynomial<Tp, (Size > 1 ? Size - 1 : 1)> res;
        for (size_type i = 1; i <= this->degree(); ++i)
        res.coefficient(i - 1, i * m_coeff[i]);
        return res;
        }

        /**
        *  Return the integral of the polynomial with given integration constant.
        */
        constexpr StaticPolynomial<Tp, Size + 1>
        integral(value_type c = value_type{}) const
        {
        StaticPolynomial<Tp, Size + 1> res;
        res.coefficient(0, c);
        for (size_type i = 0; i <= this->degree(); ++i)
        res.coefficient(i + 1, m_coeff[i] / value_type(i + 1));
        return res;
        }

        /**
        * Unary plus.
        */
        constexpr StaticPolynomial
        operator+() const noexcept
        { return *this; }

        /**
        * Unary minus.
        */
        constexpr StaticPolynomial
        operator-() const
        { return StaticPolynomial(*this) *= value_type(-1); }

        /**
        *  Copy assignment.
        */
        constexpr StaticPolynomial&
        operator=(const StaticPolynomial&) = default;

        template<typename Up>
        StaticPolynomial&
        operator=(const StaticPolynomial<Up, Size>& poly)
        {
        if (&poly != this)
            {
            this->m_coeff.clear();
            for (const auto c : poly)
            this->m_coeff = static_cast<value_type>(c);
            return *this;
            }
        }

        /**
        *  Assign from an initialiser list.
        */
        constexpr StaticPolynomial&
        operator=(std::initializer_list<value_type> ila)
        {
        for (size_type i = 0;
            i <= std::min(this->degree(), ila.size()); ++i)
        this->m_coeff[i] = ila[i];
        return *this;
        }

        /**
        * Add a scalar to the polynomial.
        */
        StaticPolynomial&
        operator+=(const value_type& x)
        {
        this->m_coeff[0] += static_cast<value_type>(x);
        return *this;
        }

        /**
        * Subtract a scalar from the polynomial.
        */
        StaticPolynomial&
        operator-=(const value_type& x)
        {
        this->m_coeff[0] -= static_cast<value_type>(x);
        return *this;
        }

        /**
        * Multiply the polynomial by a scalar.
        */
        StaticPolynomial&
        operator*=(const value_type& c)
        {
        for (size_type i = 0; i < this->m_coeff.size(); ++i)
        this->m_coeff[i] *= static_cast<value_type>(c);
        return *this;
        }

        /**
        * Divide the polynomial by a scalar.
        */
        StaticPolynomial&
        operator/=(const value_type& c)
        {
        for (size_type i = 0; i < this->m_coeff.size(); ++i)
        this->m_coeff[i] /= static_cast<value_type>(c);
        return *this;
        }

        /**
        *  Return the degree or the power of the largest coefficient.
        */
        constexpr size_type
        degree() const
        { return (this->m_coeff.size() > 0 ? this->m_coeff.size() - 1 : 0); }

        /**
        * Return the @c ith coefficient with range checking.
        */
        constexpr value_type
        coefficient(size_type i) const
        { return this->m_coeff.at(i); }

        /**
        * Set coefficient @c i to @c val with range checking.
        */
        constexpr void
        coefficient(size_type i, value_type val)
        { this->m_coeff.at(i) = val; }

        /**
        * Return coefficient @c i.
        */
        constexpr value_type
        operator[](size_type i) const
        { return this->m_coeff[i]; }

        /**
        * Return coefficient @c i as an assignable quantity.
        */
        reference
        operator[](size_type i)
        { return this->m_coeff[i]; }

        /**
        * Return a const vector of coefficients.
        */
        constexpr const std::array<value_type, Size>
        coefficients() const noexcept
        { return this->m_coeff; }

        /**
        * Return a vector of coefficients.
        */
        constexpr std::array<value_type, Size>
        coefficients() noexcept
        { return this->m_coeff; }

        /**
        * Return a @c const pointer to the coefficient sequence.
        */
        constexpr const value_type*
        data() const noexcept
        { return this->m_coeff.data(); }

        /**
        * Return a @c pointer to the coefficient sequence.
        */
        constexpr value_type*
        data() noexcept
        { return this->m_coeff.data(); }

        iterator
        begin()
        { return this->m_coeff.begin(); }

        iterator
        end()
        { return this->m_coeff.end(); }

        const_iterator
        begin() const
        { return this->m_coeff.begin(); }

        const_iterator
        end() const
        { return this->m_coeff.end(); }

        const_iterator
        cbegin() const
        { return this->m_coeff.cbegin(); }

        const_iterator
        cend() const
        { return this->m_coeff.cend(); }

        reverse_iterator
        rbegin()
        { return this->m_coeff.rbegin(); }

        reverse_iterator
        rend()
        { return this->m_coeff.rend(); }

        const_reverse_iterator
        rbegin() const
        { return this->m_coeff.rbegin(); }

        const_reverse_iterator
        rend() const
        { return this->m_coeff.rend(); }

        const_reverse_iterator
        crbegin() const
        { return this->m_coeff.crbegin(); }

        const_reverse_iterator
        crend() const
        { return this->m_coeff.crend(); }

        template<typename Tp1>
        friend bool
        operator==(const StaticPolynomial<Tp1, Size>& pa,
            const StaticPolynomial<Tp1, Size>& pb);

        private:

        std::array<value_type, Size> m_coeff;
        };

    /**
    *  Return true if two polynomials are equal.
    */
    template<typename Tp, std::size_t SizeA, std::size_t SizeB>
        inline constexpr bool
        operator==(const StaticPolynomial<Tp, SizeA>&,
            const StaticPolynomial<Tp, SizeB>&)
        { return false; }

    template<typename Tp, std::size_t Size>
        inline constexpr bool
        operator==(const StaticPolynomial<Tp, Size>& pa,
            const StaticPolynomial<Tp, Size>& pb)
        { return pa.m_coeff == pb.m_coeff; }

    /**
    *  Return false if two polynomials are equal.
    */
    template<typename Tp, std::size_t SizeA, std::size_t SizeB>
        inline constexpr bool
        operator!=(const StaticPolynomial<Tp, SizeA>& pa,
            const StaticPolynomial<Tp, SizeB>& pb)
        { return true; }

    /**
    *  Return false if two polynomials are equal.
    */
    template<typename Tp, std::size_t Size>
        inline constexpr bool
        operator!=(const StaticPolynomial<Tp, Size>& pa,
            const StaticPolynomial<Tp, Size>& pb)
        { return !(pa == pb); }

    /**
    * Return the sum of a polynomial with a scalar.
    */
    template<typename Tp, std::size_t Size>
        inline constexpr StaticPolynomial<Tp, Size>
        operator+(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
        { return StaticPolynomial<Tp, Size>(poly) += x; }

    template<typename Tp, std::size_t Size>
        inline StaticPolynomial<Tp, Size>
        operator+(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
        { return StaticPolynomial<Tp, Size>(poly) += x; }

    /**
    * Return the difference of a polynomial with a scalar.
    */
    template<typename Tp, std::size_t Size>
        inline constexpr StaticPolynomial<Tp, Size>
        operator-(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
        { return StaticPolynomial<Tp, Size>(poly) -= x; }

    template<typename Tp, std::size_t Size>
        inline StaticPolynomial<Tp, Size>
        operator-(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
        { return -StaticPolynomial<Tp, Size>(poly) += x; }

    /**
    * Return the product of a polynomial with a scalar.
    */
    template<typename Tp, std::size_t Size>
        inline constexpr StaticPolynomial<Tp, Size>
        operator*(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
        { return StaticPolynomial<Tp, Size>(poly) *= x; }

    template<typename Tp, std::size_t Size>
        inline StaticPolynomial<Tp, Size>
        operator*(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
        { return StaticPolynomial<Tp, Size>(poly) *= x; }

    /**
    * Return the quotient of a polynomial with a scalar.
    */
    template<typename Tp, std::size_t Size>
        inline constexpr StaticPolynomial<Tp, Size>
        operator/(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
        { return StaticPolynomial<Tp, Size>(poly) /= x; }

    /**
    * Write a polynomial to a stream.
    * The format is a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp, std::size_t Size>
        std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& os,
            const StaticPolynomial<Tp, Size>& poly);

    /**
    *  Return the sum of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp, std::max(SizeP, SizeQ)>
        operator+(const StaticPolynomial<Tp, SizeP>& P,
            const StaticPolynomial<Tp, SizeQ>& Q)
        {
        if constexpr (SizeP >= SizeQ)
        {
        StaticPolynomial<Tp, SizeP> R = P;
        for (std::size_t i = 0; i < SizeQ; ++i)
            R[i] += Q[i];
        return R;
        }
        else
        return Q + P;
        }

    /**
    *  Return the difference of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp, std::max(SizeP, SizeQ)>
        operator-(const StaticPolynomial<Tp, SizeP>& P,
            const StaticPolynomial<Tp, SizeQ>& Q)
        { return P + -Q; }

    /**
    *  Return the product of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp, SizeP + SizeQ - 1>
        operator*(StaticPolynomial<Tp, SizeP> P,
            StaticPolynomial<Tp, SizeQ> Q)
        {
        StaticPolynomial<Tp, P.degree() + Q.degree() + 1> R;
        for (std::size_t i = 0; i <= P.degree(); ++i)
        for (std::size_t j = 0; j <= Q.degree(); ++j)
        R[i + j] = P[i] * Q[j];
        return R;
        }

    /**
    * Return the product of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp, SizeP + SizeQ - 1>
        operator*(StaticPolynomial<Tp, SizeP> P,
            StaticPolynomial<Tp, SizeQ> Q);

    /**
    * Return type for divmod.
    */
    template<typename Tp, std::size_t SizeN, std::size_t SizeD>
        struct divmod_t
        {
        static constexpr std::size_t
        SizeQuo = (SizeD <= SizeN) ? SizeN - SizeD + 1 : 1;
        static constexpr std::size_t
        SizeRem = (SizeD > 1) ? SizeD - 1 : 1;

        StaticPolynomial<Tp, SizeQuo> quo;
        StaticPolynomial<Tp, SizeRem> rem;
        };

    /**
    * Divide two polynomials returning the quotient and remainder.
    */
    template<typename Tp, std::size_t SizeN, std::size_t SizeD>
        constexpr divmod_t<Tp, SizeN, SizeD>
        divmod(StaticPolynomial<Tp, SizeN> num,
        StaticPolynomial<Tp, SizeD> den);

    /**
    * Return the quotient of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp,
                divmod_t<Tp, SizeP, SizeQ>::SizeQuo>
        operator/(StaticPolynomial<Tp, SizeP> P,
            StaticPolynomial<Tp, SizeQ> Q)
        { return divmod(P, Q).quo; }

    /**
    * Return the remainder of two polynomials.
    */
    template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
        inline constexpr StaticPolynomial<Tp,
                divmod_t<Tp, SizeP, SizeQ>::SizeRem>
        operator%(StaticPolynomial<Tp, SizeP> P,
            StaticPolynomial<Tp, SizeQ> Q)
        { return divmod(P, Q).rem; }


    /**
    * Write a polynomial to a stream.
    * The format is a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp, std::size_t Size>
        std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& os,
            const StaticPolynomial<Tp, Size>& poly)
        {
        int old_prec = os.precision(std::numeric_limits<Tp>::max_digits10);
        os << "(";
        for (size_t i = 0; i < poly.degree(); ++i)
            os << poly.coefficient(i) << ",";
        os << poly.coefficient(poly.degree());
        os << ")";
        os.precision(old_prec);
        return os;
        }

    /**
    * Divide two polynomials returning the quotient and remainder.
    */
    template<typename Tp, std::size_t SizeN, std::size_t SizeD>
        constexpr divmod_t<Tp, SizeN, SizeD>
        divmod(StaticPolynomial<Tp, SizeN> num,
        StaticPolynomial<Tp, SizeD> den)
        {
        constexpr auto DegN = num.degree();
        constexpr auto DegD = den.degree();
        auto rem = num;
        auto quo = StaticPolynomial<Tp, SizeN>{};
        if (DegD <= DegN)
        {
        for (std::ptrdiff_t k = DegN - DegD; k >= 0; --k)
            {
            quo.coefficient(k, rem.coefficient(DegD + k)
                    / den.coefficient(DegD));
            for (int j = DegD + k - 1; j >= k; --j)
            rem.coefficient(j, rem.coefficient(j)
                        - quo.coefficient(k)
                        * den.coefficient(j - k));
            }
        }
        divmod_t<Tp, SizeN, SizeD> ret;
        for (std::size_t i = 0ULL; i < divmod_t<Tp, SizeN, SizeD>::SizeQuo; ++i)
            ret.quo[i] = quo[i];
        for (std::size_t i = 0ULL; i < divmod_t<Tp, SizeN, SizeD>::SizeRem; ++i)
            ret.rem[i] = rem[i];
        return ret;
        }
    /**
    * 
    OTOH, What does it mean to get roots of integer polynomials - the roots will
    often be real or even complex.
    What does it even mean to get roots of polynomial polynomials.
    I think you'd want to pick a parameter for the inner coefficient polynomials.
    Wait! This has nothing to do with the coefficient type.  It's the evaluation type
    that determines this. It may be sane to make real_type be double if it comes in
    as integral.
    OTOOH, scaling down an integer polynomial will either screw up the polynomial
    or it will make a polynomial of a different type.
    This scaling thing can only apply to real or complex polynomials.
    template<typename Tp>
        std::enable_if_t<std::is_integral_v<Tp>>
        Polynomial<Tp>::m_set_scale()
        { }
    */

    /**
    * 
    template<typename Tp>
        //std::enable_if_t<std::is_floating_point_v<Tp>>
        void
        Polynomial<Tp>::m_set_scale()
        {
        constexpr real_type s_eps = std::numeric_limits<real_type>::epsilon();
        constexpr real_type
        s_base = real_type{std::numeric_limits<real_type>::radix};
        constexpr real_type s_tiny = s_eps * s_eps * s_eps; // ~ 1.0e-50;
        constexpr real_type s_huge = std::numeric_limits<real_type>::max();
        constexpr real_type s_low = s_tiny / s_eps;
        // Find largest and smallest moduli of coefficients.
        auto a_max = real_type{0};
        auto a_min = s_huge;
        for (int i = 0; i <= this->degree(); ++i)
        {
        const auto x = get_scale(this->m_coeff[i]);
        if (x > a_max)
            a_max = x;
        if (x != real_type{0} && x < a_min)
            a_min = x;
        }
        // Scale if there are large or very tiny coefficients.
        // Computes a scale factor to multiply the coefficients
        // of the polynomial. The scaling is done to avoid overflow
        // and to avoid undetected underflow interfering
        // with the convergence criterion.
        // The factor is a power of the base.
        auto scale = s_low / a_min;
        bool rescale = true;
        if (scale > real_type{1} && s_huge / scale < a_max)
        rescale = false;
        if (scale <= real_type{1} && a_max < real_type{10})
        rescale = false;
        this->m_scale = real_type{1};
        if (rescale)
        {
        // Scale polynomial.
        if (scale == real_type{0})
            scale = s_tiny;
        const auto lg = std::ilogb(scale);
        this->m_scale = std::pow(s_base, lg);
        //if (this->m_scale != real_type{1})
        //  for (int i = 0; i <= this->degree(); ++i)
        //    this->m_coeff[i] *= this->m_scale;
        }
        }
    */

    /**
    * Scaling specialization for polynomial coefficient polynomials.
    template<typename Tp>
        void
        Polynomial<Polynomial<Tp>>::m_set_scale()
        {
        for (int i = 0; i <= this->degree(); ++i)
        this->m_coeff[i].m_set_scale();
        }
    */

    /**
    * Evaluate the polynomial using a modification of Horner's rule which
    * exploits the fact that the polynomial coefficients are all real.
    *
    * The algorithm is discussed in detail in:
    * Knuth, D. E., The Art of Computer Programming: Seminumerical
    * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
    *
    * If n is the degree of the polynomial,
    * n - 3 multiplies and 4 * n - 6 additions are saved.
    */
    template<typename Tp>
        template<typename Up>
        auto
        Polynomial<Tp>::operator()(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
        {
        const auto r = Tp{2} * std::real(z);
        const auto s = std::norm(z);
        size_type n = this->degree();
        auto aa = this->coefficient(n);
        if (n > 0)
        {
            auto bb = this->coefficient(n - 1);
            for (size_type j = 2; j <= n; ++j)
            bb = this->coefficient(n - j)
            - s * std::exchange(aa, bb + r * aa);
            return aa * z + bb;
        }
        else
        return aa;
        };

    //  Could/should this be done by output iterator range?
    template<typename Tp>
        template<typename Polynomial<Tp>::size_type N>
        void
        Polynomial<Tp>::eval(typename Polynomial<Tp>::value_type x,
                    std::array<Polynomial<Tp>::value_type, N>& arr)
        {
        if (arr.size() > 0)
        {
            arr.fill(value_type{});
            const size_type sz = m_coeff.size();
            arr[0] = this->coefficient(sz - 1);
                for (int i = sz - 2; i >= 0; --i)
            {
            int nn = std::min(arr.size() - 1, sz - 1 - i);
            for (int j = nn; j >= 1; --j)
            arr[j] = std::fma(arr[j], x, arr[j - 1]);
            arr[0] = std::fma(arr[0], x, this->coefficient(i));
            }
            //  Now put in the factorials.
            value_type fact = value_type(1);
            for (size_type n = arr.size(), i = 2; i < n; ++i)
            {
            fact *= value_type(i);
            arr[i] *= fact;
            }
        }
        }

    /**
    * Evaluate the polynomial and its derivatives at the point x.
    * The values are placed in the output range starting with the
    * polynomial value and continuing through higher derivatives.
    */
    template<typename Tp>
        template<typename OutIter>
        void
        Polynomial<Tp>::eval(typename Polynomial<Tp>::value_type x,
                    OutIter b, OutIter e)
        {
        if(b != e)
        {
            std::fill(b, e, value_type{});
            const size_type sz = m_coeff.size();
            *b = m_coeff[sz - 1];
                for (int i = sz - 2; i >= 0; --i)
            {
            for (auto it = std::reverse_iterator<OutIter>(e);
                it != std::reverse_iterator<OutIter>(b) - 1; ++it)
            *it = std::fma(*it, x, *(it + 1));
            *b = std::fma(*b, x, m_coeff[i]);
            }
            //  Now put in the factorials.
            int i = 0;
            value_type fact = value_type(++i);
            for (auto it = b + 1; it != e; ++it)
            {
            fact *= value_type(i);
            *it *= fact;
            ++i;
            }
        }
        }

    /**
    * Evaluate the even part of the polynomial at the input point.
    */
    template<typename Tp>
        typename Polynomial<Tp>::value_type
        Polynomial<Tp>::eval_even(typename Polynomial<Tp>::value_type x) const
        {
        if (this->degree() > 0)
        {
        const auto odd = this->degree() % 2;
        const auto xx = x * x;
        auto poly(this->coefficient(this->degree() - odd));
        for (int i = this->degree() - odd - 2; i >= 0; i -= 2)
            poly = std::fma(xx, poly, this->coefficient(i));
        return poly;
        }
        else
        return value_type{};
        }

    /**
    * Evaluate the odd part of the polynomial at the input point.
    */
    template<typename Tp>
        typename Polynomial<Tp>::value_type
        Polynomial<Tp>::eval_odd(typename Polynomial<Tp>::value_type x) const
        {
        if (this->degree() > 0)
        {
        const auto even = (this->degree() % 2 == 0 ? 1 : 0);
        const auto xx = x * x;
        auto poly(this->coefficient(this->degree() - even));
        for (int i = this->degree() - even - 2; i >= 0; i -= 2)
            poly = std::fma(xx, poly, this->coefficient(i));
        return x * poly;
        }
        else
        return value_type{};
        }

    /**
    * Evaluate the even part of the polynomial using a modification
    * of Horner's rule which exploits the fact that the polynomial
    * coefficients are all real.
    *
    * The algorithm is discussed in detail in:
    * Knuth, D. E., The Art of Computer Programming: Seminumerical
    * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
    *
    * If n is the degree of the polynomial,
    * n - 3 multiplies and 4 * n - 6 additions are saved.
    */
    template<typename Tp>
        template<typename Up>
        auto
        Polynomial<Tp>::eval_even(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
        {
        using real_t = std::decay_t<decltype(value_type{} * Up{})>;
        using cmplx_t = std::complex<real_t>;
        if (this->degree() > 0)
        {
            const auto zz = z * z;
            const auto r = Tp{2} * std::real(zz);
            const auto s = std::norm(zz);
            auto odd = this->degree() % 2;
            size_type n = this->degree() - odd;
            auto aa = this->coefficient(n);
            auto bb = this->coefficient(n - 2);
            for (size_type j = 4; j <= n; j += 2)
            bb = std::fma(-s, std::exchange(aa, bb + r * aa),
                    this->coefficient(n - j));
            return emsr::fma(cmplx_t(aa), cmplx_t(zz), cmplx_t(bb));
        }
        else
        return cmplx_t{};
        };

    /**
    * Evaluate the odd part of the polynomial using a modification
    * of Horner's rule which exploits the fact that the polynomial
    * coefficients are all real.
    *
    * The algorithm is discussed in detail in:
    * Knuth, D. E., The Art of Computer Programming: Seminumerical
    * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
    *
    * If n is the degree of the polynomial,
    * n - 3 multiplies and 4 * n - 6 additions are saved.
    */
    template<typename Tp>
        template<typename Up>
        auto
        Polynomial<Tp>::eval_odd(const std::complex<Up>& z) const
        -> std::enable_if_t<!has_imag_v<Tp>,
                std::complex<std::decay_t<
            decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
        {
        using real_t = std::decay_t<decltype(value_type{} * Up{})>;
        using cmplx_t = std::complex<real_t>;
        if (this->degree() > 0)
        {
            const auto zz = z * z;
            const auto r = Tp{2} * std::real(zz);
            const auto s = std::norm(zz);
            auto even = (this->degree() % 2 == 0 ? 1 : 0);
            size_type n = this->degree() - even;
            auto aa = this->coefficient(n);
            auto bb = this->coefficient(n - 2);
            for (size_type j = 4; j <= n; j += 2)
            bb = std::fma(-s, std::exchange(aa, bb + r * aa),
                    this->coefficient(n - j));
            return z
            * emsr::fma(cmplx_t(aa), cmplx_t(zz), cmplx_t(bb));
        }
        else
        return cmplx_t{};
        };

        /**
        * Multiply the polynomial by another polynomial.
        */
    template<typename Tp>
        template<typename Up>
        Polynomial<Tp>&
        Polynomial<Tp>::operator*=(const Polynomial<Up>& poly)
        {
        //  Test for zero size polys and do special processing?
        const size_type m = this->degree();
        const size_type n = poly.degree();
        std::vector<value_type> new_coeff(m + n + 1);
        for (size_type i = 0; i <= m; ++i)
        for (size_type j = 0; j <= n; ++j)
            new_coeff[i + j] += this->m_coeff[i]
                    * static_cast<value_type>(poly.m_coeff[j]);
        this->m_coeff = new_coeff;
        return *this;
        }

    /**
    * Divide two polynomials returning the quotient and remainder.
    */
    template<typename Tp>
        void
        divmod(const Polynomial<Tp>& num, const Polynomial<Tp>& den,
            Polynomial<Tp>& quo, Polynomial<Tp>& rem)
        {
        rem = num;
        quo = Polynomial<Tp>(Tp{}, num.degree());
        const std::size_t d_num = num.degree();
        const std::size_t d_den = den.degree();
        if (d_den <= d_num)
        {
        for (int k = d_num - d_den; k >= 0; --k)
            {
            quo.coefficient(k, rem.coefficient(d_den + k)
                    / den.coefficient(d_den));
            for (int j = d_den + k - 1; j >= k; --j)
            rem.coefficient(j, rem.coefficient(j)
                        - quo.coefficient(k)
                        * den.coefficient(j - k));
            }
        quo.degree(d_num - d_den);
        rem.degree(d_den > 0 ? d_den - 1 : 0);
        }
        else
        quo.degree(0);
        }

    /**
    * Write a polynomial to a stream.
    * The format is a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& os,
            const Polynomial<Tp>& poly)
        {
        int old_prec = os.precision(std::numeric_limits<Tp>::max_digits10);
        os << "(";
        for (size_t i = 0; i < poly.degree(); ++i)
            os << poly.coefficient(i) << ",";
        os << poly.coefficient(poly.degree());
        os << ")";
        os.precision(old_prec);
        return os;
        }

    /**
    * Read a polynomial from a stream.
    * The input format can be a plain scalar (zero degree polynomial)
    * or a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& is,
            Polynomial<Tp>& poly)
        {
        Tp x;
        CharT ch;
        is >> ch;
        if (ch == '(')
        {
        do
            {
            is >> x >> ch;
            poly.m_coeff.push_back(x);
            }
        while (ch == ',');
        if (ch != ')')
            is.setstate(std::ios_base::failbit);
        }
        else
        {
        is.putback(ch);
        is >> x;
        poly = x;
        }
        // No null polynomial.
        if (poly.size() == 0)
        poly.m_coeff.resize(1, Tp{});

        return is;
        }


    /**
    * Perform compile-time evaluation of a constant zero-order polynomial.
    */
    template<typename ArgT, typename Coef0>
    constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
    horner(ArgT, Coef0 c0)
    {
        using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
                        double, ArgT>;
        return arg_t(c0);
    }

    /**
    * Perform compile-time evaluation of a constant polynomial.
    * The polynomial coefficients are lowest-order first.
    */
    template<typename ArgT, typename Coef0, typename... Coef>
    constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
    horner(ArgT x, Coef0 c0, Coef... c)
    {
        using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
                        double, ArgT>;
        return arg_t(c0) + x * horner(x, c...);
    }


    /**
    * Perform compile-time evaluation of a constant zero-order polynomial.
    * The polynomial coefficients are highest-order first.
    */
    template<typename ArgT, typename Coef0>
    constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
    horner_big_end(ArgT, Coef0 c0)
    {
        using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
                        double, ArgT>;
        return arg_t(c0);
    }

    /**
    * Perform compile-time evaluation of a constant first-order polynomial.
    * The polynomial coefficients are highest-order first.
    */
    template<typename ArgT, typename Coef1, typename Coef0>
    constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
    horner_big_end(ArgT x, Coef1 c1, Coef0 c0)
    {
        using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
                        double, ArgT>;
        return horner_big_end(x, x * arg_t(c1) + arg_t(c0));
    }

    /**
    * Perform compile-time evaluation of a constant polynomial.
    * The polynomial coefficients are highest-order first.
    */
    template<typename ArgT, typename CoefN, typename CoefNm1, typename... Coef>
    constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
    horner_big_end(ArgT x, CoefN cn, CoefNm1 cnm1, Coef... c)
    {
        using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
                        double, ArgT>;
        return horner_big_end(x, x * arg_t(cn) + arg_t(cnm1), c...);
    }

        /**
    *
    */
    template<typename Tp>
    class RationalPolynomial
    {
    public:
        /**
        * Typedefs.
        */
        using polynomial_type = Polynomial<Tp>;
        using value_type = typename polynomial_type::value_type;
    // The above should be Polynomial<Tp>to be consistent with rational.h
    /* These might not make sense.
        using reference = typename polynomial_type::reference;
        using const_reference = typename polynomial_type::const_reference;
        using pointer = typename polynomial_type::pointer;
        using const_pointer = typename polynomial_type::const_pointer;
        using iterator = typename polynomial_type::iterator;
        using const_iterator = typename polynomial_type::const_iterator;
        using reverse_iterator = typename polynomial_type::reverse_iterator;
        using const_reverse_iterator = typename polynomial_type::const_reverse_iterator;
    */
        using size_type = typename polynomial_type::size_type;
        using difference_type = typename polynomial_type::difference_type;

        /**
        * Create a zero degree polynomial with value zero.
        */
        RationalPolynomial()
        : m_num(), m_den()
        { }

        /**
        * Copy ctor.
        */
        RationalPolynomial(const RationalPolynomial&) = default;

        /**
        * Move ctor.
        */
        RationalPolynomial(RationalPolynomial&&) = default;

        RationalPolynomial(const Polynomial<Tp>& num,
                    const Polynomial<Tp>& den)
        : m_num(num), m_den(den)
        { }

        explicit RationalPolynomial(const Polynomial<Tp>& num)
        : m_num(num), m_den(Polynomial<Tp>(Tp{1}))
        { }

        /**
        * Evaluate the polynomial at the input point.
        */
        value_type
        operator()(value_type x) const
        { return this->m_num(x) / this->m_den(x); }

        /**
        * Unary plus.
        */
        RationalPolynomial
        operator+() const
        { return *this; }

        /**
        * Unary minus.
        */
        RationalPolynomial
        operator-() const
        { return RationalPolynomial(*this) *= value_type(-1); }

        /**
        * Copy assignment.
        */
        RationalPolynomial&
        operator=(const RationalPolynomial&) = default;

        /**
        * Move assignment.
        */
        RationalPolynomial&
        operator=(RationalPolynomial&&) = default;

        /**
        * Add a rational polynomial to this rational polynomial.
        */
        RationalPolynomial&
        operator+=(const RationalPolynomial& x)
        {
            this->numer() = this->numer() * x.denom() + this->denom() * x.numer();
            this->denom() *= x.denom();
        return *this;
        }

        /**
        * Subtract a rational polynomial from this rational polynomial.
        */
        RationalPolynomial&
        operator-=(const RationalPolynomial& x)
        {
            this->numer() = this->numer() * x.denom() - this->denom() * x.numer();
            this->denom() *= x.denom();
        return *this;
        }

        /**
        * Multiply this rational polynomial by a rational polynomial.
        */
        RationalPolynomial&
        operator*=(const RationalPolynomial& x)
        {
        this->numer() *= x.numer();
        this->denom() *= x.denom();
        return *this;
        }

        /**
        * Divide this rational polynomial by a rational polynomial.
        */
        RationalPolynomial&
        operator/=(const RationalPolynomial& x)
        {
        this->numer() *= x.denom();
        this->denom() *= x.numer();
        return *this;
        }

        const Polynomial<value_type>&
        numer() const
        { return this->m_num; }

        Polynomial<value_type>&
        numer()
        { return this->m_num; }

        value_type
        numer(value_type x) const
        { return this->m_num(x); }

        const Polynomial<value_type>&
        denom() const
        { return this->m_den; }

        Polynomial<value_type>&
        denom()
        { return this->m_den; }

        value_type
        denom(value_type x) const
        { return this->m_den(x); }

        private:

        Polynomial<value_type> m_num;
        Polynomial<value_type> m_den;
        };

    /**
    * Write a polynomial to a stream.
    * The format is a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& os, const RationalPolynomial<Tp>& poly)
        {
        os << poly.numer() << "/" << poly.denom();
        return os;
        }

    /**
    * Read a polynomial from a stream.
    * The input format can be a plain scalar (zero degree polynomial)
    * or a parenthesized comma-delimited list of coefficients.
    */
    template<typename CharT, typename Traits, typename Tp>
        std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& is, RationalPolynomial<Tp>& poly)
        {
        Polynomial<Tp> numer, denom;
        is >> numer;
        if (!is.fail())
        {
        CharT ch;
        is >> ch;
        if (ch != '/')
            is.setstate(std::ios_base::failbit);
        else
            {
            is >> denom;
            if (!is.fail())
            {
            poly.numer() = numer;
            poly.denom() = denom;
            }
            }
        }
        return is;
        }

    template<typename Real>
        struct Solution
        : public std::variant<std::monostate, Real, std::complex<Real>>
        {
        using Base = std::variant<std::monostate, Real, std::complex<Real>>;

        Solution() noexcept
        : Base()
        {}

        Solution(Real x) noexcept
        : Base(x)
        {}

        Solution(const std::complex<Real>& x) noexcept
        : Base(x)
        {}
        };

    template<typename Real>
        constexpr bool
        is_valid(const Solution<Real>& x)
        { return x.index() != 0; }

    template<typename Real>
        constexpr Real
        real(const Solution<Real>& x)
        {
        if (x.index() == 0)
        return std::numeric_limits<Real>::quiet_NaN();
        else if (x.index() == 1)
        return std::get<1>(x);
        else
        return std::real(std::get<2>(x));
        }

    template<typename Real>
        constexpr Real
        imag(const Solution<Real>& x)
        {
        if (x.index() == 0)
        return std::numeric_limits<Real>::quiet_NaN();
        else if (x.index() == 1)
        return Real{0};
        else
        return std::imag(std::get<2>(x));
        }

    template<typename Real>
        constexpr Real
        abs(const Solution<Real>& x)
        {
        if (x.index() == 0)
        return std::numeric_limits<Real>::quiet_NaN();
        else if (x.index() == 1)
        return std::abs(std::get<1>(x));
        else
        return std::abs(std::get<2>(x));
        }

    /**
    * Return the solution as a complex number or NaN.
    */
    template<typename Real>
        constexpr Solution<Real>
        to_complex(const Solution<Real>& x)
        {
        if (x.index() == 0)
        return Solution<Real>();
        else if (x.index() == 1)
        return Solution<Real>(std::complex<Real>(std::get<1>(x)));
        else
        return x;
        }

    /**
    * Unary +/-
    */
    template<typename Real>
        constexpr Solution<Real>
        operator+(const Solution<Real>& x)
        { return x; }

    template<typename Real>
        constexpr Solution<Real>
        operator-(const Solution<Real>& x)
        {
        if (x.index() == 0)
        return x;
        else if (x.index() == 1)
        return Solution<Real>(-std::get<1>(x));
        else
        return Solution<Real>(-std::get<2>(x));
        }

    /**
    * Addition operators...
    */
    template<typename Real>
        constexpr Solution<Real>
        operator+(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0)
        return x;
        if (y.index() == 0)
        return y;
        else if (x.index() == 1)
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<1>(x) + std::get<1>(y));
        else
            return Solution<Real>(std::get<1>(x) + std::get<2>(y));
        }
        else
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<2>(x) + std::get<1>(y));
        else
            return Solution<Real>(std::get<2>(x) + std::get<2>(y));
        }
        }

    template<typename Real>
        constexpr Solution<Real>
        operator+(const Solution<Real>& x, Real y)
        { return operator+(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator+(Real x, const Solution<Real>& y)
        { return operator+(Solution<Real>(x), y); }

    template<typename Real>
        constexpr Solution<Real>
        operator+(const Solution<Real>& x, const std::complex<Real>& y)
        { return operator+(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator+(const std::complex<Real>& x, const Solution<Real>& y)
        { return operator+(Solution<Real>(x), y); }

    /**
    * Subtraction operators...
    */
    template<typename Real>
        constexpr Solution<Real>
        operator-(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0)
        return x;
        if (y.index() == 0)
        return y;
        else if (x.index() == 1)
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<1>(x) - std::get<1>(y));
        else
            return Solution<Real>(std::get<1>(x) - std::get<2>(y));
        }
        else
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<2>(x) - std::get<1>(y));
        else
            return Solution<Real>(std::get<2>(x) - std::get<2>(y));
        }
        }

    template<typename Real>
        constexpr Solution<Real>
        operator-(const Solution<Real>& x, Real y)
        { return operator-(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator-(Real x, const Solution<Real>& y)
        { return operator-(Solution<Real>(x), y); }

    template<typename Real>
        constexpr Solution<Real>
        operator-(const Solution<Real>& x, const std::complex<Real>& y)
        { return operator-(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator-(const std::complex<Real>& x, const Solution<Real>& y)
        { return operator-(Solution<Real>(x), y); }

    /**
    * Multiplication operators...
    */
    template<typename Real>
        constexpr Solution<Real>
        operator*(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0)
        return x;
        if (y.index() == 0)
        return y;
        else if (x.index() == 1)
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<1>(x) * std::get<1>(y));
        else
            return Solution<Real>(std::get<1>(x) * std::get<2>(y));
        }
        else
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<2>(x) * std::get<1>(y));
        else
            return Solution<Real>(std::get<2>(x) * std::get<2>(y));
        }
        }

    template<typename Real>
        constexpr Solution<Real>
        operator*(const Solution<Real>& x, Real y)
        { return operator*(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator*(Real x, const Solution<Real>& y)
        { return operator*(Solution<Real>(x), y); }

    template<typename Real>
        constexpr Solution<Real>
        operator*(const Solution<Real>& x, const std::complex<Real>& y)
        { return operator*(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator*(const std::complex<Real>& x, const Solution<Real>& y)
        { return operator*(Solution<Real>(x), y); }

    /**
    * division operators...
    */
    template<typename Real>
        constexpr Solution<Real>
        operator/(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0)
        return x;
        if (y.index() == 0)
        return y;
        else if (x.index() == 1)
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<1>(x) / std::get<1>(y));
        else
            return Solution<Real>(std::get<1>(x) / std::get<2>(y));
        }
        else
        {
        if (y.index() == 1)
            return Solution<Real>(std::get<2>(x) / std::get<1>(y));
        else
            return Solution<Real>(std::get<2>(x) / std::get<2>(y));
        }
        }

    template<typename Real>
        constexpr Solution<Real>
        operator/(const Solution<Real>& x, Real y)
        { return operator/(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator/(Real x, const Solution<Real>& y)
        { return operator/(Solution<Real>(x), y); }

    template<typename Real>
        constexpr Solution<Real>
        operator/(const Solution<Real>& x, const std::complex<Real>& y)
        { return operator/(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr Solution<Real>
        operator/(const std::complex<Real>& x, const Solution<Real>& y)
        { return operator/(Solution<Real>(x), y); }

    /**
    * Test for equality and inequality.
    */
    template<typename Real>
        constexpr bool
        operator==(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0 || y.index() == 0)
        return false;
        else if (x.index() == y.index())
        {
        if (x.index() == 1)
            return std::get<1>(x) == std::get<1>(y);
        else
            return std::get<2>(x) == std::get<2>(y);
        }
        else
        return false;
        }

    template<typename Real>
        bool
        operator==(const Solution<Real>& x, Real y)
        { return x == Solution<Real>(y); }

    template<typename Real>
        constexpr bool
        operator==(Real x, const Solution<Real>& y)
        { return Solution<Real>(x) == y; }

    template<typename Real>
        constexpr bool
        operator==(const Solution<Real>& x, const std::complex<Real>& y)
        { return x == Solution<Real>(y); }

    template<typename Real>
        constexpr bool
        operator==(const std::complex<Real>& x, const Solution<Real>& y)
        { return Solution<Real>(x) == y; }

    template<typename Real>
        constexpr bool
        operator!=(const Solution<Real>& x, const Solution<Real>& y)
        { return !(x == y); }

    template<typename Real>
        bool
        operator!=(const Solution<Real>& x, Real y)
        { return !(x == y); }

    template<typename Real>
        constexpr bool
        operator!=(Real x, const Solution<Real>& y)
        { return !(x == y); }

    template<typename Real>
        constexpr bool
        operator!=(const Solution<Real>& x, const std::complex<Real>& y)
        { return !(x == y); }

    template<typename Real>
        constexpr bool
        operator!=(const std::complex<Real>& x, const Solution<Real>& y)
        { return !(x == y); }

    /**
    * Lexicographic order of solutions as complex numbers.
    * Null solutions compare as less than except to another null solution.
    *
    * A tribool might be a good thing for this when either
    * of the solutions is null.
    */
    template<typename Real>
        constexpr bool
        operator<(const Solution<Real>& x, const Solution<Real>& y)
        {
        if (x.index() == 0 && y.index() == 0)
        return false;
        else if (x.index() == 0)
        return true;
        else if (y.index() == 0)
        return false;
        else
        {
        const auto rex = emsr::real(x);
        const auto rey = emsr::real(y);
        if (rex < rey)
            return true;
        else if (rex == rey)
            return emsr::imag(x) < emsr::imag(y);
        else
            return false;
        }
        }

    template<typename Real>
        constexpr bool
        operator<(const Solution<Real>& x, Real y)
        { return operator<(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr bool
        operator<(Real x, const Solution<Real>& y)
        { return operator<(Solution<Real>(x), y); }

    template<typename Real>
        constexpr bool
        operator<(const Solution<Real>& x, const std::complex<Real>& y)
        { return operator<(x, Solution<Real>(y)); }

    template<typename Real>
        constexpr bool
        operator<(const std::complex<Real>& x, const Solution<Real>& y)
        { return operator<(Solution<Real>(x), y); }

    /**
    * Output a solution to a stream.
    */
    template<typename Char, typename Real>
        std::basic_ostream<Char>&
        operator<<(std::basic_ostream<Char>& out,
            const emsr::Solution<Real>& sln)
        {
        const auto idx = sln.index();
        if (idx == 0)
        out << "null";
        else if (idx == 1)
        out << std::get<1>(sln);
        else if (idx == 2)
        out << std::get<2>(sln);
        return out;
        }

    template<std::size_t Dim, typename Iter, typename NumTp>
        NumTp
        refine_solution_newton(NumTp z, const Iter& CC);

    template<std::size_t Dim, typename Iter, typename NumTp>
        NumTp
        refine_solution_halley(NumTp z, const Iter& CC);

    template<typename Real, typename Iter>
        std::array<Solution<Real>, 2>
        quadratic(const Iter& coef);

    template<typename Real>
        inline std::array<Solution<Real>, 2>
        quadratic(Real c0, Real c1, Real c2)
        {
        return quadratic<Real>(std::array<Real, 3>{c0, c1, c2});
        }

    template<typename Real, typename Iter>
        std::array<Solution<Real>, 3>
        cubic(const Iter& coef);

    template<typename Real>
        inline std::array<Solution<Real>, 3>
        cubic(Real c0, Real c1, Real c2, Real c3)
        {
        return cubic<Real>(std::array<Real, 4>{c0, c1, c2, c3});
        }

    template<typename Real, typename Iter>
        std::array<Solution<Real>, 4>
        quartic(const Iter& coef);

    template<typename Real>
        inline std::array<Solution<Real>, 4>
        quartic(Real c0, Real c1, Real c2, Real c3, Real c4)
        {
        return quartic<Real>(std::array<Real, 5>{c0, c1, c2, c3, c4});
        }

        /**
    * Refine a solution using the Newton method:
    * @f[
    *    x_{n+1} - x_n = -\frac{P(x_n)}{P'(x_n)}
    * @f]
    */
    template<std::size_t _Dim, typename _Iter, typename _NumTp>
        _NumTp
        refine_solution_newton(_NumTp z, const _Iter& _CC)
        {
        for (int i = 0; i < 3; ++i)
        {
        auto f = _NumTp(_CC[_Dim - 1]);
        for (std::size_t j = _Dim - 1; j > 0; --j)
            f = _NumTp(_CC[j - 1]) + z * f;

        auto df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
        for (std::size_t j = _Dim - 1; j > 1; --j)
            df = _NumTp((j - 1) * _CC[j - 1]) + z * df;

        const auto del = f / df;
        z -= del;
        }
        return z;
        }

    /**
    * Refine a solution using the Halley method:
    * @f[
    *    x_{n+1} - x_n = -\frac{2 P(x_n) P'(x_n)}
    *                    {2 [P'(x_n)]^2 - P(x_n) P''(x_n)}
    * @f]
    * This form indicates the close relationship to the Newton method:
    * @f[
    *    x_{n+1} - x_n = -\frac{P'(x_n)}
    *                    {P'(x_n) - [P(x_n) P''(x_n)]/[2P'(x_n)]}
    * @f]
    */
    template<std::size_t _Dim, typename _Iter, typename _NumTp>
        _NumTp
        refine_solution_halley(_NumTp z, const _Iter& _CC)
        {
        for (int i = 0; i < 3; ++i)
        {
        auto f = _NumTp(_CC[_Dim - 1]);
        for (std::size_t j = _Dim - 1; j > 0; --j)
            f = _NumTp(_CC[j - 1]) + z * f;

        auto df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
        for (std::size_t j = _Dim - 1; j > 1; --j)
            df = _NumTp((j - 1) * _CC[j - 1]) + z * df;

        auto d2f = _NumTp((_Dim - 2) * (_Dim - 1) * _CC[_Dim - 1]);
        for (std::size_t j = _Dim - 1; j > 2; --j)
            d2f = _NumTp((j - 2) * (j - 1) * _CC[j - 1]) + z * d2f;

        const auto del = _NumTp{2} * f * df
                / (_NumTp{2} * df * df - f * d2f);

        z -= del;
        }
        return z;
        }

    template<std::size_t _Dim, typename _Iter, typename Real>
        void
        refine_solutions(std::array<Solution<Real>, _Dim - 1>& _ZZ, const _Iter& _CC)
        {
        for (std::size_t i = 0; i < _Dim - 1; ++i)
        {
        if (_ZZ[i].index() == 0)
            continue;
        else if (_ZZ[i].index() == 1)
            _ZZ[i] = refine_solution_newton<_Dim>(std::get<1>(_ZZ[i]), _CC);
        else if (_ZZ[i].index() == 2)
            _ZZ[i] = refine_solution_newton<_Dim>(std::get<2>(_ZZ[i]), _CC);
        }
        }

    /**
    * @brief Finds the roots of a quadratic equation of the form:
    * @f[
    *    a_2 x^2 + a_1 x + a_0 = 0
    * @f]
    * for real coefficients @f$ a_k @f$.
    *
    * For non-degenerate coefficients two roots are returned:
    * Either the roots are real or the roots are a complex conjugate pair.
    *
    * If the quadratic coefficient @f$ a_2 @f$ is zero (degenerate case)
    * at most one valid root is returned.
    * If the linear coefficient @f$ a_1 @f$ is also zero
    * no valid root is returned.
    *
    * @param[in] _CC Array that contains the three coefficients
    *                  of the quadratic equation.
    */
    template<typename Real, typename _Iter>
        std::array<Solution<Real>, 2>
        quadratic(const _Iter& _CC)
        {
        std::array<Solution<Real>, 2> _ZZ;

        if (_CC[2] == Real{0})
        {
        // Equation is linear (or completely degenerate).
        if (_CC[1] == Real{0})
            return _ZZ;
        else
            {
            _ZZ[0] = -_CC[0] / _CC[1];
            return _ZZ;
            }
        }
        else if (_CC[0] == Real{0})
        {
        _ZZ[0] = Real{0};
        if (_CC[2] == Real{0})
            return _ZZ;
        else
            {
            _ZZ[1] = -_CC[1] / _CC[2];
            return _ZZ;
            }
        }
        else
        {
        // The discriminant of a quadratic equation
        const auto _QQ = _CC[1] * _CC[1] - Real{4} * _CC[2] * _CC[0];

        if (_QQ < Real{0})
            {
            // The roots are complex conjugates.
            const auto _ReZZ = -_CC[1] / (Real{2} * _CC[2]);
            const auto _ImZZ = std::sqrt(std::abs(_QQ)) / (Real{2} * _CC[2]);
            _ZZ[0] = std::complex<Real>(_ReZZ, -_ImZZ);
            _ZZ[1] = std::complex<Real>(_ReZZ, _ImZZ);
            }
        else
            {
            // The roots are real.
            Real temp = -(_CC[1]
                + std::copysign(std::sqrt(_QQ), _CC[1])) / Real{2};
            _ZZ[0] = temp / _CC[2];
            _ZZ[1] = _CC[0] / temp;
            }
        }

        return _ZZ;
        }


    /**
    * @brief Finds the roots of a cubic equation of the form:
    * @f[
    *    a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
    * @f]
    * for real coefficients @f$ a_k @f$.
    *
    * In the non-degenerate case there are three roots:
    * - All three roots are real
    * - One root is real and the other two are a complex conjugate pair
    *
    * If the cubic coefficient @f$ a_3 @f$ is zero (degenerate case)
    * the problem is referred to the quadratic solver to return, at most,
    * two valid roots.
    *
    * @param[in] _CC Array that contains the four coefficients
    *                  of the cubic equation
    */
    template<typename Real, typename Iter>
        std::array<Solution<Real>, 3>
        cubic(const Iter& CC)
        {
        std::array<Solution<Real>, 3> ZZ;

        if (CC[3] == Real{0})
        {
        // Last root is null, remaining equation is quadratic.
        const auto ZZ2 = quadratic<Real>(CC);
        ZZ[0] = ZZ2[0];
        ZZ[1] = ZZ2[1];
        }
        else if (CC[0] == Real{0})
        {
        // First root is zero, remaining equation is quadratic.
        ZZ[0] = Real{0};
        const auto ZZ2 = quadratic<Real>(CC[1], CC[2], CC[3]);
        ZZ[1] = ZZ2[0];
        ZZ[2] = ZZ2[1];
        }
        else
        {
        // Normalize cubic equation coefficients.
        std::array<Real, 4> AA3;
        AA3[3] = Real{1};
        AA3[2] = CC[2] / CC[3];
        AA3[1] = CC[1] / CC[3];
        AA3[0] = CC[0] / CC[3];

        const auto S_2pi
            = 2 * Real{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
        const auto PP = AA3[2] / Real{3};
        const auto QQ = (AA3[2] * AA3[2] - Real{3} * AA3[1])
                / Real{9};
        const auto QQp3 = QQ * QQ * QQ;
        const auto RR = (Real{2} * AA3[2] * AA3[2] * AA3[2]
                - Real{9} * AA3[2] * AA3[1]
                + Real{27} * AA3[0]) / Real{54};
        const auto RRp2 = RR * RR;

        if (QQp3 - RRp2 > Real{0})
            {
            // Calculate the three real roots.
            const auto phi = std::acos(RR / std::sqrt(QQp3));
            const auto fact = -Real{2} * std::sqrt(QQ);
            for (int i = 0; i < 3; ++i)
            ZZ[i] = fact * std::cos((phi + i * S_2pi) / Real{3}) - PP;
            }
        else
            {
            // Calculate the single real root.
            const auto fact = std::cbrt(std::abs(RR)
                        + std::sqrt(RRp2 - QQp3));
            const auto BB = -std::copysign(fact + QQ / fact, RR);
            ZZ[0] = BB - PP;

            // Find the other two roots which are complex conjugates.
            std::array<Real, 3> AA2;
            AA2[2] = Real{1};
            AA2[1] = BB;
            AA2[0] = BB * BB - Real{3} * QQ;
            const auto ZZ2 = quadratic<Real>(AA2);
            ZZ[1] = std::get<2>(ZZ2[0]) - PP;
            ZZ[2] = std::get<2>(ZZ2[1]) - PP;
            }
        }

        return ZZ;
        }


    /**
    * @brief Finds the roots a quartic equation of the form:
    * @f[
    * 	 a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
    * @f]
    * for real coefficients @f$ a_k @f$.
    *
    * In the non-degenerate case there are four roots:
    * - All four roots are real
    * - Two roots real and two complex roots are a complex conjugate pair
    * - Four complex roots in two complex conjugate pairs
    *
    * If the qartic coefficient @f$ a_4 @f$ is zero (degenerate case)
    * the problem is referred to the cubic solver to return, at most,
    * three valid roots.
    *
    * @param[in] CC Array that contains the five(5) coefficients
    *                  of the quartic equation.
    */
    template<typename Real, typename Iter>
        std::array<Solution<Real>, 4>
        quartic(const Iter& CC)
        {
        std::array<Solution<Real>, 4> ZZ;

        if (CC[4] == Real{0})
        {
        const auto ZZ3 = cubic<Real>(CC);
        ZZ[0] = ZZ3[0];
        ZZ[1] = ZZ3[1];
        ZZ[2] = ZZ3[2];
        }
        else if (CC[0] == Real{0})
        {
        ZZ[0] = Real{0};
        const auto ZZ3 = cubic<Real>(CC[1], CC[2], CC[3], CC[4]);
        ZZ[1] = ZZ3[0];
        ZZ[2] = ZZ3[1];
        }
        else if (CC[3] == Real{0} && CC[1] == Real{0})
        {
        // Solve the biquadratic equation.
        std::array<Real, 3> AA2{{CC[0], CC[2], CC[4]}};
        const auto ZZ2 = quadratic<Real>(AA2);
        auto sqrt = [](Solution<Real> z) -> Solution<Real>
                {
                const auto idx = z.index();
                if (idx == 0)
                    return z;
                else if (idx == 1)
                    {
                    auto zz = std::get<1>(z);
                    return zz < Real{0}
                    ? Solution<Real>(std::sqrt(std::complex<Real>(zz)))
                    : Solution<Real>(std::sqrt(zz));
                    }
                else
                    return Solution<Real>(std::sqrt(std::get<2>(z)));
                };
        ZZ[0] = sqrt(ZZ2[0]);
        ZZ[1] = sqrt(ZZ2[1]);
        ZZ[2] = -ZZ[0];
        ZZ[3] = -ZZ[1];
        }
        else
        {
        // Normalize quartic equation coefficients.
        std::array<Real, 5> AA4;
        AA4[4] = Real{1};
        AA4[3] = CC[3] / CC[4];
        AA4[2] = CC[2] / CC[4];
        AA4[1] = CC[1] / CC[4];
        AA4[0] = CC[0] / CC[4];

        // Calculate the coefficients of the resolvent cubic equation.
        std::array<Real, 4> AA3;
        AA3[3] = Real{1};
        AA3[2] = -AA4[2];
        AA3[1] = AA4[3] * AA4[1] - Real{4} * AA4[0];
        AA3[0] = AA4[0] * (Real{4} * AA4[2] - AA4[3] * AA4[3])
            - AA4[1] * AA4[1];

        // Find the algebraically largest real root of the cubic equation
        // Note: A cubic equation has either three real roots or one
        //       real root and two complex roots that are complex
        //       conjugates. If there is only a single real root then
        //       subroutine cubic always returns that single real root
        //       (and therefore the algebraically largest real root of
        //       the cubic equation) as root[0].
        Real Z3max;
        auto ZZ3 = cubic<Real>(AA3);
        if (ZZ3[1].index() == 1 && ZZ3[2].index() == 1)
                {
            // There is some horrible bug with swap and this variant.
            // They may need to hold the same type.
            if (ZZ3[0] < ZZ3[1])
            //std::swap(ZZ3[0], ZZ3[1]);
            {
            const auto tmp = ZZ3[0];
            ZZ3[0] = ZZ3[1];
            ZZ3[1] = tmp;
            }
            if (ZZ3[0] < ZZ3[2])
            //std::swap(ZZ3[0], ZZ3[2]);
            {
            const auto tmp = ZZ3[0];
            ZZ3[0] = ZZ3[2];
            ZZ3[2] = tmp;
            }
            Z3max = std::get<1>(ZZ3[0]);
                }
        else
            Z3max = std::get<1>(ZZ3[0]);

        // Calculate the coefficients for the two quadratic equations
        const auto capa = Real{0.5L} * AA4[3];
        const auto capb = Real{0.5L} * Z3max;
        const auto capc = std::sqrt(capa * capa - AA4[2] + Z3max);
        const auto capd = std::sqrt(capb * capb - AA4[0]);
        const auto cp = capa + capc;
        const auto cm = capa - capc;
        auto dp = capb + capd;
        auto dm = capb - capd;
        const auto t1 = cp * dm + cm * dp;
        const auto t2 = cp * dp + cm * dm;
        if (std::abs(t2 - AA4[1]) < std::abs(t1 - AA4[1]))
            std::swap(dp, dm);

        // Coefficients for the first quadratic equation and find the roots.
        std::array<Real, 3> AA2;
        AA2[2] = Real{1};
        AA2[1] = cp;
        AA2[0] = dp;
        const auto ZZ2p = quadratic<Real>(AA2);
        ZZ[0] = ZZ2p[0];
        ZZ[1] = ZZ2p[1];

        // Coefficients for the second quadratic equation and find the roots.
        AA2[2] = Real{1};
        AA2[1] = cm;
        AA2[0] = dm;
        const auto ZZ2m = quadratic<Real>(AA2);
        ZZ[2] = ZZ2m[0];
        ZZ[3] = ZZ2m[1];
        }

        return ZZ;
        }


    //
    /// @todo: If you *don't* reverse the input array you solve for 1/z.
    /// @todo Take a max_error.  If this->m_eps grows larger than max_error throw.
    ///
    template<typename Real>
    class BairstowSolver
    {
    public:

        BairstowSolver(const std::vector<Real>& coeff,
                unsigned int seed = std::random_device()())
        : m_coeff{coeff.rbegin(), coeff.rend()},
        m_b(coeff.size()), m_c(coeff.size()),
        m_order(coeff.size() - 1),
        m_urng(seed), m_pdf(Real{0}, Real{2})
        {
        if (this->m_coeff.size() == 0)
        throw std::domain_error("BairstowSolver: Coefficient size must be nonzero.");

        if (this->m_coeff[0] == Real{0})
        throw std::domain_error("BairstowSolver: Leading-order coefficient must be nonzero.");

        const auto scale = this->m_coeff[0];
        for (int i = 0; i <= this->m_order; ++i)
        this->m_coeff[i] /= scale;

        this->m_zero.reserve(this->m_coeff.size());
        }

        std::vector<Solution<Real>> solve();
        std::vector<Real> equations() const;

    private:

        void m_iterate();

        template<int Index>
        static void s_refine_quadratic_eqn(Real& dr, Real& r,
                        Real& ds, Real& s,
                        std::array<Solution<Real>, 2>& w);

        void
        m_add_zero(Solution<Real> z)
        {
        this->m_zero.push_back(z);
        --this->m_order;
        }

        static constexpr Real s_ratio
        = Real{std::numeric_limits<Real>::digits}
        / std::numeric_limits<double>::digits;
        static constexpr int s_max_rand_iter = 200 * s_ratio;
        static constexpr int s_max_error_iter = 500 * s_ratio;
        static constexpr auto s_eps_factor = Real{10} * s_ratio;
        static constexpr auto s_eps
        = s_eps_factor * std::numeric_limits<Real>::epsilon();

        std::vector<Real> m_coeff;
        std::vector<Real> m_b;
        std::vector<Real> m_c;
        std::vector<Solution<Real>> m_zero;
        Real m_eps = s_eps;
        int m_order;
        bool m_precision_error = false;
        std::mt19937 m_urng;
        std::uniform_real_distribution<Real> m_pdf;
    };

    /**
    * Attempt to find a quadratic factor of the polynomial by Bairstow's method.
    */
    template<typename Real>
    void
    BairstowSolver<Real>::m_iterate()
    {
        auto r = Real{0};
        auto s = Real{0};
        auto dr = Real{1};
        auto ds = Real{0};

        int iter = 1;
        this->m_precision_error = false;

        while (std::abs(dr) + std::abs(ds) > this->m_eps)
        {
        if (iter % s_max_rand_iter == 0)
        r = this->m_pdf(this->m_urng);
        if (iter % s_max_error_iter == 0)
        {
            this->m_eps *= s_eps_factor;
            this->m_precision_error = true;
            std::cout << "Loss of precision: " << this->m_eps << '\n';
        }

        this->m_b[1] = this->m_coeff[1] - r;
        this->m_c[1] = this->m_b[1] - r;
        for (int i = 2; i <= this->m_order; ++i)
        {
            this->m_b[i] = this->m_coeff[i]
                - r * this->m_b[i - 1] - s * this->m_b[i - 2];
            this->m_c[i] = this->m_b[i]
                - r * this->m_c[i - 1] - s * this->m_c[i - 2];
        }

        auto dn = this->m_c[this->m_order - 1]
            * this->m_c[this->m_order - 3]
            - this->m_c[this->m_order - 2]
            * this->m_c[this->m_order - 2];
        if (std::abs(dn) < std::numeric_limits<Real>::epsilon())
        {
            dr = Real{1};
            ds = Real{1};
        }
        else
        {
            auto drn = this->m_b[this->m_order]
                * this->m_c[this->m_order - 3]
                - this->m_b[this->m_order - 1]
                * this->m_c[this->m_order - 2];
            auto dsn = this->m_b[this->m_order - 1]
                * this->m_c[this->m_order - 1]
                - this->m_b[this->m_order]
                * this->m_c[this->m_order - 2];
            dr = drn / dn;
            ds = dsn / dn;
        }

        r += dr;
        s += ds;
        ++iter;
        if (std::abs(dr) + std::abs(ds) <= this->m_eps)
        {
            // Before exiting, solve the quadratic, refine the roots,
            // multiply out the resulting polynomial, extract
            // the new coefficients, get new dr, r, ds, s.
            auto w = quadratic(s, r, Real{1});
            if (w[0].index() == 1 && w[1].index() == 1)
            s_refine_quadratic_eqn<1>(dr, r, ds, s, w);
            else if (w[0].index() == 2 && w[1].index() == 2)
            s_refine_quadratic_eqn<2>(dr, r, ds, s, w);
        }
        }
        for (int i = 0; i < this->m_order - 1; ++i)
        this->m_coeff[i] = this->m_b[i];

        std::array<Solution<Real>, 2> z2 = quadratic(s, r, Real{1});

        this->m_add_zero(z2[0]);
        this->m_add_zero(z2[1]);
    }

    template<typename Real>
    std::vector<Solution<Real>>
    BairstowSolver<Real>::solve()
    {
        this->m_eps = s_eps;

        this->m_b[0] = Real{1};
        this->m_c[0] = Real{1};

        while (this->m_order > 2)
        this->m_iterate();
        if (this->m_order == 1)
        this->m_add_zero(-this->m_coeff[1]);

        return this->m_zero;
    }

    /**
    * Return the equations fouble by Bairstow's method.
    * There will be order/2 quadratic equations of the form
    * a[k] + a[k+1]x + x^2 = 0
    * and, if there is an odd term, a linear equation
    * a[order] + t = 0.
    */
    template<typename Real>
    std::vector<Real>
    BairstowSolver<Real>::equations() const
    {
        std::vector<Real> eqs(this->m_coeff.rbegin(), this->m_coeff.rend());
        eqs.erase(eqs.begin() + eqs.size() - 1, eqs.end());
        eqs[eqs.size() - 1] *= Real{-1};
        return eqs;
    }

    template<typename Real>
    template<int Index>
        void
        BairstowSolver<Real>::s_refine_quadratic_eqn(Real& dr, Real& r,
                        Real& ds, Real& s,
                        std::array<Solution<Real>, 2>& w)
        {
        const auto p = std::array<Real, 3>{s, r, Real{1}};
        w[0] = refine_solution_newton<3>(std::get<Index>(w[0]), p);
        w[1] = refine_solution_newton<3>(std::get<Index>(w[1]), p);
        const auto sp = std::get<Index>(w[0]) * std::get<Index>(w[1]);
        const auto rp = -(std::get<Index>(w[0]) + std::get<Index>(w[1]));
        dr = std::real(rp - r);
        r = std::real(rp);
        ds = std::real(sp - s);
        s = std::real(sp);
        }

    /**
    * A solver for real-coefficient polynomials due to Jenkins and Traub.
    */
    template<typename Real>
    class JenkinsTraubSolver
    {
    public:

        JenkinsTraubSolver(const std::vector<Real>& op);
        JenkinsTraubSolver(std::vector<Real>&& op);

        std::vector<Solution<Real>> solve();

    private:

        enum NormalizationType
        {
        none,
        divide_by_c,
        divide_by_d,
        near_h_root
        };

        void quadratic(Real a, Real b, Real c,
            Solution<Real> &z_small, Solution<Real> &z_large);
        int fxshfr(int l2);
        int iter_quadratic(Real uu, Real vv);
        int iter_real(Real sss, int& iflag);
        NormalizationType init_next_h_poly();
        void next_h_poly(NormalizationType type);
        std::pair<Real, Real> quadratic_coefficients(NormalizationType type);
        void remquo_quadratic(int n, Real u, Real v,
                std::vector<Real>& poly, std::vector<Real>& quot,
                Real& a, Real& b);

        static constexpr auto s_eps = std::numeric_limits<Real>::epsilon();
        static constexpr auto s_base = Real{std::numeric_limits<Real>::radix};
        static constexpr auto s_tiny = s_eps * s_eps * s_eps; // 1.0e-50; //std::numeric_limits<Real>::min();
        static constexpr auto s_huge = std::numeric_limits<Real>::max();
        //static constexpr auto s_low = s_tiny / s_eps;
        const Real s_low = s_tiny / s_eps;
        static constexpr auto s_pi = Real{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
        static constexpr auto s_rotation = Real{94} * s_pi / Real{180};

        int m_max_iter_quadratic = 20;
        Real m_min_log_deriv = Real{0.005L};
        int m_max_iter_real = 10;
        // Epsilon parameter.
        Real m_are = s_eps;
        // Epsilon parameter.
        Real m_mre = s_eps;

        std::vector<Real> m_P;
        std::vector<Real> m_P_quot;
        std::vector<Real> m_H, m_H_quot, m_H_save;
        Real m_sr, m_si;
        Real m_u, m_v;
        Real m_a;
        Real m_b;
        Real m_c;
        Real m_d;
        Real m_e;
        Real m_f;
        Real m_g;
        Real m_h;
        Real m_a1;
        Real m_a2;
        Real m_a3;
        Real m_a7;
        Solution<Real> m_z_small;
        Solution<Real> m_z_large;
        int m_order;
        bool m_zerok;
        int m_num_iters = 0;
    };


    /**
    * Constructor from input polynomial.
    */
    template<typename Real>
    JenkinsTraubSolver<Real>::
    JenkinsTraubSolver(const std::vector<Real>& op)
    : m_P(op)
    {
        if (this->m_P.size() == 0)
        throw std::domain_error("Polynomial degree must be at least 1.");

        // Algorithm fails if the leading coefficient is zero.
        // We could erase leading-order zero coefficients.
        if (this->m_P[0] == Real{0})
        throw std::domain_error("Leading coefficient must be nonzero.");

        const auto degree = this->m_P.size() - 1;
        this->m_order = degree;
        this->m_P_quot.resize(degree + 1);
        this->m_H.resize(degree + 1);
        this->m_H_quot.resize(degree + 1);
        this->m_H_save.resize(degree + 1);
    }

    /**
    *
    */
    template<typename Real>
    std::vector<Solution<Real>>
    JenkinsTraubSolver<Real>::solve()
    {
        // Initialization of constants for shift rotation.
        auto xx = 1 / Real{1.4142'13562'37309'50488'01688'72420'96980'78569e+0L};
        auto yy = -xx;
        const auto cosr = std::cos(s_rotation);
        const auto sinr = std::sin(s_rotation);

        std::vector<Solution<Real>> zero;
        zero.reserve(this->m_P.size());

        // Remove the zeros at the origin, if any.
        while (this->m_P[this->m_order] == Real{0})
        {
        zero.push_back(Real{0});
        --this->m_order;
        }
        if (this->m_order < 1)
        return zero;

        std::vector<Real> pt(this->m_order + 1);
        std::vector<Real> _H_temp(this->m_order + 1);

        while (true)
        {
        // Start the algorithm for one zero.
        this->m_num_iters = 0;
        if (this->m_order == 1)
        {
            zero.push_back(-this->m_P[1] / this->m_P[0]);
            --this->m_order;
            return zero;
        }
        // Calculate the final zero or pair of zeros.
        if (this->m_order == 2)
        {
            Solution<Real> z_small, z_large;
            this->quadratic(this->m_P[0], this->m_P[1], this->m_P[2],
                    z_small, z_large);
            if (z_small.index() != 0)
            {
            zero.push_back(z_small);
            --this->m_order;
            }
            if (z_large.index() != 0)
            {
            zero.push_back(z_large);
            --this->m_order;
            }
            return zero;
        }

        // Find largest and smallest moduli of coefficients.
        auto a_max = Real{0};
        auto a_min = s_huge;
        for (int i = 0; i <= this->m_order; ++i)
        {
            const auto x = std::abs(this->m_P[i]);
            if (x > a_max)
            a_max = x;
            if (x != Real{0} && x < a_min)
            a_min = x;
        }
        // Scale if there are large or very tiny coefficients.
        // Computes a scale factor to multiply the coefficients
        // of the polynomial. The scaling is done to avoid overflow
        // and to avoid undetected underflow interfering
        // with the convergence criterion.
        // The factor is a power of the base.
        auto scale = s_low / a_min;
        bool rescale = true;
        if (scale > Real{1} && s_huge / scale < a_max)
        rescale = false;
        if (scale <= Real{1} && a_max < Real{10})
        rescale = false;

        if (rescale)
        {
            // Scale polynomial.
            if (scale == Real{0})
            scale = s_tiny;
            const auto l = std::ilogb(scale);
            const auto factor = std::pow(s_base, l);
            if (factor != Real{1})
            for (int i = 0; i <= this->m_order; ++i)
            this->m_P[i] *= factor;
        }

        // Compute lower bound on moduli of roots.
        for (int i = 0; i <= this->m_order; ++i)
        pt[i] = std::abs(this->m_P[i]);
        pt[this->m_order] = -pt[this->m_order];
        // Compute upper estimate of bound.
        auto x = std::exp((std::log(-pt[this->m_order])
                    - std::log(pt[0])) / Real(this->m_order));
        // If Newton step at the origin is better, use it.	
        if (pt[this->m_order - 1] != Real{0})
        {
            const auto xm = -pt[this->m_order] / pt[this->m_order - 1];
            if (xm < x)
            x = xm;
        }
        // Chop the interval (0,x) until ff <= 0.
        while (true)
        {
            auto xm = x * Real{0.1L};
            auto ff = pt[0];
            for (int i = 1; i <= this->m_order; ++i)
            ff = ff * xm + pt[i];
            if (ff <= Real{0})
            break;
            x = xm;
        }
        // Do Newton interation until x converges to two decimal places.
        auto dx = x;
        while (std::abs(dx / x) > this->m_min_log_deriv)
        {
            auto ff = pt[0];
            auto df = ff;
            for (int i = 1; i < this->m_order; ++i)
            {
            ff = ff * x + pt[i];
            df = df * x + ff;
            }
            ff = ff * x + pt[this->m_order];
            dx = ff / df;
            x -= dx;
            ++this->m_num_iters;
        }
        const auto bound = x;
        // Compute the derivative as the initial _H polynomial
        // and do 5 steps with no shift.
        const auto nm1 = this->m_order - 1;
        for (int i = 1; i < this->m_order; ++i)
        this->m_H[i] = Real(this->m_order - i) * this->m_P[i]
                / Real(this->m_order);
        this->m_H[0] = this->m_P[0];
        const auto aa = this->m_P[this->m_order];
        const auto bb = this->m_P[this->m_order - 1];
        this->m_zerok = (this->m_H[this->m_order - 1] == Real{0});
        for(int jj = 0; jj < 5; ++jj)
        {
            ++this->m_num_iters;
            auto cc = this->m_H[this->m_order - 1];
            if (!this->m_zerok)
            {
            // Use a scaled form of recurrence if value of H at 0
            // is nonzero.	
            const auto t = -aa / cc;
            for (int i = 0; i < nm1; ++i)
            {
                const auto j = this->m_order - i - 1;
                this->m_H[j] = t * this->m_H[j - 1] + this->m_P[j];
            }
            this->m_H[0] = this->m_P[0];
            this->m_zerok = (std::abs(this->m_H[this->m_order - 1])
                    <= Real{10} * s_eps * std::abs(bb));
            }
            else
            {
            // Use unscaled form of recurrence.
            for (int i = 0; i < nm1; ++i)
            {
                const auto j = this->m_order - i - 1;
                this->m_H[j] = this->m_H[j - 1];
            }
            this->m_H[0] = Real{0};
            this->m_zerok = (this->m_H[this->m_order - 1] == Real{0});
            }
        }
        // Save H for restarts with new shifts.
        _H_temp = this->m_H;

        // Loop to select the quadratic corresponding to each new shift.
        for (int count = 0; count < 20; ++count)
        {
            /*  Quadratic corresponds to a Real shift to a	
            *  non-real point and its complex conjugate. The point
            *  has modulus bound and amplitude rotated by 94 degrees
            *  from the previous shift.
            */
            const auto xxx = cosr * xx - sinr * yy;
            yy = sinr * xx + cosr * yy;
            auto xx = xxx;
            this->m_sr = bound * xx;
            this->m_si = bound * yy;
            this->m_u = -Real{2} * this->m_sr;
            this->m_v = bound;
            auto num_zeros = this->fxshfr(20 * (count + 1));
            bool cycle = false;
            if (num_zeros != 0)
            {
            /*  The second stage jumps directly to one of the third
            *  stage iterations and returns here if successful.
            *  Deflate the polynomial, store the zero or zeros
            *  and return to the main algorithm.
            */
            zero.push_back(this->m_z_small);
            this->m_order -= num_zeros;
            this->m_P = this->m_P_quot;
            if (num_zeros != 1)
            zero.push_back(this->m_z_large);
            cycle = true;
            break;
            }
            if (cycle)
            continue;

            // If the iteration is unsuccessful another quadratic
            // is chosen after restoring H.
            this->m_H = _H_temp;
        }
        }
    }


    /**
    * Computes up to L2 fixed shift H-polynomials, testing for convergence
    * in the linear or quadratic case.
    * Initiates one of the variable shift iterations and returns
    * the number of zeros found.
    */
    template<typename Real>
    int
    JenkinsTraubSolver<Real>::fxshfr(int l2)
    {
        Real ts_old{}, tv_old{};
        int iflag;

        int num_zeros = 0;

        auto betav = Real{0.25};
        auto betas = Real{0.25};
        auto ss_old = this->m_sr;
        auto vv_old = this->m_v;
        // Evaluate polynomial by synthetic division.
        this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
                this->m_P, this->m_P_quot,
                this->m_a, this->m_b);
        auto type = this->init_next_h_poly();
        for (int j = 0; j < l2; ++j)
        {
        // Calculate next H polynomial and estimate v.
        this->next_h_poly(type);
        type = this->init_next_h_poly();
        auto [ui, vi] = this->quadratic_coefficients(type);
        auto vv = vi;
        // Estimate s.
        auto ss = Real{0};
        if (this->m_H[this->m_order - 1] != Real{0})
        ss = -this->m_P[this->m_order] / this->m_H[this->m_order - 1];
        auto tv = Real{1};
        auto ts = Real{1};
        if (j == 0 || type == near_h_root)
        {
            vv_old = vv;
            ss_old = ss;
            tv_old = tv;
            ts_old = ts;
            continue;
        }

        // Compute relative measures of convergence of s and v sequences.
        if (vv != Real{0})
        tv = std::abs((vv - vv_old) / vv);
        if (ss != Real{0})
        ts = std::abs((ss - ss_old) / ss);

        // If decreasing, multiply two most recent convergence measures.
        const auto tvv = tv < tv_old ? tv * tv_old : Real{1};
        const auto tss = ts < ts_old ? ts * ts_old : Real{1};

        // Compare with convergence criteria.
        const auto vpass = tvv < betav;
        const auto spass = tss < betas;
        if (!(spass || vpass))
        {
            vv_old = vv;
            ss_old = ss;
            tv_old = tv;
            ts_old = ts;
            continue;
        }

        // At least one sequence has passed the convergence test.
        // Store variables before iterating.
        const auto u_save = this->m_u;
        const auto v_save = this->m_v;
        this->m_H_save = this->m_H;
        const auto s = ss;

        // Choose iteration according to the fastest converging sequence.
        auto vtry = false;
        auto stry = false;
        if ((spass && !vpass) || tss < tvv)
        goto TRY_LINEAR;

    TRY_QUADRATIC:
        num_zeros = this->iter_quadratic(ui, vi);
        if (num_zeros > 0)
        return num_zeros;
        // Quadratic iteration has failed. Flag that it has
        // been tried and decrease the convergence criterion.
        vtry = true;
        betav *= Real{0.25};
        // Try linear iteration if it has not been tried and
        // the S sequence is converging.
        if (stry || !spass)
        goto RESTORE_VARS;
        this->m_H = this->m_H_save;

    TRY_LINEAR:
        num_zeros = this->iter_real(s, iflag);
        if (num_zeros > 0)
        return num_zeros;
        // Linear iteration has failed. Flag that it has been
        // tried and decrease the convergence criterion.
        stry = true;
        betas *= Real{0.25};
        if (iflag == 0)
        goto RESTORE_VARS;
        // If linear iteration signals an almost real
        // zero attempt quadratic iteration.
        ui = -Real{2} * s;
        vi = s * s;
        goto TRY_QUADRATIC;

    RESTORE_VARS:
        // Restore variables.
        this->m_u = u_save;
        this->m_v = v_save;
        this->m_H = this->m_H_save;

        // Try quadratic iteration if it has not been tried
        // and the V sequence is converging.
        if (!vtry && vpass)
        goto TRY_QUADRATIC;

        // Recompute polynomial quotient and remainder
            // to continue the second stage.
        this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
                    this->m_P, this->m_P_quot,
                    this->m_a, this->m_b);
        type = this->init_next_h_poly();
        }
        return num_zeros;
    }


    /**
    * Variable-shift H-polynomial iteration for a quadratic factor
    * converges only if the zeros are equimodular or nearly so.
    * @param uu The linear coefficient of the starting quadratic equation.
    * @param vv The constant  coefficient of the starting quadratic equation.
    * @return The number of zeros found.
    */
    template<typename Real>
    int
    JenkinsTraubSolver<Real>::iter_quadratic(Real uu, Real vv)
    {
        Real mp, mp_old{}, ee, relstp = std::sqrt(s_eps), t, zm;
        NormalizationType type;

        int num_zeros = 0;
        bool tried = false;
        this->m_u = uu;
        this->m_v = vv;
        int j = 0;

        while (true)
        {
        ++this->m_num_iters;
        this->quadratic(Real{1}, this->m_u, this->m_v,
                this->m_z_small, this->m_z_large);
        // Return if roots of the quadratic are real and not
        // close to multiple or nearly equal and of opposite sign.
        if (std::abs(std::abs(real(this->m_z_small))
            - std::abs(real(this->m_z_large)))
            > Real{0.01L} * std::abs(real(this->m_z_large)))
        return num_zeros;
        // Evaluate polynomial by quadratic synthetic division.
        this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
                    this->m_P, this->m_P_quot, this->m_a, this->m_b);
        mp = std::abs(this->m_a - real(this->m_z_small) * this->m_b)
        + std::abs(imag(this->m_z_small) * this->m_b);
        // Compute a rigorous bound on the rounding error in evaluating _P.
        zm = std::sqrt(std::abs(this->m_v));
        ee = Real{2} * std::abs(this->m_P_quot[0]);
        t = -real(this->m_z_small) * this->m_b;
        for (int i = 1; i < this->m_order; ++i)
        ee = ee * zm + std::abs(this->m_P_quot[i]);
        ee = ee * zm + std::abs(this->m_a + t);
        ee *= (Real{5} * this->m_mre + Real{4} * this->m_are);
        ee -= (Real{5} * this->m_mre + Real{2} * this->m_are)
            * (std::abs(this->m_a + t) + std::abs(this->m_b) * zm);
        ee += Real{2} * this->m_are * std::abs(t);
        // Iteration has converged sufficiently if the
        // polynomial value is less than 20 times this bound.
        if (mp <= Real{20} * ee)
        {
            num_zeros = 2;
            return num_zeros;
        }
        ++j;
        // Stop iteration after 20 steps.
        if (j > this->m_max_iter_quadratic)
        return num_zeros;
        if (j < 2 || relstp > Real{0.01L} || mp < mp_old || tried)
        {
            mp_old = mp;
            // Calculate next H polynomial and new u and v.
            type = this->init_next_h_poly();
            this->next_h_poly(type);
            type = this->init_next_h_poly();
            const auto [ui, vi] = this->quadratic_coefficients(type);
            // If vi is zero the iteration is not converging.
            if (vi == Real{0})
            return num_zeros;
            relstp = std::abs((vi - this->m_v) / vi);
            this->m_u = ui;
            this->m_v = vi;
            continue;
        }
        // A cluster appears to be stalling the convergence.
        // Five fixed shift steps are taken with a u, v close to the cluster.
        if (relstp < s_eps)
        relstp = s_eps;
        relstp = std::sqrt(relstp);
        this->m_u -= this->m_u * relstp;
        this->m_v += this->m_v * relstp;
        this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
                    this->m_P, this->m_P_quot,
                    this->m_a, this->m_b);
        for (int i = 0; i < 5; ++i)
        {
            type = this->init_next_h_poly();
            this->next_h_poly(type);
        }
        tried = true;
        j = 0;
        }
    }


    /**
    * Variable-shift H polynomial iteration for a real zero.
    * @param sss Starting iterate.
    * @param iflag Flag to indicate a pair of zeros near real axis.
    * @return The number of zeros found.
    */
    template<typename Real>
    int
    JenkinsTraubSolver<Real>::iter_real(Real sss, int& iflag)
    {
        auto t = Real{0};
        decltype(std::abs(this->m_P[0])) mp_old{};

        int num_zeros = 0;
        auto s = sss;
        iflag = 0;
        int i_real = 0;

        while (true)
        {
        ++this->m_num_iters;
        auto pval = this->m_P[0];
        // Evaluate P at s.
        this->m_P_quot[0] = pval;
        for (int i = 1; i <= this->m_order; ++i)
        {
            pval = pval * s + this->m_P[i];
            this->m_P_quot[i] = pval;
        }
        auto mp = std::abs(pval);
        // Compute a rigorous bound on the error in evaluating P.
        const auto ms = std::abs(s);
        auto ee = (this->m_mre / (this->m_are + this->m_mre))
            * std::abs(this->m_P_quot[0]);
        for (int i = 1; i <= this->m_order; ++i)
        ee = ee * ms + std::abs(this->m_P_quot[i]);
        // Iteration has converged sufficiently if the polynomial
        // value is less than 20 times this bound.
        if (mp <= Real{20}
            * ((this->m_are + this->m_mre) * ee - this->m_mre * mp))
        {
            num_zeros = 1;
            this->m_z_small = s;
            return num_zeros;
        }
        ++i_real;
        // Stop iteration after max_iter_real steps.
        if (i_real > this->m_max_iter_real)
        return num_zeros;
        else if (i_real < 2
        || std::abs(t) > Real{0.001L} * std::abs(s - t)
        || mp < mp_old)
        {
            // Return if the polynomial value has increased significantly.
            mp_old = mp;

            // Compute t, the next polynomial, and the new iterate.
            auto hval = this->m_H[0];
            this->m_H_quot[0] = hval;
            for (int i = 1; i < this->m_order; ++i)
            {
            hval = hval * s + this->m_H[i];
            this->m_H_quot[i] = hval;
            }

            if (std::abs(hval)
            <= std::abs(this->m_H[this->m_order - 1]) * Real{10} * s_eps)
            {
            // Use unscaled form.
            this->m_H[0] = Real{0};
            for (int i = 1; i < this->m_order; ++i)
            this->m_H[i] = this->m_H_quot[i-1];
            }
            else
            {
            // Use the scaled form of the recurrence if the value
            // of H at s is nonzero.
            t = -pval / hval;
            this->m_H[0] = this->m_P_quot[0];
            for (int i = 1; i < this->m_order; ++i)
            this->m_H[i] = t * this->m_H_quot[i - 1]
                    + this->m_P_quot[i];
            }

            hval = this->m_H[0];
            for (int i = 1; i < this->m_order; ++i)
            hval = hval * s + this->m_H[i];
            auto t = Real{0};
            if (std::abs(hval)
            > std::abs(this->m_H[this->m_order - 1] * Real{10} * s_eps))
            t = -pval / hval;
            s += t;
        }
        else
        {
            // A cluster of zeros near the real axis has been encountered.
            // Return with iflag set to initiate a quadratic iteration.
            iflag = 1;
            sss = s;
            return num_zeros;
        }
        }
    }

    /**
    * This routine calculates scalar quantities used to compute
    * the next H-polynomial and new estimates of the quadratic coefficients.
    *
    * @return Flag indicating how the calculations are normalized
    *         to avoid overflow.
    */
    template<typename Real>
    typename JenkinsTraubSolver<Real>::NormalizationType
    JenkinsTraubSolver<Real>::init_next_h_poly()
    {
        const auto eps = Real{100} * s_eps;
        // Synthetic division of H by the quadratic 1, u, v
        NormalizationType type = none;
        this->remquo_quadratic(this->m_order - 1, this->m_u, this->m_v,
                this->m_H, this->m_H_quot, this->m_c, this->m_d);
        if (std::abs(this->m_c) > eps * std::abs(this->m_H[this->m_order - 1])
        || std::abs(this->m_d) > eps * std::abs(this->m_H[this->m_order - 2]))
        {
        if (std::abs(this->m_d) < std::abs(this->m_c))
        {
            type = divide_by_c;
            this->m_e = this->m_a / this->m_c;
            this->m_f = this->m_d / this->m_c;
            this->m_g = this->m_u * this->m_e;
            this->m_h = this->m_v * this->m_b;
            this->m_a3 = this->m_a * this->m_e
                + (this->m_h / this->m_c + this->m_g) * this->m_b;
            this->m_a1 = this->m_b - this->m_a * (this->m_d / this->m_c);
            this->m_a7 = this->m_a
                + this->m_g * this->m_d
                + this->m_h * this->m_f;
            return type;
        }
        else
        {
            type = divide_by_d;
            this->m_e = this->m_a / this->m_d;
            this->m_f = this->m_c / this->m_d;
            this->m_g = this->m_u * this->m_b;
            this->m_h = this->m_v * this->m_b;
            this->m_a3 = (this->m_a + this->m_g) * this->m_e
                + this->m_h * (this->m_b / this->m_d);
            this->m_a1 = this->m_b * this->m_f - this->m_a;
            this->m_a7 = (this->m_f + this->m_u) * this->m_a + this->m_h;
            return type;
        }
        }
        else
        {
        type = near_h_root;
        return type;
        }
    }


    /**
    * Computes the next H polynomials using scalars computed in init_next_h_poly.
    */
    template<typename Real>
    void
    JenkinsTraubSolver<Real>::next_h_poly(NormalizationType type)
    {
        if (type == near_h_root)
        {
        // Use unscaled form of the recurrence if type is 3.
        this->m_H[0] = Real{0};
        this->m_H[1] = Real{0};
        for (int i = 2; i < this->m_order; ++i)
        this->m_H[i] = this->m_H_quot[i-2];
        return;
        }
        auto ab_temp = this->m_a;
        if (type == divide_by_c)
        ab_temp = this->m_b;
        if (std::abs(this->m_a1) <= std::abs(ab_temp) * s_eps * Real{10})
        {
        // If a1 is nearly zero then use a special form of the recurrence.
        this->m_H[0] = Real{0};
        this->m_H[1] = -this->m_a7 * this->m_P_quot[0];
        for(int i = 2; i < this->m_order; ++i)
        this->m_H[i] = this->m_a3 * this->m_H_quot[i-2]
                - this->m_a7 * this->m_P_quot[i-1];
        }
        else
        {
        // Use scaled form of the recurrence.
        this->m_a7 /= this->m_a1;
        this->m_a3 /= this->m_a1;
        this->m_H[0] = this->m_P_quot[0];
        this->m_H[1] = this->m_P_quot[1] - this->m_a7 * this->m_P_quot[0];
        for (int i = 2; i < this->m_order; ++i)
        this->m_H[i] = this->m_a3 * this->m_H_quot[i-2]
                - this->m_a7 * this->m_P_quot[i-1]
                + this->m_P_quot[i];
        }
    }

    /**
    * Compute new estimates of the quadratic coefficients
    * using the scalars computed in init_next_h_poly.
    */
    template<typename Real>
    std::pair<Real, Real>
    JenkinsTraubSolver<Real>::quadratic_coefficients(NormalizationType type)
    {
        if (type == near_h_root)
        return std::make_pair(Real{0}, Real{0});

        Real a4, a5;
        if (type == divide_by_d)
        {
        a4 = (this->m_a + this->m_g) * this->m_f + this->m_h;
        a5 = (this->m_f + this->m_u) * this->m_c + this->m_v * this->m_d;
        }
        else
        {
        a4 = this->m_a + this->m_u * this->m_b + this->m_h * this->m_f;
        a5 = this->m_c + (this->m_u + this->m_v * this->m_f) * this->m_d;
        }

        // Evaluate new quadratic coefficients.
        const auto n = this->m_order;
        const auto b1 = -this->m_H[n - 1] / this->m_P[n];
        const auto b2 = -(this->m_H[n - 2] + b1 * this->m_P[n - 1]) / this->m_P[n];
        const auto c1 = this->m_v * b2 * this->m_a1;
        const auto c2 = b1 * this->m_a7;
        const auto c3 = b1 * b1 * this->m_a3;
        const auto c4 = c1 - c2 - c3;
        if (auto temp = a5 + b1 * a4 - c4; temp == Real{0})
        return std::make_pair(Real{0}, Real{0});
        else
        {
        auto uu = this->m_u - (this->m_u * (c3 + c2)
            + this->m_v * (b1 * this->m_a1 + b2 * this->m_a7)) / temp;
        auto vv = this->m_v * (Real{1} + c4 / temp);
        return std::make_pair(uu, vv);
        }
    }

    /**
    * Divides the polynomial P by the quadratic 1x^2 + ux + v
    * placing the quotient in q and the remainder in a, b.
    */
    template<typename Real>
    void
    JenkinsTraubSolver<Real>::remquo_quadratic(int nn, Real u, Real v,
                            std::vector<Real>& poly,
                            std::vector<Real>& quot,
                            Real& a, Real& b)
    {
        b = poly[0];
        quot[0] = b;
        a = poly[1] - b * u;
        quot[1] = a;
        for (int i = 2; i <= nn; ++i)
        {
        auto c = poly[i] - a * u - b * v;
        quot[i] = c;
        b = a;
        a = c;
        }	
    }


    /**
    * Calculate the zeros of the quadratic az^2 + bz + c.
    * The quadratic formula, modified to avoid overflow, is used
    * to find the larger zero if the zeros are real and both
    * are complex. The smaller real zero is found directly from
    * the product of the zeros c/a.
    */
    template<typename Real>
    void
    JenkinsTraubSolver<Real>::quadratic(Real a, Real b, Real c,
                        Solution<Real>& z_small,
                        Solution<Real>& z_large)
    {
        z_small = {};
        z_large = {};
        if (a == Real{0})
        { // Less than two roots.
        if (b != Real{0})
        z_small = -c / b;
        return;
        }

        if (c == Real{0})
        { // one real root, one zero root.
        z_small = Real{0};
        z_large = -b / a;
        return;
        }

        // Compute discriminant avoiding overflow.
        auto b2 = b / Real{2};

        Real d, e;
        if (std::abs(b2) < std::abs(c))
        {
        e = std::copysign(a, c);
        e = b2 * (b2 / std::abs(c)) - e;
        d = std::sqrt(std::abs(e)) * std::sqrt(std::abs(c));
        }
        else
        {
        e = Real{1} - (a / b2) * (c / b2);
        d = std::sqrt(std::abs(e)) * std::abs(b2);
        }

        if (e < Real{0})
        { // complex conjugate zeros.
        z_small = std::complex<Real>(-b2 / a, +std::abs(d / a));
        z_large = std::complex<Real>(-b2 / a, -std::abs(d / a));
        }
        else
        {
        if (b2 >= Real{0})
        d = -d; // Real zeros.
        z_large = (-b2 + d) / a;
        z_small = Real{0};
        if (z_large != Real{0})
        z_small = (c / z_large) / a;
        }
    }

    /**
    * A solver for complex-coefficient polynomials due to Laguerre.
    */
    template<typename Real>
        class LaguerreSolver
        {
        public:

        LaguerreSolver(Polynomial<std::complex<Real>>& P)
        : m_poly(P), m_num_iters{0}
        { }

        std::vector<std::complex<Real>> solve();

        std::complex<Real>
        step()
        {
        using cmplx = std::complex<Real>;

        const auto z0 = this->m_root_laguerre();

        Polynomial<cmplx> zpoly({-z0, cmplx{1}});
        this->m_poly /= zpoly;

        return z0;
        }

        int
        num_iters() const
        { return this->m_num_iters; }

        int
        max_num_iters() const
        { return this->m_max_iter(); }

        int
        num_steps_per_frac() const
        { return this->m_steps_per_frac; }

        const Polynomial<std::complex<Real>>&
        polynomial() const
        { return this->m_poly; }

        LaguerreSolver&
        num_steps_per_frac(int num)
        {
        this->m_steps_per_frac = num;
        return *this;
        }

        private:

        // Estimated fractional roundoff error.
        static constexpr Real s_eps = std::numeric_limits<Real>::epsilon();

        // Number of fractional values.
        static constexpr int s_num_fracs = 8;
        // Fractions used to break a limit cycle (in a heap).
        static constexpr Real
        s_frac[s_num_fracs + 1]
        {0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 1.0};

        // Number of steps taken before trying a new fraction.
        int m_steps_per_frac = 10;

        int
        m_max_iter()
        { return this->m_steps_per_frac * s_num_fracs; }

        std::complex<Real> m_root_laguerre();

        Polynomial<std::complex<Real>> m_poly;

        int m_num_iters = 0;
        };

        /**
    * Find a root of a complex-coefficient polynomial by Laguerre's method.
    * This routine can be iterated by dividing the original polynomial
    * by the root factor (z - x) where x is the found root and finding
    * the next root.
    */
    template<typename Real>
    std::complex<Real>
    LaguerreSolver<Real>::m_root_laguerre()
    {
        using cmplx = std::complex<Real>;

        const auto m = this->m_poly.degree();
        const int max_iter = this->m_max_iter();

        this->m_num_iters = 0;

        //if (m == 1)
        //return -this->m_poly.cefficient(0) / this->m_poly.cefficient(1);

        cmplx x{};
        for (int iter = 1; iter <= max_iter; ++iter)
        {
        ++this->m_num_iters;

        // Efficient computation of the polynomial
        // and its first two derivatives. F stores P''(x)/2.
        auto b = this->m_poly[m];
        auto err = std::abs(b);
        const auto abx = std::abs(x);
        cmplx d{}, f{};
        for (int j = m - 1; j >= 0; --j)
            {
            f = x * f + d;
            d = x * d + b;
            b = x * b + this->m_poly[j];
            err = abx * err + std::abs(b);
            }
        err *= s_eps;
        // Estimate of roundoff error in evaluating polynomial.
        if (std::abs(b) <= err) // We have the root.
            return x;

        // Use Laguerre's formula.
        const auto g = d / b;
        const auto g2 = g * g;
        const auto h = g2 - Real{2} * f / b;
        const auto sq = std::sqrt(Real(m - 1)
                    * (Real(m) * h - g2));
        auto gp = g + sq;
        const auto gm = g - sq;
        const auto abp = std::abs(gp);
        const auto abm = std::abs(gm);
        if (abp < abm)
            gp = gm;
        const auto dx = std::max(abp, abm) > Real{0}
                ? Real(m) / gp
                : std::polar(Real{1} + abx, Real(iter));
        const auto x1 = x - dx;
        if (x == x1)
            return x;
        if (iter % this->m_steps_per_frac != 0)
            x = x1;
        else
            x -= s_frac[iter / this->m_steps_per_frac] * dx;
        }

        throw std::runtime_error("m_root_laguerre: Maximum number of iterations exceeded");
    }    

    /**
 * Return the L1 sum of absolute values or Manhattan metric of a complex number.
 */
template<typename Real>
  inline Real
  norm_l1(std::complex<Real> z)
  {
    return std::abs(std::real(z)) + std::abs(std::imag(z));
  }

/**
 * Return the L2 modulus of the complex number (this is std::norm).
 */
template<typename Real>
  inline Real
  norm_l2(std::complex<Real> z)
  {
    return std::norm(z);
  }

template<typename Real>
  class SolverMadsenReid
  {
  public:

    using Cmplx = std::complex<Real>;

    /**
     * Constructor.
     *
     * @param a_in Coefficient of the polynomial in "big-endian" - largest degree coefficient first - order.
     */
    SolverMadsenReid(const std::vector<Cmplx>& a_in)
    : poly(a_in),
      poly_work(a_in.size())
    {}

    /**
     * Solve the polynomial.
     */
    std::vector<Cmplx>
    solve()
    {
        std::vector<Cmplx> root;
        if (poly.size() <= 1)
            return root;

        int degree = poly.size() - 1;
        root.resize(degree);
        solve(poly, degree, root, poly_work);
        return root;
    }

  private:

    static constexpr Real DIGITS = std::numeric_limits<Real>::max_digits10;
    static constexpr Real BIG = std::numeric_limits<Real>::max(); // Overflow limit
    static constexpr Real SMALL = std::numeric_limits<Real>::min(); // Underflow limit.
    static constexpr Real BASE = std::numeric_limits<Real>::radix;
    static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

    /// Big-endian polynomial.
    std::vector<Cmplx> poly;

    /// Big-endian working polynomial.
    std::vector<Cmplx> poly_work;

    /**
     * Evaluate polynomial at z, set fz, return squared modulus.
     *
     * @param  z     Argument of the polynomial.
     * @param  fz    Polynomial value at the given argument.
     * @param  size  Size (degree + 1) of the polynomial polynomial.
     * @param  a     Polynomial coefficients.
     * @return  Squared modulus of the function value.
     */
    Real
    eval(Cmplx z, Cmplx& fz, int size, const std::vector<Cmplx>& a)
    {
        auto deg = size - 1;
        auto p = a[0];
        for (int i = 0; i < deg; ++i)
        {
            p = p * z + a[i + 1];
        }

        fz = p;

        return norm_l2(p);
    }

    /**
     * Store the root.
     *
     * @param  a1  Working polynomial.
     * @param  root  Roots of the polynomial.
     * @param  z  New root of the polynomial.
     */
    void
    push_root(std::vector<Cmplx>& a1, std::vector<Cmplx>& root, int& n, Cmplx z)
    {
        a1[n - 1] = root[n - 1];
        root[n - 1] = z;
        --n;
    }

    /**
     * Deflate the polynomial.
     *
     * @param  a  Polynomial
     * @param  n  New degree.
     * @param  z  New root of the polynomial.
     */
    void
    deflate(std::vector<Cmplx>& a, int n, Cmplx z)
    {
        for (int k = 1; k < n; ++k)
        {
            a[k] = a[k - 1] * z + a[k];
        }
    }

    /**
     * Root search...
     *
     * @param  a     Input polynomial
     * @param  n     Degree.
     * @param  root  Roots of the polynomial.
     * @param  a     Working polynomial.
     */
    void solve(std::vector<Cmplx>& a1, int m, std::vector<Cmplx>& root, std::vector<Cmplx>& a);
  };

    /**
    * Find all solutions of the  polynomial by stepping and 
    */
    template<typename Real>
    std::vector<std::complex<Real>>
    LaguerreSolver<Real>::solve()
    {
        using cmplx = std::complex<Real>;

        std::vector<cmplx> roots;
        const auto deg = this->m_poly.degree();
        roots.reserve(deg);
        for (unsigned i = 0; i < deg; ++i)
        {
        const auto z0 = this->m_root_laguerre();

        Polynomial<cmplx> zpoly({-z0, cmplx{1}});
        this->m_poly /= zpoly;

        roots.push_back(z0);
        }
        return roots;
        }

   
        /**
    * Root search...
    */
    template<typename Real>
    void
    SolverMadsenReid<Real>::
    solve(std::vector<std::complex<Real>>& a1, int m,
        std::vector<std::complex<Real>>& root, std::vector<std::complex<Real>>& a)
        {
            using Cmplx = std::complex<Real>;

            Cmplx z0, f0z, z, dz, f1z, fz, w, fw, dzk;
            Real f0, ff, f, fa, fmin, f2;
            bool stage1, div2;

            Real r0, u0, r, r1;

            const Real SSMALL = std::sqrt(SMALL);
            const Real ALOGB = std::log(BASE);

            const auto THETA = std::atan(Real{3} / Real{4});
            const auto PHASE = std::polar(Real{1}, -THETA);

            int n = m;

            // Store original polynomial in a and in root.
            int j = m;
            a[j] = a1[0];
            for (int i = 0; i < m; ++i)
            {
                root[i] = a1[i];
                a[i] = a1[j];
                --j;
            }

            // Test for zeros at infinity.
            while (norm_l1(a[0]) <= 0 && n > 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    a[i] = a[i + 1];
                }
                root[n] = BIG;
                --n;
            }

            while (n >= 0)
            {
                if (n == 1)
                {
                    z = -a[1] / a[0];
                    a1[n-1] = root[n-1];
                    root[n-1] = z;
                    return;
                }

                // Scale the coefficients.
                auto u1 = Real {0};
                auto u2 = BIG;
                for(int k = 0; k <= n; ++k)
                {
                    auto u = norm_l1(a[k]);
                    if (u <= Real{0})
                        continue;
                    if (u > u1)
                        u1 = u;
                    if (u < u2)
                        u2 = u;
                }
                auto u = std::sqrt(u1) * std::sqrt(u2);
                int i = -std::log(u) / ALOGB; // ilogb?
                u = std::pow(BASE, Real(i));
                for (int k = 0; k < n; ++k)
                {
                    a[k] = u * a[k];
                    a1[k] = a[k] * Real(n - k);
                }
                a[n] = u * a[n];

                // Test for zeros at (0, 0)
                z = Cmplx(0, 0);
                if (norm_l1(a[n]) <= SSMALL)
                {
                    push_root(a1, root, n, z);
                    continue;
                }
                z0 = Cmplx(0, 0);
                f0 = norm_l2(a[n]);
                fmin = f0 * std::pow(Real(n) * DIGITS * EPS, 2);

                // z is the current point
                // f = |f(z)|^2
                // z0 is the last point
                // f0z = f'(z0)
                // f0 = |f(z0)|^2
                // r0 = 3|z - z0|
                // dz is the last tentative step if the last step was successful or is the required next step.
                ff = f0;
                u0 = f0;
                auto t = BIG;
                for (int k = 0; k < n; ++k)
                {
                    u = norm_l2(a[k]);
                    if (u == Real{0})
                        continue;
                    u = std::log(u0 / u) / Real(2 *(n - k));
                    if (u < t)
                        t = u;
                }
                t = std::exp(t);
                f0z = a[n - 1];
                z = Cmplx(1, 0);
                if (norm_l1(f0z) > Real{0})
                {
                    z = -a[n] / a[n - 1];
                }
                u = 0.5 * t / norm_l1(z);
                z = u * z;
                dz = z;
                f =  eval(z, fz, n + 1, a);
                r0 = 0.5 * t;

            _120:
                // Calculate tentative step.
                u = eval(z, f1z, n, a1);
                if (u == Real{0})
                {
                    dz *= Real{3} * PHASE;
                    stage1 = true;
                }
                else
                {
                    dz = -fz / f1z;
                    f2 = norm_l2(f0z - f1z) / norm_l2(z0 - z);
                    stage1 = (f * f2 / u > 0.25 * u) || (f != ff);
                    r = norm_l1(dz);
                    if (r > Real{3} * r0)
                    {
                        dz *= Real{3} * PHASE * (r0 / r);
                    }
                }

                f0z = f1z;

            _160:
                // Find the next point in the iteration.
                // This is where iteration starts if the previous one was unsuccessful.
                z0 = z;
                f0 = f;
                dzk = dz;
                z = z0 + dz;
                // If either part of z is small replace by zero to avoid underflows.
                if (std::abs(std::real(z)) < EPS * std::abs(std::imag(z)))
                    z = Cmplx(Real{0}, std::imag(z));
                if (std::abs(std::imag(z)) < EPS * std::abs(std::real(z)))
                    z = Cmplx(std::real(z), Real{0});
                w = z;
                f = eval(z, fz, n + 1, a);
                ff = f;
                if (stage1)
                {
                    int j = 1;
                    div2 = f >= f0;

                    do
                    {
                        if (div2)
                        {
                            dz *= Real{0.5L};
                            w = z0 + dz;
                        }
                        else
                        {
                            w += dz;
                        }

                        fa = eval(w, fw, n + 1, a);
                        if (fa >= f)
                        {
                            break;
                        }
                        f = fa;
                        fz = fw;
                        z = w;
                        ++j;
                        if (div2 && j == 3)
                        {
                            dz *= PHASE;
                            z = z0 + dz;
                            f = eval(z, fz, n + 1, a);
                            break;
                        }
                    }
                    while (j <= n);
                }

                r0 = norm_l1(z0 - z);

                // Convergence test.
                if (f >= f0)
                {
                    z = z0;
                }

                r1 = norm_l1(z);
                if (r0 < EPS * r1)
                {
                    deflate(a, n, z);
                    push_root(a1, root, n, z);
                    continue;
                }

                if (f < f0)
                {
                    goto _120;
                }

                f = f0;
                if (f < fmin)
                {
                    deflate(a, n, z);
                    push_root(a1, root, n, z);
                    continue;
                }

                dz = -Real{0.5} * PHASE * dzk;
                stage1 = true;
                goto _160;
            }

            return;
        }

        /**
    * A solver for complex-coefficient polynomials due to Laguerre.
    */
    template<typename Real>
        class QuadraticSolver
        {
        public:

        QuadraticSolver(Polynomial<std::complex<Real>>& P)
        : m_poly(P), m_num_iters{0}
        { }

        //std::vector<Solution<Real>> solve();
        std::vector<std::complex<Real>> solve();

        Polynomial<std::complex<Real>>
        step()
        {
        const auto q = this->m_root_quadratic();
        this->m_poly.deflate(q, Real{10} * s_eps);
        return q;
        }

        int
        num_iters() const
        { return this->m_num_iters; }

        int
        max_num_iters() const
        { return this->m_max_num_iters; }

        const Polynomial<std::complex<Real>>&
        polynomial() const
        { return this->m_poly; }

        private:

        // Estimated fractional roundoff error.
        static constexpr Real s_eps = std::numeric_limits<Real>::epsilon();
        static constexpr Real s_tiny = Real{10} * s_eps;

        // Fractional roundoff error.
        Real m_eps = Real{100} * s_eps;

        int m_max_iter = 50;

        Polynomial<std::complex<Real>> m_root_quadratic();

        Polynomial<std::complex<Real>> m_poly;

        int m_num_iters = 0;
        };


   /**
    * I think this is trying to factor out a quadratic
    * from a complex-coefficient polynomial.
    */
    template<typename Tp>
    Polynomial<std::complex<Tp>>
    QuadraticSolver<Tp>::m_root_quadratic()
    {
        using Cmplx = std::complex<Tp>;
        using Poly = Polynomial<Cmplx>;

        if (this->m_poly.degree() <= 2)
        return this->m_poly;

        this->m_num_iters = 0;

        Cmplx c, b;
        Poly q, qq, rem;
        for (int iter = 0; iter < this->m_max_iter; ++iter)
        {
        ++this->m_num_iters;

        Poly d({c, b, Cmplx{1}});

        // First division: r, s.
        divmod(this->m_poly, d, q, rem);
        const auto s = rem[0];
        const auto r = rem[1];

        // Second division: partial r, s with respect to c.
        divmod(q, d, qq, rem);
        const auto sc = -rem[0];
        const auto rc = -rem[1];
        const auto sb = -c * rc;
        const auto rb = -b * rc + sc;

        // Solve 2x2 equation.
        const auto dv = Tp{1} / (sb * rc - sc * rb);
        const auto delb = ( r * sc - s * rc) * dv;
        b += delb;
        const auto delc = (-r * sb + s * rb) * dv;
        c += delc;
        if ((std::abs(delb) <= this->m_eps * std::abs(b)
            || std::abs(b) < s_tiny)
            && (std::abs(delc) <= this->m_eps * std::abs(c)
            || std::abs(c) < s_tiny))
            return Poly({c, b, Cmplx{1}});
        }
        throw std::runtime_error("m_root_quadratic: Maximum number of iterations exceeded");
    }

} // namespace emsr