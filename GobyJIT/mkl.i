%module mkl
%{
#include "Mkl.hpp"
using namespace MKL;
%}

%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_complex.i"
//%include "lua_fnptr.i"

%inline %{
typedef float MklType;

template<typename T>
std::vector<T> interleave(std::vector<T> & left, std::vector<T> & right) {
    assert(left.size() == right.size());
    std::vector<T> r(left.size()+right.size());
    size_t x = 0;
    for(size_t i = 0; i < left.size(); i++) {
        r[x++] = left[i];
        r[x++] = right[i];
    }
    return r;
}
template<typename T>
void deinterleave(std::vector<T> & chans, std::vector<T> & left, std::vector<T> & right) {    
    left.resize(chans.size()/2);
    right.resize(chans.size()/2);
    for(size_t i = 0; i < chans.size(); i+=2)
    {
        left[i] = chans[i];
        right[i] = chans[i+1];
    }
}

template<typename T>
MKL::Vector<T> interleave(MKL::Vector<T> & left, MKL::Vector<T> & right) 
{
    assert(left.size() == right.size());
    std::vector<T> r(left.size()+right.size());
    size_t x = 0;
    for(size_t i = 0; i < left.size(); i++) {
        r[x++] = left[i];
        r[x++] = right[i];
    }
    return r;
}

template<typename T>
void deinterleave(MKL::Vector<T> & chans, MKL::Vector<T> & left, MKL::Vector<T> & right) {    
    left.resize(chans.size()/2);
    right.resize(chans.size()/2);
    for(size_t i = 0; i < chans.size(); i+=2)
    {
        left[i] = chans[i];
        right[i] = chans[i+1];
    }
}

%}

%template(stdvector) std::vector<MklType>;
//%template(luavector) std::vector<SWIGLUA_REF>;
%template(stringmap) std::map<std::string,std::string>;
//%template(luamap) std::map<std::string,SWIGLUA_REF>;

%include "Mkl.hpp"

%extend MKL::Vector
{
    T    __getitem__(size_t i) { return $self->vector[i]; }
    void __setitem__(size_t i, T x) { $self->vector[i] = x; }

    void println() {
        std::cout << "vector[" << $self->size() << "]=";
        for(size_t i = 0; i < $self->size()-1; i++) 
            std::cout << $self->vector[i] << ",";
        std::cout << $self->vector[$self->size()-1] << std::endl;
    }
}

%extend MKL::Matrix
{
    MatrixView<T> __getitem__(size_t row) { return (*$self)[row]; }

    void println() {
        std::cout << "matrix[" << $self->rows() << "," << $self->cols() << "]=";
        for(size_t j = 0; j < $self->rows(); j++)
        {
            for(size_t i = 0; i < $self->cols()-1; i++) 
                std::cout << (*$self)(j,i) << ",";
            std::cout << (*$self)(j,$self->cols()-1)<< std::endl;
        }
    }
}

%template(Vector) MKL::Vector<MklType>;
%template(Matrix) MKL::Matrix<MklType>;
%template(MatrixView) MKL::MatrixView<MklType>;

%template(flt_complex) std::complex<float>;
%template(dbl_omplex) std::complex<double>;

%template(matmul) MKL::matmul<MklType>;
%template(sqr) MKL::sqr<MklType>;
%template(abs) MKL::abs<MklType>;
%template(inc) MKL::inv<MklType>;
%template(sqrt) MKL::sqrt<MklType>;
%template(rsqrt) MKL::rsqrt<MklType>;
%template(cbrt) MKL::cbrt<MklType>;
%template(rcbrt) MKL::rcbrt<MklType>;
%template(pow) MKL::pow<MklType>;
%template(pow2o3) MKL::pow2o3<MklType>;
%template(pow3o2) MKL::pow3o2<MklType>;
%template(hypot) MKL::hypot<MklType>;
%template(exp) MKL::exp<MklType>;
%template(exp2) MKL::exp2<MklType>;
%template(exp10) MKL::exp10<MklType>;
%template(expm1) MKL::expm1<MklType>;
%template(ln) MKL::ln<MklType>;
%template(log10) MKL::log10<MklType>;
%template(log2) MKL::log2<MklType>;
%template(logb) MKL::logb<MklType>;
%template(log1p) MKL::log1p<MklType>;
%template(cos) MKL::cos<MklType>;
%template(sin) MKL::sin<MklType>;
%template(tan) MKL::tan<MklType>;
%template(acos) MKL::acos<MklType>;
%template(asin) MKL::asin<MklType>;
%template(atan) MKL::atan<MklType>;
%template(atan2) MKL::atan2<MklType>;
%template(cosh) MKL::cosh<MklType>;
%template(sinh) MKL::sinh<MklType>;
%template(tanh) MKL::tanh<MklType>;
%template(acosh) MKL::acosh<MklType>;
%template(asinh) MKL::asinh<MklType>;
%template(atanh) MKL::atanh<MklType>;
%template(sincos) MKL::sincos<MklType>;
%template(erf) MKL::erf<MklType>;
%template(erfinv) MKL::erfinv<MklType>;
%template(erfc) MKL::erfc<MklType>;
%template(cdfnorm) MKL::cdfnorm<MklType>;
%template(cdfnorminv) MKL::cdfnorminv<MklType>;
%template(floor) MKL::floor<MklType>;
%template(ceil) MKL::ceil<MklType>;
%template(trunc) MKL::trunc<MklType>;
%template(round) MKL::round<MklType>;
%template(nearbyint) MKL::nearbyint<MklType>;
%template(rint) MKL::rint<MklType>;
%template(fmod) MKL::fmod<MklType>;
/*
%template(mulbyconj) MKL::mulbyconj<MklType>;
%template(conj) MKL::conj<MklType>;
%template(arg) MKL::arg<MklType>;
%template(CIS) MKL::CIS<MklType>;
*/
%template(sinpi) MKL::sinpi<MklType>;
%template(cospi) MKL::cospi<MklType>;
%template(tanpi) MKL::tanpi<MklType>;
%template(asinpi) MKL::asinpi<MklType>;
%template(acospi) MKL::acospi<MklType>;
%template(atanpi) MKL::atanpi<MklType>;
%template(atan2pi) MKL::atan2pi<MklType>;
%template(cosd) MKL::cosd<MklType>;
%template(sind) MKL::sind<MklType>;
%template(tand) MKL::tand<MklType>;
%template(lgamma) MKL::lgamma<MklType>;
%template(tgamma) MKL::tgamma<MklType>;
%template(expint1) MKL::expint1<MklType>;
%template(copy) MKL::copy<MklType>;
%template(add) MKL::add<MklType>;
%template(sub) MKL::sub<MklType>;
%template(dot) MKL::dot<MklType>;
%template(nrm2) MKL::nrm2<MklType>;
%template(scale) MKL::scale<MklType>;
%template(min_index) MKL::min_index<MklType>;
%template(max_index) MKL::max_index<MklType>;
%template(linspace) MKL::linspace<MklType>;

%template(interleave) interleave<MklType>;
%template(dinterleave) deinterleave<MklType>;

