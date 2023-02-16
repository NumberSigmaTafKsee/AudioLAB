%module se
%{
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include "SimpleEigen/Eigen.h"
#include "SimpleEigen/functions.h"

%}

//%include "stdint.i"
typedef signed char int8_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long int int64_t;
typedef unsigned long int uint64_t;
typedef float f32_t;
typedef double f64_t;
typedef long double f80_t;

%include "std_common.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "stl.i"
%include "std_common.i"
%include "lua_fnptr.i"

%template(float_vector) std::vector<float>;
%template(float_matrix) std::vector<std::vector<float>>;
%template(double_vector) std::vector<double>;
%template(double_matrix) std::vector<std::vector<double>>;
%template(ldouble_vector) std::vector<long double>;
%template(ldouble_matrix) std::vector<std::vector<long double>>;
%template(char_vector) std::vector<char>;
%template(byte_matrix) std::vector<std::vector<unsigned char>>;
%template(short_vector) std::vector<short>;
%template(ushort_matrix) std::vector<std::vector<unsigned short>>;
%template(int_vector) std::vector<int>;
%template(uint_matrix) std::vector<std::vector<unsigned int>>;
%template(long_vector) std::vector<long>;
%template(ulong_matrix) std::vector<std::vector<unsigned long>>;
%template(ll_vector) std::vector<long long>;
%template(ull_matrix) std::vector<std::vector<unsigned long long>>;

%include "Eigen.h"
%include "Base.h"
%include "Complex.h"
%include "Scalar.h"
%include "RowVector.h"
%include "ColVector.h"
%include "Matrix.h"
%include "ColMatrix.h"
%include "Array.h"
%include "Array2D.h"
%include "SparseRowVector.h"
%include "SparseColVector.h"
%include "SparseMatrix.h"



%template(ComplexFloat)     SimpleEigen::Complex<float>;
%template(ComplexDouble)    SimpleEigen::Complex<double>;

%template(FloatScalar)  SimpleEigen::Scalar<float>;
%template(DoubleScalar) SimpleEigen::Scalar<double>;
%template(LongDoubleScalar) SimpleEigen::Scalar<long double>;
%template(Int8Scalar)  SimpleEigen::Scalar<int8_t>;
%template(UInt8Scalar) SimpleEigen::Scalar<uint8_t>;
%template(Int16Scalar)  SimpleEigen::Scalar<int16_t>;
%template(UInt16Scalar) SimpleEigen::Scalar<uint16_t>;
%template(Int32Scalar)  SimpleEigen::Scalar<int32_t>;
%template(UInt32Scalar) SimpleEigen::Scalar<uint32_t>;
%template(LongScalar)  SimpleEigen::Scalar<long>;
%template(ULongScalar) SimpleEigen::Scalar<unsigned long>;
%template(LLongScalar)  SimpleEigen::Scalar<long long>;
%template(ULLongScalar) SimpleEigen::Scalar<unsigned long long>;

%template(FloatRowVector)  SimpleEigen::RowVector<float>;
%template(DoubleRowVector) SimpleEigen::RowVector<double>;

%template(FloatColVector)  SimpleEigen::ColVector<float>;
%template(DoubleColVector) SimpleEigen::ColVector<double>;

%template(FloatMatrixView)  SimpleEigen::MatrixView<float>;
%template(DoubleMatrixView) SimpleEigen::MatrixView<double>;

%template(FloatMatrix)  SimpleEigen::Matrix<float>;
%template(DoubleMatrix) SimpleEigen::Matrix<double>;

%template(FloatColMatrixView)  SimpleEigen::ColMatrixView<float>;
%template(DoubleCOlMatrixView) SimpleEigen::ColMatrixView<double>;

%template(FloatColMatrix)  SimpleEigen::ColMatrix<float>;
%template(DoubleColMatrix) SimpleEigen::ColMatrix<double>;

%template(FloatArray)  SimpleEigen::Array<float>;
%template(DoubleArray) SimpleEigen::Array<double>;

%template(FloatArray2D)  SimpleEigen::Array2D<float>;
%template(DoubleArray2D) SimpleEigen::Array2D<double>;



%constant double Huge =  HUGE_VAL; 
%constant float Hugef =  HUGE_VALF; 
%constant double Infinity =  INFINITY;
%constant double NaN =  NAN;  

%include "functions.h"

%template(hadamard_float)  SimpleEigen::hadamard<float>;
%template(hadamard_double) SimpleEigen::hadamard<double>;

%template(sigmoid_float) SimpleEigen::sigmoid<float>;
%template(sigmoid_double) SimpleEigen::sigmoid<double>;

%template(sigmoid_deriv_float) SimpleEigen::sigmoidd<float>;
%template(sigmoid_deriv_double) SimpleEigen::sigmoidd<double>;

%template(tanh_float) SimpleEigen::tanH<float>;
%template(tanh_double) SimpleEigen::tanH<double>;

%template(tanh_deriv_float) SimpleEigen::tanHd<float>;
%template(tanh_deriv_double) SimpleEigen::tanHd<double>;


%template(relu_float) SimpleEigen::relu<float>;
%template(relu_double) SimpleEigen::relu<double>;

%template(relu_deriv_float) SimpleEigen::relud<float>;
%template(relu_deriv_double) SimpleEigen::relud<double>;


%template(softmax_float) SimpleEigen::softmax<float>;
%template(softmax_double) SimpleEigen::softmax<double>;

%template(noalias_float) SimpleEigen::noalias<float>;
%template(noalias_double) SimpleEigen::noalias<double>;

%template(mish_activate_float) SimpleEigen::mish_activate<float>;
%template(mish_activate_double) SimpleEigen::mish_activate<double>;

%template(mish_apply_jacobian_float) SimpleEigen::mish_apply_jacobian<float>;
%template(mish_apply_jacobian_double) SimpleEigen::mish_apply_jacobian<double>;


%inline %{
double SwigRefTableGet(SWIGLUA_REF  r, size_t index){
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    lua_rawgeti(r.L,-1,index);
    return lua_tonumber(r.L,-1);
}
double SwigRefTableGet2(SWIGLUA_REF  r, size_t row, size_t col){
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    lua_rawgeti(r.L,-1,row);
    assert(lua_istable(r.L,-1));
    lua_rawgeti(r.L,-1,col);
    return lua_tonumber(r.L,-1);
}
void SwigRefTableSet(SWIGLUA_REF  r, size_t index, double val){
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    lua_pushnumber(r.L,val);
    lua_rawseti(r.L,-1,index);    
}
void SwigRefTableSet2(SWIGLUA_REF  r, size_t row, size_t col, double val){
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    lua_rawgeti(r.L,-1,row);
    assert(lua_istable(r.L,-1));
    lua_rawseti(r.L,-1,val);
}
SimpleEigen::Matrix<double>& Lua_TablesToMatrix(SWIGLUA_REF r, SimpleEigen::Matrix<double> & d){
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t rows = lua_objlen(r.L,-1);
    lua_rawgeti(r.L,-1,1);
    assert(lua_istable(r.L,-1));
    size_t cols = lua_objlen(r.L,-1);
    d.resize(rows,cols);
    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < cols; j++)
            d.matrix(i,j) = SwigRefTableGet2(r,i,j);
    return d;
}

template<typename T>
SimpleEigen::RowVector<T> Lua_RowVector(SWIGLUA_REF r){
    SimpleEigen::RowVector<T> vector;
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t n = lua_objlen(r.L,-1);
    assert(n > 0);    
    vector.resize(n);
    for(size_t i = 0; i < n; i++)    {
        lua_rawgeti(r.L,-1,i+1);
        vector.vector(i) = lua_tonumber(r.L,-1);
        lua_pop(r.L,1);
    }
    return vector;
}
SimpleEigen::RowVector<float> CreateRowVectorFloat(SWIGLUA_REF r){
    return Lua_RowVector<float>(r);
}
SimpleEigen::RowVector<double> CreateRowVectorDouble(SWIGLUA_REF r){
    return Lua_RowVector<double>(r);
}

template<typename T>
SimpleEigen::ColVector<T> Lua_ColVector(SWIGLUA_REF r){
    SimpleEigen::ColVector<T> vector;
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t n = lua_objlen(r.L,-1);
    assert(n > 0);
    vector.resize(n);
    for(size_t i = 0; i < n; i++)    {
        lua_rawgeti(r.L,-1,i+1);
        vector.vector(i) = lua_tonumber(r.L,-1);
        lua_pop(r.L,1);
    }
    return vector;
}
SimpleEigen::ColVector<float> CreateColVectorFloat(SWIGLUA_REF r){
    return Lua_ColVector<float>(r);
}
SimpleEigen::ColVector<double> CreateColVectorDouble(SWIGLUA_REF r){
    return Lua_ColVector<double>(r);
}

template<typename T>
SimpleEigen::Array<T> Lua_CreateArray(SWIGLUA_REF r){
    SimpleEigen::Array<T> vector;
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t n = lua_objlen(r.L,-1);
    assert(n > 0);
    vector.resize(n);
    for(size_t i = 0; i < n; i++)    {
        lua_rawgeti(r.L,-1,i+1);
        vector.array(i) = lua_tonumber(r.L,-1);
        lua_pop(r.L,1);
    }
    return vector;
}
SimpleEigen::Array<float> CreateArrayFloat(SWIGLUA_REF r){
    return Lua_CreateArray<float>(r);
}
SimpleEigen::Array<double> CreateArrayDouble(SWIGLUA_REF r){
    return Lua_CreateArray<double>(r);
}

template<typename T>
SimpleEigen::Matrix<T> Lua_CreateMatrix(SWIGLUA_REF r){
    SimpleEigen::Matrix<T> matrix;
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t row = lua_objlen(r.L,-1);
    assert(row > 0);
    lua_rawgeti(r.L,-1,1);
    assert(lua_istable(r.L,-1));
    size_t col = lua_objlen(r.L,-1);
    matrix.resize(row,col);
    lua_pop(r.L,1);
    for(size_t i = 0; i < row; i++)    {
        lua_rawgeti(r.L,-1,i+1);
        for(size_t j = 0; j < col; j++)
        {
            lua_rawgeti(r.L,-1,j+1);
            matrix.matrix(i,j) = lua_tonumber(r.L,-1);
            lua_pop(r.L,1);
        }
        lua_pop(r.L,1);
    }
    return matrix;
}

SimpleEigen::Matrix<float> CreateMatrixFloat(SWIGLUA_REF r){
    return Lua_CreateMatrix<float>(r);
}
SimpleEigen::Matrix<double> CreateMatrixDouble(SWIGLUA_REF r){
    return Lua_CreateMatrix<double>(r);
}

template<typename T>
SimpleEigen::Array2D<T> Lua_CreateArray2D(SWIGLUA_REF r){
    SimpleEigen::Array2D<T> matrix;
    swiglua_ref_get(&r);
    assert(lua_istable(r.L,-1));
    size_t row = lua_objlen(r.L,-1);
    assert(row > 0);
    lua_rawgeti(r.L,-1,1);
    assert(lua_istable(r.L,-1));
    size_t col = lua_objlen(r.L,-1);
    matrix.resize(row,col);
    lua_pop(r.L,1);
    for(size_t i = 0; i < row; i++)    {
        lua_rawgeti(r.L,-1,i+1);
        for(size_t j = 0; j < col; j++)
        {
            lua_rawgeti(r.L,-1,j+1);
            matrix.array(i,j) = lua_tonumber(r.L,-1);
            lua_pop(r.L,1);
        }
        lua_pop(r.L,1);
    }
    return matrix;
}    

SimpleEigen::Array2D<float> CreateArray2DFloat(SWIGLUA_REF r){
    return Lua_CreateArray2D<float>(r);
}
SimpleEigen::Array2D<double> CreateArray2DDouble(SWIGLUA_REF r){
    return Lua_CreateArray2D<double>(r);
}


%}

%template(absf) SimpleEigen::abs<float>;
%template(inversef) SimpleEigen::inverse<float>;
%template(expf) SimpleEigen::exp<float>;
%template(logf) SimpleEigen::log<float>;
%template(log1pf) SimpleEigen::log1p<float>;
%template(log10f) SimpleEigen::log10<float>;
%template(powf) SimpleEigen::pow<float>;
%template(sqrtf) SimpleEigen::sqrt<float>;
%template(rsqrtf) SimpleEigen::rsqrt<float>;
%template(squaref) SimpleEigen::square<float>;
%template(cubef) SimpleEigen::cube<float>;
%template(abs2f) SimpleEigen::abs2<float>;
%template(sinf) SimpleEigen::sin<float>;
%template(cosf) SimpleEigen::cos<float>;
%template(tanf) SimpleEigen::tan<float>;
%template(asinf) SimpleEigen::asin<float>;
%template(acosf) SimpleEigen::acos<float>;
%template(atanf) SimpleEigen::atan<float>;
%template(sinhf) SimpleEigen::sinh<float>;
%template(coshf) SimpleEigen::cosh<float>;
%template(tanhf) SimpleEigen::tanh<float>;
%template(asinhf) SimpleEigen::asinh<float>;
%template(acoshf) SimpleEigen::acosh<float>;
%template(atanhf) SimpleEigen::atanh<float>;
%template(atan2f) SimpleEigen::atan2<float>;
%template(ceilf) SimpleEigen::ceil<float>;
%template(floorf) SimpleEigen::floor<float>;
%template(roundf) SimpleEigen::round<float>;
%template(rintf) SimpleEigen::rint<float>;

%template(cbrtf) SimpleEigen::cbrt<float>;
%template(copysignf) SimpleEigen::copysign<float>;
%template(erff) SimpleEigen::erf<float>;
%template(erfcf) SimpleEigen::erfc<float>;
%template(exp2f) SimpleEigen::exp2<float>;
%template(expm1f) SimpleEigen::expm1<float>;
%template(fdimf) SimpleEigen::fdim<float>;
%template(fmaf) SimpleEigen::fma<float>;
%template(fmaxf) SimpleEigen::fmax<float>;
%template(fminf) SimpleEigen::fmin<float>;
%template(fmodf) SimpleEigen::fmod<float>;
%template(fpclassifyf) SimpleEigen::fpclassify<float>;
%template(hypot) SimpleEigen::hypot<float>;
%template(ilogbf) SimpleEigen::ilogb<float>;
%template(isfinitef) SimpleEigen::isfinite<float>;
%template(isgreaterf) SimpleEigen::isgreater<float>;
%template(isgreaterequalf) SimpleEigen::isgreaterequal<float>;
%template(isinff) SimpleEigen::isinf<float>;
%template(islessf) SimpleEigen::isless<float>;
%template(islessequalf) SimpleEigen::islessequal<float>;
%template(islessgreaterf) SimpleEigen::islessgreater<float>;
%template(isnanf) SimpleEigen::isnan<float>;
%template(isnormalf) SimpleEigen::isnormal<float>;
%template(isunorderedf) SimpleEigen::isunordered<float>;
%template(ldexpf) SimpleEigen::ldexp<float>;
%template(lgammaf) SimpleEigen::lgamma<float>;
%template(llrintf) SimpleEigen::llrint<float>;
%template(llroundf) SimpleEigen::llround<float>;
%template(log2f) SimpleEigen::log2<float>;
%template(logbf) SimpleEigen::logb<float>;
%template(lrintf) SimpleEigen::lrint<float>;
%template(lroundf) SimpleEigen::lround<float>;
%template(nanf) SimpleEigen::nan<float>;
%template(nanff) SimpleEigen::nanf<float>;
%template(nanlf) SimpleEigen::nanl<float>;
%template(nearbyintf) SimpleEigen::nearbyint<float>;
%template(nextafterf) SimpleEigen::nextafter<float>;
%template(nexttowardf) SimpleEigen::nexttoward<float>;
%template(scalblnf) SimpleEigen::scalbln<float>;
%template(scalbnf) SimpleEigen::scalbn<float>;
%template(tgammaf) SimpleEigen::tgamma<float>;
%template(truncf) SimpleEigen::trunc<float>;


%inline %{
SimpleEigen::Scalar<float>   SF32(const float val) { return SimpleEigen::Scalar<float>(val); }
SimpleEigen::Scalar<double>  SF64(const float val) { return SimpleEigen::Scalar<double>(val); }
SimpleEigen::Scalar<int8_t>  SFI8(const int8_t val) { return SimpleEigen::Scalar<int8_t>(val); }
SimpleEigen::Scalar<uint8_t> SFUI8(const uint8_t val) { return SimpleEigen::Scalar<uint8_t>(val); }
SimpleEigen::Scalar<int16_t>  SFI16(const int16_t val) { return SimpleEigen::Scalar<int16_t>(val); }
SimpleEigen::Scalar<uint16_t> SFUI16(const uint16_t val) { return SimpleEigen::Scalar<uint16_t>(val); }
SimpleEigen::Scalar<int32_t>  SFI32(const int32_t val) { return SimpleEigen::Scalar<int32_t>(val); }
SimpleEigen::Scalar<uint32_t> SFUI32(const uint32_t val) { return SimpleEigen::Scalar<uint32_t>(val); }
SimpleEigen::Scalar<int64_t>  SFI64(const int64_t val) { return SimpleEigen::Scalar<int64_t>(val); }
SimpleEigen::Scalar<uint64_t> SFUI64(const uint64_t val) { return SimpleEigen::Scalar<uint64_t>(val); }
%}
