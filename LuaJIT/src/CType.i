%module CType
%{
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include "expr.h"
%}

%include "stdint.i"
%include "expr.h"

%inline %{

    template<typename T>
    struct CType
    {
        T data;

        CType()    { data = T(); }
        CType(T x) { data = x;   }

        CType(const CType& c) {
            data = c.data;
        }

        CType<T>& operator = (const CType & c) {
            data = c.data;
            return *this;
        }

        CType<T> operator + (const CType & c) {
            return CType<T>(data + c.data);
        }
        CType<T> operator - (const CType & c) {
            return CType<T>(data  - c.data);
        }
        CType<T> operator * (const CType & c) {
            return CType<T>(data * c.data);
        }
        CType<T> operator / (const CType & c) {
            return CType<T>(data / c.data);
        }

        CType<T> __unm__() { return CType<T>(-data); }
        CType<T> __pow__(const CType<T> & x) { return CType<T>(std::pow((double)data,(double)x.data)); }
        CType<T> __mod__(const CType<T> & x) { return CType<T>(std::fmod((double)data,(double)x.data)); }      
        CType<T> __pow__(const T & x) { return CType<T>(std::pow((double)data,(double)x)); }
        CType<T> __mod__(const T & x) { return CType<T>(std::fmod((double)data,(double)x)); }      

        void print() const { std::cout << data << std::endl; }        

        CType<float> cast_float() { return CType<float>(data); }
        CType<double> cast_double() { return CType<double>(data); }

        CType<int8_t> cast_int8() { return CType<int8_t>(data); }
        CType<uint8_t> cast_uint8() { return CType<uint8_t>(data); }

        CType<int16_t> cast_int16() { return CType<int16_t>(data); }
        CType<uint16_t> cast_uint16() { return CType<uint16_t>(data); }

        CType<int32_t> cast_int32() { return CType<int32_t>(data); }
        CType<uint32_t> cast_uint32() { return CType<uint32_t>(data); }

        CType<long> cast_int64() { return CType<long>(data); }
        CType<unsigned long> cast_uint64() { return CType<unsigned long>(data); }
    };

    template<typename T>
    CType<T> band(const CType<T> & a, const CType<T> & b) {
        return CType<int32_t>(a.data & b.data);
    }
    template<typename T>
    CType<T> bor(const CType<T> & a, const CType<T> & b) {
        return CType<int32_t>(a.data | b.data);    
    }
    template<typename T>
    CType<T> bxor(const CType<T> & a, const CType<T> & b) {
        return CType<int32_t>(a.data ^ b.data);
    }

    template<typename T>
    CType<T> shiftl(const CType<T> & a, const CType<T> & b) {
        return CType<int32_t>(a.data << b.data);
    }
    template<typename T>
    CType<T> shiftr(const CType<T> & a, const CType<T> & b) {
        return CType<int32_t>(a.data >> b.data);
    }

    template<typename T>
    CType<T> bnot(const CType<T> & a) {
        return CType<int32_t>(~a.data);
    }
%}

%template(float_ctype) CType<float>;
%template(double_ctype) CType<double>;

%template(int8_ctype) CType<int8_t>;
%template(uint8_ctype) CType<uint8_t>;

%template(int16_ctype) CType<int16_t>;
%template(uint16_ctype) CType<uint16_t>;

%template(int32_ctype) CType<int32_t>;
%template(uint32_ctype) CType<uint32_t>;

%template(int64_ctype) CType<long int>;
%template(uint64_ctype) CType<unsigned long int>;






