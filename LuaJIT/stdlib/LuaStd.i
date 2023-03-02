%module LuaStd
%{
#include "std.h"
#include <cctype>
#include <cassert>
#include <cstdio>
#include <cstdlib>
using namespace Std;
%}
%include "lua_fnptr.i"
//%include "stdint.i"

%inline %{
typedef signed char i8;
typedef unsigned char u8;
typedef signed short i16;
typedef unsigned short u16;
typedef signed int i32;
typedef unsigned int u32;
typedef signed long long i64;
typedef unsigned long long u64;
%}

%include "std_vector.i"
%include "std_string.i"
%include "std.h"

%template( float_vector) std::vector<float>;
%template( double_vector) std::vector<double>;
%template( i8_vector) std::vector<i8>;
%template( u8_vector) std::vector<u8>;
%template( i16_vector) std::vector<i16>;
%template( u16_vector) std::vector<u16>;
%template( i32_vector) std::vector<i32>;
%template( u32_vector) std::vector<u32>;
%template( i64_vector) std::vector<i64>;
%template( u64_vector) std::vector<u64>;

%template( StdVector ) Std::StdVector<SWIGLUA_REF>;
%template( StdMatrix ) Std::StdMatrix<SWIGLUA_REF>;

%template (FloatVector)  Std::StdVector<float>;
%template (DoubleVector) Std::StdVector<double>;
%template( I8Vector) Std::StdVector<i8>;
%template( U8Vector) Std::StdVector<u8>;
%template( I16Vector) Std::StdVector<i16>;
%template( U16Vector) Std::StdVector<u16>;
%template( I32Vector) Std::StdVector<i32>;
%template( U32Vector) Std::StdVector<u32>;
%template( I64Vector) Std::StdVector<i64>;
%template( U64Vector) Std::StdVector<u64>;

%template (FloatArray) Std::StdArray<float>;
%template (DoubleArray) Std::StdArray<double>;
%template( I8Array)  Std::StdArray<i8>;
%template( U8Array)  Std::StdArray<u8>;
%template( I16Array) Std::StdArray<i16>;
%template( U16Array) Std::StdArray<u16>;
%template( I32Array) Std::StdArray<i32>;
%template( U32Array) Std::StdArray<u32>;
%template( I64Array) Std::StdArray<i64>;
%template( U64Array) Std::StdArray<u64>;

int isalnum(int ch);
int isalpha(int ch);
int islower(int ch);
int isupper(int ch);
int isdigit(int ch);
int isxdigit(int ch);
int iscntrl(int ch);
int isgraph(int ch);
int isspace(int ch);
int isblank(int ch);
int isprint(int ch);
int ispunct(int ch);

%inline %{
    struct StdBinaryFileReader
    {
        FILE * f;

        StdBinaryFileReader() {
            f = NULL;
        }
        StdBinaryFileReader(const char * filename) {
            f = fopen(filename,"rb");
            assert(f != NULL);            
        }
        ~StdBinaryFileReader() {
            if(f) fclose(f);
        }

        void open(const char * filename) {
            if(f) fclose(f);
            f = fopen(filename,"rb");
        }
        void close() {
            if(f) fclose(f);
            f = NULL;
        }

        
        template<typename T>
        void read(T & value) {
            fread(&value,sizeof(T),1,f);
        }
        void read(StdString & s) {
            size_t len;
            fread(&len,sizeof(size_t),1,f);
            s.resize(len);
            fread(s.data(),sizeof(char),len,f);
        }
        template<typename T>
        void read(StdVector<T> & values) {
            size_t len;
            fread(&len,sizeof(size_t),1,f);
            values.resize(len);
            fread(values.data(),sizeof(T),len,f);
        }

        template<typename T>
        void read(StdMatrix<T> & values) {
            size_t M,N;
            fread(&M,sizeof(size_t),1,f);
            fread(&N,sizeof(size_t),1,f);
            values.resize(M,N);
            fread(values.data(),sizeof(T),M*N,f);
        }
    };    
%}

%template( readf )   StdBinaryFileReader::read<float>;
%template( readd )   StdBinaryFileReader::read<double>;
%template( readi8 )  StdBinaryFileReader::read<i8>;
%template( readu8 )  StdBinaryFileReader::read<u8>;
%template( readi16 ) StdBinaryFileReader::read<i16>;
%template( readu16 ) StdBinaryFileReader::read<u16>;
%template( readi32 ) StdBinaryFileReader::read<i32>;
%template( readu32 ) StdBinaryFileReader::read<u32>;
%template( readi64 ) StdBinaryFileReader::read<i64>;
%template( readu64 ) StdBinaryFileReader::read<u64>;
