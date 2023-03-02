%module Std 
%{
extern "C" {
#include "gnuplot_i.h"
}
#include "Std/StdObject.h"
#include "Std/Std.h"   
#include "Std/datalogic.h"
#include "Std/StdCsv.h" 
#include "Std/StdJson.h"
#include "Std/StdFileSystem.h"
#include "Std/StdPosix.h"
#include "Std/StdMath.h"
#include "Std/StdPlot.h"
#include "Std/StdAssocMap.h"
#include "Std/StdCBuffer.h"
#include "Std/StdComplex.h"
#include "Std/StdScalar.h"
#include "Std/StdDeque.h"
#include "Std/StdForwardList.h"
#include "Std/StdList.h"
#include "Std/StdRandom.h"
#include "Std/StdQueue.h"
#include "Std/StdSet.h"
#include "Std/StdStack.h"
#include "Std/StdDir.h"
#include "Std/StdFile.h"
#include "Std/StdValArray.h"
#include "Std/StdVecarray.h"

//#include "SQLite++/sqlite3pp.h"

using namespace Std;
using namespace Json;
//using namespace std;
%}

%inline %{
    typedef signed char int8_t;
    typedef unsigned char uint8_t;
    typedef signed short int16_t;
    typedef unsigned short uint16_t;
    typedef signed int int32_t;
    typedef unsigned int uint32_t;
    typedef signed long long i64_t;
    typedef unsigned long long u64_t;
%}

%include "stdint.i"
%include "std_vector.i"
%include "std_string.i"
%include "lua_fnptr.i"
%include "Std/StdObject.h"
%include "Std/Std.h"   
%include "Std/datalogic.h"
%include "Std/StdCsv.h" 
%include "Std/StdJson.h"
%include "Std/StdFileSystem.h"
%include "Std/StdPosix.h"
%include "Std/StdMath.h"
%include "Std/StdPlot.h"
%include "Std/StdAssocMap.h"
%include "Std/StdCBuffer.h"
%include "Std/StdComplex.h"
%include "Std/StdScalar.h"
%include "Std/StdDeque.h"
%include "Std/StdForwardList.h"
%include "Std/StdList.h"
%include "Std/StdRandom.h"
%include "Std/StdQueue.h"
%include "Std/StdSet.h"
%include "Std/StdStack.h"
%include "Std/StdDir.h"
%include "Std/StdFile.h"
%include "Std/StdValArray.h"
%include "Std/StdVecarray.h"


%ignore Std::aligned_allocator;

%template(lua_vector)     Std::StdVector<SWIGLUA_REF>;  
%template(string_vector)  Std::StdVector<Std::StdString>;
%template(float_vector)   Std::StdVector<float>;
%template(double_vector)  Std::StdVector<double>;
%template(char_vector)    Std::StdVector<signed char>;
%template(uchar_vector)   Std::StdVector<unsigned char>;
%template(short_vector)   Std::StdVector<signed short>;
%template(ushort_vector)  Std::StdVector<unsigned short>;
%template(int_vector)     Std::StdVector<signed int>;
%template(uint_vector)    Std::StdVector<unsigned int>;
%template(long_vector)    Std::StdVector<signed long>;
%template(ulong_vector)   Std::StdVector<unsigned long>;
%template(llong_vector)   Std::StdVector<signed long long int>;
%template(ullong_vector)  Std::StdVector<unsigned long long int>;

%inline %{
    void srand() { std::srand(::time(NULL)); }
    float randomf(float min = 0.0f, float max = 1.0f) { return min + (max-min)* ((float)rand() / (float)RAND_MAX); }
    int system(const char * cmd) { return std::system(cmd); }
    void abort() { std::abort(); }
    void exit(int exit_code=-1) { std::exit(exit_code); }
    char* getenv(const char * var) { return std::getenv(var); }
    int   setenv(const char * name, const char * val, int overwrite) { return ::setenv(name,val,overwrite); }
    int unsetenv(const char * name) { return ::unsetenv(name); }
    int putenv(char * string) { return ::putenv(string); }
    int raise(int sig) { return std::raise(sig); }
   
%}