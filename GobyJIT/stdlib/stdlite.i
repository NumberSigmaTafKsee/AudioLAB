// StdLite
%module stdlite
%{
#include "stdlite.h"
#include "ctype.h"
#include <cctype>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>
using namespace Std;
%}

%inline %{
typedef signed char i8;
typedef unsigned char u8;
typedef signed short i16;
typedef unsigned short u16;
typedef signed int i32;
typedef unsigned int u32;
typedef signed long long i64;
typedef unsigned long long u64;

struct Void
{
    void * ptr;
    void ** pptr;

    Void(void * p) {
        ptr = p;
        pptr = &ptr;
    }
    Void(char * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(unsigned char * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(short * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(unsigned short * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(int * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(unsigned * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(long * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(unsigned long * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(long long * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(unsigned long long * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(float * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(double * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
    Void(const char * p) {
        ptr = (void*)p;
        pptr = &ptr;
    }
};
%}

%include "stdint.i"
%include "std_common.i"
%include "std_vector.i"
%include "std_string.i"
%include "lua_fnptr.i"
%include "ctype.h"

%template (lua_matrix) std::vector<std::vector<SWIGLUA_REF>>;
%template (float_matrix) std::vector<std::vector<float>>;
%template (double_matrix) std::vector<std::vector<double>>;
%template (ldouble_matrix) std::vector<std::vector<long double>>;
%template (char_matrix) std::vector<std::vector<signed char>>;
%template (uchar_matrix) std::vector<std::vector<unsigned char>>;
%template (short_matrix) std::vector<std::vector<signed short>>;
%template (ushort_matrix) std::vector<std::vector<unsigned short>>;
%template (int_matrix) std::vector<std::vector<signed int>>;
%template (uint_matrix) std::vector<std::vector<unsigned int>>;
%template (long_matrix) std::vector<std::vector<signed long>>;
%template (ulong_matrix) std::vector<std::vector<unsigned long>>;
%template (llong_matrix) std::vector<std::vector<signed long long>>;
%template (ullong_matrix) std::vector<std::vector<unsigned long long>>;


%template(lua_vector) std::vector<SWIGLUA_REF>;
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(ldouble_vector) std::vector<long double>;
%template(char_vector) std::vector<signed char>;
%template(uchar_vector) std::vector<unsigned char>;
%template(short_vector) std::vector<signed short>;
%template(ushort_vector) std::vector<unsigned short>;
%template(int_vector) std::vector<signed int>;
%template(uint_vector) std::vector<unsigned int>;
%template(long_vector) std::vector<signed long>;
%template(ulong_vector) std::vector<unsigned long>;
%template(llong_vector) std::vector<signed long long int>;
%template(ullong_vector) std::vector<unsigned long long int>;

%include "stdlite.h"


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

%template( FloatArray ) Std::StdArray<float>;
%template( DoubleArray) Std::StdArray<double>;
%template( I8Array )  Std::StdArray<i8>;
%template( U8Array )  Std::StdArray<u8>;
%template( I16Array) Std::StdArray<i16>;
%template( U16Array) Std::StdArray<u16>;
%template( I32Array) Std::StdArray<i32>;
%template( U32Array) Std::StdArray<u32>;
%template( I64Array) Std::StdArray<i64>;
%template( U64Array) Std::StdArray<u64>;


void     clearerr(FILE *);
char    *ctermid(char *);
//int      dprintf(int, const char *, ...)
int      fclose(FILE *);
FILE    *fdopen(int, const char *);
int      feof(FILE *);
int      ferror(FILE *);
int      fflush(FILE *);
int      fgetc(FILE *);
int      fgetpos(FILE *, fpos_t *);
char    *fgets(char *, int, FILE *);
int      fileno(FILE *);
void     flockfile(FILE *);
FILE    *fmemopen(void *, size_t, const char *);
FILE    *fopen(const char *, const char *);
//int      fprintf(FILE *, const char *, ...);
int      fputc(int, FILE *);
int      fputs(const char *, FILE *);
size_t   fread(void *, size_t, size_t, FILE *);
FILE    *freopen(const char *, const char *, FILE *);
//int      fscanf(FILE *, const char *, ...);
int      fseek(FILE *, long, int);
int      fseeko(FILE *, off_t, int);
int      fsetpos(FILE *, const fpos_t *);
long     ftell(FILE *);
off_t    ftello(FILE *);
int      ftrylockfile(FILE *);
void     funlockfile(FILE *);
size_t   fwrite(const void *, size_t, size_t, FILE *);
int      getc(FILE *);
int      getchar(void);
int      getc_unlocked(FILE *);
int      getchar_unlocked(void);
ssize_t  getdelim(char **, size_t *, int, FILE *);
ssize_t  getline(char **, size_t *, FILE *);
//char    *gets(char *);
FILE    *open_memstream(char **, size_t *);
int      pclose(FILE *);
void     perror(const char *);
FILE    *popen(const char *, const char *);
//int      printf(const char *, ...);
int      putc(int, FILE *);
int      putchar(int);
int      putc_unlocked(int, FILE *);
int      putchar_unlocked(int);
int      puts(const char *);
int      remove(const char *);
int      rename(const char *, const char *);
int      renameat(int, const char *, int, const char *);
void     rewind(FILE *);
//int      scanf(const char *, ...);
void     setbuf(FILE *, char *);
int      setvbuf(FILE *, char *, int, size_t);
//int      snprintf(char *, size_t, const char *, ...);
//int      sprintf(char *, const char *, ...);
//int      sscanf(const char *, const char *, ...);
char    *tempnam(const char *, const char *);
FILE    *tmpfile(void);
char    *tmpnam(char *);
int      ungetc(int, FILE *);
/*
int      vdprintf(int, const char *, va_list);
int      vfprintf(FILE *, const char *, va_list);
int      vfscanf(FILE *, const char *, va_list);
int      vprintf(const char *, va_list);
int      vscanf(const char *, va_list);
int      vsnprintf(char *, size_t, const char *,
            va_list);
int      vsprintf(char *, const char *, va_list);
int      vsscanf(const char *, const char *, va_list);
*/

int   isalnum(int);
int   isalpha(int);
int   isascii(int);
int   isblank(int);
int   iscntrl(int);
int   isdigit(int);
int   isgraph(int);
int   islower(int);
int   isprint(int);
int   ispunct(int);
int   isspace(int);
int   isupper(int);
int   isxdigit(int);
int   toascii(int);
int   tolower(int);
int   toupper(int);
int   _toupper(int);
int   _tolower(int);

/*
%inline %{ 
    struct tm {
        int tm_sec;    
        int tm_min;    
        int tm_hour;   
        int tm_mday;   
        int tm_mon;    
        int tm_year;   
        int tm_wday;   
        int tm_yday;   
        int tm_isdst;  
    };
%}
*/
char *asctime(const struct tm *tm);
char *asctime_r(const struct tm * tm, char * buf);

char *ctime(const time_t *timep);
char *ctime_r(const time_t * timep, char * buf);

struct tm *gmtime(const time_t *timep);
struct tm *gmtime_r(const time_t * timep,
                    struct tm * result);

struct tm *localtime(const time_t *timep);
struct tm *localtime_r(const time_t * timep,
                    struct tm * result);

time_t mktime(struct tm *tm);

%constant int char_bit = CHAR_BIT;
%constant int schar_min = SCHAR_MIN;
%constant int schar_max = SCHAR_MAX;
%constant int uchar_max = UCHAR_MAX;
%constant int char_min = CHAR_MIN;
%constant int char_max = CHAR_MAX;
%constant int mb_len_max = MB_LEN_MAX;
%constant int shrt_min = SHRT_MIN;
%constant int shrt_max = SHRT_MAX;
%constant int ushrt_max = USHRT_MAX;
%constant int int_min = INT_MIN;
%constant int int_max = INT_MAX;
%constant int uint_max = UINT_MAX;
%constant int long_min = LONG_MIN;
%constant int long_max = LONG_MAX;
%constant int ulong_max = ULONG_MAX;
%constant int llong_min = LLONG_MIN;
%constant int llong_max = LLONG_MAX;
%constant int ullong_max = ULLONG_MAX;



// string.h
void    *memccpy(void *, const void *, int, size_t);
void    *memchr(const void *, int, size_t);
int      memcmp(const void *, const void *, size_t);
void    *memcpy(void *, const void *, size_t);
void    *memmove(void *, const void *, size_t);
void    *memset(void *, int, size_t);
char    *stpcpy(char *, const char *);
char    *stpncpy(char *, const char *, size_t);
char    *strcat(char *, const char *);
char    *strchr(const char *, int);
int      strcmp(const char *, const char *);
int      strcoll(const char *, const char *);
char    *strcpy(char *, const char *);
size_t   strcspn(const char *, const char *);
char    *strdup(const char *);
char    *strerror(int);
int      strerror_r(int, char *, size_t);
size_t   strlen(const char *);
char    *strncat(char *, const char *, size_t);
int      strncmp(const char *, const char *, size_t);
char    *strncpy(char *, const char *, size_t);
char    *strndup(const char *, size_t);
size_t   strnlen(const char *, size_t);
char    *strpbrk(const char *, const char *);
char    *strrchr(const char *, int);
char    *strsignal(int);
size_t   strspn(const char *, const char *);
char    *strstr(const char *, const char *);
char    *strtok(char *, const char *);
char    *strtok_r(char *, const char *, char **);
size_t   strxfrm(char *, const char *, size_t);

// stdlib.h
void          _Exit(int);
long          a64l(const char *);
void          abort(void);
int           abs(int);
int           atexit(void (*)(void));
double        atof(const char *);
int           atoi(const char *);
long          atol(const char *);
long long     atoll(const char *);
void         *bsearch(const void *, const void *, size_t, size_t, int (*)(const void *, const void *));
void         *calloc(size_t, size_t);
div_t         div(int, int);
double        drand48(void);
//char         *ecvt(double, int, int *, int *); (LEGACY )
double        erand48(unsigned short[3]);
void          exit(int);
//char         *fcvt(double, int, int *, int *); (LEGACY )
void          free(void *);
//char         *gcvt(double, int, char *); (LEGACY )
char         *getenv(const char *);
int           getsubopt(char **, char *const *, char **);
int           grantpt(int);
char         *initstate(unsigned, char *, size_t);
long          jrand48(unsigned short[3]);
char         *l64a(long);
long          labs(long);
void          lcong48(unsigned short[7]);
ldiv_t        ldiv(long, long);
long long     llabs(long long);
lldiv_t       lldiv(long long, long long);
long          lrand48(void);
void         *malloc(size_t);
int           mblen(const char *, size_t);
size_t        mbstowcs(wchar_t *, const char *, size_t);
int           mbtowc(wchar_t *, const char *, size_t);
//char         *mktemp(char *); (LEGACY )
int           mkstemp(char *);
long          mrand48(void);
long          nrand48(unsigned short[3]);
int           posix_memalign(void **, size_t, size_t);
int           posix_openpt(int);
char         *ptsname(int);
int           putenv(char *);
void          qsort(void *, size_t, size_t, int (*)(const void *,const void *));
int           rand(void);
int           rand_r(unsigned *);
long          random(void);
void         *realloc(void *, size_t);
char         *realpath(const char *, char *);
unsigned short seed48(unsigned short[3]);
int           setenv(const char *, const char *, int);
//void          setkey(const char *);
char         *setstate(const char *);
void          srand(unsigned);
void          srand48(long);
void          srandom(unsigned);
double        strtod(const char *, char **);
float         strtof(const char *, char **);
long          strtol(const char *, char **, int);
long double   strtold(const char *, char **);
long long     strtoll(const char *, char **, int);
unsigned long strtoul(const char *, char **, int);
unsigned long long strtoull(const char *, char **, int);
int           system(const char *);
int           unlockpt(int);
int           unsetenv(const char *);



%include "ctype.h"