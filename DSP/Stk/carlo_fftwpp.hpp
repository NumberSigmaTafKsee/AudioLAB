#define __FFTWPP_H_VERSION__ 2.09

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <fftw3.h>
#include <cerrno>
#include <map>
#include <iostream>
#include <sstream>
#include <climits>
#include <cstdlib>
#include <cerrno>
#include <complex>
#include <sys/time.h>
#include <cmath>

#ifndef _OPENMP
#ifndef FFTWPP_SINGLE_THREAD
#define FFTWPP_SINGLE_THREAD
#endif
#endif

#ifndef FFTWPP_SINGLE_THREAD
#include <omp.h>
#endif


inline int get_thread_num()
{
#ifdef FFTWPP_SINGLE_THREAD
  return 0;
#else
  return omp_get_thread_num();
#endif
}

inline int get_max_threads()
{
#ifdef FFTWPP_SINGLE_THREAD
  return 1;
#else
  return omp_get_max_threads();
#endif
}

#ifndef FFTWPP_SINGLE_THREAD
#define PARALLEL(code)                                  \
  if(threads > 1) {                                     \
    _Pragma("omp parallel for num_threads(threads)")    \
      code                                              \
      } else {                                          \
    code                                                \
      }
#else
#define PARALLEL(code)                          \
  {                                             \
    code                                        \
      }
#endif


using std::istream;
using std::ostream;
using std::ws;

class Complex
{
public:

  double re;
  double im;

public:

  Complex() {}
  Complex(double r, double i=0) : re(r), im(i) {}
  Complex(const Complex& y) : re(y.re), im(y.im) {}

  ~Complex() {}

  double real() const {return re;}
  double imag() const {return im;}

  const Complex& operator = (const Complex& y);

  const Complex& operator += (const Complex& y);
  const Complex& operator += (double y);
  const Complex& operator -= (const Complex& y);
  const Complex& operator -= (double y);
  const Complex& operator *= (const Complex& y);
  const Complex& operator *= (double y);
  const Complex& operator /= (const Complex& y);
  const Complex& operator /= (double y);

  void error(char* msg) const;
};

// inline members

inline const Complex& Complex::operator = (const Complex& y)

{
  re = y.re; im = y.im; return *this;
}

inline const Complex& Complex::operator += (const Complex& y)
{
  re += y.re;  im += y.im; return *this;
}

inline const Complex& Complex::operator += (double y)
{
  re += y; return *this;
}

inline const Complex& Complex::operator -= (const Complex& y)
{
  re -= y.re;  im -= y.im; return *this;
}

inline const Complex& Complex::operator -= (double y)
{
  re -= y; return *this;
}

inline const Complex& Complex::operator *= (const Complex& y)
{
  double r = re * y.re - im * y.im;
  im = re * y.im + im * y.re;
  re = r;
  return *this;
}

inline const Complex& Complex::operator *= (double y)
{
  re *= y; im *= y; return *this;
}

inline const Complex& Complex::operator /= (const Complex& y)
{
  double t1,t2,t3;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im; t3=re;
  re *= t1; re += im*t2;
  im *= t1; im -= t3*t2;
  return *this;
}

inline const Complex& Complex::operator /= (double y)
{
  re /= y;
  im /= y;
  return *this;
}

//      functions

inline int operator == (const Complex& x, const Complex& y)
{
  return x.re == y.re && x.im == y.im;
}

inline int operator == (const Complex& x, double y)
{
  return x.im == 0.0 && x.re == y;
}

inline int operator != (const Complex& x, const Complex& y)
{
  return x.re != y.re || x.im != y.im;
}

inline int operator != (const Complex& x, double y)
{
  return x.im != 0.0 || x.re != y;
}

inline Complex operator - (const Complex& x)
{
  return Complex(-x.re, -x.im);
}

inline Complex conj(const Complex& x)
{
  return Complex(x.re, -x.im);
}

inline Complex operator + (const Complex& x, const Complex& y)
{
  return Complex(x.re+y.re, x.im+y.im);
}

inline Complex operator + (const Complex& x, double y)
{
  return Complex(x.re+y, x.im);
}

inline Complex operator + (double x, const Complex& y)
{
  return Complex(x+y.re, y.im);
}

inline Complex operator - (const Complex& x, const Complex& y)
{
  return Complex(x.re-y.re, x.im-y.im);
}

inline Complex operator - (const Complex& x, double y)
{
  return Complex(x.re-y, x.im);
}

inline Complex operator - (double x, const Complex& y)
{
  return Complex(x-y.re, -y.im);
}

inline Complex operator * (const Complex& x, const Complex& y)
{
  return Complex(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re);
}

inline Complex multconj(const Complex& x, const Complex& y)
{
  return Complex(x.re*y.re+x.im*y.im,x.im*y.re-x.re*y.im);
}

inline Complex operator * (const Complex& x, double y)
{
  return Complex(x.re*y, x.im*y);
}

inline Complex operator * (double x, const Complex& y)
{
  return Complex(x*y.re, x*y.im);
}

inline Complex operator / (const Complex& x, const Complex& y)
{
  double t1,t2;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im;
  return Complex(x.im*t2+x.re*t1, x.im*t1-x.re*t2);
}

inline Complex operator / (const Complex& x, double y)
{
  return Complex(x.re/y,x.im/y);
}

inline Complex operator / (double x, const Complex& y)
{
  double factor;
  factor=1.0/(y.re*y.re+y.im*y.im);
  return Complex(x*y.re*factor,-x*y.im*factor);
}

inline double real(const Complex& x)
{
  return x.re;
}

inline double imag(const Complex& x)
{
  return x.im;
}

inline double abs2(const Complex& x)
{
  return x.re*x.re+x.im*x.im;
}

inline double abs(const Complex& x)
{
  return sqrt(abs2(x));
}

inline double arg(const Complex& x)
{
  return x.im != 0.0 ? atan2(x.im, x.re) : 0.0;
}

// Return the principal branch of the square root (non-negative real part).
inline Complex sqrt(const Complex& x)
{
  double mag=abs(x);
  if(mag == 0.0) return Complex(0.0,0.0);
  else if(x.re > 0) {
    double re=sqrt(0.5*(mag+x.re));
    return Complex(re,0.5*x.im/re);
  } else {
    double im=sqrt(0.5*(mag-x.re));
    if(x.im < 0) im=-im;
    return Complex(0.5*x.im/im,im);
  }
}

inline Complex polar(double r, double t)
{
  return Complex(r*cos(t), r*sin(t));
}

// Complex exponentiation
inline Complex pow(const Complex& z, const Complex& w)
{
  double u=w.re;
  double v=w.im;
  if(z == 0.0) return w == 0.0 ? 1.0 : 0.0;
  double logr=0.5*log(abs2(z));
  double th=arg(z);
  double phi=logr*v+th*u;
  return exp(logr*u-th*v)*Complex(cos(phi),sin(phi));
}

inline Complex pow(const Complex& z, double u)
{
  if(z == 0.0) return u == 0.0 ? 1.0 : 0.0;
  double logr=0.5*log(abs2(z));
  double theta=u*arg(z);
  return exp(logr*u)*Complex(cos(theta),sin(theta));
}

inline istream& operator >> (istream& s, Complex& y)
{
  char c;
  s >> ws >> c;
  if(c == '(') {
    s >> y.re >> c;
    if(c == ',') s >> y.im >> c;
    else y.im=0.0;
  } else {
    s.putback(c);
    s >> y.re; y.im=0.0;
  }
  return s;
}

inline ostream& operator << (ostream& s, const Complex& y)
{
  s << "(" << y.re << "," << y.im << ")";
  return s;
}

inline bool isfinite(const Complex& z)
{
#ifdef _WIN32
  return _finite(z.re) && _finite(z.im);
#else
  return !(std::isinf(z.re) || std::isnan(z.re) || std::isinf(z.im) || std::isnan(z.im));
#endif
}



#define __ARRAY_H_VERSION__ 1.55

// Defining NDEBUG improves optimization but disables argument checking.
// Defining __NOARRAY2OPT inhibits special optimization of Array2[].


#ifdef NDEBUG
#define __check(i,n,dim,m)
#define __checkSize()
#define __checkEqual(a,b,dim,m)
#define __checkActivate(i,align) this->Activate(align)
#else
#define __check(i,n,dim,m) this->Check(i,n,dim,m)
#define __checkSize() this->CheckSize()
#define __checkEqual(a,b,dim,m) this->CheckEqual(a,b,dim,m)
#define __checkActivate(i,align) this->CheckActivate(i,align)
#ifndef __NOARRAY2OPT
#define __NOARRAY2OPT
#endif
#endif

#ifndef HAVE_POSIX_MEMALIGN

#ifdef __GLIBC_PREREQ
#if __GLIBC_PREREQ(2,3)
#define HAVE_POSIX_MEMALIGN
#endif
#else
#ifdef _POSIX_SOURCE
#define HAVE_POSIX_MEMALIGN
#endif
#endif

#else

#ifdef _AIX
extern "C" int posix_memalign(void **memptr, size_t alignment, size_t size);
#endif

#endif

namespace Array {
inline std::ostream& _newl(std::ostream& s) {s << '\n'; return s;}

inline void ArrayExit(const char *x);

#ifndef __ExternalArrayExit
inline void ArrayExit(const char *x)
{
  std::cerr << _newl << "ERROR: " << x << "." << std::endl;
  exit(1);
}
#endif

#ifndef __fftwpp_h__

// Adapted from FFTW aligned malloc/free.  Assumes that malloc is at least
// sizeof(void*)-aligned. Allocated memory must be freed with free0.
inline int posix_memalign0(void **memptr, size_t alignment, size_t size)
{
  if(alignment % sizeof (void *) != 0 || (alignment & (alignment - 1)) != 0)
    return EINVAL;
  void *p0=malloc(size+alignment);
  if(!p0) return ENOMEM;
  void *p=(void *)(((size_t) p0+alignment)&~(alignment-1));
  *((void **) p-1)=p0;
  *memptr=p;
  return 0;
}

inline void free0(void *p)
{
  if(p) free(*((void **) p-1));
}

template<class T>
inline void newAlign(T *&v, size_t len, size_t align)
{
  void *mem=NULL;
  const char *invalid="Invalid alignment requested";
  const char *nomem="Memory limits exceeded";
#ifdef HAVE_POSIX_MEMALIGN
  int rc=posix_memalign(&mem,align,len*sizeof(T));
#else
  int rc=posix_memalign0(&mem,align,len*sizeof(T));
#endif
  if(rc == EINVAL) Array::ArrayExit(invalid);
  if(rc == ENOMEM) Array::ArrayExit(nomem);
  v=(T *) mem;
  for(size_t i=0; i < len; i++) new(v+i) T;
}

template<class T>
inline void deleteAlign(T *v, size_t len)
{
  for(size_t i=len-1; i > 0; i--) v[i].~T();
  v[0].~T();
#ifdef HAVE_POSIX_MEMALIGN
  free(v);
#else
  free0(v);
#endif
}

#endif

template<class T>
class array1 {
protected:
  T *v;
  unsigned int size;
  mutable int state;
public:
  enum alloc_state {unallocated=0, allocated=1, temporary=2, aligned=4};
  virtual unsigned int Size() const {return size;}
  void CheckSize() const {
    if(!test(allocated) && size == 0)
      ArrayExit("Operation attempted on unallocated array");
  }
  void CheckEqual(int a, int b, unsigned int dim, unsigned int m) const {
    if(a != b) {
      std::ostringstream buf;
      buf << "Array" << dim << " index ";
      if(m) buf << m << " ";
      buf << "is incompatible in assignment (" << a << " != " << b << ")";
      const std::string& s=buf.str();
      ArrayExit(s.c_str());
    }
  }

  int test(int flag) const {return state & flag;}
  void clear(int flag) const {state &= ~flag;}
  void set(int flag) const {state |= flag;}
  void Activate(size_t align=0) {
    if(align) {
      newAlign(v,size,align);
      set(allocated | aligned);
    } else {
      v=new T[size];
      set(allocated);
    }
  }
  void CheckActivate(int dim, size_t align=0) {
    Deallocate();
    Activate(align);
  }
  void Deallocate() const {
    if(test(allocated)) {
      if(test(aligned)) deleteAlign(v,size);
      else delete [] v;
      state=unallocated;
    }
  }
  virtual void Dimension(unsigned int nx0) {size=nx0;}
  void Dimension(unsigned int nx0, T *v0) {
    Dimension(nx0); v=v0; clear(allocated);
  }
  void Dimension(const array1<T>& A) {
    Dimension(A.size,A.v); state=A.test(temporary);
  }

  void CheckActivate(size_t align=0) {
    __checkActivate(1,align);
  }

  void Allocate(unsigned int nx0, size_t align=0) {
    Dimension(nx0);
    CheckActivate(align);
  }

  void Reallocate(unsigned int nx0, size_t align=0) {
    Deallocate();
    Allocate(nx0,align);
  }

  array1() : v(NULL), size(0), state(unallocated) {}
  array1(const void *) : size(0), state(unallocated) {}
  array1(unsigned int nx0, size_t align=0) : state(unallocated) {
    Allocate(nx0,align);
  }
  array1(unsigned int nx0, T *v0) : state(unallocated) {Dimension(nx0,v0);}
  array1(T *v0) : state(unallocated) {Dimension(INT_MAX,v0);}
  array1(const array1<T>& A) : v(A.v), size(A.size),
                               state(A.test(temporary)) {}

  virtual ~array1() {Deallocate();}

  void Freeze() {state=unallocated;}
  void Hold() {if(test(allocated)) {state=temporary;}}
  void Purge() const {if(test(temporary)) {Deallocate(); state=unallocated;}}

  virtual void Check(int i, int n, unsigned int dim, unsigned int m,
                     int o=0) const {
    if(i < 0 || i >= n) {
      std::ostringstream buf;
      buf << "Array" << dim << " index ";
      if(m) buf << m << " ";
      buf << "is out of bounds (" << i+o;
      if(n == 0) buf << " index given to empty array";
      else {
        if(i < 0) buf << " < " << o;
        else buf << " > " << n+o-1;
      }
      buf << ")";
      const std::string& s=buf.str();
      ArrayExit(s.c_str());
    }
  }

  unsigned int Nx() const {return size;}

#ifdef NDEBUG
  typedef T *opt;
#else
  typedef array1<T> opt;
#endif

  T& operator [] (int ix) const {__check(ix,size,1,1); return v[ix];}
  T& operator () (int ix) const {__check(ix,size,1,1); return v[ix];}
  T* operator () () const {return v;}
  operator T* () const {return v;}

  array1<T> operator + (int i) const {return array1<T>(size-i,v+i);}

  void Load(T a) const {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i]=a;
  }
  void Load(const T *a) const {
    for(unsigned int i=0; i < size; i++) v[i]=a[i];
  }
  void Store(T *a) const {
    for(unsigned int i=0; i < size; i++) a[i]=v[i];
  }
  void Set(T *a) {v=a; clear(allocated);}
  T Min() {
    if(size == 0)
      ArrayExit("Cannot take minimum of empty array");
    T min=v[0];
    for(unsigned int i=1; i < size; i++) if(v[i] < min) min=v[i];
    return min;
  }
  T Max() {
    if(size == 0)
      ArrayExit("Cannot take maximum of empty array");
    T max=v[0];
    for(unsigned int i=1; i < size; i++) if(v[i] > max) max=v[i];
    return max;
  }

  std::istream& Input (std::istream &s) const {
    __checkSize();
    for(unsigned int i=0; i < size; i++) s >> v[i];
    return s;
  }

  array1<T>& operator = (T a) {Load(a); return *this;}
  array1<T>& operator = (const T *a) {Load(a); return *this;}
  array1<T>& operator = (const array1<T>& A) {
    if(size != A.Size()) {
      Deallocate();
      Allocate(A.Size());
    }
    Load(A());
    A.Purge();
    return *this;
  }

  array1<T>& operator += (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array1<T>& operator -= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
  array1<T>& operator *= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= A(i);
    return *this;
  }
  array1<T>& operator /= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] /= A(i);
    return *this;
  }

  array1<T>& operator += (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += a;
    return *this;
  }
  array1<T>& operator -= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= a;
    return *this;
  }
  array1<T>& operator *= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= a;
    return *this;
  }
  array1<T>& operator /= (T a) {
    __checkSize();
    T ainv=1.0/a;
    for(unsigned int i=0; i < size; i++) v[i] *= ainv;
    return *this;
  }

  double L1() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) norm += abs(v[i]);
    return norm;
  }
#ifdef __ArrayExtensions
  double Abs2() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) norm += abs2(v[i]);
    return norm;
  }
  double L2() const {
    return sqrt(Abs2());
  }
  double LInfinity() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) {
      T a=abs(v[i]);
      if(a > norm) norm=a;
    }
    return norm;
  }
  double LMinusInfinity() const {
    __checkSize();
    double norm=DBL_MAX;
    for(unsigned int i=0; i < size; i++) {
      T a=abs(v[i]);
      if(a < norm) norm=a;
    }
    return norm;
  }
#endif
};

template<class T>
void swaparray(T& A, T& B)
{
  T C;
  C.Dimension(A);
  A.Dimension(B);
  B.Dimension(C);
}

template<class T>
void leftshiftarray(T& A, T& B, T& C)
{
  T D;
  D.Dimension(A);
  A.Dimension(B);
  B.Dimension(C);
  C.Dimension(D);
}

template<class T>
void rightshiftarray(T& A, T& B, T& C)
{
  T D;
  D.Dimension(C);
  C.Dimension(B);
  B.Dimension(A);
  A.Dimension(D);
}

template<class T>
std::ostream& operator << (std::ostream& s, const array1<T>& A)
{
  T *p=A();
  for(unsigned int i=0; i < A.Nx(); i++) {
    s << *(p++) << " ";
  }
  return s;
}

template<class T>
std::istream& operator >> (std::istream& s, const array1<T>& A)
{
  return A.Input(s);
}

template<class T>
class array2 : public array1<T> {
protected:
  unsigned int nx;
  unsigned int ny;
public:
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0) {
    nx=nx0; ny=ny0;
    this->size=nx*ny;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, T *v0) {
    Dimension(nx0,ny0);
    this->v=v0;
    this->clear(this->allocated);
  }
  void Dimension(const array1<T> &A) {ArrayExit("Operation not implemented");}

  void Allocate(unsigned int nx0, unsigned int ny0, size_t align=0) {
    Dimension(nx0,ny0);
    __checkActivate(2,align);
  }

  array2() : nx(0), ny(0) {}
  array2(unsigned int nx0, unsigned int ny0, size_t align=0) {
    Allocate(nx0,ny0,align);
  }
  array2(unsigned int nx0, unsigned int ny0, T *v0) {Dimension(nx0,ny0,v0);}

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}

#ifndef __NOARRAY2OPT
  T *operator [] (int ix) const {
    return this->v+ix*ny;
  }
#else
  array1<T> operator [] (int ix) const {
    __check(ix,nx,2,1);
    return array1<T>(ny,this->v+ix*ny);
  }
#endif
  T& operator () (int ix, int iy) const {
    __check(ix,nx,2,1);
    __check(iy,ny,2,2);
    return this->v[ix*ny+iy];
  }
  T& operator () (int i) const {
    __check(i,this->size,2,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}

  array2<T>& operator = (T a) {this->Load(a); return *this;}
  array2<T>& operator = (T *a) {this->Load(a); return *this;}
  array2<T>& operator = (const array2<T>& A) {
    __checkEqual(nx,A.Nx(),2,1);
    __checkEqual(ny,A.Ny(),2,2);
    this->Load(A());
    A.Purge();
    return *this;
  }

  array2<T>& operator += (const array2<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] += A(i);
    return *this;
  }
  array2<T>& operator -= (const array2<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] -= A(i);
    return *this;
  }
  array2<T>& operator *= (const array2<T>& A);

  array2<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] += a;
    return *this;
  }
  array2<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] -= a;
    return *this;
  }
  array2<T>& operator *= (T a) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] *= a;
    return *this;
  }

  void Identity() {
    this->Load((T) 0);
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i]=(T) 1;
  }
};

template<class T>
std::ostream& operator << (std::ostream& s, const array2<T>& A)
{
  T *p=A();
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      s << *(p++) << " ";
    }
    s << _newl;
  }
  s << std::flush;
  return s;
}

template<class T>
std::istream& operator >> (std::istream& s, const array2<T>& A)
{
  return A.Input(s);
}

template<class T>
class array3 : public array1<T> {
protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nyz;
public:
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0) {
    nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
    this->size=nx*nyz;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) {
    Dimension(nx0,ny0,nz0);
    this->v=v0;
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                size_t align=0) {
    Dimension(nx0,ny0,nz0);
    __checkActivate(3,align);
  }

  array3() : nx(0), ny(0), nz(0), nyz(0) {}
  array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         size_t align=0) {
    Allocate(nx0,ny0,nz0,align);
  }
  array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) {
    Dimension(nx0,ny0,nz0,v0);
  }

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return nz;}

  array2<T> operator [] (int ix) const {
    __check(ix,nx,3,1);
    return array2<T>(ny,nz,this->v+ix*nyz);
  }
  T& operator () (int ix, int iy, int iz) const {
    __check(ix,nx,3,1);
    __check(iy,ny,3,2);
    __check(iz,nz,3,3);
    return this->v[ix*nyz+iy*nz+iz];
  }
  T& operator () (int i) const {
    __check(i,this->size,3,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}

  array3<T>& operator = (T a) {this->Load(a); return *this;}
  array3<T>& operator = (T *a) {this->Load(a); return *this;}
  array3<T>& operator = (const array3<T>& A) {
    __checkEqual(nx,A.Nx(),3,1);
    __checkEqual(ny,A.Ny(),3,2);
    __checkEqual(nz,A.Nz(),3,3);
    this->Load(A());
    A.Purge();
    return *this;
  }

  array3<T>& operator += (array3<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] += A(i);
    return *this;
  }
  array3<T>& operator -= (array3<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] -= A(i);
    return *this;
  }

  array3<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyz+nz+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] += a;
    return *this;
  }
  array3<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyz+nz+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] -= a;
    return *this;
  }
};

template<class T>
std::ostream& operator << (std::ostream& s, const array3<T>& A)
{
  T *p=A();
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      for(unsigned int k=0; k < A.Nz(); k++) {
        s << *(p++) << " ";
      }
      s << _newl;
    }
    s << _newl;
  }
  s << std::flush;
  return s;
}

template<class T>
std::istream& operator >> (std::istream& s, const array3<T>& A)
{
  return A.Input(s);
}

template<class T>
class array4 : public array1<T> {
protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nw;
  unsigned int nyz;
  unsigned int nzw;
  unsigned int nyzw;
public:
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
    this->size=nx*nyzw;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0);
    this->v=v0;
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                unsigned int nw0, size_t align=0) {
    Dimension(nx0,ny0,nz0,nw0);
    __checkActivate(4,align);
  }

  array4() : nx(0), ny(0), nz(0), nw(0), nyz(0), nzw(0), nyzw(0) {}
  array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, size_t align=0) {Allocate(nx0,ny0,nz0,nw0,align);}
  array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,v0);
  }

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return nz;}
  unsigned int N4() const {return nw;}

  array3<T> operator [] (int ix) const {
    __check(ix,nx,3,1);
    return array3<T>(ny,nz,nw,this->v+ix*nyzw);
  }
  T& operator () (int ix, int iy, int iz, int iw) const {
    __check(ix,nx,4,1);
    __check(iy,ny,4,2);
    __check(iz,nz,4,3);
    __check(iw,nw,4,4);
    return this->v[ix*nyzw+iy*nzw+iz*nw+iw];
  }
  T& operator () (int i) const {
    __check(i,this->size,4,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}

  array4<T>& operator = (T a) {this->Load(a); return *this;}
  array4<T>& operator = (T *a) {this->Load(a); return *this;}
  array4<T>& operator = (const array4<T>& A) {
    __checkEqual(nx,A.Nx(),4,1);
    __checkEqual(ny,A.Ny(),4,2);
    __checkEqual(nz,A.Nz(),4,3);
    __checkEqual(nw,A.N4(),4,4);
    this->Load(A());
    A.Purge();
    return *this;
  }

  array4<T>& operator += (array4<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] += A(i);
    return *this;
  }
  array4<T>& operator -= (array4<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] -= A(i);
    return *this;
  }

  array4<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyzw+nzw+nw+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] += a;
    return *this;
  }
  array4<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyzw+nzw+nw+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] -= a;
    return *this;
  }
};

template<class T>
std::ostream& operator << (std::ostream& s, const array4<T>& A)
{
  T *p=A;
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      for(unsigned int k=0; k < A.Nz(); k++) {
        for(unsigned int l=0; l < A.N4(); l++) {
          s << *(p++) << " ";
        }
        s << _newl;
      }
      s << _newl;
    }
    s << _newl;
  }
  s << std::flush;
  return s;
}

template<class T>
std::istream& operator >> (std::istream& s, const array4<T>& A)
{
  return A.Input(s);
}

template<class T>
class array5 : public array1<T> {
protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nw;
  unsigned int nv;
  unsigned int nwv;
  unsigned int nzwv;
  unsigned int nyzwv;
public:
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0, unsigned int nv0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nv=nv0; nwv=nw*nv; nzwv=nz*nwv;
    nyzwv=ny*nzwv;
    this->size=nx*nyzwv;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0, unsigned int nv0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,nv0);
    this->v=v0;
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                unsigned int nw0, unsigned int nv0, size_t align=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0);
    __checkActivate(5,align);
  }

  array5() : nx(0), ny(0), nz(0), nw(0), nv(0), nwv(0), nzwv(0), nyzwv(0) {}
  array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, unsigned int nv0, size_t align=0) {
    Allocate(nx0,ny0,nz0,nw0,nv0,align);
  }
  array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, unsigned int nv0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,nv0);
  }

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return nz;}
  unsigned int N4() const {return nw;}
  unsigned int N5() const {return nv;}

  array4<T> operator [] (int ix) const {
    __check(ix,nx,4,1);
    return array4<T>(ny,nz,nw,nv,this->v+ix*nyzwv);
  }
  T& operator () (int ix, int iy, int iz, int iw, int iv) const {
    __check(ix,nx,5,1);
    __check(iy,ny,5,2);
    __check(iz,nz,5,3);
    __check(iw,nw,5,4);
    __check(iv,nv,5,5);
    return this->v[ix*nyzwv+iy*nzwv+iz*nwv+iw*nv+iv];
  }
  T& operator () (int i) const {
    __check(i,this->size,5,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}

  array5<T>& operator = (T a) {this->Load(a); return *this;}
  array5<T>& operator = (T *a) {this->Load(a); return *this;}
  array5<T>& operator = (const array5<T>& A) {
    __checkEqual(nx,A.Nx(),5,1);
    __checkEqual(ny,A.Ny(),5,2);
    __checkEqual(nz,A.Nz(),5,3);
    __checkEqual(nw,A.N4(),5,4);
    __checkEqual(nv,A.N5(),5,5);
    this->Load(A());
    A.Purge();
    return *this;
  }

  array5<T>& operator += (array5<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] += A(i);
    return *this;
  }
  array5<T>& operator -= (array5<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < this->size; i++) this->v[i] -= A(i);
    return *this;
  }

  array5<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyzwv+nzwv+nwv+nv+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] += a;
    return *this;
  }
  array5<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyzwv+nzwv+nwv+nv+1;
    for(unsigned int i=0; i < this->size; i += inc) this->v[i] -= a;
    return *this;
  }
};

template<class T>
std::ostream& operator << (std::ostream& s, const array5<T>& A)
{
  T *p=A;
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      for(unsigned int k=0; k < A.Nz(); k++) {
        for(unsigned int l=0; l < A.N4(); l++) {
          for(unsigned int l=0; l < A.N5(); l++) {
            s << *(p++) << " ";
          }
          s << _newl;
        }
        s << _newl;
      }
      s << _newl;
    }
    s << _newl;
  }
  s << std::flush;
  return s;
}

template<class T>
std::istream& operator >> (std::istream& s, const array5<T>& A)
{
  return A.Input(s);
}

#undef __check

#ifdef NDEBUG
#define __check(i,n,o,dim,m)
#else
#define __check(i,n,o,dim,m) this->Check(i-o,n,dim,m,o)
#endif

template<class T>
class Array1 : public array1<T> {
protected:
  T *voff; // Offset pointer to memory block
  int ox;
public:
  void Offsets() {
    voff=this->v-ox;
  }
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, int ox0=0) {
    this->size=nx0;
    ox=ox0;
    Offsets();
  }
  void Dimension(unsigned int nx0, T *v0, int ox0=0) {
    this->v=v0;
    Dimension(nx0,ox0);
    this->clear(this->allocated);
  }
  void Dimension(const Array1<T>& A) {
    Dimension(A.size,A.v,A.ox); this->state=A.test(this->temporary);
  }

  void Allocate(unsigned int nx0, int ox0=0, size_t align=0) {
    Dimension(nx0,ox0);
    __checkActivate(1,align);
    Offsets();
  }

  void Reallocate(unsigned int nx0, int ox0=0, size_t align=0) {
    this->Deallocate();
    Allocate(nx0,ox0,align);
  }

  Array1() : ox(0) {}
  Array1(unsigned int nx0, int ox0=0, size_t align=0) {
    Allocate(nx0,ox0,align);
  }
  Array1(unsigned int nx0, T *v0, int ox0=0) {
    Dimension(nx0,v0,ox0);
  }
  Array1(T *v0, int ox0=0) {
    Dimension(INT_MAX,v0,ox0);
  }

#ifdef NDEBUG
  typedef T *opt;
#else
  typedef Array1<T> opt;
#endif

  T& operator [] (int ix) const {__check(ix,this->size,ox,1,1); return voff[ix];}
  T& operator () (int i) const {__check(i,this->size,0,1,1); return this->v[i];}
  T* operator () () const {return this->v;}
  operator T* () const {return this->v;}

  Array1<T> operator + (int i) const {return Array1<T>(this->size-i,this->v+i,ox);}
  void Set(T *a) {this->v=a; Offsets(); this->clear(this->allocated);}

  Array1<T>& operator = (T a) {this->Load(a); return *this;}
  Array1<T>& operator = (const T *a) {this->Load(a); return *this;}
  Array1<T>& operator = (const Array1<T>& A) {
    __checkEqual(this->size,A.Size(),1,1);
    __checkEqual(ox,A.Ox(),1,1);
    this->Load(A());
    A.Purge();
    return *this;
  }
  Array1<T>& operator = (const array1<T>& A) {
    __checkEqual(this->size,A.Size(),1,1);
    __checkEqual(ox,0,1,1);
    this->Load(A());
    A.Purge();
    return *this;
  }

  int Ox() const {return ox;}
};

template<class T>
class Array2 : public array2<T> {
protected:
  T *voff,*vtemp;
  int ox,oy;
public:
  void Offsets() {
    vtemp=this->v-ox*(int) this->ny;
    voff=vtemp-oy;
  }
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
    this->nx=nx0; this->ny=ny0;
    this->size=this->nx*this->ny;
    ox=ox0; oy=oy0;
    Offsets();
  }
  void Dimension(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0,
                 int oy0=0) {
    this->v=v0;
    Dimension(nx0,ny0,ox0,oy0);
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0,
                size_t align=0) {
    Dimension(nx0,ny0,ox0,oy0);
    __checkActivate(2,align);
    Offsets();
  }

  Array2() : ox(0), oy(0) {}
  Array2(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0,
         size_t align=0) {
    Allocate(nx0,ny0,ox0,oy0,align);
  }
  Array2(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0, int oy0=0) {
    Dimension(nx0,ny0,v0,ox0,oy0);
  }

#ifndef __NOARRAY2OPT
  T *operator [] (int ix) const {
    return voff+ix*(int) this->ny;
  }
#else
  Array1<T> operator [] (int ix) const {
    __check(ix,this->nx,ox,2,1);
    return Array1<T>(this->ny,vtemp+ix*(int) this->ny,oy);
  }
#endif

  T& operator () (int ix, int iy) const {
    __check(ix,this->nx,ox,2,1);
    __check(iy,this->ny,oy,2,2);
    return voff[ix*(int) this->ny+iy];
  }
  T& operator () (int i) const {
    __check(i,this->size,0,2,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}
  void Set(T *a) {this->v=a; Offsets(); this->clear(this->allocated);}

  Array2<T>& operator = (T a) {this->Load(a); return *this;}
  Array2<T>& operator = (T *a) {this->Load(a); return *this;}
  Array2<T>& operator = (const Array2<T>& A) {
    __checkEqual(this->nx,A.Nx(),2,1);
    __checkEqual(this->ny,A.Ny(),2,2);
    __checkEqual(ox,A.Ox(),2,1);
    __checkEqual(oy,A.Oy(),2,2);
    this->Load(A());
    A.Purge();
    return *this;
  }
  Array2<T>& operator = (const array2<T>& A) {
    __checkEqual(this->nx,A.Nx(),2,1);
    __checkEqual(this->ny,A.Ny(),2,2);
    __checkEqual(ox,0,2,1);
    __checkEqual(oy,0,2,2);
    this->Load(A());
    A.Purge();
    return *this;
  }

  int Ox() const {return ox;}
  int Oy() const {return oy;}

};

template<class T>
class Array3 : public array3<T> {
protected:
  T *voff,*vtemp;
  int ox,oy,oz;
public:
  void Offsets() {
    vtemp=this->v-ox*(int) this->nyz;
    voff=vtemp-oy*(int) this->nz-oz;
  }
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 int ox0=0, int oy0=0, int oz0=0) {
    this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nyz=this->ny*this->nz;
    this->size=this->nx*this->nyz;
    ox=ox0; oy=oy0; oz=oz0;
    Offsets();
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 T *v0, int ox0=0, int oy0=0, int oz0=0) {
    this->v=v0;
    Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                int ox0=0, int oy0=0, int oz0=0, size_t align=0) {
    Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
    __checkActivate(3,align);
    Offsets();
  }

  Array3() : ox(0), oy(0), oz(0) {}
  Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         int ox0=0, int oy0=0, int oz0=0, size_t align=0) {
    Allocate(nx0,ny0,nz0,ox0,oy0,oz0,align);
  }
  Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0,
         int ox0=0, int oy0=0, int oz0=0) {
    Dimension(nx0,ny0,nz0,v0,ox0,oy0,oz0);
  }

  Array2<T> operator [] (int ix) const {
    __check(ix,this->nx,ox,3,1);
    return Array2<T>(this->ny,this->nz,vtemp+ix*(int) this->nyz,oy,oz);
  }
  T& operator () (int ix, int iy, int iz) const {
    __check(ix,this->nx,ox,3,1);
    __check(iy,this->ny,oy,3,2);
    __check(iz,this->nz,oz,3,3);
    return voff[ix*(int) this->nyz+iy*(int) this->nz+iz];
  }
  T& operator () (int i) const {
    __check(i,this->size,0,3,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}
  void Set(T *a) {this->v=a; Offsets(); this->clear(this->allocated);}

  Array3<T>& operator = (T a) {this->Load(a); return *this;}
  Array3<T>& operator = (T *a) {this->Load(a); return *this;}
  Array3<T>& operator = (const Array3<T>& A) {
    __checkEqual(this->nx,A.Nx(),3,1);
    __checkEqual(this->ny,A.Ny(),3,2);
    __checkEqual(this->nz,A.Nz(),3,3);
    __checkEqual(ox,A.Ox(),3,1);
    __checkEqual(oy,A.Oy(),3,2);
    __checkEqual(oz,A.Oz(),3,3);
    this->Load(A());
    A.Purge();
    return *this;
  }
  Array3<T>& operator = (const array3<T>& A) {
    __checkEqual(this->nx,A.Nx(),3,1);
    __checkEqual(this->ny,A.Ny(),3,2);
    __checkEqual(this->nz,A.Nz(),3,3);
    __checkEqual(ox,0,3,1);
    __checkEqual(oy,0,3,2);
    __checkEqual(oz,0,3,3);
    this->Load(A());
    A.Purge();
    return *this;
  }

  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}

};

template<class T>
class Array4 : public array4<T> {
protected:
  T *voff,*vtemp;
  int ox,oy,oz,ow;
public:
  void Offsets() {
    vtemp=this->v-ox*(int) this->nyzw;
    voff=vtemp-oy*(int) this->nzw-oz*(int) this->nw-ow;
  }
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0,
                 int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nw=nw0;
    this->nzw=this->nz*this->nw; this->nyzw=this->ny*this->nzw;
    this->size=this->nx*this->nyzw;
    ox=ox0; oy=oy0; oz=oz0; ow=ow0;
    Offsets();
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0, T *v0,
                 int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    this->v=v0;
    Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                unsigned int nw0,
                int ox0=0, int oy0=0, int oz0=0, int ow0=0, size_t align=0) {
    Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
    __checkActivate(4,align);
    Offsets();
  }

  Array4() : ox(0), oy(0), oz(0), ow(0) {}
  Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0,
         int ox0=0, int oy0=0, int oz0=0, int ow0=0, size_t align=0) {
    Allocate(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0,align);
  }
  Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, T *v0,
         int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    Dimension(nx0,ny0,nz0,nw0,v0,ox0,oy0,oz0,ow0);
  }

  Array3<T> operator [] (int ix) const {
    __check(ix,this->nx,ox,3,1);
    return Array3<T>(this->ny,this->nz,this->nw,vtemp+ix*(int) this->nyzw,
                     oy,oz,ow);
  }
  T& operator () (int ix, int iy, int iz, int iw) const {
    __check(ix,this->nx,ox,4,1);
    __check(iy,this->ny,oy,4,2);
    __check(iz,this->nz,oz,4,3);
    __check(iw,this->nw,ow,4,4);
    return voff[ix*(int) this->nyzw+iy*(int) this->nzw+iz*(int) this->nw+iw];
  }
  T& operator () (int i) const {
    __check(i,this->size,0,4,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}
  void Set(T *a) {this->v=a; Offsets(); this->clear(this->allocated);}

  Array4<T>& operator = (T a) {this->Load(a); return *this;}
  Array4<T>& operator = (T *a) {this->Load(a); return *this;}

  Array4<T>& operator = (const Array4<T>& A) {
    __checkEqual(this->nx,A.Nx(),4,1);
    __checkEqual(this->ny,A.Ny(),4,2);
    __checkEqual(this->nz,A.Nz(),4,3);
    __checkEqual(this->nw,A.N4(),4,4);
    __checkEqual(ox,A.Ox(),4,1);
    __checkEqual(oy,A.Oy(),4,2);
    __checkEqual(oz,A.Oz(),4,3);
    __checkEqual(ow,A.O4(),4,4);
    this->Load(A());
    A.Purge();
    return *this;
  }
  Array4<T>& operator = (const array4<T>& A) {
    __checkEqual(this->nx,A.Nx(),4,1);
    __checkEqual(this->ny,A.Ny(),4,2);
    __checkEqual(this->nz,A.Nz(),4,3);
    __checkEqual(this->nw,A.N4(),4,4);
    __checkEqual(this->nx,A.Nx(),4,1);
    __checkEqual(this->ny,A.Nx(),4,2);
    __checkEqual(this->nz,A.Nx(),4,3);
    __checkEqual(this->nw,A.Nx(),4,4);
    __checkEqual(ox,0,4,1);
    __checkEqual(oy,0,4,2);
    __checkEqual(oz,0,4,3);
    __checkEqual(ow,0,4,4);
    this->Load(A());
    A.Purge();
    return *this;
  }

  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}
  int O4() const {return ow;}
};

template<class T>
class Array5 : public array5<T> {
protected:
  T *voff,*vtemp;
  int ox,oy,oz,ow,ov;
public:
  void Offsets() {
    vtemp=this->v-ox*(int) this->nyzwv;
    voff=vtemp-oy*(int) this->nzwv-oz*(int) this->nwv-ow*(int) this->nv-ov;
  }
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0,  unsigned int nv0,
                 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nw=nw0; this->nv=nv0;
    this->nwv=this->nw*this->nv; this->nzwv=this->nz*this->nwv;
    this->nyzwv=this->ny*this->nzwv;
    this->size=this->nx*this->nyzwv;
    ox=ox0; oy=oy0; oz=oz0; ow=ow0; ov=ov0;
    Offsets();
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                 unsigned int nw0, unsigned int nv0, T *v0,
                 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    this->v=v0;
    Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
    this->clear(this->allocated);
  }

  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                unsigned int nw0, unsigned int nv0,
                int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0,
                size_t align=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
    __checkActivate(5,align);
    Offsets();
  }

  Array5() : ox(0), oy(0), oz(0), ow(0), ov(0) {}
  Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, unsigned int nv0, int ox0=0, int oy0=0,
         int oz0=0, int ow0=0, int ov0=0, size_t align=0) {
    Allocate(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0,align);
  }
  Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
         unsigned int nw0, unsigned int nv0, T *v0,
         int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,v0,ox0,oy0,oz0,ow0,ov0);
  }

  Array4<T> operator [] (int ix) const {
    __check(ix,this->nx,ox,4,1);
    return Array4<T>(this->ny,this->nz,this->nw,this->nv,
                     vtemp+ix*(int) this->nyzwv,oy,oz,ow,ov);
  }
  T& operator () (int ix, int iy, int iz, int iw, int iv) const {
    __check(ix,this->nx,ox,5,1);
    __check(iy,this->ny,oy,5,2);
    __check(iz,this->nz,oz,5,3);
    __check(iw,this->nw,ow,5,4);
    __check(iv,this->nv,ov,5,5);
    return voff[ix*(int) this->nyzwv+iy*(int) this->nzwv+iz*(int) this->nwv
                +iw*(int) this->nv+iv];
  }
  T& operator () (int i) const {
    __check(i,this->size,0,5,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}
  void Set(T *a) {this->v=a; Offsets(); this->clear(this->allocated);}

  Array5<T>& operator = (T a) {this->Load(a); return *this;}
  Array5<T>& operator = (T *a) {this->Load(a); return *this;}

  Array5<T>& operator = (const Array5<T>& A) {
    __checkEqual(this->nx,A.Nx(),5,1);
    __checkEqual(this->ny,A.Ny(),5,2);
    __checkEqual(this->nz,A.Nz(),5,3);
    __checkEqual(this->nw,A.N4(),5,4);
    __checkEqual(this->nv,A.N5(),5,5);
    __checkEqual(ox,A.Ox(),5,1);
    __checkEqual(oy,A.Oy(),5,2);
    __checkEqual(oz,A.Oz(),5,3);
    __checkEqual(ow,A.O4(),5,4);
    __checkEqual(ov,A.O5(),5,5);
    this->Load(A());
    A.Purge();
    return *this;
  }
  Array5<T>& operator = (const array5<T>& A) {
    __checkEqual(this->nx,A.Nx(),5,1);
    __checkEqual(this->ny,A.Ny(),5,2);
    __checkEqual(this->nz,A.Nz(),5,3);
    __checkEqual(this->nw,A.N4(),5,4);
    __checkEqual(this->nv,A.N5(),5,5);
    __checkEqual(ox,0,5,1);
    __checkEqual(oy,0,5,2);
    __checkEqual(oz,0,5,3);
    __checkEqual(ow,0,5,4);
    __checkEqual(ov,0,5,5);
    this->Load(A());
    A.Purge();
    return *this;
  }
  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}
  int O4() const {return ow;}
  int O5() const {return ov;}
};

template<class T>
inline bool Active(array1<T>& A)
{
  return A.Size();
}

template<class T>
inline bool Active(T *A)
{
  return A;
}

template<class T>
inline void Set(T *&A, T *v)
{
  A=v;
}

template<class T>
inline void Set(array1<T>& A, T *v)
{
  A.Set(v);
}

template<class T>
inline void Set(array1<T>& A, const array1<T>& B)
{
  A.Set(B());
}

template<class T>
inline void Set(Array1<T>& A, T *v)
{
  A.Set(v);
}

template<class T>
inline void Set(Array1<T>& A, const array1<T>& B)
{
  A.Set(B());
}

template<class T>
inline void Null(T *&A)
{
  A=NULL;
}

template<class T>
inline void Null(array1<T>& A)
{
  A.Dimension(0);
}

template<class T>
inline void Dimension(T *&, unsigned int)
{
}

template<class T>
inline void Dimension(array1<T> &A, unsigned int n)
{
  A.Dimension(n);
}

template<class T>
inline void Dimension(T *&A, unsigned int, T *v)
{
  A=v;
}

template<class T>
inline void Dimension(array1<T>& A, unsigned int n, T *v)
{
  A.Dimension(n,v);
}

template<class T>
inline void Dimension(Array1<T>& A, unsigned int n, T *v)
{
  A.Dimension(n,v,0);
}

template<class T>
inline void Dimension(T *&A, T *v)
{
  A=v;
}

template<class T>
inline void Dimension(array1<T>& A, const array1<T>& B)
{
  A.Dimension(B);
}

template<class T>
inline void Dimension(Array1<T>& A, const Array1<T>& B)
{
  A.Dimension(B);
}

template<class T>
inline void Dimension(Array1<T>& A, const array1<T>& B)
{
  A.Dimension(B);
}

template<class T>
inline void Dimension(array1<T>& A, unsigned int n, const array1<T>& B)
{
  A.Dimension(n,B);
}

template<class T>
inline void Dimension(Array1<T>& A, unsigned int n, const array1<T>& B, int o)
{
  A.Dimension(n,B,o);
}

template<class T>
inline void Dimension(Array1<T>& A, unsigned int n, T *v, int o)
{
  A.Dimension(n,v,o);
}

template<class T>
inline void Dimension(T *&A, unsigned int, T *v, int o)
{
  A=v-o;
}

template<class T>
inline void Allocate(T *&A, unsigned int n, size_t align=0)
{
  if(align) newAlign(A,n,align);
  else A=new T[n];
}

template<class T>
inline void Allocate(array1<T>& A, unsigned int n, size_t align=0)
{
  A.Allocate(n,align);
}

template<class T>
inline void Allocate(Array1<T>& A, unsigned int n, size_t align=0)
{
  A.Allocate(n,align);
}

template<class T>
inline void Allocate(T *&A, unsigned int n, int o, size_t align=0)
{
  Allocate(A,n,align);
  A -= o;
}

template<class T>
inline void Allocate(Array1<T>& A, unsigned int n, int o, size_t align=0)
{
  A.Allocate(n,o,align);
}

template<class T>
inline void Deallocate(T *A)
{
  if(A) delete [] A;
}

template<class T>
inline void Deallocate(array1<T>& A)
{
  A.Deallocate();
}

template<class T>
inline void Deallocate(Array1<T>& A)
{
  A.Deallocate();
}

template<class T>
inline void Deallocate(T *A, int o)
{
  if(A) delete [] (A+o);
}

template<class T>
inline void Deallocate(Array1<T>& A, int)
{
  A.Deallocate();
}

template<class T>
inline void Reallocate(T *&A, unsigned int n, size_t align=0)
{
  if(A) delete [] A;
  Allocate(A,n,align);
}

template<class T>
inline void Reallocate(array1<T>& A, unsigned int n)
{
  A.Reallocate(n);
}

template<class T>
inline void Reallocate(Array1<T>& A, unsigned int n)
{
  A.Reallocate(n);
}

template<class T>
inline void Reallocate(T *&A, unsigned int n, int o, size_t align=0)
{
  if(A) delete [] A;
  Allocate(A,n,align);
  A -= o;
}

template<class T>
inline void Reallocate(Array1<T>& A, unsigned int n, int o, size_t align=0)
{
  A.Reallocate(n,o,align);
}

template<class T>
inline void CheckReallocate(T& A, unsigned int n, unsigned int& old,
                            size_t align=0)
{
  if(n > old) {A.Reallocate(n,align); old=n;}
}

template<class T>
inline void CheckReallocate(T& A, unsigned int n, int o, unsigned int& old,
                            size_t align=0)
{
  if(n > old) {A.Reallocate(n,o,align); old=n;}
}

}

#undef __check
#undef __checkSize
#undef __checkActivate

namespace fftwpp::utils
{
    inline double totalseconds()
    {
    timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec+tv.tv_usec/1000000.0;
    }

    inline double seconds()
    {
    static double lastseconds=totalseconds();
    double t=totalseconds();
    double seconds=t-lastseconds;
    lastseconds=t;
    return seconds;
    }

    class statistics {
        unsigned int N;
        double A;
        double varL;
        double varH;
        public:
        statistics() {clear();}
        void clear() {N=0; A=varL=varH=0.0;}
        double count() {return N;}
        double mean() {return A;}
        void add(double t) {
            ++N;
            double diff=t-A;
            A += diff/N;
            double v=diff*(t-A);
            if(diff < 0.0)
            varL += v;
            else
            varH += v;
        }
        double stdev(double var, double f) {
            double factor=N > f ? f/(N-f) : 0.0;
            return sqrt(var*factor);
        }
        double stdev() {
            return stdev(varL+varH,1.0);
        }
        double stdevL() {
            return stdev(varL,2.0);
        }
        double stdevH() {
            return stdev(varH,2.0);
        }
        void output(const char *text, unsigned int m) {
            std::cout << text << ":\n"
                    << m << "\t"
                    << A << "\t"
                    << stdevL() << "\t"
                    << stdevH() << std::endl;
        }
    };
}



namespace fftwpp::utils {

    inline unsigned int ceilquotient(unsigned int a, unsigned int b)
    {
    return (a+b-1)/b;
    }

    inline Complex *ComplexAlign(size_t size)
    {
    if(size == 0) return NULL;
    Complex *v;
    Array::newAlign(v,size,sizeof(Complex));
    return v;
    }

    inline double *doubleAlign(size_t size)
    {
    double *v;
    Array::newAlign(v,size,sizeof(Complex));
    return v;
    }

    template<class T>
    inline void deleteAlign(T *p)
    {
    #ifdef HAVE_POSIX_MEMALIGN
    free(p);
    #else
    Array::free0(p);
    #endif
    }
}


namespace fftwpp {

// Obsolete names:
#define FFTWComplex ComplexAlign
#define FFTWdouble doubleAlign
#define FFTWdelete deleteAlign

// Return the memory alignment used by FFTW.
// Use of this function requires applying patches/fftw-3.3.8-alignment.patch
// to the FFTW source, recompiling, and reinstalling the FFW library.
extern "C" size_t fftw_alignment();

class fftw;

extern "C" fftw_plan Planner(fftw *F, Complex *in, Complex *out);
void LoadWisdom();
void SaveWisdom();

extern const char *inout;

struct threaddata {
  unsigned int threads;
  double mean;
  double stdev;
  threaddata() : threads(0), mean(0.0), stdev(0.0) {}
  threaddata(unsigned int threads, double mean, double stdev) :
    threads(threads), mean(mean), stdev(stdev) {}
};

class fftw;

class ThreadBase
{
protected:
  unsigned int threads;
  unsigned int innerthreads;
public:
  ThreadBase();
  ThreadBase(unsigned int threads) : threads(threads) {}
  void Threads(unsigned int nthreads) {threads=nthreads;}
  unsigned int Threads() {return threads;}

  void multithread(unsigned int nx) {
    if(nx >= threads) {
      innerthreads=1;
    } else {
      innerthreads=threads;
      threads=1;
    }
  }
};

inline unsigned int realsize(unsigned int n, Complex *in, Complex *out=NULL)
{
  return (!out || in == out) ? 2*(n/2+1) : n;
}

inline unsigned int realsize(unsigned int n, Complex *in, double *out)
{
  return realsize(n,in,(Complex *) out);
}

inline unsigned int realsize(unsigned int n, double *in, Complex *out)
{
  return realsize(n,(Complex *) in,out);
}

// Base clase for fft routines
//
class fftw : public ThreadBase {
protected:
  unsigned int doubles; // number of double precision values in dataset
  int sign;
  unsigned int threads;
  double norm;

  fftw_plan plan;
  bool inplace;

  unsigned int Dist(unsigned int n, size_t stride, size_t dist) {
    return dist ? dist : ((stride == 1) ? n : 1);
  }

  static const double twopi;

public:
  static unsigned int effort;
  static unsigned int maxthreads;
  static double testseconds;
  static const char *WisdomName;
  static fftw_plan (*planner)(fftw *f, Complex *in, Complex *out);

  virtual unsigned int Threads() {return threads;}

  static const char *oddshift;

  // In-place shift of Fourier origin to (nx/2,0) for even nx.
  static void Shift(Complex *data, unsigned int nx, unsigned int ny,
                    unsigned int threads) {
    unsigned int nyp=ny/2+1;
    unsigned int stop=nx*nyp;
    if(nx % 2 == 0) {
      unsigned int inc=2*nyp;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=nyp; i < stop; i += inc) {
        Complex *p=data+i;
        for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
      }
    } else {
      std::cerr << oddshift << std::endl;
      exit(1);
    }
  }

  // Out-of-place shift of Fourier origin to (nx/2,0) for even nx.
  static void Shift(double *data, unsigned int nx, unsigned int ny,
                    unsigned int threads) {
    if(nx % 2 == 0) {
      unsigned int stop=nx*ny;
      unsigned int inc=2*ny;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=ny; i < stop; i += inc) {
        double *p=data+i;
        for(unsigned int j=0; j < ny; j++) p[j]=-p[j];
      }
    } else {
      std::cerr << oddshift << std::endl;
      exit(1);
    }
  }

  // In-place shift of Fourier origin to (nx/2,ny/2,0) for even nx and ny.
  static void Shift(Complex *data, unsigned int nx, unsigned int ny,
                    unsigned int nz, unsigned int threads) {
    unsigned int nzp=nz/2+1;
    unsigned int nyzp=ny*nzp;
    if(nx % 2 == 0 && ny % 2 == 0) {
      unsigned int pinc=2*nzp;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; i++) {
        Complex *pstart=data+i*nyzp;
        Complex *pstop=pstart+nyzp;
        for(Complex *p=pstart+(1-(i % 2))*nzp; p < pstop; p += pinc) {
          for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
        }
      }
    } else {
      std::cerr << oddshift << " or odd ny" << std::endl;
      exit(1);
    }
  }

  // Out-of-place shift of Fourier origin to (nx/2,ny/2,0) for even nx and ny.
  static void Shift(double *data, unsigned int nx, unsigned int ny,
                    unsigned int nz, unsigned int threads) {
    unsigned int nyz=ny*nz;
    if(nx % 2 == 0 && ny % 2 == 0) {
      unsigned int pinc=2*nz;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; i++) {
        double *pstart=data+i*nyz;
        double *pstop=pstart+nyz;
        for(double *p=pstart+(1-(i % 2))*nz; p < pstop; p += pinc) {
          for(unsigned int k=0; k < nz; k++) p[k]=-p[k];
        }
      }
    } else {
      std::cerr << oddshift << " or odd ny" << std::endl;
      exit(1);
    }
  }

  fftw() : plan(NULL) {}
  fftw(unsigned int doubles, int sign, unsigned int threads,
       unsigned int n=0) :
    doubles(doubles), sign(sign), threads(threads),
    norm(1.0/(n ? n : doubles/2)), plan(NULL) {
#ifndef FFTWPP_SINGLE_THREAD
    fftw_init_threads();
#endif
  }

  virtual ~fftw() {
    if(plan) fftw_destroy_plan(plan);
  }

  virtual fftw_plan Plan(Complex *in, Complex *out) {return NULL;};

  inline void CheckAlign(Complex *p, const char *s) {
    if((size_t) p % sizeof(Complex) == 0) return;
    std::cerr << "WARNING: " << s << " array is not " << sizeof(Complex)
              << "-byte aligned: address " << p << std::endl;
  }

  void noplan() {
    std::cerr << "Unable to construct FFTW plan" << std::endl;
    exit(1);
  }

  static void planThreads(unsigned int threads) {
#ifndef FFTWPP_SINGLE_THREAD
    omp_set_num_threads(threads);
    fftw_plan_with_nthreads(threads);
#endif
  }

  threaddata time(fftw_plan plan1, fftw_plan planT, Complex *in, Complex *out,
                  unsigned int Threads) {
    utils::statistics S,ST;
    double stop=utils::totalseconds()+testseconds;
    threads=1;
    plan=plan1;
    fft(in,out);
    threads=Threads;
    plan=planT;
    fft(in,out);
    unsigned int N=1;
    unsigned int ndoubles=doubles/2;
    for(;;) {
      double t0=utils::totalseconds();
      threads=1;
      plan=plan1;
      for(unsigned int i=0; i < N; ++i) {
        for(unsigned int i=0; i < ndoubles; ++i) out[i]=i;
        fft(in,out);
      }
      double t1=utils::totalseconds();
      threads=Threads;
      plan=planT;
      for(unsigned int i=0; i < N; ++i) {
        for(unsigned int i=0; i < ndoubles; ++i) out[i]=i;
        fft(in,out);
      }
      double t=utils::totalseconds();
      S.add(t1-t0);
      ST.add(t-t1);
      if(S.mean() < 100.0/CLOCKS_PER_SEC) {
        N *= 2;
        S.clear();
        ST.clear();
      }
      if(S.count() >= 10) {
        double error=S.stdev();
        double diff=ST.mean()-S.mean();
        if(diff >= 0.0 || t > stop) {
          threads=1;
          plan=plan1;
          fftw_destroy_plan(planT);
          break;
        }
        if(diff < -error) {
          threads=Threads;
          fftw_destroy_plan(plan1);
          break;
        }
      }
    }
    return threaddata(threads,S.mean(),S.stdev());
  }

  virtual threaddata lookup(bool inplace, unsigned int threads) {
    return threaddata();
  }
  virtual void store(bool inplace, const threaddata& data) {}

  inline Complex *CheckAlign(Complex *in, Complex *out, bool constructor=true)
  {
#ifndef NO_CHECK_ALIGN
    CheckAlign(in,constructor ? "constructor input" : "input");
    if(out) CheckAlign(out,constructor ? "constructor output" : "output");
    else out=in;
#else
    if(!out) out=in;
#endif
    return out;
  }

  threaddata Setup(Complex *in, Complex *out=NULL) {
    bool alloc=!in;
    if(alloc) in=utils::ComplexAlign((doubles+1)/2);
    out=CheckAlign(in,out);
    inplace=(out==in);

    threaddata data;
    unsigned int Threads=threads;

    if(threads > 1) data=lookup(inplace,threads);
    else data=threaddata(1,0.0,0.0);

    threads=data.threads > 0 ? data.threads : 1;
    planThreads(threads);
    plan=(*planner)(this,in,out);
    if(!plan) noplan();

    fftw_plan planT;
    if(fftw::maxthreads > 1) {
      threads=Threads;
      planThreads(threads);
      planT=(*planner)(this,in,out);

      if(data.threads == 0) {
        if(planT)
          data=time(plan,planT,in,out,threads);
        else noplan();
        store(inplace,threaddata(threads,data.mean,data.stdev));
      }
    }

    if(alloc) Array::deleteAlign(in,(doubles+1)/2);
    return data;
  }

  threaddata Setup(Complex *in, double *out) {
    return Setup(in,(Complex *) out);
  }

  threaddata Setup(double *in, Complex *out=NULL) {
    return Setup((Complex *) in,out);
  }

  virtual void Execute(Complex *in, Complex *out, bool=false) {
    fftw_execute_dft(plan,(fftw_complex *) in,(fftw_complex *) out);
  }

  Complex *Setout(Complex *in, Complex *out) {
    out=CheckAlign(in,out,false);
    if(inplace ^ (out == in)) {
      std::cerr << "ERROR: fft " << inout << std::endl;
      exit(1);
    }
    return out;
  }

  void fft(Complex *in, Complex *out=NULL) {
    out=Setout(in,out);
    Execute(in,out);
  }

  void fft(double *in, Complex *out=NULL) {
    fft((Complex *) in,out);
  }

  void fft(Complex *in, double *out) {
    fft(in,(Complex *) out);
  }

  void fft0(Complex *in, Complex *out=NULL) {
    out=Setout(in,out);
    Execute(in,out,true);
  }

  void fft0(double *in, Complex *out=NULL) {
    fft0((Complex *) in,out);
  }

  void fft0(Complex *in, double *out) {
    fft0(in,(Complex *) out);
  }

  void Normalize(Complex *out) {
    unsigned int stop=doubles/2;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < stop; i++) out[i] *= norm;
  }

  void Normalize(double *out) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < doubles; i++) out[i] *= norm;
  }

  virtual void fftNormalized(Complex *in, Complex *out=NULL, bool shift=false)
  {
    out=Setout(in,out);
    Execute(in,out,shift);
    Normalize(out);
  }

  virtual void fftNormalized(Complex *in, double *out, bool shift=false) {
    out=(double *) Setout(in,(Complex *) out);
    Execute(in,(Complex *) out,shift);
    Normalize(out);
  }

  virtual void fftNormalized(double *in, Complex *out, bool shift=false) {
    fftNormalized((Complex *) in,out,shift);
  }

  template<class I, class O>
  void fft0Normalized(I in, O out) {
    fftNormalized(in,out,true);
  }

  template<class O>
  void Normalize(unsigned int nx, unsigned int M, size_t ostride,
                 size_t odist, O *out) {
    unsigned int stop=nx*ostride;
    O *outMdist=out+M*odist;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < stop; i += ostride) {
      O *pstop=outMdist+i;
      for(O *p=out+i; p < pstop; p += odist) {
        *p *= norm;
      }
    }
  }

  template<class I, class O>
  void fftNormalized(unsigned int nx, unsigned int M, size_t ostride,
                     size_t odist, I *in, O *out=NULL, bool shift=false) {
    out=(O *) Setout((Complex *) in,(Complex *) out);
    Execute((Complex *) in,(Complex *) out,shift);
    Normalize(nx,M,ostride,odist,out);
  }

}; // class fftw

class Transpose {
  fftw_plan plan;
  bool inplace;
public:
  template<class T>
  Transpose(unsigned int rows, unsigned int cols, unsigned int length,
            T *in, T *out=NULL, unsigned int threads=fftw::maxthreads) {
    unsigned int size=sizeof(T);
    if(size % sizeof(double) != 0) {
      std::cerr << "ERROR: Transpose is not implemented for type of size "
                << size;
      exit(1);
    }
    plan=NULL;
    if(!out) out=in;
    inplace=(out==in);
    if(rows == 0 || cols == 0) return;
    size /= sizeof(double);
    length *= size;

    fftw::planThreads(threads);

    fftw_iodim dims[3];

    dims[0].n=rows;
    dims[0].is=cols*length;
    dims[0].os=length;

    dims[1].n=cols;
    dims[1].is=length;
    dims[1].os=rows*length;

    dims[2].n=length;
    dims[2].is=1;
    dims[2].os=1;

    // A plan with rank=0 is a transpose.
    plan=fftw_plan_guru_r2r(0,NULL,3,dims,(double *) in,(double *) out,
                            NULL,fftw::effort);
  }

  ~Transpose() {
    if(plan) fftw_destroy_plan(plan);
  }

  template<class T>
  void transpose(T *in, T *out=NULL) {
    if(!plan) return;
    if(!out) out=in;
    if(inplace ^ (out == in)) {
      std::cerr << "ERROR: Transpose " << inout << std::endl;
      exit(1);
    }
    fftw_execute_r2r(plan,(double *) in,(double*) out);
  }
};

template<class T, class L>
class Threadtable {
public:
  typedef std::map<T,threaddata,L> Table;

  threaddata Lookup(Table& table, T key) {
    typename Table::iterator p=table.find(key);
    return p == table.end() ? threaddata() : p->second;
  }

  void Store(Table& threadtable, T key, const threaddata& data) {
    threadtable[key]=data;
  }
};

struct keytype1 {
  unsigned int nx;
  unsigned int threads;
  bool inplace;
  keytype1(unsigned int nx, unsigned int threads, bool inplace) :
    nx(nx), threads(threads), inplace(inplace) {}
};

struct keyless1 {
  bool operator()(const keytype1& a, const keytype1& b) const {
    return a.nx < b.nx || (a.nx == b.nx &&
                           (a.threads < b.threads || (a.threads == b.threads &&
                                                      a.inplace < b.inplace)));
  }
};

struct keytype2 {
  unsigned int nx;
  unsigned int ny;
  unsigned int threads;
  bool inplace;
  keytype2(unsigned int nx, unsigned int ny, unsigned int threads,
           bool inplace) :
    nx(nx), ny(ny), threads(threads), inplace(inplace) {}
};

struct keyless2 {
  bool operator()(const keytype2& a, const keytype2& b) const {
    return a.nx < b.nx || (a.nx == b.nx &&
                           (a.ny < b.ny || (a.ny == b.ny &&
                                            (a.threads < b.threads ||
                                             (a.threads == b.threads &&
                                              a.inplace < b.inplace)))));
  }
};

struct keytype3 {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int threads;
  bool inplace;
  keytype3(unsigned int nx, unsigned int ny, unsigned int nz,
           unsigned int threads, bool inplace) :
    nx(nx), ny(ny), nz(nz), threads(threads), inplace(inplace) {}
};

struct keyless3 {
  bool operator()(const keytype3& a, const keytype3& b) const {
    return a.nx < b.nx || (a.nx == b.nx &&
                           (a.ny < b.ny || (a.ny == b.ny &&
                                            (a.nz < b.nz ||
                                             (a.nz == b.nz &&
                                              (a.threads < b.threads ||
                                               (a.threads == b.threads &&
                                                a.inplace < b.inplace)))))));
  }
};

// Compute the complex Fourier transform of n complex values.
// Before calling fft(), the arrays in and out (which may coincide) must be
// allocated as Complex[n].
//
// Out-of-place usage:
//
//   fft1d Forward(n,-1,in,out);
//   Forward.fft(in,out);
//
//   fft1d Backward(n,1,in,out);
//   Backward.fft(in,out);
//
//   fft1d Backward(n,1,in,out);
//   Backward.fftNormalized(in,out); // True inverse of Forward.fft(out,in);
//
// In-place usage:
//
//   fft1d Forward(n,-1);
//   Forward.fft(in);
//
//   fft1d Backward(n,1);
//   Backward.fft(in);
//
class fft1d : public fftw, public Threadtable<keytype1,keyless1> {
  unsigned int nx;
  static Table threadtable;
public:
  fft1d(unsigned int nx, int sign, Complex *in=NULL, Complex *out=NULL,
        unsigned int threads=maxthreads)
    : fftw(2*nx,sign,threads), nx(nx) {Setup(in,out);}

#ifdef __Array_h__
  fft1d(int sign, const Array::array1<Complex>& in,
        const Array::array1<Complex>& out=Array::NULL1,
        unsigned int threads=maxthreads)
    : fftw(2*in.Nx(),sign,threads), nx(in.Nx()) {Setup(in,out);}
#endif

  threaddata lookup(bool inplace, unsigned int threads) {
    return this->Lookup(threadtable,keytype1(nx,threads,inplace));
  }
  void store(bool inplace, const threaddata& data) {
    this->Store(threadtable,keytype1(nx,data.threads,inplace),data);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_1d(nx,(fftw_complex *) in,(fftw_complex *) out,
                            sign,effort);
  }
};

template<class I, class O>
inline bool Hermitian(I in, O out) {
  return false;
}

inline bool Hermitian(double in, fftw_complex out) {
  return true;
}

inline bool Hermitian(fftw_complex in, double out) {
  return true;
}

template<class I, class O>
class fftwblock : public virtual fftw {
public:
  int nx;
  unsigned int M;
  size_t istride,ostride;
  size_t idist,odist;
  fftw_plan plan1,plan2;
  unsigned int T,Q,R;
  fftwblock(unsigned int nx, unsigned int M,
            size_t istride, size_t ostride, size_t idist, size_t odist,
            Complex *in, Complex *out, unsigned int Threads)
    : fftw(), nx(nx), M(M), istride(istride), ostride(ostride),
      idist(Dist(nx,istride,idist)), odist(Dist(nx,ostride,odist)),
      plan1(NULL), plan2(NULL) {
    T=1;
    Q=M;
    R=0;

    threaddata S1=Setup(in,out);
    fftw_plan planT1=plan;
    threads=S1.threads;
    I input;
    O output;
    bool hermitian=Hermitian(input,output);

    if(fftw::maxthreads > 1 && (!hermitian || ostride*(nx/2+1) < idist)) {
      if(Threads > 1) {
        T=std::min(M,Threads);
        Q=T > 0 ? M/T : 0;
        R=M-Q*T;

        threads=Threads;
        threaddata ST=Setup(in,out);

        if(R > 0 && threads == 1 && plan1 != plan2) {
          fftw_destroy_plan(plan2);
          plan2=plan1;
        }

        if(ST.mean > S1.mean-S1.stdev) { // Use FFTW's multi-threading
          fftw_destroy_plan(plan);
          if(R > 0) {
            fftw_destroy_plan(plan2);
            plan2=NULL;
          }
          T=1;
          Q=M;
          R=0;
          plan=planT1;
        } else {                         // Do the multi-threading ourselves
          fftw_destroy_plan(planT1);
          threads=ST.threads;
        }
      } else
        Setup(in,out); // Synchronize wisdom
    }
  }

  fftw_plan Plan(int Q, fftw_complex *in, fftw_complex *out) {
    return fftw_plan_many_dft(1,&nx,Q,in,NULL,istride,idist,
                              out,NULL,ostride,odist,sign,effort);
  }

  fftw_plan Plan(int Q, double *in, fftw_complex *out) {
    return fftw_plan_many_dft_r2c(1,&nx,Q,in,NULL,istride,idist,
                                  out,NULL,ostride,odist,effort);
  }

  fftw_plan Plan(int Q, fftw_complex *in, double *out) {
    return fftw_plan_many_dft_c2r(1,&nx,Q,in,NULL,istride,idist,
                                  out,NULL,ostride,odist,effort);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    if(R > 0) {
      plan2=Plan(Q+1,(I *) in,(O *) out);
      if(!plan2) return NULL;
      if(threads == 1) plan1=plan2;
    }
    return Plan(Q,(I *) in,(O *) out);
  }

  void Execute(fftw_plan plan, fftw_complex *in, fftw_complex *out) {
    fftw_execute_dft(plan,in,out);
  }

  void Execute(fftw_plan plan, double *in, fftw_complex *out) {
    fftw_execute_dft_r2c(plan,in,out);
  }

  void Execute(fftw_plan plan, fftw_complex *in, double *out) {
    fftw_execute_dft_c2r(plan,in,out);
  }

  void Execute(Complex *in, Complex *out, bool=false) {
    if(T == 1)
      Execute(plan,(I *) in,(O *) out);
    else {
      unsigned int extra=T-R;
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(T)
#endif
      for(unsigned int i=0; i < T; ++i) {
        unsigned int iQ=i*Q;
        if(i < extra)
          Execute(plan,(I *) in+iQ*idist,(O *) out+iQ*odist);
        else {
          unsigned int offset=iQ+i-extra;
          Execute(plan2,(I *) in+offset*idist,(O *) out+offset*odist);
        }
      }
    }
  }

  unsigned int Threads() {return std::max(T,threads);}

  ~fftwblock() {
    if(plan2) fftw_destroy_plan(plan2);
  }
};

// Compute the complex Fourier transform of M complex vectors, each of
// length n.
// Before calling fft(), the arrays in and out (which may coincide) must be
// allocated as Complex[M*n].
//
// Out-of-place usage:
//
//   mfft1d Forward(n,-1,M,stride,dist,in,out);
//   Forward.fft(in,out);
//
//   mfft1d Forward(n,-1,M,istride,ostride,idist,odist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mfft1d Forward(n,-1,M,stride,dist);
//   Forward.fft(in);
//
//
//
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors.
//
//
class mfft1d : public fftwblock<fftw_complex,fftw_complex>,
               public Threadtable<keytype3,keyless3> {
  static Table threadtable;
public:
  mfft1d(unsigned int nx, int sign, unsigned int M=1,
         Complex *in=NULL, Complex *out=NULL,
         unsigned int threads=maxthreads) :
    fftw(2*((nx-1)+(M-1)*nx+1),sign,threads,nx),
    fftwblock<fftw_complex,fftw_complex>
    (nx,M,1,1,nx,nx,in,out,threads) {}

  mfft1d(unsigned int nx, int sign, unsigned int M, size_t stride=1,
         size_t dist=0, Complex *in=NULL, Complex *out=NULL,
         unsigned int threads=maxthreads) :
    fftw(2*((nx-1)*stride+(M-1)*Dist(nx,stride,dist)+1),sign,threads,nx),
    fftwblock<fftw_complex,fftw_complex>
    (nx,M,stride,stride,dist,dist,in,out,threads) {}

  mfft1d(unsigned int nx, int sign, unsigned int M,
         size_t istride, size_t ostride, size_t idist, size_t odist,
         Complex *in, Complex *out, unsigned int threads=maxthreads):
    fftw(std::max(2*((nx-1)*istride+(M-1)*Dist(nx,istride,idist)+1),
                  2*((nx-1)*ostride+(M-1)*Dist(nx,ostride,odist)+1)),sign,
         threads, nx),
    fftwblock<fftw_complex,fftw_complex>(nx,M,istride,ostride,idist,odist,in,
                                         out,threads) {}

  threaddata lookup(bool inplace, unsigned int threads) {
    return Lookup(threadtable,keytype3(nx,Q,R,threads,inplace));
  }
  void store(bool inplace, const threaddata& data) {
    Store(threadtable,keytype3(nx,Q,R,data.threads,inplace),data);
  }
};

// Compute the complex Fourier transform of n real values, using phase sign -1.
// Before calling fft(), the array in must be allocated as double[n] and
// the array out must be allocated as Complex[n/2+1]. The arrays in and out
// may coincide, allocated as Complex[n/2+1].
//
// Out-of-place usage:
//
//   rcfft1d Forward(n,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   rcfft1d Forward(n);
//   Forward.fft(out);
//
// Notes:
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class rcfft1d : public fftw, public Threadtable<keytype1,keyless1> {
  unsigned int nx;
  static Table threadtable;
public:
  rcfft1d(unsigned int nx, Complex *out=NULL, unsigned int threads=maxthreads)
    : fftw(2*(nx/2+1),-1,threads,nx), nx(nx) {Setup(out,(double*) NULL);}

  rcfft1d(unsigned int nx, double *in, Complex *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(2*(nx/2+1),-1,threads,nx), nx(nx) {Setup(in,out);}

  threaddata lookup(bool inplace, unsigned int threads) {
    return Lookup(threadtable,keytype1(nx,threads,inplace));
  }
  void store(bool inplace, const threaddata& data) {
    Store(threadtable,keytype1(nx,data.threads,inplace),data);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_1d(nx,(double *) in,(fftw_complex *) out, effort);
  }

  void Execute(Complex *in, Complex *out, bool=false) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};

// Compute the real inverse Fourier transform of the n/2+1 Complex values
// corresponding to the non-negative part of the frequency spectrum, using
// phase sign +1.
// Before calling fft(), the array in must be allocated as Complex[n/2+1]
// and the array out must be allocated as double[n]. The arrays in and out
// may coincide, allocated as Complex[n/2+1].
//
// Out-of-place usage (input destroyed):
//
//   crfft1d Backward(n,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   crfft1d Backward(n);
//   Backward.fft(in);
//
// Notes:
//   in contains the first n/2+1 Complex Fourier values.
//   out contains the n real values stored as a Complex array;
//
class crfft1d : public fftw, public Threadtable<keytype1,keyless1> {
  unsigned int nx;
  static Table threadtable;
public:
  crfft1d(unsigned int nx, double *out=NULL, unsigned int threads=maxthreads)
    : fftw(2*(nx/2+1),1,threads,nx), nx(nx) {Setup(out);}

  crfft1d(unsigned int nx, Complex *in, double *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(realsize(nx,in,out),1,threads,nx), nx(nx) {Setup(in,out);}

  threaddata lookup(bool inplace, unsigned int threads) {
    return Lookup(threadtable,keytype1(nx,threads,inplace));
  }
  void store(bool inplace, const threaddata& data) {
    Store(threadtable,keytype1(nx,data.threads,inplace),data);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_1d(nx,(fftw_complex *) in,(double *) out,effort);
  }

  void Execute(Complex *in, Complex *out, bool=false) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};

// Compute the real Fourier transform of M real vectors, each of length n,
// using phase sign -1. Before calling fft(), the array in must be
// allocated as double[M*n] and the array out must be allocated as
// Complex[M*(n/2+1)]. The arrays in and out may coincide,
// allocated as Complex[M*(n/2+1)].
//
// Out-of-place usage:
//
//   mrcfft1d Forward(n,M,istride,ostride,idist,odist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mrcfft1d Forward(n,M,istride,ostride,idist,odist);
//   Forward.fft(out);
//
// Notes:
//   istride is the spacing between the elements of each real vector;
//   ostride is the spacing between the elements of each Complex vector;
//   idist is the spacing between the first elements of the real vectors;
//   odist is the spacing between the first elements of the Complex vectors;
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class mrcfft1d : public fftwblock<double,fftw_complex>,
                 public Threadtable<keytype3,keyless3> {
  static Table threadtable;
public:
  mrcfft1d(unsigned int nx, unsigned int M,
           size_t istride, size_t ostride,
           size_t idist, size_t odist,
           double *in=NULL, Complex *out=NULL,
           unsigned int threads=maxthreads)
    : fftw(std::max((realsize(nx,in,out)-2)*istride+(M-1)*idist+2,
                    2*(nx/2*ostride+(M-1)*odist+1)),-1,threads,nx),
      fftwblock<double,fftw_complex>
    (nx,M,istride,ostride,idist,odist,(Complex *) in,out,threads) {}

  threaddata lookup(bool inplace, unsigned int threads) {
    return Lookup(threadtable,keytype3(nx,Q,R,threads,inplace));
  }

  void store(bool inplace, const threaddata& data) {
    Store(threadtable,keytype3(nx,Q,R,data.threads,inplace),data);
  }

  void Normalize(Complex *out) {
    fftw::Normalize<Complex>(nx/2+1,M,ostride,odist,out);
  }

  void fftNormalized(double *in, Complex *out=NULL, bool shift=false) {
    fftw::fftNormalized<double,Complex>(nx/2+1,M,ostride,odist,in,out,false);
  }

  void fft0Normalized(double *in, Complex *out=NULL) {
    fftw::fftNormalized<double,Complex>(nx/2+1,M,ostride,odist,in,out,true);
  }
};

// Compute the real inverse Fourier transform of M complex vectors, each of
// length n/2+1, corresponding to the non-negative parts of the frequency
// spectra, using phase sign +1. Before calling fft(), the array in must be
// allocated as Complex[M*(n/2+1)] and the array out must be allocated as
// double[M*n]. The arrays in and out may coincide,
// allocated as Complex[M*(n/2+1)].
//
// Out-of-place usage (input destroyed):
//
//   mcrfft1d Backward(n,M,istride,ostride,idist,odist,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   mcrfft1d Backward(n,M,istride,ostride,idist,odist);
//   Backward.fft(out);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//   in contains the first n/2+1 Complex Fourier values;
//   out contains the n real values stored as a Complex array.
//
class mcrfft1d : public fftwblock<fftw_complex,double>,
                 public Threadtable<keytype3,keyless3> {
  static Table threadtable;
public:
  mcrfft1d(unsigned int nx, unsigned int M, size_t istride, size_t ostride,
           size_t idist, size_t odist, Complex *in=NULL, double *out=NULL,
           unsigned int threads=maxthreads)
    : fftw(std::max(2*(nx/2*istride+(M-1)*idist+1),
                    (realsize(nx,in,out)-2)*ostride+(M-1)*odist+2),1,threads,nx),
      fftwblock<fftw_complex,double>
    (nx,M,istride,ostride,idist,odist,in,(Complex *) out,threads) {}

  threaddata lookup(bool inplace, unsigned int threads) {
    return Lookup(threadtable,keytype3(nx,Q,R,threads,inplace));
  }

  void store(bool inplace, const threaddata& data) {
    Store(threadtable,keytype3(nx,Q,R,data.threads,inplace),data);
  }

  void Normalize(double *out) {
    fftw::Normalize<double>(nx,M,ostride,odist,out);
  }

  void fftNormalized(Complex *in, double *out=NULL, bool shift=false) {
    fftw::fftNormalized<Complex,double>(nx,M,ostride,odist,in,out,false);
  }

  void fft0Normalized(Complex *in, double *out=NULL) {
    fftw::fftNormalized<Complex,double>(nx,M,ostride,odist,in,out,true);
  }
};

// Compute the complex two-dimensional Fourier transform of nx times ny
// complex values. Before calling fft(), the arrays in and out (which may
// coincide) must be allocated as Complex[nx*ny].
//
// Out-of-place usage:
//
//   fft2d Forward(nx,ny,-1,in,out);
//   Forward.fft(in,out);
//
//   fft2d Backward(nx,ny,1,in,out);
//   Backward.fft(in,out);
//
//   fft2d Backward(nx,ny,1,in,out);
//   Backward.fftNormalized(in,out); // True inverse of Forward.fft(out,in);
//
// In-place usage:
//
//   fft2d Forward(nx,ny,-1);
//   Forward.fft(in);
//
//   fft2d Backward(nx,ny,1);
//   Backward.fft(in);
//
// Note:
//   in[ny*i+j] contains the ny Complex values for each i=0,...,nx-1.
//
class fft2d : public fftw, public Threadtable<keytype2,keyless2> {
  unsigned int nx;
  unsigned int ny;
  static Table threadtable;
public:
  fft2d(unsigned int nx, unsigned int ny, int sign, Complex *in=NULL,
        Complex *out=NULL, unsigned int threads=maxthreads)
    : fftw(2*nx*ny,sign,threads), nx(nx), ny(ny) {Setup(in,out);}

#ifdef __Array_h__
  fft2d(int sign, const Array::array2<Complex>& in,
        const Array::array2<Complex>& out=Array::NULL2,
        unsigned int threads=maxthreads)
    : fftw(2*in.Size(),sign,threads), nx(in.Nx()), ny(in.Ny()) {
    Setup(in,out);
  }
#endif

  threaddata lookup(bool inplace, unsigned int threads) {
    return this->Lookup(threadtable,keytype2(nx,ny,threads,inplace));
  }
  void store(bool inplace, const threaddata& data) {
    this->Store(threadtable,keytype2(nx,ny,data.threads,inplace),data);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_2d(nx,ny,(fftw_complex *) in,(fftw_complex *) out,
                            sign,effort);
  }

  void Execute(Complex *in, Complex *out, bool=false) {
    fftw_execute_dft(plan,(fftw_complex *) in,(fftw_complex *) out);
  }
};

// Compute the complex two-dimensional Fourier transform of nx times ny real
// values, using phase sign -1.
// Before calling fft(), the array in must be allocated as double[nx*ny] and
// the array out must be allocated as Complex[nx*(ny/2+1)]. The arrays in
// and out may coincide, allocated as Complex[nx*(ny/2+1)].
//
// Out-of-place usage:
//
//   rcfft2d Forward(nx,ny,in,out);
//   Forward.fft(in,out);       // Origin of Fourier domain at (0,0)
//   Forward.fft0(in,out);      // Origin of Fourier domain at (nx/2,0);
//                                 input destroyed.
//
// In-place usage:
//
//   rcfft2d Forward(nx,ny);
//   Forward.fft(in);           // Origin of Fourier domain at (0,0)
//   Forward.fft0(in);          // Origin of Fourier domain at (nx/2,0)
//
// Notes:
//   in contains the nx*ny real values stored as a Complex array;
//   out contains the upper-half portion (ky >= 0) of the Complex transform.
//
class rcfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:
  rcfft2d(unsigned int nx, unsigned int ny, Complex *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(2*nx*(ny/2+1),-1,threads,nx*ny), nx(nx), ny(ny) {Setup(out);}

  rcfft2d(unsigned int nx, unsigned int ny, double *in, Complex *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(2*nx*(ny/2+1),-1,threads,nx*ny), nx(nx), ny(ny) {
    Setup(in,out);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_2d(nx,ny,(double *) in,(fftw_complex *) out,
                                effort);
  }

  void Execute(Complex *in, Complex *out, bool shift=false) {
    if(shift) {
      if(inplace) Shift(in,nx,ny,threads);
      else Shift((double *) in,nx,ny,threads);
    }
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }

  // Set Nyquist modes of even shifted transforms to zero.
  void deNyquist(Complex *f) {
    unsigned int nyp=ny/2+1;
    if(nx % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int j=0; j < nyp; ++j)
        f[j]=0.0;
    if(ny % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i)
        f[(i+1)*nyp-1]=0.0;
  }
};

// Compute the real two-dimensional inverse Fourier transform of the
// nx*(ny/2+1) Complex values corresponding to the spectral values in the
// half-plane ky >= 0, using phase sign +1.
// Before calling fft(), the array in must be allocated as
// Complex[nx*(ny/2+1)] and the array out must be allocated as
// double[nx*ny]. The arrays in and out may coincide,
// allocated as Complex[nx*(ny/2+1)].
//
// Out-of-place usage (input destroyed):
//
//   crfft2d Backward(nx,ny,in,out);
//   Backward.fft(in,out);      // Origin of Fourier domain at (0,0)
//   Backward.fft0(in,out);     // Origin of Fourier domain at (nx/2,0)
//
// In-place usage:
//
//   crfft2d Backward(nx,ny);
//   Backward.fft(in);          // Origin of Fourier domain at (0,0)
//   Backward.fft0(in);         // Origin of Fourier domain at (nx/2,0)
//
// Notes:
//   in contains the upper-half portion (ky >= 0) of the Complex transform;
//   out contains the nx*ny real values stored as a Complex array.
//
class crfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:
  crfft2d(unsigned int nx, unsigned int ny, double *out=NULL,
          unsigned int threads=maxthreads) :
    fftw(2*nx*(ny/2+1),1,threads,nx*ny), nx(nx), ny(ny) {Setup(out);}

  crfft2d(unsigned int nx, unsigned int ny, Complex *in, double *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(nx*realsize(ny,in,out),1,threads,nx*ny), nx(nx), ny(ny) {
    Setup(in,out);
  }

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_2d(nx,ny,(fftw_complex *) in,(double *) out,
                                effort);
  }

  void Execute(Complex *in, Complex *out, bool shift=false) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) {
      if(inplace) Shift(out,nx,ny,threads);
      else Shift((double *) out,nx,ny,threads);
    }
  }

  // Set Nyquist modes of even shifted transforms to zero.
  void deNyquist(Complex *f) {
    unsigned int nyp=ny/2+1;
    if(nx % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int j=0; j < nyp; ++j)
        f[j]=0.0;
    if(ny % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i)
        f[(i+1)*nyp-1]=0.0;
  }
};

// Compute the complex three-dimensional Fourier transform of
// nx times ny times nz complex values. Before calling fft(), the arrays in
// and out (which may coincide) must be allocated as Complex[nx*ny*nz].
//
// Out-of-place usage:
//
//   fft3d Forward(nx,ny,nz,-1,in,out);
//   Forward.fft(in,out);
//
//   fft3d Backward(nx,ny,nz,1,in,out);
//   Backward.fft(in,out);
//
//   fft3d Backward(nx,ny,nz,1,in,out);
//   Backward.fftNormalized(in,out); // True inverse of Forward.fft(out,in);
//
// In-place usage:
//
//   fft3d Forward(nx,ny,nz,-1);
//   Forward.fft(in);
//
//   fft3d Backward(nx,ny,nz,1);
//   Backward.fft(in);
//
// Note:
//   in[nz*(ny*i+j)+k] contains the (i,j,k)th Complex value,
//   indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
//
class fft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:
  fft3d(unsigned int nx, unsigned int ny, unsigned int nz,
        int sign, Complex *in=NULL, Complex *out=NULL,
        unsigned int threads=maxthreads)
    : fftw(2*nx*ny*nz,sign,threads), nx(nx), ny(ny), nz(nz) {Setup(in,out);}

#ifdef __Array_h__
  fft3d(int sign, const Array::array3<Complex>& in,
        const Array::array3<Complex>& out=Array::NULL3,
        unsigned int threads=maxthreads)
    : fftw(2*in.Size(),sign,threads), nx(in.Nx()), ny(in.Ny()), nz(in.Nz())
  {Setup(in,out);}
#endif

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *) in,
                            (fftw_complex *) out, sign, effort);
  }
};

// Compute the complex two-dimensional Fourier transform of
// nx times ny times nz real values, using phase sign -1.
// Before calling fft(), the array in must be allocated as double[nx*ny*nz]
// and the array out must be allocated as Complex[nx*ny*(nz/2+1)]. The
// arrays in and out may coincide, allocated as Complex[nx*ny*(nz/2+1)].
//
// Out-of-place usage:
//
//   rcfft3d Forward(nx,ny,nz,in,out);
//   Forward.fft(in,out);       // Origin of Fourier domain at (0,0)
//   Forward.fft0(in,out);      // Origin of Fourier domain at (nx/2,ny/2,0);
//                                 input destroyed
// In-place usage:
//
//   rcfft3d Forward(nx,ny,nz);
//   Forward.fft(in);           // Origin of Fourier domain at (0,0)
//   Forward.fft0(in);          // Origin of Fourier domain at (nx/2,ny/2,0)
//
// Notes:
//   in contains the nx*ny*nz real values stored as a Complex array;
//   out contains the upper-half portion (kz >= 0) of the Complex transform.
//
class rcfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:
  rcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(2*nx*ny*(nz/2+1),-1,threads,nx*ny*nz), nx(nx), ny(ny), nz(nz) {
    Setup(out);
  }

  rcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, double *in,
          Complex *out=NULL, unsigned int threads=maxthreads)
    : fftw(2*nx*ny*(nz/2+1),-1,threads,nx*ny*nz),
      nx(nx), ny(ny), nz(nz) {Setup(in,out);}

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_3d(nx,ny,nz,(double *) in,(fftw_complex *) out,
                                effort);
  }

  void Execute(Complex *in, Complex *out, bool shift=false) {
    if(shift) {
      if(inplace) Shift(in,nx,ny,nz,threads);
      else Shift((double *) in,nx,ny,nz,threads);
    }
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }

  // Set Nyquist modes of even shifted transforms to zero.
  void deNyquist(Complex *f) {
    unsigned int nzp=nz/2+1;
    unsigned int yz=ny*nzp;
    if(nx % 2 == 0) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int k=0; k < yz; ++k)
        f[k]=0.0;
    }

    if(ny % 2 == 0) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i) {
        unsigned int iyz=i*yz;
        for(unsigned int k=0; k < nzp; ++k)
          f[iyz+k]=0.0;
      }
    }

    if(nz % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i)
        for(unsigned int j=0; j < ny; ++j)
          f[i*yz+(j+1)*nzp-1]=0.0;
  }
};

// Compute the real two-dimensional inverse Fourier transform of the
// nx*ny*(nz/2+1) Complex values corresponding to the spectral values in the
// half-plane kz >= 0, using phase sign +1.
// Before calling fft(), the array in must be allocated as
// Complex[nx*ny*(nz+1)/2] and the array out must be allocated as
// double[nx*ny*nz]. The arrays in and out may coincide,
// allocated as Complex[nx*ny*(nz/2+1)].
//
// Out-of-place usage (input destroyed):
//
//   crfft3d Backward(nx,ny,nz,in,out);
//   Backward.fft(in,out);      // Origin of Fourier domain at (0,0)
//   Backward.fft0(in,out);     // Origin of Fourier domain at (nx/2,ny/2,0)
//
// In-place usage:
//
//   crfft3d Backward(nx,ny,nz);
//   Backward.fft(in);          // Origin of Fourier domain at (0,0)
//   Backward.fft0(in);         // Origin of Fourier domain at (nx/2,ny/2,0)
//
// Notes:
//   in contains the upper-half portion (kz >= 0) of the Complex transform;
//   out contains the nx*ny*nz real values stored as a Complex array.
//
class crfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:
  crfft3d(unsigned int nx, unsigned int ny, unsigned int nz, double *out=NULL,
          unsigned int threads=maxthreads)
    : fftw(2*nx*ny*(nz/2+1),1,threads,nx*ny*nz), nx(nx), ny(ny), nz(nz)
  {Setup(out);}

  crfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in,
          double *out=NULL, unsigned int threads=maxthreads)
    : fftw(nx*ny*(realsize(nz,in,out)),1,threads,nx*ny*nz), nx(nx), ny(ny),
      nz(nz) {Setup(in,out);}

  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_3d(nx,ny,nz,(fftw_complex *) in,(double *) out,
                                effort);
  }

  void Execute(Complex *in, Complex *out, bool shift=false) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) {
      if(inplace) Shift(out,nx,ny,nz,threads);
      else Shift((double *) out,nx,ny,nz,threads);
    }
  }

  // Set Nyquist modes of even shifted transforms to zero.
  void deNyquist(Complex *f) {
    unsigned int nzp=nz/2+1;
    unsigned int yz=ny*nzp;
    if(nx % 2 == 0) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int k=0; k < yz; ++k)
        f[k]=0.0;
    }

    if(ny % 2 == 0) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i) {
        unsigned int iyz=i*yz;
        for(unsigned int k=0; k < nzp; ++k)
          f[iyz+k]=0.0;
      }
    }

    if(nz % 2 == 0)
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < nx; ++i)
        for(unsigned int j=0; j < ny; ++j)
          f[i*yz+(j+1)*nzp-1]=0.0;
  }
};

namespace fftwpp {

#ifdef __SSE2__

#include <emmintrin.h>

typedef __m128d Vec;

union uvec {
  unsigned u[4];
  Vec v;
};

extern const union uvec sse2_pm;
extern const union uvec sse2_mm;

#if defined(__INTEL_COMPILER) || !defined(__GNUC__)
static inline Vec operator -(const Vec& a)
{
  return _mm_xor_pd(sse2_mm.v,a);
}

static inline Vec operator +(const Vec& a, const Vec& b)
{
  return _mm_add_pd(a,b);
}

static inline Vec operator -(const Vec& a, const Vec& b)
{
  return _mm_sub_pd(a,b);
}

static inline Vec operator *(const Vec& a, const Vec& b)
{
  return _mm_mul_pd(a,b);
}

static inline void operator +=(Vec& a, const Vec& b)
{
  a=_mm_add_pd(a,b);
}

static inline void operator -=(Vec& a, const Vec& b)
{
  a=_mm_sub_pd(a,b);
}

static inline void operator *=(Vec& a, const Vec& b)
{
  a=_mm_mul_pd(a,b);
}
#endif

// Return (z.x,w.x)
static inline Vec UNPACKL(const Vec& z, const Vec& w)
{
  return _mm_unpacklo_pd(z,w);
}

// Return (z.y,w.y)
static inline Vec UNPACKH(const Vec& z, const Vec& w)
{
  return _mm_unpackhi_pd(z,w);
}

// Return (z.y,z.x)
static inline Vec FLIP(const Vec& z)
{
  return _mm_shuffle_pd(z,z,1);
}

// Return (z.x,-z.y)
static inline Vec CONJ(const Vec& z)
{
  return _mm_xor_pd(sse2_pm.v,z);
}

static inline Vec LOAD(double x)
{
  return _mm_load1_pd(&x);
}

static inline Vec SQRT(const Vec& z)
{
  return _mm_sqrt_pd(z);
}

#else

class Vec {
public:
  double x;
  double y;

  Vec() {};
  Vec(double x, double y) : x(x), y(y) {};
  Vec(const Vec &v) : x(v.x), y(v.y) {};
  Vec(const Complex &z) : x(z.re), y(z.im) {};

  const Vec& operator += (const Vec& v) {
    x += v.x;
    y += v.y;
    return *this;
  }

  const Vec& operator -= (const Vec& v) {
    x -= v.x;
    y -= v.y;
    return *this;
  }

  const Vec& operator *= (const Vec& v) {
    x *= v.x;
    y *= v.y;
    return *this;
  }
};

static inline Vec operator -(const Vec& a)
{
  return Vec(-a.x,-a.y);
}

static inline Vec operator +(const Vec& a, const Vec& b)
{
  return Vec(a.x+b.x,a.y+b.y);
}

static inline Vec operator -(const Vec& a, const Vec& b)
{
  return Vec(a.x-b.x,a.y-b.y);
}

static inline Vec operator *(const Vec& a, const Vec& b)
{
  return Vec(a.x*b.x,a.y*b.y);
}

static inline Vec UNPACKL(const Vec& z, const Vec& w)
{
  return Vec(z.x,w.x);
}

static inline Vec UNPACKH(const Vec& z, const Vec& w)
{
  return Vec(z.y,w.y);
}

static inline Vec FLIP(const Vec& z)
{
  return Vec(z.y,z.x);
}

static inline Vec CONJ(const Vec& z)
{
  return Vec(z.x,-z.y);
}

static inline Vec LOAD(double x)
{
  return Vec(x,x);
}

static inline Vec SQRT(const Vec& z)
{
  return Vec(sqrt(z.x),sqrt(z.y));
}

#endif

static inline Vec LOAD(const Complex *z)
{
  return *(const Vec *) z;
}

static inline void STORE(Complex *z, const Vec& v)
{
  *(Vec *) z = v;
}

static inline Vec LOAD(const double *z)
{
  return *(const Vec *) z;
}

static inline void STORE(double *z, const Vec& v)
{
  *(Vec *) z = v;
}

// Return I*z.
static inline Vec ZMULTI(const Vec& z)
{
  return FLIP(CONJ(z));
}

// Return the complex product of z and w.
static inline Vec ZMULT(const Vec& z, const Vec& w)
{
  return w*UNPACKL(z,z)+UNPACKH(z,z)*ZMULTI(w);
}

// Return the complex product of CONJ(z) and w.
static inline Vec ZMULTC(const Vec& z, const Vec& w)
{
  return w*UNPACKL(z,z)-UNPACKH(z,z)*ZMULTI(w);
}

// Return the complex product of z and I*w.
static inline Vec ZMULTI(const Vec& z, const Vec& w)
{
  return ZMULTI(w)*UNPACKL(z,z)-UNPACKH(z,z)*w;
}

// Return the complex product of CONJ(z) and I*w.
static inline Vec ZMULTIC(const Vec& z, const Vec& w)
{
  return ZMULTI(w)*UNPACKL(z,z)+UNPACKH(z,z)*w;
}

static inline Vec ZMULT(const Vec& x, const Vec& y, const Vec& w)
{
  return x*w+y*FLIP(w);
}

static inline Vec ZMULTI(const Vec& x, const Vec& y, const Vec& w)
{
  Vec z=CONJ(w);
  return x*FLIP(z)+y*z;
}

}

namespace utils {
extern unsigned int defaultmpithreads;

struct mpiOptions {
  int a; // Block divisor: -1=sqrt(size), 0=Tune
  int alltoall; // -1=Tune, 0=Optimized, 1=MPI, 2=Inplace
  unsigned int threads;
  unsigned int verbose;
  mpiOptions(int a=0, int alltoall=-1,
             unsigned int threads=defaultmpithreads,
             unsigned int verbose=0) :
    a(a), alltoall(alltoall), threads(threads), verbose(verbose) {}
};

}

namespace fftwpp {

extern const double sqrt3;
extern const double hsqrt3;

extern const Complex hSqrt3;
extern const Complex mhsqrt3;
extern const Complex mhalf;
extern const Complex zeta3;
extern const double twopi;

inline unsigned int min(unsigned int a, unsigned int b)
{
  return (a < b) ? a : b;
}

inline unsigned int max(unsigned int a, unsigned int b)
{
  return (a > b) ? a : b;
}

// Build the factored zeta tables.
unsigned int BuildZeta(double arg, unsigned int m,
                       Complex *&ZetaH, Complex *&ZetaL,
                       unsigned int threads=1, unsigned int s=0);

unsigned int BuildZeta(unsigned int n, unsigned int m,
                       Complex *&ZetaH, Complex *&ZetaL,
                       unsigned int threads=1, unsigned int s=0);

struct convolveOptions {
  unsigned int nx,ny,nz;           // |
  unsigned int stride2,stride3;    // | Used internally by the MPI interface.
  utils::mpiOptions mpi;           // |
  bool toplevel;

  convolveOptions(unsigned int nx, unsigned int ny, unsigned int nz,
                  unsigned int stride2, unsigned int stride3) :
    nx(nx), ny(ny), nz(nz), stride2(stride2), stride3(stride3),
    toplevel(true) {}

  convolveOptions(unsigned int nx, unsigned int ny, unsigned int stride2,
                  utils::mpiOptions mpi, bool toplevel=true) :
    nx(nx), ny(ny), stride2(stride2), mpi(mpi), toplevel(toplevel) {}

  convolveOptions(unsigned int ny, unsigned int nz,
                  unsigned int stride2, unsigned int stride3,
                  utils::mpiOptions mpi, bool toplevel=true) :
    ny(ny), nz(nz), stride2(stride2), stride3(stride3), mpi(mpi),
    toplevel(toplevel) {}

  convolveOptions(bool toplevel=true) : nx(0), ny(0), nz(0),
                                        toplevel(toplevel) {}
};

static const convolveOptions defaultconvolveOptions;

typedef void multiplier(Complex **, unsigned int m,
                        const unsigned int indexsize,
                        const unsigned int *index,
                        unsigned int r, unsigned int threads);
typedef void realmultiplier(double **, unsigned int m,
                            const unsigned int indexsize,
                            const unsigned int *index,
                            unsigned int r, unsigned int threads);

// Multipliers for binary convolutions.

multiplier multautoconvolution;
multiplier multautocorrelation;
multiplier multbinary;
multiplier multcorrelation;
multiplier multbinary2;
multiplier multbinary3;
multiplier multbinary4;
multiplier multbinary8;

realmultiplier multbinary;
realmultiplier multbinary2;
realmultiplier multadvection2;

struct general {};
struct pretransform1 {};
struct pretransform2 {};
struct pretransform3 {};
struct pretransform4 {};

// In-place implicitly dealiased 1D complex convolution using
// function pointers for multiplication
class ImplicitConvolution : public ThreadBase {
private:
  unsigned int m;
  Complex **U;
  unsigned int A;
  unsigned int B;
  Complex *u;
  unsigned int s;
  Complex *ZetaH, *ZetaL;
  fft1d *BackwardsO,*ForwardsO;
  fft1d *Backwards,*Forwards;
  bool pointers;
  bool allocated;
  unsigned int indexsize;
public:
  unsigned int *index;

  void initpointers(Complex **&U, Complex *u) {
    unsigned int C=max(A,B);
    U=new Complex *[C];
    for(unsigned int a=0; a < C; ++a)
      U[a]=u+a*m;
    pointers=true;
  }

  void deletepointers(Complex **&U) {
    delete [] U;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
  }

  void init() {
    indexsize=0;

    Complex* U0=U[0];
    Complex* U1=A == 1 ? utils::ComplexAlign(m) : U[1];

    BackwardsO=new fft1d(m,1,U0,U1);
    ForwardsO=new fft1d(m,-1,U0,U1);
    threads=std::min(threads,max(BackwardsO->Threads(),ForwardsO->Threads()));

    if(A == B) {
      Backwards=new fft1d(m,1,U0);
      threads=std::min(threads,Backwards->Threads());
    }
    if(A <= B) {
      Forwards=new fft1d(m,-1,U0);
      threads=std::min(threads,Forwards->Threads());
    }

    if(A == 1) utils::deleteAlign(U1);

    s=BuildZeta(2*m,m,ZetaH,ZetaL,threads);
  }

  // m is the number of Complex data values.
  // U is an array of C distinct work arrays each of size m, where C=max(A,B)
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitConvolution(unsigned int m, Complex **U, unsigned int A=2,
                      unsigned int B=1,
                      unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), U(U), A(A), B(B), pointers(false),
      allocated(false) {
    init();
  }

  // m is the number of Complex data values.
  // u is a work array of C*m Complex values.
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitConvolution(unsigned int m, Complex *u,
                      unsigned int A=2, unsigned int B=1,
                      unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), A(A), B(B), u(u), allocated(false) {
    initpointers(U,u);
    init();
  }

  // m is the number of Complex data values.
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitConvolution(unsigned int m,
                      unsigned int A=2, unsigned int B=1,
                      unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), A(A), B(B), allocated(true) {
    u=utils::ComplexAlign(max(A,B)*m);
    initpointers(U,u);
    init();
  }

  ~ImplicitConvolution() {
    utils::deleteAlign(ZetaH);
    utils::deleteAlign(ZetaL);

    if(pointers) deletepointers(U);
    if(allocated) utils::deleteAlign(u);

    if(A == B)
      delete Backwards;
    if(A <= B)
      delete Forwards;

    delete ForwardsO;
    delete BackwardsO;
  }

  // F is an array of C pointers to distinct data blocks each of
  // size m, shifted by offset (contents not preserved).
  void convolve(Complex **F, multiplier *pmult, unsigned int i=0,
                unsigned int offset=0);

  void autoconvolve(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautoconvolution);
  }

  void autocorrelate(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautocorrelation);
  }

  // Binary convolution:
  void convolve(Complex *f, Complex *g) {
    Complex *F[]={f,g};
    convolve(F,multbinary);
  }

  // Binary correlation:
  void correlate(Complex *f, Complex *g) {
    Complex *F[]={f,g};
    convolve(F,multcorrelation);
  }

  template<class T>
  inline void pretransform(Complex **F, unsigned int k, Vec& Zetak);

  template<class T>
  void pretransform(Complex **F);

  void posttransform(Complex *f, Complex *u);
};

// In-place implicitly dealiased 1D Hermitian convolution.
class ImplicitHConvolution : public ThreadBase {
protected:
  unsigned int m;
  unsigned int c;
  bool compact;
  Complex **U;
  unsigned int A;
  unsigned int B;
  Complex *u;
  unsigned int s;
  Complex *ZetaH,*ZetaL;
  rcfft1d *rc,*rco,*rcO;
  crfft1d *cr,*cro,*crO;
  Complex *w; // Work array of size max(A,B) to hold f[c] in even case.
  bool pointers;
  bool allocated;
  bool even;
  unsigned int indexsize;
public:
  unsigned int *index;

  void initpointers(Complex **&U, Complex *u) {
    unsigned int C=max(A,B);
    U=new Complex *[C];
    unsigned stride=c+1;
    for(unsigned int a=0; a < C; ++a)
      U[a]=u+a*stride;
    pointers=true;
  }

  void deletepointers(Complex **&U) {
    delete [] U;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
  }

  void init() {
    even=m == 2*c;
    indexsize=0;
    Complex* U0=U[0];

    rc=new rcfft1d(m,U0);
    cr=new crfft1d(m,U0);

    Complex* U1=A == 1 ? utils::ComplexAlign(m) : U[1];
    rco=new rcfft1d(m,(double *) U0,U1);
    cro=new crfft1d(m,U1,(double *) U0);
    if(A == 1) utils::deleteAlign(U1);

    if(A != B) {
      rcO=rco;
      crO=cro;
    } else {
      rcO=rc;
      crO=cr;
    }

    threads=std::min(threads,std::max(rco->Threads(),cro->Threads()));
    s=BuildZeta(3*m,c+2,ZetaH,ZetaL,threads);
    w=even ? utils::ComplexAlign(max(A,B)) : u;
  }

  // m is the number of independent data values
  // U is an array of max(A,B) distinct work arrays of size c+1, where c=m/2
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitHConvolution(unsigned int m, Complex **U, unsigned int A=2,
                       unsigned int B=1, unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), c(m/2), compact(true), U(U), A(A), B(B),
      pointers(false), allocated(false) {
    init();
  }

  ImplicitHConvolution(unsigned int m, bool compact, Complex **U,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), c(m/2), compact(compact), U(U), A(A), B(B),
      pointers(false), allocated(false) {
    init();
  }

  // m is the number of independent data values
  // u is a work array of max(A,B)*(c+1) Complex values, where c=m/2
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitHConvolution(unsigned int m, Complex *u,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), c(m/2), compact(true), A(A), B(B), u(u),
      allocated(false) {
    initpointers(U,u);
    init();
  }

  ImplicitHConvolution(unsigned int m, bool compact, Complex *u,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), c(m/2), compact(compact), A(A), B(B), u(u),
      allocated(false) {
    initpointers(U,u);
    init();
  }

  // m is the number of independent data values
  // u is a work array of max(A,B)*(c+1) Complex values, where c=m/2
  // A is the number of inputs.
  // B is the number of outputs.
  ImplicitHConvolution(unsigned int m, bool compact=true, unsigned int A=2,
                       unsigned int B=1, unsigned int threads=fftw::maxthreads)
    : ThreadBase(threads), m(m), c(m/2), compact(compact), A(A), B(B),
      u(utils::ComplexAlign(max(A,B)*(c+1))), allocated(true) {
    initpointers(U,u);
    init();
  }

  virtual ~ImplicitHConvolution() {
    if(even) utils::deleteAlign(w);
    utils::deleteAlign(ZetaH);
    utils::deleteAlign(ZetaL);

    if(pointers) deletepointers(U);
    if(allocated) utils::deleteAlign(u);

    if(A != B) {
      delete cro;
      delete rco;
    }

    delete cr;
    delete rc;
  }

  // F is an array of A pointers to distinct data blocks each of size m,
  // shifted by offset (contents not preserved).
  void convolve(Complex **F, realmultiplier *pmult, unsigned int i=0,
                unsigned int offset=0);

  void pretransform(Complex *F, Complex *f1c, Complex *U);
  void posttransform(Complex *F, const Complex& f1c, Complex *U);

  // Binary convolution:
  void convolve(Complex *f, Complex *g) {
    Complex *F[]={f,g};
    convolve(F,multbinary);
  }
};


// Compute the scrambled implicitly m-padded complex Fourier transform of M
// complex vectors, each of length m.
// The arrays in and out (which may coincide), along with the array u, must
// be allocated as Complex[M*m].
//
//   fftpad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fftpad {
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int s;
  Complex *ZetaH, *ZetaL;
  unsigned int threads;
public:
  mfft1d *Backwards;
  mfft1d *Forwards;

  fftpad(unsigned int m, unsigned int M,
         unsigned int stride, Complex *u=NULL,
         unsigned int Threads=fftw::maxthreads)
    : m(m), M(M), stride(stride), threads(Threads) {
    Backwards=new mfft1d(m,1,M,stride,1,u,NULL,threads);
    Forwards=new mfft1d(m,-1,M,stride,1,u,NULL,threads);

    threads=std::max(Backwards->Threads(),Forwards->Threads());

    s=BuildZeta(2*m,m,ZetaH,ZetaL,threads);
  }

  ~fftpad() {
    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }

  void expand(Complex *f, Complex *u);
  void reduce(Complex *f, Complex *u);

  void backwards(Complex *f, Complex *u);
  void forwards(Complex *f, Complex *u);
};

// Compute the scrambled implicitly m-padded complex Fourier transform of M
// complex vectors, each of length 2m-1 with the origin at index m-1,
// containing physical data for wavenumbers -m+1 to m-1.
// The arrays in and out (which may coincide) must be allocated as
// Complex[M*(2m-1)]. The array u must be allocated as Complex[M*(m+1)].
//
//   fft0pad fft(m,M,stride,u);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft0pad {
protected:
  unsigned int m;
  unsigned int M;
  unsigned int s;
  unsigned int stride;
  Complex *ZetaH, *ZetaL;
  unsigned int threads;
public:
  mfft1d *Forwards;
  mfft1d *Backwards;

  fft0pad(unsigned int m, unsigned int M, unsigned int stride, Complex *u=NULL,
          unsigned int Threads=fftw::maxthreads)
    : m(m), M(M), stride(stride), threads(Threads) {
    Backwards=new mfft1d(m,1,M,stride,1,u,NULL,threads);
    Forwards=new mfft1d(m,-1,M,stride,1,u,NULL,threads);

    s=BuildZeta(3*m,m,ZetaH,ZetaL);
  }

  virtual ~fft0pad() {
    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }

  // Unscramble indices, returning spatial index stored at position i
  inline static unsigned findex(unsigned i, unsigned int m) {
    return i < m-1 ? 3*i : 3*i+4-3*m; // for i >= m-1: j=3*(i-(m-1))+1
  }

  inline static unsigned uindex(unsigned i, unsigned int m) {
    return i > 0 ? (i < m ? 3*i-1 : 3*m-3) : 3*m-1;
  }

  virtual void expand(Complex *f, Complex *u);
  virtual void reduce(Complex *f, Complex *u);

  void backwards(Complex *f, Complex *u);
  virtual void forwards(Complex *f, Complex *u);

  virtual void Backwards1(Complex *f, Complex *u);
  virtual void Forwards0(Complex *f);
  virtual void Forwards1(Complex *f, Complex *u);
};

// Compute the scrambled implicitly m-padded complex Fourier transform of M
// complex vectors, each of length 2m with the origin at index m,
// corresponding to wavenumbers -m to m-1.
// The arrays in and out (which may coincide) must be allocated as
// Complex[M*2m]. The array u must be allocated as Complex[M*m].
//
//   fft1pad fft(m,M,stride,u);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft1pad : public fft0pad {
public:
  fft1pad(unsigned int m, unsigned int M, unsigned int stride,
          Complex *u=NULL, unsigned int threads=fftw::maxthreads) :
    fft0pad(m,M,stride,u,threads) {}

  // Unscramble indices, returning spatial index stored at position i
  inline static unsigned findex(unsigned i, unsigned int m) {
    return i < m ? 3*i : 3*(i-m)+1;
  }

  inline static unsigned uindex(unsigned i, unsigned int m) {
    return i > 0 ? 3*i-1 : 3*m-1;
  }

  void expand(Complex *f, Complex *u);
  void reduce(Complex *f, Complex *u);

  void forwards(Complex *f, Complex *u);

  void Backwards1(Complex *f, Complex *u);
  void Forwards0(Complex *f);
  void Forwards1(Complex *f, Complex *u);
};

// In-place implicitly dealiased 2D complex convolution.
class ImplicitConvolution2 : public ThreadBase {
protected:
  unsigned int mx,my;
  Complex *u1;
  Complex *u2;
  unsigned int A,B;
  fftpad *xfftpad;
  ImplicitConvolution **yconvolve;
  Complex **U2;
  bool allocated;
  unsigned int indexsize;
  bool toplevel;
public:
  unsigned int *index;

  void initpointers2(Complex **&U2, Complex *u2, unsigned int stride) {
    U2=new Complex *[A];
    for(unsigned int a=0; a < A; ++a)
      U2[a]=u2+a*stride;

    if(toplevel) allocateindex(1,new unsigned int[1]);
  }

  void deletepointers2(Complex **&U2) {
    if(toplevel) {
      delete [] index;

      for(unsigned int t=1; t < threads; ++t)
        delete [] yconvolve[t]->index;
    }

    delete [] U2;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
    yconvolve[0]->allocateindex(n,i);
    for(unsigned int t=1; t < threads; ++t)
      yconvolve[t]->allocateindex(n,new unsigned int[n]);
  }

  void init(const convolveOptions& options) {
    toplevel=options.toplevel;
    xfftpad=new fftpad(mx,options.ny,options.ny,u2,threads);
    unsigned int C=max(A,B);
    yconvolve=new ImplicitConvolution*[threads];
    for(unsigned int t=0; t < threads; ++t)
      yconvolve[t]=new ImplicitConvolution(my,u1+t*my*C,A,B,innerthreads);
    initpointers2(U2,u2,options.stride2);
  }

  void set(convolveOptions& options) {
    if(options.nx == 0) options.nx=mx;
    if(options.ny == 0) {
      options.ny=my;
      options.stride2=mx*my;
    }
  }

  // u1 is a temporary array of size my*C*threads.
  // u2 is a temporary array of size mx*my*C.
  // A is the number of inputs.
  // B is the number of outputs.
  // Here C=max(A,B).
  ImplicitConvolution2(unsigned int mx, unsigned int my,
                       Complex *u1, Complex *u2,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads,
                       convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), u1(u1), u2(u2), A(A), B(B),
    allocated(false) {
    set(options);
    multithread(options.nx);
    init(options);
  }

  ImplicitConvolution2(unsigned int mx, unsigned int my,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads,
                       convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), A(A), B(B), allocated(true) {
    set(options);
    multithread(options.nx);
    unsigned int C=max(A,B);
    u1=utils::ComplexAlign(my*C*threads);
    u2=utils::ComplexAlign(options.stride2*C);
    init(options);
  }

  virtual ~ImplicitConvolution2() {
    deletepointers2(U2);

    for(unsigned int t=0; t < threads; ++t)
      delete yconvolve[t];
    delete [] yconvolve;

    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(u2);
      utils::deleteAlign(u1);
    }
  }

  void backwards(Complex **F, Complex **U2, unsigned int offset) {
    for(unsigned int a=0; a < A; ++a)
      xfftpad->backwards(F[a]+offset,U2[a]);
  }

  void subconvolution(Complex **F, multiplier *pmult,
                      unsigned int r, unsigned int M, unsigned int stride,
                      unsigned int offset=0) {
    if(threads > 1) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < M; ++i)
        yconvolve[get_thread_num()]->convolve(F,pmult,2*i+r,offset+i*stride);
    } else {
      ImplicitConvolution *yconvolve0=yconvolve[0];
      for(unsigned int i=0; i < M; ++i)
        yconvolve0->convolve(F,pmult,2*i+r,offset+i*stride);
    }
  }

  void forwards(Complex **F, Complex **U2, unsigned int offset) {
    for(unsigned int b=0; b < B; ++b)
      xfftpad->forwards(F[b]+offset,U2[b]);
  }

  // F is a pointer to A distinct data blocks each of size mx*my,
  // shifted by offset (contents not preserved).
  virtual void convolve(Complex **F, multiplier *pmult, unsigned int i=0,
                        unsigned int offset=0) {
    if(!toplevel) {
      index[indexsize-2]=i;
      if(threads > 1) {
        for(unsigned int t=1; t < threads; ++t) {
          unsigned int *Index=yconvolve[t]->index;
          for(unsigned int i=0; i < indexsize; ++i)
            Index[i]=index[i];
        }
      }
    }
    backwards(F,U2,offset);
    subconvolution(F,pmult,0,mx,my,offset);
    subconvolution(U2,pmult,1,mx,my);
    forwards(F,U2,offset);
  }

  // Binary convolution:
  void convolve(Complex *f, Complex *g) {
    Complex *F[]={f,g};
    convolve(F,multbinary);
  }

  // Binary correlation:
  void correlate(Complex *f, Complex *g) {
    Complex *F[]={f, g};
    convolve(F,multcorrelation);
  }

  void autoconvolve(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautoconvolution);
  }

  void autocorrelate(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautocorrelation);
  }
};

inline void HermitianSymmetrizeX(unsigned int mx, unsigned int my,
                                 unsigned int xorigin, Complex *f)
{
  unsigned int offset=xorigin*my;
  unsigned int stop=mx*my;
  f[offset].im=0.0;
  for(unsigned int i=my; i < stop; i += my)
    f[offset-i]=conj(f[offset+i]);
}

// Enforce 3D Hermiticity using specified (x,y > 0,z=0) and (x >= 0,y=0,z=0)
// data.
inline void HermitianSymmetrizeXY(unsigned int mx, unsigned int my,
                                  unsigned int mz, unsigned int xorigin,
                                  unsigned int yorigin, Complex *f,
                                  unsigned int threads=fftw::maxthreads)
{
  int stride=(yorigin+my)*mz;
  int mxstride=mx*stride;
  unsigned int myz=my*mz;
  unsigned int origin=xorigin*stride+yorigin*mz;

  f[origin].im=0.0;

  for(int i=stride; i < mxstride; i += stride)
    f[origin-i]=conj(f[origin+i]);

  PARALLEL(
    for(int i=stride-mxstride; i < mxstride; i += stride) {
      int stop=i+myz;
      for(int j=i+mz; j < stop; j += mz) {
        f[origin-j]=conj(f[origin+j]);
      }
    }
    );
}

typedef unsigned int IndexFunction(unsigned int, unsigned int m);

class ImplicitHConvolution2 : public ThreadBase {
protected:
  unsigned int mx,my;
  bool xcompact,ycompact;
  Complex *u1;
  Complex *u2;
  unsigned int A,B;
  fft0pad *xfftpad;
  ImplicitHConvolution **yconvolve;
  Complex **U2;
  bool allocated;
  unsigned int indexsize;
  bool toplevel;
public:
  unsigned int *index;

  void initpointers2(Complex **&U2, Complex *u2, unsigned int stride)
  {
    unsigned int C=max(A,B);
    U2=new Complex *[C];
    for(unsigned int a=0; a < C; ++a)
      U2[a]=u2+a*stride;

    if(toplevel) allocateindex(1,new unsigned int[1]);
  }

  void deletepointers2(Complex **&U2) {
    if(toplevel) {
      delete [] index;

      for(unsigned int t=1; t < threads; ++t)
        delete [] yconvolve[t]->index;
    }

    delete [] U2;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
    yconvolve[0]->allocateindex(n,i);
    for(unsigned int t=1; t < threads; ++t)
      yconvolve[t]->allocateindex(n,new unsigned int[n]);
  }

  void init(const convolveOptions& options) {
    unsigned int C=max(A,B);
    toplevel=options.toplevel;
    xfftpad=xcompact ? new fft0pad(mx,options.ny,options.ny,u2) :
      new fft1pad(mx,options.ny,options.ny,u2);

    yconvolve=new ImplicitHConvolution*[threads];
    for(unsigned int t=0; t < threads; ++t)
      yconvolve[t]=new ImplicitHConvolution(my,ycompact,u1+t*(my/2+1)*C,A,B,
                                            innerthreads);
    initpointers2(U2,u2,options.stride2);
  }

  void set(convolveOptions& options) {
    if(options.nx == 0) options.nx=mx;
    if(options.ny == 0) {
      options.ny=my+!ycompact;
      options.stride2=(mx+xcompact)*options.ny;
    }
  }

  // u1 is a temporary array of size (my/2+1)*C*threads.
  // u2 is a temporary array of size (mx+xcompact)*(my+!ycompact)*C;
  // A is the number of inputs.
  // B is the number of outputs.
  // Here C=max(A,B).
  ImplicitHConvolution2(unsigned int mx, unsigned int my,
                        Complex *u1, Complex *u2,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), xcompact(true), ycompact(true),
    u1(u1), u2(u2), A(A), B(B), allocated(false) {
    set(options);
    multithread(options.nx);
    init(options);
  }

  ImplicitHConvolution2(unsigned int mx, unsigned int my,
                        bool xcompact, bool ycompact,
                        Complex *u1, Complex *u2,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my),
    xcompact(xcompact), ycompact(ycompact), u1(u1), u2(u2), A(A), B(B),
    allocated(false) {
    set(options);
    multithread(options.nx);
    init(options);
  }

  ImplicitHConvolution2(unsigned int mx, unsigned int my,
                        bool xcompact=true, bool ycompact=true,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my),
    xcompact(xcompact), ycompact(ycompact), A(A), B(B), allocated(true) {
    set(options);
    multithread(options.nx);
    unsigned int C=max(A,B);
    u1=utils::ComplexAlign((my/2+1)*C*threads);
    u2=utils::ComplexAlign(options.stride2*C);
    init(options);
  }

  virtual ~ImplicitHConvolution2() {
    deletepointers2(U2);

    for(unsigned int t=0; t < threads; ++t)
      delete yconvolve[t];
    delete [] yconvolve;

    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(u2);
      utils::deleteAlign(u1);
    }
  }

  void backwards(Complex **F, Complex **U2, unsigned int ny,
                 bool symmetrize, unsigned int offset) {
    for(unsigned int a=0; a < A; ++a) {
      Complex *f=F[a]+offset;
      if(symmetrize)
        HermitianSymmetrizeX(mx,ny,mx-xcompact,f);
      xfftpad->backwards(f,U2[a]);
    }
  }

  void subconvolution(Complex **F, realmultiplier *pmult,
                      IndexFunction indexfunction,
                      unsigned int M, unsigned int stride,
                      unsigned int offset=0) {
    if(threads > 1) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < M; ++i)
        yconvolve[get_thread_num()]->convolve(F,pmult,indexfunction(i,mx),
                                              offset+i*stride);
    } else {
      ImplicitHConvolution *yconvolve0=yconvolve[0];
      for(unsigned int i=0; i < M; ++i)
        yconvolve0->convolve(F,pmult,indexfunction(i,mx),offset+i*stride);
    }
  }

  void forwards(Complex **F, Complex **U2, unsigned int offset) {
    for(unsigned int b=0; b < B; ++b)
      xfftpad->forwards(F[b]+offset,U2[b]);
  }

  // F is a pointer to A distinct data blocks each of size
  // (2mx-xcompact)*(my+!ycompact), shifted by offset (contents not preserved).
  virtual void convolve(Complex **F, realmultiplier *pmult,
                        bool symmetrize=true, unsigned int i=0,
                        unsigned int offset=0) {
    if(!toplevel) {
      index[indexsize-2]=i;
      if(threads > 1) {
        for(unsigned int t=1; t < threads; ++t) {
          unsigned int *Index=yconvolve[t]->index;
          for(unsigned int i=0; i < indexsize; ++i)
            Index[i]=index[i];
        }
      }
    }
    unsigned stride=my+!ycompact;
    backwards(F,U2,stride,symmetrize,offset);
    subconvolution(F,pmult,xfftpad->findex,2*mx-xcompact,stride,offset);
    subconvolution(U2,pmult,xfftpad->uindex,mx+xcompact,stride);
    forwards(F,U2,offset);
  }

  // Binary convolution:
  void convolve(Complex *f, Complex *g, bool symmetrize=true) {
    Complex *F[]={f,g};
    convolve(F,multbinary,symmetrize);
  }
};

// In-place implicitly dealiased 3D complex convolution.
class ImplicitConvolution3 : public ThreadBase {
protected:
  unsigned int mx,my,mz;
  Complex *u1;
  Complex *u2;
  Complex *u3;
  unsigned int A,B;
  fftpad *xfftpad;
  ImplicitConvolution2 **yzconvolve;
  Complex **U3;
  bool allocated;
  unsigned int indexsize;
  bool toplevel;
public:
  unsigned int *index;

  void initpointers3(Complex **&U3, Complex *u3, unsigned int stride) {
    unsigned int C=max(A,B);
    U3=new Complex *[C];
    for(unsigned int a=0; a < C; ++a)
      U3[a]=u3+a*stride;

    if(toplevel) allocateindex(2,new unsigned int[2]);
  }

  void deletepointers3(Complex **&U3) {
    if(toplevel) {
      delete [] index;

      for(unsigned int t=1; t < threads; ++t)
        delete [] yzconvolve[t]->index;
    }

    delete [] U3;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
    yzconvolve[0]->allocateindex(n,i);
    for(unsigned int t=1; t < threads; ++t)
      yzconvolve[t]->allocateindex(n,new unsigned int[n]);
  }

  void init(const convolveOptions& options) {
    toplevel=options.toplevel;
    unsigned int nyz=options.ny*options.nz;
    xfftpad=new fftpad(mx,nyz,nyz,u3,threads);

    if(options.nz == mz) {
      unsigned int C=max(A,B);
      yzconvolve=new ImplicitConvolution2*[threads];
      for(unsigned int t=0; t < threads; ++t)
        yzconvolve[t]=new ImplicitConvolution2(my,mz,u1+t*mz*C*innerthreads,
                                               u2+t*options.stride2*C,A,B,
                                               innerthreads,false);
      initpointers3(U3,u3,options.stride3);
    } else yzconvolve=NULL;
  }

  void set(convolveOptions &options)
  {
    if(options.ny == 0) {
      options.ny=my;
      options.nz=mz;
      options.stride2=my*mz;
      options.stride3=mx*my*mz;
    }
  }

  // u1 is a temporary array of size mz*C*threads.
  // u2 is a temporary array of size my*mz*C*threads.
  // u3 is a temporary array of size mx*my*mz*C.
  // A is the number of inputs.
  // B is the number of outputs.
  // Here C=max(A,B).
  ImplicitConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                       Complex *u1, Complex *u2, Complex *u3,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads,
                       convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), mz(mz),
    u1(u1), u2(u2), u3(u3), A(A), B(B), allocated(false) {
    set(options);
    multithread(mx);
    init(options);
  }

  ImplicitConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                       unsigned int A=2, unsigned int B=1,
                       unsigned int threads=fftw::maxthreads,
                       convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), mz(mz), A(A), B(B),
    allocated(true) {
    set(options);
    multithread(mx);
    unsigned int C=max(A,B);
    u1=utils::ComplexAlign(mz*C*threads*innerthreads);
    u2=utils::ComplexAlign(options.stride2*C*threads);
    u3=utils::ComplexAlign(options.stride3*C);
    init(options);
  }

  virtual ~ImplicitConvolution3() {
    if(yzconvolve) {
      deletepointers3(U3);

      for(unsigned int t=0; t < threads; ++t)
        delete yzconvolve[t];
      delete [] yzconvolve;
    }

    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(u3);
      utils::deleteAlign(u2);
      utils::deleteAlign(u1);
    }
  }

  void backwards(Complex **F, Complex **U3, unsigned int offset) {
    for(unsigned int a=0; a < A; ++a)
      xfftpad->backwards(F[a]+offset,U3[a]);
  }

  void subconvolution(Complex **F, multiplier *pmult,
                      unsigned int r, unsigned int M, unsigned int stride,
                      unsigned int offset=0) {
    if(threads > 1) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < M; ++i)
        yzconvolve[get_thread_num()]->convolve(F,pmult,2*i+r,offset+i*stride);
    } else {
      ImplicitConvolution2 *yzconvolve0=yzconvolve[0];
      for(unsigned int i=0; i < M; ++i) {
        yzconvolve0->convolve(F,pmult,2*i+r,offset+i*stride);
      }
    }
  }

  void forwards(Complex **F, Complex **U3, unsigned int offset=0) {
    for(unsigned int b=0; b < B; ++b)
      xfftpad->forwards(F[b]+offset,U3[b]);
  }

  // F is a pointer to A distinct data blocks each of size mx*my*mz,
  // shifted by offset
  virtual void convolve(Complex **F, multiplier *pmult, unsigned int i=0,
                        unsigned int offset=0)
  {
    if(!toplevel) {
      index[indexsize-3]=i;
      if(threads > 1) {
        for(unsigned int t=1; t < threads; ++t) {
          unsigned int *Index=yzconvolve[t]->index;
          for(unsigned int i=0; i < indexsize; ++i)
            Index[i]=index[i];
        }
      }
    }
    unsigned int stride=my*mz;
    backwards(F,U3,offset);
    subconvolution(F,pmult,0,mx,stride,offset);
    subconvolution(U3,pmult,1,mx,stride);
    forwards(F,U3,offset);
  }

  // Binary convolution:
  void convolve(Complex *f, Complex *g) {
    Complex *F[]={f,g};
    convolve(F,multbinary);
  }

  // Binary correlation:
  void correlate(Complex *f, Complex *g) {
    Complex *F[]={f, g};
    convolve(F,multcorrelation);
  }

  void autoconvolve(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautoconvolution);
  }

  void autocorrelate(Complex *f) {
    Complex *F[]={f};
    convolve(F,multautocorrelation);
  }
};

// In-place implicitly dealiased 3D Hermitian convolution.
class ImplicitHConvolution3 : public ThreadBase {
protected:
  unsigned int mx,my,mz;
  bool xcompact,ycompact,zcompact;
  Complex *u1;
  Complex *u2;
  Complex *u3;
  unsigned int A,B;
  fft0pad *xfftpad;
  ImplicitHConvolution2 **yzconvolve;
  Complex **U3;
  bool allocated;
  unsigned int indexsize;
  bool toplevel;
public:
  unsigned int *index;

  void initpointers3(Complex **&U3, Complex *u3, unsigned int stride) {
    unsigned int C=max(A,B);
    U3=new Complex *[C];
    for(unsigned int a=0; a < C; ++a)
      U3[a]=u3+a*stride;

    if(toplevel) allocateindex(2,new unsigned int[2]);
  }

  void deletepointers3(Complex **&U3) {
    if(toplevel) {
      delete [] index;

      for(unsigned int t=1; t < threads; ++t)
        delete [] yzconvolve[t]->index;
    }

    delete [] U3;
  }

  void allocateindex(unsigned int n, unsigned int *i) {
    indexsize=n;
    index=i;
    yzconvolve[0]->allocateindex(n,i);
    for(unsigned int t=1; t < threads; ++t)
      yzconvolve[t]->allocateindex(n,new unsigned int[n]);
  }

  void init(const convolveOptions& options) {
    toplevel=options.toplevel;
    unsigned int nyz=options.ny*options.nz;
    xfftpad=xcompact ? new fft0pad(mx,nyz,nyz,u3) :
      new fft1pad(mx,nyz,nyz,u3);

    if(options.nz == mz+!zcompact) {
      unsigned int C=max(A,B);
      yzconvolve=new ImplicitHConvolution2*[threads];
      for(unsigned int t=0; t < threads; ++t)
        yzconvolve[t]=new ImplicitHConvolution2(my,mz,
                                                ycompact,zcompact,
                                                u1+t*(mz/2+1)*C*innerthreads,
                                                u2+t*options.stride2*C,
                                                A,B,innerthreads,false);
      initpointers3(U3,u3,options.stride3);
    } else yzconvolve=NULL;
  }

  void set(convolveOptions& options) {
    if(options.ny == 0) {
      options.ny=2*my-ycompact;
      options.nz=mz+!zcompact;
      options.stride2=(my+ycompact)*options.nz;
      options.stride3=(mx+xcompact)*options.ny*options.nz;
    }
  }

  // u1 is a temporary array of size (mz/2+1)*C*threads.
  // u2 is a temporary array of size (my+ycompact)*(mz+!zcompact)*C*threads.
  // u3 is a temporary array of size
  //                             (mx+xcompact)*(2my-ycompact)*(mz+!zcompact)*C.
  // A is the number of inputs.
  // B is the number of outputs.
  // Here C=max(A,B).
  ImplicitHConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                        Complex *u1, Complex *u2, Complex *u3,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), mz(mz),
    xcompact(true), ycompact(true), zcompact(true), u1(u1), u2(u2), u3(u3),
    A(A), B(B),
    allocated(false) {
    set(options);
    multithread(mx);
    init(options);
  }

  ImplicitHConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                        bool xcompact, bool ycompact, bool zcompact,
                        Complex *u1, Complex *u2, Complex *u3,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), mz(mz),
    xcompact(xcompact), ycompact(ycompact), zcompact(zcompact),
    u1(u1), u2(u2), u3(u3), A(A), B(B), allocated(false) {
    set(options);
    multithread(mx);
    init(options);
  }

  ImplicitHConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                        bool xcompact=true, bool ycompact=true,
                        bool zcompact=true,
                        unsigned int A=2, unsigned int B=1,
                        unsigned int threads=fftw::maxthreads,
                        convolveOptions options=defaultconvolveOptions) :
    ThreadBase(threads), mx(mx), my(my), mz(mz),
    xcompact(xcompact), ycompact(ycompact), zcompact(zcompact), A(A), B(B),
    allocated(true) {
    set(options);
    multithread(mx);
    unsigned int C=max(A,B);
    u1=utils::ComplexAlign((mz/2+1)*C*threads*innerthreads);
    u2=utils::ComplexAlign(options.stride2*C*threads);
    u3=utils::ComplexAlign(options.stride3*C);
    init(options);
  }

  virtual ~ImplicitHConvolution3() {
    if(yzconvolve) {
      deletepointers3(U3);

      for(unsigned int t=0; t < threads; ++t)
        delete yzconvolve[t];
      delete [] yzconvolve;
    }

    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(u3);
      utils::deleteAlign(u2);
      utils::deleteAlign(u1);
    }
  }

  virtual void HermitianSymmetrize(Complex *f, Complex *u)
  {
    HermitianSymmetrizeXY(mx,my,mz+!zcompact,mx-xcompact,my-ycompact,f,
                          threads);
  }

  void backwards(Complex **F, Complex **U3, bool symmetrize,
                 unsigned int offset) {
    for(unsigned int a=0; a < A; ++a) {
      Complex *f=F[a]+offset;
      Complex *u=U3[a];
      if(symmetrize)
        HermitianSymmetrize(f,u);
      xfftpad->backwards(f,u);
    }
  }

  void subconvolution(Complex **F, realmultiplier *pmult,
                      IndexFunction indexfunction,
                      unsigned int M, unsigned int stride,
                      unsigned int offset=0) {
    if(threads > 1) {
#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
      for(unsigned int i=0; i < M; ++i)
        yzconvolve[get_thread_num()]->convolve(F,pmult,false,
                                               indexfunction(i,mx),
                                               offset+i*stride);
    } else {
      ImplicitHConvolution2 *yzconvolve0=yzconvolve[0];
      for(unsigned int i=0; i < M; ++i)
        yzconvolve0->convolve(F,pmult,false,indexfunction(i,mx),
                              offset+i*stride);
    }
  }

  void forwards(Complex **F, Complex **U3, unsigned int offset=0) {
    for(unsigned int b=0; b < B; ++b)
      xfftpad->forwards(F[b]+offset,U3[b]);
  }

  // F is a pointer to A distinct data blocks each of size
  // (2mx-compact)*(2my-ycompact)*(mz+!zcompact), shifted by offset
  // (contents not preserved).
  virtual void convolve(Complex **F, realmultiplier *pmult,
                        bool symmetrize=true, unsigned int i=0,
                        unsigned int offset=0) {
    if(!toplevel) {
      index[indexsize-3]=i;
      if(threads > 1) {
        for(unsigned int t=1; t < threads; ++t) {
          unsigned int *Index=yzconvolve[t]->index;
          for(unsigned int i=0; i < indexsize; ++i)
            Index[i]=index[i];
        }
      }
    }
    unsigned int stride=(2*my-ycompact)*(mz+!zcompact);
    backwards(F,U3,symmetrize,offset);
    subconvolution(F,pmult,xfftpad->findex,2*mx-xcompact,stride,offset);
    subconvolution(U3,pmult,xfftpad->uindex,mx+xcompact,stride);
    forwards(F,U3,offset);
  }

  // Binary convolution:
  void convolve(Complex *f, Complex *g, bool symmetrize=true) {
    Complex *F[]={f,g};
    convolve(F,multbinary,symmetrize);
  }
};

// In-place implicitly dealiased Hermitian ternary convolution.
class ImplicitHTConvolution : public ThreadBase {
protected:
  unsigned int m;
  Complex *u,*v,*w;
  unsigned int M;
  unsigned int s;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  Complex *ZetaH, *ZetaL;
  Complex **W;
  bool allocated;
  unsigned int twom;
  unsigned int stride;
public:
  void initpointers(Complex **&W, Complex *w) {
    W=new Complex *[M];
    unsigned int m1=m+1;
    for(unsigned int s=0; s < M; ++s)
      W[s]=w+s*m1;
  }

  void deletepointers(Complex **&W) {
    delete [] W;
  }

  void init() {
    twom=2*m;
    stride=twom+2;

    rc=new rcfft1d(twom,u);
    cr=new crfft1d(twom,u);

    rco=new rcfft1d(twom,(double *) u,v);
    cro=new crfft1d(twom,v,(double *) u);

    threads=std::min(threads,std::max(rco->Threads(),cro->Threads()));

    s=BuildZeta(4*m,m,ZetaH,ZetaL,threads);

    initpointers(W,w);
  }

  // u, v, and w are distinct temporary arrays each of size (m+1)*M.
  ImplicitHTConvolution(unsigned int m, Complex *u, Complex *v,
                        Complex *w, unsigned int M=1) :
    m(m), u(u), v(v), w(w), M(M), allocated(false) {
    init();
  }

  ImplicitHTConvolution(unsigned int m, unsigned int M=1) :
    m(m), u(utils::ComplexAlign(m*M+M)), v(utils::ComplexAlign(m*M+M)),
    w(utils::ComplexAlign(m*M+M)), M(M), allocated(true) {
    init();
  }

  ~ImplicitHTConvolution() {
    deletepointers(W);

    if(allocated) {
      utils::deleteAlign(w);
      utils::deleteAlign(v);
      utils::deleteAlign(u);
    }
    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete cro;
    delete rco;
    delete cr;
    delete rc;
  }

  void mult(double *a, double *b, double **C, unsigned int offset=0);

  void convolve(Complex **F, Complex **G, Complex **H,
                Complex *u, Complex *v, Complex **W,
                unsigned int offset=0);

  // F, G, and H are distinct pointers to M distinct data blocks each of size
  // m+1, shifted by offset (contents not preserved).
  // The output is returned in F[0].
  void convolve(Complex **F, Complex **G, Complex **H, unsigned int offset=0) {
    convolve(F,G,H,u,v,W,offset);
  }

  // Constructor for special case M=1:
  void convolve(Complex *f, Complex *g, Complex *h) {
    convolve(&f,&g,&h);
  }
};

// In-place implicitly dealiased Hermitian ternary convolution.
// Special case G=H, M=1.
class ImplicitHFGGConvolution : public ThreadBase {
protected:
  unsigned int m;
  Complex *u,*v;
  unsigned int s;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  Complex *ZetaH, *ZetaL;
  bool allocated;
  unsigned int twom;
  unsigned int stride;
public:
  void init() {
    twom=2*m;
    stride=twom+2;

    rc=new rcfft1d(twom,u);
    cr=new crfft1d(twom,u);

    rco=new rcfft1d(twom,(double *) u,v);
    cro=new crfft1d(twom,v,(double *) u);

    threads=std::min(threads,std::max(rco->Threads(),cro->Threads()));

    s=BuildZeta(4*m,m,ZetaH,ZetaL,threads);
  }

  // u and v are distinct temporary arrays each of size m+1.
  ImplicitHFGGConvolution(unsigned int m, Complex *u, Complex *v) :
    m(m), u(u), v(v), allocated(false) {
    init();
  }

  ImplicitHFGGConvolution(unsigned int m) :
    m(m), u(utils::ComplexAlign(m+1)), v(utils::ComplexAlign(m+1)),
    allocated(true) {
    init();
  }

  ~ImplicitHFGGConvolution() {
    if(allocated) {
      utils::deleteAlign(v);
      utils::deleteAlign(u);
    }
    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete cro;
    delete rco;
    delete cr;
    delete rc;
  }

  void mult(double *a, double *b);

  void convolve(Complex *f, Complex *g, Complex *u, Complex *v);

  // f and g are distinct pointers to data of size m+1 (contents not
  // preserved). The output is returned in f.
  void convolve(Complex *f, Complex *g) {
    convolve(f,g,u,v);
  }
};

// In-place implicitly dealiased Hermitian ternary convolution.
// Special case F=G=H, M=1.
class ImplicitHFFFConvolution : public ThreadBase {
protected:
  unsigned int m;
  Complex *u;
  unsigned int s;
  rcfft1d *rc;
  crfft1d *cr;
  Complex *ZetaH, *ZetaL;
  bool allocated;
  unsigned int twom;
  unsigned int stride;
public:
  void mult(double *a);

  void init() {
    twom=2*m;
    stride=twom+2;

    rc=new rcfft1d(twom,u);
    cr=new crfft1d(twom,u);

    threads=std::min(threads,std::max(rc->Threads(),cr->Threads()));

    s=BuildZeta(4*m,m,ZetaH,ZetaL,threads);
  }

  // u is a distinct temporary array of size m+1.
  ImplicitHFFFConvolution(unsigned int m, Complex *u) :
    m(m), u(u), allocated(false) {
    init();
  }

  ImplicitHFFFConvolution(unsigned int m) :
    m(m), u(utils::ComplexAlign(m+1)), allocated(true) {
    init();
  }

  ~ImplicitHFFFConvolution() {
    if(allocated)
      utils::deleteAlign(u);

    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete cr;
    delete rc;
  }

  void convolve(Complex *f, Complex *u);
  // f is a pointer to data of size m+1 (contents not preserved).
  // The output is returned in f.
  void convolve(Complex *f) {
    convolve(f,u);
  }
};

// Compute the scrambled implicitly 2m-padded complex Fourier transform of M
// complex vectors, each of length 2m with the Fourier origin at index m.
// The arrays in and out (which may coincide), along
// with the array u, must be allocated as Complex[M*2m].
//
//   fft0bipad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft0bipad {
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int s;
  mfft1d *Backwards;
  mfft1d *Forwards;
  Complex *ZetaH, *ZetaL;
  unsigned int threads;
public:
  fft0bipad(unsigned int m, unsigned int M, unsigned int stride,
            Complex *f, unsigned int Threads=fftw::maxthreads) :
    m(m), M(M), stride(stride), threads(Threads) {
    unsigned int twom=2*m;
    Backwards=new mfft1d(twom,1,M,stride,1,f,NULL,threads);
    Forwards=new mfft1d(twom,-1,M,stride,1,f,NULL,threads);

    threads=std::min(threads,
                     std::max(Backwards->Threads(),Forwards->Threads()));

    s=BuildZeta(4*m,twom,ZetaH,ZetaL,threads);
  }

  ~fft0bipad() {
    utils::deleteAlign(ZetaL);
    utils::deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }

  void backwards(Complex *f, Complex *u);
  void forwards(Complex *f, Complex *u);
};

// In-place implicitly dealiased 2D Hermitian ternary convolution.
class ImplicitHTConvolution2 : public ThreadBase {
protected:
  unsigned int mx,my;
  Complex *u1,*v1,*w1;
  Complex *u2,*v2,*w2;
  unsigned int M;
  fft0bipad *xfftpad;
  ImplicitHTConvolution *yconvolve;
  Complex **U2,**V2,**W2;
  bool allocated;
  Complex **u,**v;
  Complex ***W;
public:
  void initpointers(Complex **&u, Complex **&v, Complex ***&W,
                    unsigned int threads) {
    u=new Complex *[threads];
    v=new Complex *[threads];
    W=new Complex **[threads];
    unsigned int my1M=(my+1)*M;
    for(unsigned int i=0; i < threads; ++i) {
      unsigned int imy1M=i*my1M;
      u[i]=u1+imy1M;
      v[i]=v1+imy1M;
      Complex *wi=w1+imy1M;
      yconvolve->initpointers(W[i],wi);
    }
  }

  void deletepointers(Complex **&u, Complex **&v, Complex ***&W,
                      unsigned int threads) {
    for(unsigned int i=0; i < threads; ++i)
      yconvolve->deletepointers(W[i]);
    delete [] W;
    delete [] v;
    delete [] u;
  }

  void initpointers(Complex **&U2, Complex **&V2, Complex **&W2,
                    Complex *u2, Complex *v2, Complex *w2) {
    U2=new Complex *[M];
    V2=new Complex *[M];
    W2=new Complex *[M];
    unsigned int mu=2*mx*(my+1);
    for(unsigned int s=0; s < M; ++s) {
      unsigned int smu=s*mu;
      U2[s]=u2+smu;
      V2[s]=v2+smu;
      W2[s]=w2+smu;
    }
  }

  void deletepointers(Complex **&U2, Complex **&V2, Complex **&W2) {
    delete [] W2;
    delete [] V2;
    delete [] U2;
  }

  void init() {
    xfftpad=new fft0bipad(mx,my,my+1,u2,threads);

    yconvolve=new ImplicitHTConvolution(my,u1,v1,w1,M);
    yconvolve->Threads(1);

    initpointers(u,v,W,threads);
    initpointers(U2,V2,W2,u2,v2,w2);
  }

  // u1, v1, and w1 are temporary arrays of size (my+1)*M*threads;
  // u2, v2, and w2 are temporary arrays of size 2mx*(my+1)*M.
  // M is the number of data blocks (each corresponding to a dot product term).
  // threads is the number of threads to use in the outer subconvolution loop.
  ImplicitHTConvolution2(unsigned int mx, unsigned int my,
                         Complex *u1, Complex *v1, Complex *w1,
                         Complex *u2, Complex *v2, Complex *w2,
                         unsigned int M=1,
                         unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my), u1(u1), v1(v1), w1(w1),
    u2(u2), v2(v2), w2(w2), M(M), allocated(false) {
    init();
  }

  ImplicitHTConvolution2(unsigned int mx, unsigned int my,
                         unsigned int M=1,
                         unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my),
    u1(utils::ComplexAlign((my+1)*M*threads)),
    v1(utils::ComplexAlign((my+1)*M*threads)),
    w1(utils::ComplexAlign((my+1)*M*threads)),
    u2(utils::ComplexAlign(2*mx*(my+1)*M)),
    v2(utils::ComplexAlign(2*mx*(my+1)*M)),
    w2(utils::ComplexAlign(2*mx*(my+1)*M)),
    M(M), allocated(true) {
    init();
  }

  ~ImplicitHTConvolution2() {
    deletepointers(U2,V2,W2);
    deletepointers(u,v,W,threads);

    delete yconvolve;
    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(w2);
      utils::deleteAlign(v2);
      utils::deleteAlign(u2);
      utils::deleteAlign(w1);
      utils::deleteAlign(v1);
      utils::deleteAlign(u1);
    }
  }

  void convolve(Complex **F, Complex **G, Complex **H,
                Complex **u, Complex **v, Complex ***W,
                Complex **U2, Complex **V2, Complex **W2,
                bool symmetrize=true, unsigned int offset=0) {
    Complex *u2=U2[0];
    Complex *v2=V2[0];
    Complex *w2=W2[0];

    unsigned int my1=my+1;
    unsigned int mu=2*mx*my1;

    for(unsigned int s=0; s < M; ++s) {
      Complex *f=F[s]+offset;
      if(symmetrize)
        HermitianSymmetrizeX(mx,my1,mx,f);
      xfftpad->backwards(f,u2+s*mu);
    }

    for(unsigned int s=0; s < M; ++s) {
      Complex *g=G[s]+offset;
      if(symmetrize)
        HermitianSymmetrizeX(mx,my1,mx,g);
      xfftpad->backwards(g,v2+s*mu);
    }

    for(unsigned int s=0; s < M; ++s) {
      Complex *h=H[s]+offset;
      if(symmetrize)
        HermitianSymmetrizeX(mx,my1,mx,h);
      xfftpad->backwards(h,w2+s*mu);
    }

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1) {
      unsigned int thread=get_thread_num();
      yconvolve->convolve(F,G,H,u[thread],v[thread],W[thread],i+offset);
    }

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1) {
      unsigned int thread=get_thread_num();
      yconvolve->convolve(U2,V2,W2,u[thread],v[thread],W[thread],i+offset);
    }

    xfftpad->forwards(F[0]+offset,u2);
  }

  // F, G, and H are distinct pointers to M distinct data blocks each of size
  // 2mx*(my+1), shifted by offset (contents not preserved).
  // The output is returned in F[0].
  void convolve(Complex **F, Complex **G, Complex **H, bool symmetrize=true,
                unsigned int offset=0) {
    convolve(F,G,H,u,v,W,U2,V2,W2,symmetrize,offset);
  }

  // Constructor for special case M=1:
  void convolve(Complex *f, Complex *g, Complex *h, bool symmetrize=true) {
    convolve(&f,&g,&h,symmetrize);
  }
};

// In-place implicitly dealiased 2D Hermitian ternary convolution.
// Special case G=H, M=1.
class ImplicitHFGGConvolution2 : public ThreadBase {
protected:
  unsigned int mx,my;
  Complex *u1,*v1;
  Complex *u2,*v2;
  fft0bipad *xfftpad;
  ImplicitHFGGConvolution *yconvolve;
  bool allocated;
  Complex **u,**v;
public:
  void initpointers(Complex **&u, Complex **&v, unsigned int threads) {
    u=new Complex *[threads];
    v=new Complex *[threads];
    unsigned int my1=my+1;
    for(unsigned int i=0; i < threads; ++i) {
      unsigned int imy1=i*my1;
      u[i]=u1+imy1;
      v[i]=v1+imy1;
    }
  }

  void deletepointers(Complex **&u, Complex **&v) {
    delete [] v;
    delete [] u;
  }

  void init() {
    xfftpad=new fft0bipad(mx,my,my+1,u2,threads);

    yconvolve=new ImplicitHFGGConvolution(my,u1,v1);
    yconvolve->Threads(1);

    initpointers(u,v,threads);
  }

  // u1 and v1 are temporary arrays of size (my+1)*threads.
  // u2 and v2 are temporary arrays of size 2mx*(my+1).
  // threads is the number of threads to use in the outer subconvolution loop.
  ImplicitHFGGConvolution2(unsigned int mx, unsigned int my,
                           Complex *u1, Complex *v1,
                           Complex *u2, Complex *v2,
                           unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my), u1(u1), v1(v1), u2(u2), v2(v2),
    allocated(false) {
    init();
  }

  ImplicitHFGGConvolution2(unsigned int mx, unsigned int my,
                           unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my),
    u1(utils::ComplexAlign((my+1)*threads)),
    v1(utils::ComplexAlign((my+1)*threads)),
    u2(utils::ComplexAlign(2*mx*(my+1))),
    v2(utils::ComplexAlign(2*mx*(my+1))),
    allocated(true) {
    init();
  }

  ~ImplicitHFGGConvolution2() {
    deletepointers(u,v);

    delete yconvolve;
    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(v2);
      utils::deleteAlign(u2);
      utils::deleteAlign(v1);
      utils::deleteAlign(u1);
    }
  }

  void convolve(Complex *f, Complex *g,
                Complex **u, Complex **v,
                Complex *u2, Complex *v2, bool symmetrize=true) {
    unsigned int my1=my+1;
    unsigned int mu=2*mx*my1;

    if(symmetrize)
      HermitianSymmetrizeX(mx,my1,mx,f);
    xfftpad->backwards(f,u2);

    if(symmetrize)
      HermitianSymmetrizeX(mx,my1,mx,g);
    xfftpad->backwards(g,v2);

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1) {
      unsigned int thread=get_thread_num();
      yconvolve->convolve(f+i,g+i,u[thread],v[thread]);
    }

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1) {
      unsigned int thread=get_thread_num();
      yconvolve->convolve(u2+i,v2+i,u[thread],v[thread]);
    }

    xfftpad->forwards(f,u2);
  }

  void convolve(Complex *f, Complex *g, bool symmetrize=true) {
    convolve(f,g,u,v,u2,v2,symmetrize);
  }
};

// In-place implicitly dealiased 2D Hermitian ternary convolution.
// Special case F=G=H, M=1.
class ImplicitHFFFConvolution2 : public ThreadBase {
protected:
  unsigned int mx,my;
  Complex *u1;
  Complex *u2;
  fft0bipad *xfftpad;
  ImplicitHFFFConvolution *yconvolve;
  bool allocated;
  Complex **u;
public:
  void initpointers(Complex **&u, unsigned int threads) {
    u=new Complex *[threads];
    unsigned int my1=my+1;
    for(unsigned int i=0; i < threads; ++i)
      u[i]=u1+i*my1;
  }

  void deletepointers(Complex **&u) {
    delete [] u;
  }

  void init() {
    xfftpad=new fft0bipad(mx,my,my+1,u2,threads);

    yconvolve=new ImplicitHFFFConvolution(my,u1);
    yconvolve->Threads(1);
    initpointers(u,threads);
  }

  // u1 is a temporary array of size (my+1)*threads.
  // u2 is a temporary array of size 2mx*(my+1).
  // threads is the number of threads to use in the outer subconvolution loop.
  ImplicitHFFFConvolution2(unsigned int mx, unsigned int my,
                           Complex *u1, Complex *u2,
                           unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my),
    u1(u1), u2(u2), allocated(false) {
    init();
  }

  ImplicitHFFFConvolution2(unsigned int mx, unsigned int my,
                           unsigned int threads=fftw::maxthreads) :
    ThreadBase(threads), mx(mx), my(my),
    u1(utils::ComplexAlign((my+1)*threads)),
    u2(utils::ComplexAlign(2*mx*(my+1))),
    allocated(true) {
    init();
  }

  ~ImplicitHFFFConvolution2() {
    deletepointers(u);
    delete yconvolve;
    delete xfftpad;

    if(allocated) {
      utils::deleteAlign(u2);
      utils::deleteAlign(u1);
    }
  }

  void convolve(Complex *f, Complex **u, Complex *u2, bool symmetrize=true) {
    unsigned int my1=my+1;
    unsigned int mu=2*mx*my1;

    if(symmetrize)
      HermitianSymmetrizeX(mx,my1,mx,f);
    xfftpad->backwards(f,u2);

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1)
      yconvolve->convolve(f+i,u[get_thread_num()]);

#ifndef FFTWPP_SINGLE_THREAD
#pragma omp parallel for num_threads(threads)
#endif
    for(unsigned int i=0; i < mu; i += my1)
      yconvolve->convolve(u2+i,u[get_thread_num()]);

    xfftpad->forwards(f,u2);
  }

  void convolve(Complex *f, bool symmetrize=true) {
    convolve(f,u,u2,symmetrize);
  }
};


} //end namespace fftwpp

}


using namespace std;

namespace fftwpp {

const double fftw::twopi=2.0*acos(-1.0);

// User settings:
unsigned int fftw::effort=FFTW_MEASURE;
const char *fftw::WisdomName="wisdom3.txt";
unsigned int fftw::maxthreads=1;
double fftw::testseconds=0.2; // Time limit for threading efficiency tests

fftw_plan (*fftw::planner)(fftw *f, Complex *in, Complex *out)=Planner;

const char *fftw::oddshift="Shift is not implemented for odd nx";
const char *inout=
  "constructor and call must be both in place or both out of place";

fft1d::Table fft1d::threadtable;
mfft1d::Table mfft1d::threadtable;
rcfft1d::Table rcfft1d::threadtable;
crfft1d::Table crfft1d::threadtable;
mrcfft1d::Table mrcfft1d::threadtable;
mcrfft1d::Table mcrfft1d::threadtable;
fft2d::Table fft2d::threadtable;

void LoadWisdom()
{
  static bool Wise=false;
  if(!Wise) {
    ifstream ifWisdom;
    ifWisdom.open(fftw::WisdomName);
    ostringstream wisdom;
    wisdom << ifWisdom.rdbuf();
    ifWisdom.close();
    const string& s=wisdom.str();
    fftw_import_wisdom_from_string(s.c_str());
    Wise=true;
  }
}

void SaveWisdom()
{
  ofstream ofWisdom;
  ofWisdom.open(fftw::WisdomName);
  char *wisdom=fftw_export_wisdom_to_string();
  ofWisdom << wisdom;
  fftw_free(wisdom);
  ofWisdom.close();
}

fftw_plan Planner(fftw *F, Complex *in, Complex *out)
{
  LoadWisdom();
  fftw::effort |= FFTW_WISDOM_ONLY;
  fftw_plan plan=F->Plan(in,out);
  fftw::effort &= ~FFTW_WISDOM_ONLY;
  if(!plan) {
    plan=F->Plan(in,out);
    SaveWisdom();
  }
  return plan;
}

ThreadBase::ThreadBase() {threads=fftw::maxthreads;}

}

namespace fftwpp::utils {
unsigned int defaultmpithreads=1;
}