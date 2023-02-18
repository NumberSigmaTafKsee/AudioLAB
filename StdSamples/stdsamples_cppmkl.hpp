#include <new>
#include <cstddef>
#include <mkl.h>

namespace cppmkl
{
    template <typename T>
    class cppmkl_allocator
    {
        public:
        typedef T                 value_type;
        typedef value_type*       pointer;
        typedef const value_type* const_pointer;
        typedef value_type&       reference;
        typedef const value_type& const_reference;
        typedef std::size_t       size_type;
        typedef std::ptrdiff_t    difference_type;
        pointer address(reference x) const {return &x;}
        const_pointer address(const_reference x) const { return &x;}
        pointer allocate(size_type n, const_pointer=0)
        {
            void *p = MKL_malloc(n*sizeof(T), 128);
            if(!p)
            {
            throw std::bad_alloc();
            }
            return static_cast<pointer>(p);
        }
        void deallocate(pointer p, size_type)
        {
            MKL_free(p);
        }
        size_type max_size() const { 
            return static_cast<size_type>(-1) / sizeof(value_type);
        }

        void construct(pointer p, const value_type& x) { 
            new(p) value_type(x); 
        }
        void destroy(pointer p) { p->~value_type(); }

        template <class U> cppmkl_allocator(const cppmkl_allocator<U>&) {}

        template <class U> struct rebind { typedef cppmkl_allocator<U> other; };

        cppmkl_allocator() {}
        cppmkl_allocator(const cppmkl_allocator&) {}
        ~cppmkl_allocator() {}

        private:
        void operator=(const cppmkl_allocator&);
    };

    template<> 
    class cppmkl_allocator<void>
    {
        typedef void        value_type;
        typedef void*       pointer;
        typedef const void* const_pointer;

        template <class U> 
        struct rebind { typedef cppmkl_allocator<U> other; };
    };

    template<class T>
    inline bool operator==(const cppmkl_allocator<T>&, 
        const cppmkl_allocator<T>&) {
    return true;
    }

    template <class T>
    inline bool operator!=(const cppmkl_allocator<T>&, 
        const cppmkl_allocator<T>&) {
        return false;
    }

    template<typename T>
    inline const T* ptr_to_first(const boost::numeric::ublas::matrix<T>& m)
    {
    const T* ptr = &(m.data()[0]);
    return ptr;
    }
    template<typename T>
    inline T* ptr_to_first(boost::numeric::ublas::matrix<T>& m)
    {
    T* ptr = &(m.data()[0]);
    return ptr;
    }

    template<typename MV_T>
    inline const typename MV_T::value_type* ptr_to_first(const MV_T& m)
    {
    const typename MV_T::value_type* ptr = &*m.begin();
    return ptr;
    }
    template<typename MV_T>
    inline typename MV_T::value_type* ptr_to_first(MV_T& m)
    {
    typename MV_T::value_type* ptr = &*m.begin();
    return ptr;
    }

    inline void vadd(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsAdd(n, a, b, r);
    }
    inline void vadd(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdAdd(n, a, b, r);
    }
    inline void vadd(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcAdd(n, a, b, r);
    }
    inline void vadd(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzAdd(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vadd(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vadd(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vsub(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsSub(n, a, b, r);
    }
    inline void vsub(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdSub(n, a, b, r);
    }
    inline void vsub(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcSub(n, a, b, r);
    }
    inline void vsub(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzSub(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vsub(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vsub(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vsqr(const MKL_INT n, const float* a, float* r)
    {
    vsSqr(n, a, r);
    }
    inline void vsqr(const MKL_INT n, const double* a, double* r)
    {
    vdSqr(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vsqr(const VECTOR_T& a, VECTOR_T& r)
    {
    vsqr(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vmul(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsMul(n, a, b, r);
    }
    inline void vmul(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdMul(n, a, b, r);
    }
    inline void vmul(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcMul(n, a, b, r);
    }
    inline void vmul(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzMul(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vmul(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vmul(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vmulbyconj(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcMulByConj(n, a, b, r);
    }
    inline void vmulbyconj(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzMulByConj(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vmulbyconj(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vmulbyconj(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vconj(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcConj(n, a, r);
    }
    inline void vconj(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzConj(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vconj(const VECTOR_T& a, VECTOR_T& r)
    {
    vconj(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vabs(const MKL_INT n, const float* a, float* r)
    {
    vsAbs(n, a, r);
    }
    inline void vabs(const MKL_INT n, const double* a, double* r)
    {
    vdAbs(n, a, r);
    }
    inline void vabs(const MKL_INT n, const MKL_Complex8* a, float* r)
    {
    vcAbs(n, a, r);
    }
    inline void vabs(const MKL_INT n, const MKL_Complex16* a, double* r)
    {
    vzAbs(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vabs(const VECTOR_T& a, VECTOR_T& r)
    {
    vabs(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vinv(const MKL_INT n, const float* a, float* r)
    {
    vsInv(n, a, r);
    }
    inline void vinv(const MKL_INT n, const double* a, double* r)
    {
    vdInv(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vinv(const VECTOR_T& a, VECTOR_T& r)
    {
    vinv(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vdiv(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsDiv(n, a, b, r);
    }
    inline void vdiv(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdDiv(n, a, b, r);
    }
    inline void vdiv(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcDiv(n, a, b, r);
    }
    inline void vdiv(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzDiv(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vdiv(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vdiv(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vsqrt(const MKL_INT n, const float* a, float* r)
    {
    vsSqrt(n, a, r);
    }
    inline void vsqrt(const MKL_INT n, const double* a, double* r)
    {
    vdSqrt(n, a, r);
    }
    inline void vsqrt(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcSqrt(n, a, r);
    }
    inline void vsqrt(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzSqrt(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vsqrt(const VECTOR_T& a, VECTOR_T& r)
    {
    vsqrt(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vinvsqrt(const MKL_INT n, const float* a, float* r)
    {
    vsInvSqrt(n, a, r);
    }
    inline void vinvsqrt(const MKL_INT n, const double* a, double* r)
    {
    vdInvSqrt(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vinvsqrt(const VECTOR_T& a, VECTOR_T& r)
    {
    vinvsqrt(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vcbrt(const MKL_INT n, const float* a, float* r)
    {
    vsCbrt(n, a, r);
    }
    inline void vcbrt(const MKL_INT n, const double* a, double* r)
    {
    vdCbrt(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vcbrt(const VECTOR_T& a, VECTOR_T& r)
    {
    vcbrt(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vinvcbrt(const MKL_INT n, const float* a, float* r)
    {
    vsInvCbrt(n, a, r);
    }
    inline void vinvcbrt(const MKL_INT n, const double* a, double* r)
    {
    vdInvCbrt(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vinvcbrt(const VECTOR_T& a, VECTOR_T& r)
    {
    vinvcbrt(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vpow2o3(const MKL_INT n, const float* a, float* r)
    {
    vsPow2o3(n, a, r);
    }
    inline void vpow2o3(const MKL_INT n, const double* a, double* r)
    {
    vdPow2o3(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vpow2o3(const VECTOR_T& a, VECTOR_T& r)
    {
    vpow2o3(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vpow3o2(const MKL_INT n, const float* a, float* r)
    {
    vsPow3o2(n, a, r);
    }
    inline void vpow3o2(const MKL_INT n, const double* a, double* r)
    {
    vdPow3o2(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vpow3o2(const VECTOR_T& a, VECTOR_T& r)
    {
    vpow3o2(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vpow(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsPow(n, a, b, r);
    }
    inline void vpow(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdPow(n, a, b, r);
    }
    inline void vpow(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8* b, MKL_Complex8* r)
    {
    vcPow(n, a, b, r);
    }
    inline void vpow(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16* b, MKL_Complex16* r)
    {
    vzPow(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vpow(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vpow(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vpowx(const MKL_INT n, const float* a, const float b, float* r)
    {
    vsPowx(n, a, b, r);
    }
    inline void vpowx(const MKL_INT n, const double* a, const double b, double* r)
    {
    vdPowx(n, a, b, r);
    }
    inline void vpowx(const MKL_INT n, const MKL_Complex8* a, const MKL_Complex8& b, MKL_Complex8* r)
    {
    vcPowx(n, a, b, r);
    }
    inline void vpowx(const MKL_INT n, const MKL_Complex16* a, const MKL_Complex16& b, MKL_Complex16* r)
    {
    vzPowx(n, a, b, r);
    }
    template <typename SCALAR_T, typename VECTOR_T>
    inline void vpowx(const VECTOR_T& a, const SCALAR_T& b, VECTOR_T& r)
    {
    vpowx(a.size(), ptr_to_first(a), b, ptr_to_first(r));  
    }
    inline void vhypot(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsHypot(n, a, b, r);
    }
    inline void vhypot(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdHypot(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vhypot(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vhypot(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vexp(const MKL_INT n, const float* a, float* r)
    {
    vsExp(n, a, r);
    }
    inline void vexp(const MKL_INT n, const double* a, double* r)
    {
    vdExp(n, a, r);
    }
    inline void vexp(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcExp(n, a, r);
    }
    inline void vexp(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzExp(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vexp(const VECTOR_T& a, VECTOR_T& r)
    {
    vexp(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vexpm1(const MKL_INT n, const float* a, float* r)
    {
    vsExpm1(n, a, r);
    }
    inline void vexpm1(const MKL_INT n, const double* a, double* r)
    {
    vdExpm1(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vexpm1(const VECTOR_T& a, VECTOR_T& r)
    {
    vexpm1(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vexp10(const MKL_INT n, const float* a, float* r)
    {
    vsExp10(n, a, r);
    }
    inline void vexp10(const MKL_INT n, const double* a, double* r)
    {
    vdExp10(n, a, r);
    }

    template <typename VECTOR_T>
    inline void vexp10(const VECTOR_T& a, VECTOR_T& r)
    {
    vexp10(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vexp2(const MKL_INT n, const float* a, float* r)
    {
    vsExp2(n, a, r);
    }
    inline void vexp2(const MKL_INT n, const double* a, double* r)
    {
    vdExp2(n, a, r);
    }

    template <typename VECTOR_T>
    inline void vexp2(const VECTOR_T& a, VECTOR_T& r)
    {
    vexp2(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vln(const MKL_INT n, const float* a, float* r)
    {
    vsLn(n, a, r);
    }
    inline void vln(const MKL_INT n, const double* a, double* r)
    {
    vdLn(n, a, r);
    }
    inline void vln(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcLn(n, a, r);
    }
    inline void vln(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzLn(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vln(const VECTOR_T& a, VECTOR_T& r)
    {
    vln(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vlog10(const MKL_INT n, const float* a, float* r)
    {
    vsLog10(n, a, r);
    }
    inline void vlog10(const MKL_INT n, const double* a, double* r)
    {
    vdLog10(n, a, r);
    }

    template <typename VECTOR_T>
    inline void vlog10(const VECTOR_T& a, VECTOR_T& r)
    {
    vlog10(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vlog1p(const MKL_INT n, const float* a, float* r)
    {
    vsLog1p(n, a, r);
    }
    inline void vlog1p(const MKL_INT n, const double* a, double* r)
    {
    vdLog1p(n, a, r);
    }
    inline void vlog2(const MKL_INT n, const float* a, float* r)
    {
    vsLog2(n, a, r);
    }
    inline void vlog2(const MKL_INT n, const double* a, double* r)
    {
    vdLog2(n, a, r);
    }

    template <typename VECTOR_T>
    inline void vlog2(const VECTOR_T& a, VECTOR_T& r)
    {
    vlog2(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    template <typename VECTOR_T>
    inline void vlog1p(const VECTOR_T& a, VECTOR_T& r)
    {
    vlog1p(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vlogb(const MKL_INT n, const float* a, float* r)
    {
    vsLogb(n, a, r);
    }
    inline void vlogb(const MKL_INT n, const double* a, double* r)
    {
    vdLogb(n, a, r);
    }

    template <typename VECTOR_T>
    inline void vlogb(const VECTOR_T& a, VECTOR_T& r)
    {
    vlogb(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vcos(const MKL_INT n, const float* a, float* r)
    {
    vsCos(n, a, r);
    }
    inline void vcos(const MKL_INT n, const double* a, double* r)
    {
    vdCos(n, a, r);
    }
    inline void vcos(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcCos(n, a, r);
    }
    inline void vcos(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzCos(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vcos(const VECTOR_T& a, VECTOR_T& r)
    {
    vcos(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vsin(const MKL_INT n, const float* a, float* r)
    {
    vsSin(n, a, r);
    }
    inline void vsin(const MKL_INT n, const double* a, double* r)
    {
    vdSin(n, a, r);
    }
    inline void vsin(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcSin(n, a, r);
    }
    inline void vsin(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzSin(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vsin(const VECTOR_T& a, VECTOR_T& r)
    {
    vsin(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vsincos(const MKL_INT n, const float* a, float* y, float* z)
    {
    vsSinCos(n, a, y, z);
    }
    inline void vsincos(const MKL_INT n, const double* a, double* y, double* z)
    {
    vdSinCos(n, a, y, z);
    }
    template <typename VECTOR_T>
    inline void vsincos(const VECTOR_T& a, VECTOR_T& y, VECTOR_T& z)
    {
    vsincos(a.size(), ptr_to_first(a), ptr_to_first(y), ptr_to_first(z));  
    }
    inline void vCIS(const MKL_INT n, const float* a, MKL_Complex8* r)
    {
    vcCIS(n, a, r);
    }
    inline void vCIS(const MKL_INT n, const double* a, MKL_Complex16* r)
    {
    vzCIS(n, a, r);
    }
    template <typename VECTOR_T_REAL, typename VECTOR_T_COMPLEX>
    inline void vCIS(const VECTOR_T_REAL& a, VECTOR_T_COMPLEX& r)
    {
    vCIS(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtan(const MKL_INT n, const float* a, float* r)
    {
    vsTan(n, a, r);
    }
    inline void vtan(const MKL_INT n, const double* a, double* r)
    {
    vdTan(n, a, r);
    }
    inline void vtan(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcTan(n, a, r);
    }
    inline void vtan(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzTan(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vtan(const VECTOR_T& a, VECTOR_T& r)
    {
    vtan(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vacos(const MKL_INT n, const float* a, float* r)
    {
    vsAcos(n, a, r);
    }
    inline void vacos(const MKL_INT n, const double* a, double* r)
    {
    vdAcos(n, a, r);
    }
    inline void vacos(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAcos(n, a, r);
    }
    inline void vacos(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAcos(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vacos(const VECTOR_T& a, VECTOR_T& r)
    {
    vacos(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vasin(const MKL_INT n, const float* a, float* r)
    {
    vsAsin(n, a, r);
    }
    inline void vasin(const MKL_INT n, const double* a, double* r)
    {
    vdAsin(n, a, r);
    }
    inline void vasin(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAsin(n, a, r);
    }
    inline void vasin(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAsin(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vasin(const VECTOR_T& a, VECTOR_T& r)
    {
    vasin(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vatan(const MKL_INT n, const float* a, float* r)
    {
    vsAtan(n, a, r);
    }
    inline void vatan(const MKL_INT n, const double* a, double* r)
    {
    vdAtan(n, a, r);
    }
    inline void vatan(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAtan(n, a, r);
    }
    inline void vatan(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAtan(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vatan(const VECTOR_T& a, VECTOR_T& r)
    {
    vatan(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vatan2(const MKL_INT n, const float* a, const float* b, float* r)
    {
    vsAtan2(n, a, b, r);
    }
    inline void vatan2(const MKL_INT n, const double* a, const double* b, double* r)
    {
    vdAtan2(n, a, b, r);
    }
    template <typename VECTOR_T>
    inline void vatan2(const VECTOR_T& a, const VECTOR_T& b, VECTOR_T& r)
    {
    vatan2(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }
    inline void vcosh(const MKL_INT n, const float* a, float* r)
    {
    vsCosh(n, a, r);
    }
    inline void vcosh(const MKL_INT n, const double* a, double* r)
    {
    vdCosh(n, a, r);
    }
    inline void vcosh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcCosh(n, a, r);
    }
    inline void vcosh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzCosh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vcosh(const VECTOR_T& a, VECTOR_T& r)
    {
    vcosh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vsinh(const MKL_INT n, const float* a, float* r)
    {
    vsSinh(n, a, r);
    }
    inline void vsinh(const MKL_INT n, const double* a, double* r)
    {
    vdSinh(n, a, r);
    }
    inline void vsinh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcSinh(n, a, r);
    }
    inline void vsinh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzSinh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vsinh(const VECTOR_T& a, VECTOR_T& r)
    {
    vsinh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtanh(const MKL_INT n, const float* a, float* r)
    {
    vsTanh(n, a, r);
    }
    inline void vtanh(const MKL_INT n, const double* a, double* r)
    {
    vdTanh(n, a, r);
    }
    inline void vtanh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcTanh(n, a, r);
    }
    inline void vtanh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzTanh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vtanh(const VECTOR_T& a, VECTOR_T& r)
    {
    vtanh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vacosh(const MKL_INT n, const float* a, float* r)
    {
    vsAcosh(n, a, r);
    }
    inline void vacosh(const MKL_INT n, const double* a, double* r)
    {
    vdAcosh(n, a, r);
    }
    inline void vacosh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAcosh(n, a, r);
    }
    inline void vacosh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAcosh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vacosh(const VECTOR_T& a, VECTOR_T& r)
    {
    vacosh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vasinh(const MKL_INT n, const float* a, float* r)
    {
    vsAsinh(n, a, r);
    }
    inline void vasinh(const MKL_INT n, const double* a, double* r)
    {
    vdAsinh(n, a, r);
    }
    inline void vasinh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAsinh(n, a, r);
    }
    inline void vasinh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAsinh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vasinh(const VECTOR_T& a, VECTOR_T& r)
    {
    vasinh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vatanh(const MKL_INT n, const float* a, float* r)
    {
    vsAtanh(n, a, r);
    }
    inline void vatanh(const MKL_INT n, const double* a, double* r)
    {
    vdAtanh(n, a, r);
    }
    inline void vatanh(const MKL_INT n, const MKL_Complex8* a, MKL_Complex8* r)
    {
    vcAtanh(n, a, r);
    }
    inline void vatanh(const MKL_INT n, const MKL_Complex16* a, MKL_Complex16* r)
    {
    vzAtanh(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vatanh(const VECTOR_T& a, VECTOR_T& r)
    {
    vatanh(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void verf(const MKL_INT n, const float* a, float* r)
    {
    vsErf(n, a, r);
    }
    inline void verf(const MKL_INT n, const double* a, double* r)
    {
    vdErf(n, a, r);
    }
    template <typename VECTOR_T>
    inline void verf(const VECTOR_T& a, VECTOR_T& r)
    {
    verf(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void verfc(const MKL_INT n, const float* a, float* r)
    {
    vsErfc(n, a, r);
    }
    inline void verfc(const MKL_INT n, const double* a, double* r)
    {
    vdErfc(n, a, r);
    }
    template <typename VECTOR_T>
    inline void verfc(const VECTOR_T& a, VECTOR_T& r)
    {
    verfc(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vcdfnorm(const MKL_INT n, const float* a, float* r)
    {
    vsCdfNorm(n, a, r);
    }
    inline void vcdfnorm(const MKL_INT n, const double* a, double* r)
    {
    vdCdfNorm(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vcdfnorm(const VECTOR_T& a, VECTOR_T& r)
    {
    vcdfnorm(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void verfinv(const MKL_INT n, const float* a, float* r)
    {
    vsErfInv(n, a, r);
    }
    inline void verfinv(const MKL_INT n, const double* a, double* r)
    {
    vdErfInv(n, a, r);
    }
    template <typename VECTOR_T>
    inline void verfinv(const VECTOR_T& a, VECTOR_T& r)
    {
    verfinv(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void verfcinv(const MKL_INT n, const float* a, float* r)
    {
    vsErfcInv(n, a, r);
    }
    inline void verfcinv(const MKL_INT n, const double* a, double* r)
    {
    vdErfcInv(n, a, r);
    }
    template <typename VECTOR_T>
    inline void verfcinv(const VECTOR_T& a, VECTOR_T& r)
    {
    verfcinv(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vcdfnorminv(const MKL_INT n, const float* a, float* r)
    {
    vsCdfNormInv(n, a, r);
    }
    inline void vcdfnorminv(const MKL_INT n, const double* a, double* r)
    {
    vdCdfNormInv(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vcdfnorminv(const VECTOR_T& a, VECTOR_T& r)
    {
    vcdfnorminv(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vfloor(const MKL_INT n, const float* a, float* r)
    {
    vsFloor(n, a, r);
    }
    inline void vfloor(const MKL_INT n, const double* a, double* r)
    {
    vdFloor(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vfloor(const VECTOR_T& a, VECTOR_T& r)
    {
    vfloor(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vceil(const MKL_INT n, const float* a, float* r)
    {
    vsCeil(n, a, r);
    }
    inline void vceil(const MKL_INT n, const double* a, double* r)
    {
    vdCeil(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vceil(const VECTOR_T& a, VECTOR_T& r)
    {
    vceil(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtrunc(const MKL_INT n, const float* a, float* r)
    {
    vsTrunc(n, a, r);
    }
    inline void vtrunc(const MKL_INT n, const double* a, double* r)
    {
    vdTrunc(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vtrunc(const VECTOR_T& a, VECTOR_T& r)
    {
    vtrunc(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vround(const MKL_INT n, const float* a, float* r)
    {
    vsRound(n, a, r);
    }
    inline void vround(const MKL_INT n, const double* a, double* r)
    {
    vdRound(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vround(const VECTOR_T& a, VECTOR_T& r)
    {
    vround(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vnearbyint(const MKL_INT n, const float* a, float* r)
    {
    vsNearbyInt(n, a, r);
    }
    inline void vnearbyint(const MKL_INT n, const double* a, double* r)
    {
    vdNearbyInt(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vnearbyint(const VECTOR_T& a, VECTOR_T& r)
    {
    vnearbyint(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vrint(const MKL_INT n, const float* a, float* r)
    {
    vsRint(n, a, r);
    }
    inline void vrint(const MKL_INT n, const double* a, double* r)
    {
    vdRint(n, a, r);
    }
    template <typename VECTOR_T>
    inline void vrint(const VECTOR_T& a, VECTOR_T& r)
    {
    vrint(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vmodf(const MKL_INT n, const float* a, float* y, float* z)
    {
    vsModf(n, a, y, z);
    }
    inline void vmodf(const MKL_INT n, const double* a, double* y, double* z)
    {
    vdModf(n, a, y, z);
    }
    template <typename VECTOR_T>
    inline void vmodf(const VECTOR_T& a, VECTOR_T& y, VECTOR_T& z)
    {
    vmodf(a.size(), ptr_to_first(a), ptr_to_first(y), ptr_to_first(z));  
    }

    inline void varg(const MKL_INT n, const MKL_Complex8* a, float * r)
    {
    vcArg(n, a, r);
    }
    inline void varg(const MKL_INT n, const MKL_Complex16* a, double* r)
    {
    vzArg(n, a, r);
    }

    inline void vcospi(const MKL_INT n, const float* a, float* r)
    {
    vsCospi(n, a, r);
    }
    inline void vcospi(const MKL_INT n, const double* a, double* r)
    {
    vdCospi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vcospi(const VECTOR_T& a, VECTOR_T& r)
    {
    vcospi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vsinpi(const MKL_INT n, const float* a, float* r)
    {
    vsSinpi(n, a, r);
    }
    inline void vsinpi(const MKL_INT n, const double* a, double* r)
    {
    vdSinpi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vsinpi(const VECTOR_T& a, VECTOR_T& r)
    {
    vsinpi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtanpi(const MKL_INT n, const float* a, float* r)
    {
    vsTanpi(n, a, r);
    }
    inline void vtanpi(const MKL_INT n, const double* a, double* r)
    {
    vdTanpi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vtanpi(const VECTOR_T& a, VECTOR_T& r)
    {
    vtanpi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }

    inline void vacospi(const MKL_INT n, const float* a, float* r)
    {
    vsAcospi(n, a, r);
    }
    inline void vacospi(const MKL_INT n, const double* a, double* r)
    {
    vdAcospi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vacospi(const VECTOR_T& a, VECTOR_T& r)
    {
    vacospi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }

    inline void vasinpi(const MKL_INT n, const float* a, float* r)
    {
    vsAsinpi(n, a, r);
    }
    inline void vasinpi(const MKL_INT n, const double* a, double* r)
    {
    vdAsinpi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vasinpi(const VECTOR_T& a, VECTOR_T& r)
    {
    vasinpi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }

    inline void vatanpi(const MKL_INT n, const float* a, float* r)
    {
    vsAtanpi(n, a, r);
    }
    inline void vatanpi(const MKL_INT n, const double* a, double* r)
    {
    vdAtanpi(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vatanpi(const VECTOR_T& a, VECTOR_T& r)
    {
    vatanpi(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }

    inline void vatan2pi(const MKL_INT n, const float* a, const float * b, float* r)
    {
    vsAtan2pi(n, a, b, r);
    }
    inline void vatan2pi(const MKL_INT n, const double* a, double * b, double* r)
    {
    vdAtan2pi(n, a, b, r);
    }  
    template <typename VECTOR_T>
    inline void vatan2pi(const VECTOR_T& a, VECTOR_T & b, VECTOR_T& r)
    {
    vatan2pi(a.size(), ptr_to_first(a), ptr_to_first(b), ptr_to_first(r));  
    }


    inline void vcosd(const MKL_INT n, const float* a, float* r)
    {
    vsCosd(n, a, r);
    }
    inline void vcosd(const MKL_INT n, const double* a, double* r)
    {
    vdCosd(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vcosD(const VECTOR_T& a, VECTOR_T& r)
    {
    vcosd(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vsind(const MKL_INT n, const float* a, float* r)
    {
    vsSind(n, a, r);
    }
    inline void vsind(const MKL_INT n, const double* a, double* r)
    {
    vdSind(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vsind(const VECTOR_T& a, VECTOR_T& r)
    {
    vsind(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtand(const MKL_INT n, const float* a, float* r)
    {
    vsTand(n, a, r);
    }
    inline void vtand(const MKL_INT n, const double* a, double* r)
    {
    vdTand(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vtand(const VECTOR_T& a, VECTOR_T& r)
    {
    vtand(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vlgamma(const MKL_INT n, const float* a, float* r)
    {
    vsLGamma(n, a, r);
    }
    inline void vlgamma(const MKL_INT n, const double* a, double* r)
    {
    vdLGamma(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vlgamma(const VECTOR_T& a, VECTOR_T& r)
    {
    vlgamma(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    inline void vtgamma(const MKL_INT n, const float* a, float* r)
    {
    vsTGamma(n, a, r);
    }
    inline void vtgamma(const MKL_INT n, const double* a, double* r)
    {
    vdTGamma(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vtgamma(const VECTOR_T& a, VECTOR_T& r)
    {
    vtgamma(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }

    inline void vexpint1(const MKL_INT n, const float* a, float* r)
    {
    vsExpInt1(n, a, r);
    }
    inline void vexpint1(const MKL_INT n, const double* a, double* r)
    {
    vdExpInt1(n, a, r);
    }  
    template <typename VECTOR_T>
    inline void vexpint1(const VECTOR_T& a, VECTOR_T& r)
    {
    vexpint1(a.size(), ptr_to_first(a), ptr_to_first(r));  
    }
    
    /**
    * Template function wrappers for BLAS routines.
    * Vector types must have a size() function. The
    * Matrix types must have size1(), size2() functions, and value_type typedef. 
    * size1() returns the number of rows, size2() the number of columns
    *
    * The template functions here assume all Matrix types are row major. 
    * This can be generalised by adding type traits.
    */
    // in order to default alpha and beta we need something that works for 
    // both real and complex, this can (hopefully) be made much simpler:
    template<typename T>
    inline T _alpha_default()
    {
        return T();
    }
    template<>
    inline double _alpha_default<double>()
    {
        return 1.0;
    }
    template<>
    inline float _alpha_default<float>()
    {
        return 1.0;
    }
    template<>
    inline MKL_Complex8 _alpha_default<MKL_Complex8>()
    {
        static MKL_Complex8 c;
        c.real = 1.0;
        c.imag = 0.0;
        return c;

    }
    template<>
    inline MKL_Complex16 _alpha_default<MKL_Complex16>()
    {
        static MKL_Complex16 c;
        c.real = 1.0;
        c.imag = 0.0;
        return c;

    }
    template<typename T>
    inline T _beta_default()
    {
        return T();
    }
    inline float cblas_asum(const MKL_INT N, const float *X, const MKL_INT incX)
    {
        return cblas_sasum(N, X, incX);
    }
    inline float cblas_asum(const MKL_INT N, const MKL_Complex8 *X, const MKL_INT incX)
    {
    return cblas_scasum(N, static_cast<const void*>(X), incX);
    }
    inline double cblas_asum(const MKL_INT N, const double *X, const MKL_INT incX)
    { 
    return cblas_dasum(N, X, incX);
    }
    inline double cblas_asum(const MKL_INT N, const MKL_Complex16 *X, const MKL_INT incX)
    {
    return cblas_dzasum(N, static_cast<const void*>(X), incX);
    }

    //always return a double, auto conversion to float in calling code
    template <typename VECTOR_T>
    double cblas_asum(const VECTOR_T& v, const MKL_INT incX=1)
    {
    return cblas_asum(v.size()/incX, ptr_to_first(v), incX);
    }

    inline void cblas_axpy(const MKL_INT N, const float alpha, const float *X,
                    const MKL_INT incX, float *Y, const MKL_INT incY)
    {
    cblas_saxpy(N, alpha, X, incX, Y, incY);
    }
    inline void cblas_axpy(const MKL_INT N, const double alpha, const double *X,
                    const MKL_INT incX, double *Y, const MKL_INT incY)
    {
    cblas_daxpy(N, alpha, X, incX, Y, incY);
    }
    inline void cblas_axpy(const MKL_INT N, const float alpha, const MKL_Complex8 *X,
                    const MKL_INT incX, MKL_Complex8 *Y, const MKL_INT incY)
    {
    cblas_caxpy(N, static_cast<const void*>(&alpha), static_cast<const void*>(X), 
                incX, static_cast<void*>(Y), incY);
    }
    inline void cblas_axpy(const MKL_INT N, const double alpha, const MKL_Complex16 *X,
                    const MKL_INT incX, MKL_Complex16 *Y, const MKL_INT incY)
    {
    cblas_zaxpy(N, static_cast<const void *>(&alpha), static_cast<const void *>(X),
                incX, static_cast<void*>(Y), incY);
    }
    template<typename VECTOR_T>
    inline void cblas_axpy(const VECTOR_T& x, VECTOR_T& y, double a=1.0, 
                            const MKL_INT incX=1, const MKL_INT incY=1)
    {
    cblas_axpy(x.size()/incX, a, ptr_to_first(x), incX, ptr_to_first(y), incY);
    }


    inline void cblas_copy(const MKL_INT N, const float* X, const MKL_INT incX, 
                            float* Y, const MKL_INT incY)
    {
    cblas_scopy(N, X, incX, Y, incY);    
    }
    inline void cblas_copy(const MKL_INT N, const double* X, const MKL_INT incX, 
                            double* Y, const MKL_INT incY)
    {
    cblas_dcopy(N, X, incX, Y, incY);    
    }
    inline void cblas_copy(const MKL_INT N, const MKL_Complex8* X, const MKL_INT incX, 
                            MKL_Complex8* Y, const MKL_INT incY)
    {
    cblas_ccopy(N, static_cast<const void*>(X), incX, static_cast<void*>(Y), incY);    
    }
    inline void cblas_copy(const MKL_INT N, const MKL_Complex16* X, const MKL_INT incX, 
                            MKL_Complex16* Y, const MKL_INT incY)
    {
    cblas_zcopy(N, static_cast<const void*>(X), incX, static_cast<void*>(Y), incY);    
    }

    template <typename VECTOR_T>
    inline void cblas_copy(const VECTOR_T& x, VECTOR_T& y, const MKL_INT incX=1, const MKL_INT incY=1)
    {
    cblas_copy(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY);
    }

    // dot product of vectors
    inline double cblas_dot(const MKL_INT N, const float* X, const MKL_INT incX, 
                            float* Y, const MKL_INT incY)
    {
    return cblas_sdot(N, X, incX, Y, incY);    
    }
    inline double cblas_dot(const MKL_INT N, const double* X, const MKL_INT incX, 
                            double* Y, const MKL_INT incY)
    {
    return cblas_ddot(N, X, incX, Y, incY);    
    }

    template <typename VECTOR_T>
    inline double cblas_dot(const VECTOR_T& x, VECTOR_T& y, const MKL_INT incX=1, const MKL_INT incY=1)
    {
    return cblas_dot(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY);
    }


    /*
    * The ?sdot routines compute the inner product of two vectors with extended precision. Both
    * routines use extended precision accumulation of the intermediate results, but the sdsdot
    * routine outputs the final result in single precision, whereas the dsdot routine outputs the double
    * precision result. The function sdsdot also adds scalar value sb to the inner product.
    *
    * not sure I understand this one, tricky to wrap, going to use dsdot for both variants
    * for now, add sb in the wrapper
    */
    inline double cblas_sdot(const MKL_INT N, const float  *X, const MKL_INT incX,
                    const float  *Y, const MKL_INT incY)
    {
    return cblas_dsdot(N, X, incX, Y, incY);    
    }

    template <typename VECTOR_T>
    inline double cblas_sdot(const VECTOR_T& x, VECTOR_T& y, const double sb=0.0,
        const MKL_INT incX=1, const MKL_INT incY=1)
    {
    return cblas_sdot(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY) + sb;
    }

    // dot product of conjugate(x) and y, complex
    inline void cblas_dotc(const MKL_INT N, const MKL_Complex8* X, const MKL_INT incX, 
                            const MKL_Complex8* Y, const MKL_INT incY, MKL_Complex8& result)
    {
    cblas_cdotc_sub(N, static_cast<const void*>(X), incX, static_cast<const void*>(Y), incY, 
                    static_cast<void*>(&result));    
    }
    inline void cblas_dotc(const MKL_INT N, const MKL_Complex16* X, const MKL_INT incX, 
                            const MKL_Complex16* Y, const MKL_INT incY, MKL_Complex16& result)
    {
    cblas_zdotc_sub(N, static_cast<const void*>(X), incX, static_cast<const void*>(Y), incY,
                            static_cast<void*>(&result));    
    }

    template <typename VECTOR_T, typename RESULT_T>
    inline void cblas_dotc(const VECTOR_T& x, const VECTOR_T& y, RESULT_T& result, 
    const MKL_INT incX=1, const MKL_INT incY=1)
    {
    return cblas_dotc(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY, result);
    }
    // dot product of x and y, complex
    inline void cblas_dotu(const MKL_INT N, const MKL_Complex8* X, const MKL_INT incX, 
                            const MKL_Complex8* Y, const MKL_INT incY, MKL_Complex8& result)
    {
    cblas_cdotu_sub(N, static_cast<const void*>(X), incX, static_cast<const void*>(Y), incY, 
                    static_cast<void*>(&result));    
    }
    inline void cblas_dotu(const MKL_INT N, const MKL_Complex16* X, const MKL_INT incX, 
                            const MKL_Complex16* Y, const MKL_INT incY, MKL_Complex16& result)
    {
    cblas_zdotu_sub(N, static_cast<const void*>(X), incX, static_cast<const void*>(Y), incY,
                            static_cast<void*>(&result));    
    }

    template <typename VECTOR_T, typename RESULT_T>
    inline void cblas_dotu(const VECTOR_T& x, const VECTOR_T& y, RESULT_T& result, 
    const MKL_INT incX=1, const MKL_INT incY=1)
    {
    return cblas_dotu(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY, result);
    }

    // vector 2-norm
    inline double cblas_nrm2(const MKL_INT N, const float *X, const MKL_INT incX)
    {
    return cblas_snrm2(N, X, incX);
    }
    inline double cblas_nrm2(const MKL_INT N, const double *X, const MKL_INT incX)
    {
    return cblas_dnrm2(N, X, incX);
    }
    inline double cblas_nrm2(const MKL_INT N, const MKL_Complex8 *X, const MKL_INT incX)
    {
    return cblas_scnrm2(N, static_cast<const void*>(X), incX);
    }
    inline double cblas_nrm2(const MKL_INT N, const MKL_Complex16 *X, const MKL_INT incX)  
    {
    return cblas_dznrm2(N, static_cast<const void *>(X), incX);
    }
    template <typename VECTOR_T>
    inline double cblas_nrm2(const VECTOR_T& x, const MKL_INT incX=1)
    {
    return cblas_nrm2(x.size() / incX, ptr_to_first(x), incX);
    }

    // rot, rotate vectors

    inline void cblas_rot(const MKL_INT N, float* X, const MKL_INT incX,
                float *Y, const MKL_INT incY, const float c, const float s)
    {
    cblas_srot(N, X, incX, Y, incY, c, s);
    }

    inline void cblas_rot(const MKL_INT N, double* X, const MKL_INT incX,
                double *Y, const MKL_INT incY, const double c, const double s)
    {
    cblas_drot(N, X, incX, Y, incY, c, s);
    }

    inline void cblas_rot(const MKL_INT N, MKL_Complex8 *X, const MKL_INT incX,
                    MKL_Complex8 *Y, const MKL_INT incY, const float c, const float s)
    {
    cblas_csrot(N, static_cast<void*>(X), incX, 
                static_cast<void *>(Y), incY, c, s);
    }
    inline void cblas_rot(const MKL_INT N, MKL_Complex16 *X, const MKL_INT incX,
                    MKL_Complex16 *Y, const MKL_INT incY, const double c, const double s)
    {
    cblas_zdrot(N, static_cast<void*>(X), incX, 
                static_cast<void *>(Y), incY, c, s);
    }

    template <typename VECTOR_T>
    inline void cblas_rot(VECTOR_T& x, VECTOR_T& y, const double c, const double s,
                        const MKL_INT incX=1, const MKL_INT incY=1)
    {
    cblas_rot(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY, c, s);
    }

    //rotg, Givens rotation params calculation, docs are wrong on parameter types here. 
    inline void cblas_rotg(float& a, float& b, float& c, float& s)
    {
    cblas_srotg(&a, &b, &c, &s);
    }
    inline void cblas_rotg(double& a, double& b, double& c, double& s)
    {
    cblas_drotg(&a, &b, &c, &s);
    }
    inline void cblas_rotg(MKL_Complex8& a, MKL_Complex8& b, float& c, MKL_Complex8& s)
    {
    cblas_crotg(static_cast<void *>(&a), static_cast<const void *>(&b), &c, static_cast<void*>(&s));
    }
    inline void cblas_rotg(MKL_Complex16& a, MKL_Complex16& b, double& c, MKL_Complex16& s)
    {
    cblas_zrotg(static_cast<void *>(&a), static_cast<const void *>(&b), &c, static_cast<void*>(&s));
    }

    //TODO rotm
    //TODO rotmg

    //?scal: multiply vector by scalar
    inline void cblas_scal(const MKL_INT N, const float alpha, float *X, const MKL_INT incX)
    {
    cblas_sscal(N, alpha, X, incX);
    }
    inline void cblas_scal(const MKL_INT N, const double alpha, double *X, const MKL_INT incX)
    {
    cblas_dscal(N, alpha, X, incX);
    }
    inline void cblas_scal(const MKL_INT N, const MKL_Complex8& alpha, MKL_Complex8 *X, const MKL_INT incX)
    {
    cblas_cscal(N, static_cast<const void *>(&alpha), static_cast<void *>(X), incX);
    }
    inline void cblas_scal(const MKL_INT N, const MKL_Complex16& alpha, MKL_Complex16 *X, const MKL_INT incX)
    {
    cblas_zscal(N, static_cast<const void *>(&alpha), static_cast<void *>(X), incX);
    }
    inline void cblas_scal(const MKL_INT N, const float alpha, MKL_Complex8* X, const MKL_INT incX)
    {
    cblas_csscal(N, alpha, static_cast<void *>(X), incX);
    }
    inline void cblas_scal(const MKL_INT N, const double alpha, MKL_Complex16 *X, const MKL_INT incX)
    {
    cblas_zdscal(N, alpha, static_cast<void *>(X), incX);
    }
    template <typename VECTOR_T, typename SCALAR_T>
    inline void cblas_scal(VECTOR_T& x, const SCALAR_T& a, const MKL_INT incX=1)
    {
    cblas_scal(x.size()/incX, a, ptr_to_first(x), incX);
    }

    // swap: swap vectors

    inline void cblas_swap(const MKL_INT N, float *X, const MKL_INT incX,
                    float *Y, const MKL_INT incY)
    {
    cblas_sswap(N, X,incX, Y, incY);
    }
    inline void cblas_swap(const MKL_INT N, double *X, const MKL_INT incX,
                    double *Y, const MKL_INT incY)
    { 
    cblas_dswap(N, X, incX, Y, incY);
    }
    inline void cblas_swap(const MKL_INT N, MKL_Complex8 *X, const MKL_INT incX,
                    MKL_Complex8 *Y, const MKL_INT incY) 
    {
    cblas_cswap(N, static_cast<void *>(X), incX, static_cast<void *>(Y), incY);
    }
    inline void cblas_swap(const MKL_INT N, MKL_Complex16 *X, const MKL_INT incX,
                    MKL_Complex16 *Y, const MKL_INT incY) 
    {
    cblas_zswap(N, static_cast<void *>(X), incX, static_cast<void *>(Y), incY);
    }

    template<typename VECTOR_T>
    inline void cblas_swap(VECTOR_T& x, VECTOR_T& y, const MKL_INT incX=1, const MKL_INT incY=1)
    {
    cblas_swap(x.size()/incX, ptr_to_first(x), incX, ptr_to_first(y), incY);
    }

    //iamax: find max element in vector
    inline CBLAS_INDEX cblas_iamax(const MKL_INT N, const float  *X, const MKL_INT incX)
    {
    return cblas_isamax(N, X, incX);
    }
    inline CBLAS_INDEX cblas_iamax(const MKL_INT N, const double *X, const MKL_INT incX)
    { 
    return cblas_idamax(N, X, incX);
    }
    inline CBLAS_INDEX cblas_iamax(const MKL_INT N, const MKL_Complex8 *X, const MKL_INT incX)
    {
    return cblas_icamax(N, static_cast<const void *>(X), incX);
    }
    inline CBLAS_INDEX cblas_iamax(const MKL_INT N, const MKL_Complex16 *X, const MKL_INT incX)
    {
    return cblas_izamax(N, static_cast<const void*>(X), incX);
    }

    template<typename VECTOR_T>
    inline CBLAS_INDEX cblas_iamax(VECTOR_T& x, const MKL_INT incX=1)
    {
    return cblas_iamax(x.size()/incX, ptr_to_first(x), incX);
    } 
    //iamin: find min element in vector
    inline CBLAS_INDEX cblas_iamin(const MKL_INT N, const float  *X, const MKL_INT incX)
    {
    return cblas_isamin(N, X, incX);
    }
    inline CBLAS_INDEX cblas_iamin(const MKL_INT N, const double *X, const MKL_INT incX)
    { 
    return cblas_idamin(N, X, incX);
    }
    inline CBLAS_INDEX cblas_iamin(const MKL_INT N, const MKL_Complex8 *X, const MKL_INT incX)
    {
    return cblas_icamin(N, static_cast<const void *>(X), incX);
    }
    inline CBLAS_INDEX cblas_iamin(const MKL_INT N, const MKL_Complex16 *X, const MKL_INT incX)
    {
    return cblas_izamin(N, static_cast<const void*>(X), incX);
    }

    template<typename VECTOR_T>
    inline CBLAS_INDEX cblas_iamin(VECTOR_T& x, const MKL_INT incX=1)
    {
    return cblas_iamin(x.size()/incX, ptr_to_first(x), incX);
    } 

    inline double cblas_cabs1(const MKL_Complex16& z)
    {
    return cblas_dcabs1(static_cast<const void*>(&z));
    }

    /**
    * CBLAS BLAS Level 2 wrappers
    */

    //gemv Matrix vector product
    inline void cblas_gemv(const  CBLAS_ORDER order,
                    const  CBLAS_TRANSPOSE TransA, const MKL_INT M, const MKL_INT N,
                    const double alpha, const float *A, const MKL_INT lda,
                    const float *X, const MKL_INT incX, const double beta,
                    float *Y, const MKL_INT incY)
    {
    cblas_sgemv(order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
    }
    inline void cblas_gemv(const  CBLAS_ORDER order,
                    const  CBLAS_TRANSPOSE TransA, const MKL_INT M, const MKL_INT N,
                    const double alpha, const double *A, const MKL_INT lda,
                    const double *X, const MKL_INT incX, const double beta,
                    double *Y, const MKL_INT incY)
    {
    cblas_dgemv(order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
    }
    inline void cblas_gemv(const  CBLAS_ORDER order,
                    const  CBLAS_TRANSPOSE TransA, const MKL_INT M, const MKL_INT N,
                    const MKL_Complex8& alpha, const MKL_Complex8 *A, const MKL_INT lda,
                    const MKL_Complex8 *X, const MKL_INT incX, const MKL_Complex8& beta,
                    MKL_Complex8 *Y, const MKL_INT incY)
    {
    cblas_cgemv(order, TransA, M, N, static_cast<const void*>(&alpha), static_cast<const void*>(A), lda, 
                static_cast<const void*>(X), incX, static_cast<const void*>(&beta), static_cast<void*>(Y), incY);
    }
    inline void cblas_gemv(const  CBLAS_ORDER order,
                    const  CBLAS_TRANSPOSE TransA, const MKL_INT M, const MKL_INT N,
                    const MKL_Complex16& alpha, const MKL_Complex16 *A, const MKL_INT lda,
                    const MKL_Complex16 *X, const MKL_INT incX, const MKL_Complex16& beta,
                    MKL_Complex16 *Y, const MKL_INT incY)
    {
    cblas_zgemv(order, TransA, M, N, static_cast<const void*>(&alpha), static_cast<const void*>(A), lda, 
                static_cast<const void*>(X), incX, static_cast<const void*>(&beta), static_cast<void*>(Y), incY);
    }

    template<typename MATRIX_T, typename VECTOR_T, typename SCALAR_T>
    inline void cblas_gemv(const MATRIX_T& A, const VECTOR_T& x, VECTOR_T& y,
                    const SCALAR_T& alpha, 
                    const SCALAR_T& beta, 
                    const CBLAS_TRANSPOSE TransA=CblasNoTrans,
                    const MKL_INT incX=1, const MKL_INT incY=1)
    {
    CBLAS_ORDER order = CblasRowMajor;
    bool doTransA = TransA != CblasNoTrans;
    const size_t opA_row_count = doTransA == false ? A.size1() : A.size2();//rows of op(A)
    const size_t opA_col_count = doTransA == false ? A.size2() : A.size1();//cols of op(A)
    const MKL_INT M = A.size1();
    const MKL_INT N = A.size2();
    const MKL_INT lda = doTransA == false ? opA_col_count : opA_row_count;
    cblas_gemv(order, TransA, M, N, alpha, ptr_to_first(A), lda, ptr_to_first(x), incX, beta, ptr_to_first(y), incY); 
    }
    //ger rank 1 update of a general matrix
    //gerc rank 1 update of a conjugated general matrix


    /** cppmkl wrapper for cblas_dgemm, BLAS level 3 matrix multiplication
    * MATRIX_T can be any row-major matrix type that has functions size1(), size2() 
    * giving row count and column count respectively 
    */
    inline void cblas_gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                    const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                    const MKL_INT K, const double alpha, const double *A,
                    const MKL_INT lda, const double *B, const MKL_INT ldb,
                    const double beta, double *C, const MKL_INT ldc)
    {
    cblas_dgemm(Order, TransA, TransB, M, N,K, alpha, A, lda, B, ldb,
                beta, C, ldc);
    }
    inline void cblas_gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                    const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                    const MKL_INT K, const float alpha, const float *A,
                    const MKL_INT lda, const float *B, const MKL_INT ldb,
                    const float beta, float *C, const MKL_INT ldc)
    {
    cblas_sgemm(Order, TransA, TransB, M, N,K, alpha, A, lda, B, ldb,
                beta, C, ldc);
    }
    inline void cblas_gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                    const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                    const MKL_INT K, const MKL_Complex8 alpha, const MKL_Complex8 *A,
                    const MKL_INT lda, const MKL_Complex8 *B, const MKL_INT ldb,
                    const MKL_Complex8 beta, MKL_Complex8 *C, const MKL_INT ldc)
    {
    cblas_cgemm(Order, TransA, TransB, M, N,K, 
        static_cast<const void*>(&alpha), 
        static_cast<const void*>(A), lda, 
        static_cast<const void*>(B), ldb,
        static_cast<const void*>(&beta), 
        static_cast<void*>(C), ldc);
    }
    inline void cblas_gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                    const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                    const MKL_INT K, const MKL_Complex16 alpha, const MKL_Complex16 *A,
                    const MKL_INT lda, const MKL_Complex16 *B, const MKL_INT ldb,
                    const MKL_Complex16 beta, MKL_Complex16 *C, const MKL_INT ldc)
    {
    cblas_zgemm(Order, TransA, TransB, M, N,K, 
        static_cast<const void*>(&alpha), 
        static_cast<const void*>(A), lda, 
        static_cast<const void*>(B), ldb,
        static_cast<const void*>(&beta), 
        static_cast<void*>(C), ldc);
    }
    template <typename MATRIX_T>
    inline void cblas_gemm(const MATRIX_T& A, const MATRIX_T& B, MATRIX_T& C,
        const CBLAS_TRANSPOSE TransA=CblasNoTrans, const CBLAS_TRANSPOSE TransB=CblasNoTrans,
        const typename MATRIX_T::value_type alpha=_alpha_default<typename MATRIX_T::value_type>(), 
        const typename MATRIX_T::value_type beta=_beta_default<typename MATRIX_T::value_type>())
    {
        CBLAS_ORDER order = CblasRowMajor;
        bool doTransA = TransA != CblasNoTrans;
        bool doTransB = TransB != CblasNoTrans;
        const size_t opA_row_count = doTransA == false ? A.size1() : A.size2();//rows of op(A)
        const size_t opB_col_count = doTransB == false ? B.size2() : B.size1();//cols of op(B) 
        const size_t opA_col_count = doTransA == false ? A.size2() : A.size1();//cols of op(A)
        //const size_t opB_row_count = doTransB == false ? B.size1() : B.size2(); //rows of op(B)
        const MKL_INT lda = doTransA == false ? opA_col_count : opA_row_count;
        const MKL_INT ldb = doTransB == false ? opB_col_count : opA_col_count;
        const MKL_INT ldc = opB_col_count; 
        //opB_row_count;
        //assert(opA_col_count == opB_row_count);
        assert(C.size1() == opA_row_count);
        assert(C.size2() == opB_col_count);
        // call the appropriate overloaded cblas_gemm function based on the type of X
        cblas_gemm(order, TransA, TransB, opA_row_count, opB_col_count, opA_col_count, alpha, ptr_to_first(A), lda, ptr_to_first(B), ldb, beta, ptr_to_first(C), ldc);     
    }
}



