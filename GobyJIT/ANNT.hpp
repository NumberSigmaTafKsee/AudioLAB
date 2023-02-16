# pragma once

#include <limits>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <cstring>
#include <random>
#include <type_traits>
#include <chrono>
#include <cassert>
#include <algorithm>
#include <thread>

using namespace std;
using namespace std::chrono;
#ifdef _MSC_VER
    #include <intrin.h>
#elif __GNUC__
    #include <x86intrin.h>
    #include <cpuid.h>
#endif

// Enable support of SSE/SSE2 instructions set
#define ANNT_USE_SSE

// Enable support of AVX instructions set
#define ANNT_USE_AVX

// Enable Open MP usage for loops parallelization
#define ANNT_USE_OMP

namespace ANNT {

// Allocate/free aligned memory
void* AlignedAlloc( std::size_t align, std::size_t size );
void AlignedFree( void* ptr );
    
// Aligned allocator for standard containers
template <typename T, std::size_t Alignment>
class XAlignedAllocator
{
public:
    // Typedefs
    typedef T               value_type;
    typedef T*              pointer;
    typedef const T*        const_pointer;
    typedef T&              reference;
    typedef const T&        const_reference;
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;

public:
    // Convert an allocator<T> to allocator<U>
    template <typename U>
    struct rebind
    {
        typedef XAlignedAllocator<U, Alignment> other;
    };

public:
    XAlignedAllocator( ) { }

    template <typename U>
    XAlignedAllocator( const XAlignedAllocator<U, Alignment> & ) { }

    // Address
    inline pointer address( reference value ) const
    {
        return std::addressof( value );
    }
    inline const_pointer address( const_reference value ) const
    {
        return std::addressof(value);
    }

    // Memory allocation
    inline pointer allocate( size_type size, const void* = nullptr )
    {
        void* p = AlignedAlloc( Alignment, sizeof( T ) * size );

        if ( ( p == nullptr ) && ( size > 0 ) )
        {
            throw std::bad_alloc( );
        }

        return static_cast<pointer>( p );
    }
    inline void deallocate( pointer ptr, size_type )
    {
        AlignedFree( ptr );
    }

    // Size
    inline size_type max_size( ) const
    {
        return ~static_cast<std::size_t>( 0 ) / sizeof( T );
    }

    // Construction/destruction
    inline void construct( pointer ptr, const_reference value )
    {
        ::new ( ptr ) value_type( value );
    }
    inline void destroy( pointer ptr )
    {
        if ( ptr )
        {
            ptr->~value_type( );
        }
    }

    inline bool operator==( const XAlignedAllocator& ) { return true; }
    inline bool operator!=( const XAlignedAllocator& rhs ) { return !operator==( rhs ); }
};

template <typename T, std::size_t Alignment>
inline bool operator==( const XAlignedAllocator<T, Alignment> &,
                        const XAlignedAllocator<T, Alignment> & )
{
    return true;
}

template <typename T, std::size_t Alignment>
inline bool operator!=( const XAlignedAllocator<T, Alignment> &,
                        const XAlignedAllocator<T, Alignment> & )
{
    return false;
}

// Numeric type used for neural network's data/callculations
// (weights, biases, errors, gradients, parameters, etc.)
#ifdef ANNT_USE_DOUBLE
typedef double float_t;
#else
typedef float  float_t;
#endif

// Vector type to use for network's input/output/error/gradient flow.
// 32 bytes aligned to enable SIMD operations on those.
typedef std::vector<float_t, XAlignedAllocator<float_t, 32>> fvector_t;

// Vector type with unsigned integers (size_t) as elements
typedef std::vector<size_t> uvector_t;

// Border handling modes for convolution and pooling
enum class BorderMode
{
    Valid,  // Output is smaller than input, since convolution is only computed
            // where input and filter fully overlap.

    Same    // Output is of the same size as input. To get this input is padded.
};

// Modes of selecting training samples into batches while running training epch.
enum class EpochSelectionMode
{
    Sequential,     // Samples are not shuffled and are chosen sequentially one after another in the provided order.

    RandomPick,     // Samples are not shuffled (order is kept), but individual items are chosed randomly into batches.

    Shuffle,        // Training samples are shuffled at the start of each epoch. Then chosen sequentially into batches.
};

// A value to represent missing connection (between inputs/outputs, neurons, layers, etc)
static const size_t ANNT_NOT_CONNECTED = std::numeric_limits<size_t>::max( );

// Macro to suppress warnings caused by unreferenced parameter
#define ANNT_UNREFERENCED_PARAMETER(param) (void)param

// Allocate aligned memory
void* AlignedAlloc( std::size_t align, std::size_t size )
{
#if defined(_MSC_VER)
    return ::_aligned_malloc( size, align );
#elif defined(__MINGW32__)
    return ::_mm_malloc( size, align );
#else  // posix assumed
    void* p;

    if ( ::posix_memalign( &p, align, size ) != 0 )
    {
        p = 0;
    }

    return p;
#endif
}

// Free aligned memory
void AlignedFree( void* ptr )
{
#if defined(_MSC_VER)
    ::_aligned_free( ptr );
#elif defined(__MINGW32__)
    ::_mm_free( ptr );
#else
    ::free( ptr );
#endif
}

// Collection of tools to encode/decode data to/from formats expected/produced by ANNs
class XDataEncodingTools
{
private:
    XDataEncodingTools( );

public:
    // Encodes single class/label using one-hot encoding - a vector of all zeros except the one element set to 1,
    // which index corresponds to the class value
    static fvector_t OneHotEncoding( size_t label, size_t labelsCount );

    // Encodes a vector of labels using one-hot encoding
    static std::vector<fvector_t> OneHotEncoding( const uvector_t& labels, size_t labelsCount );

    // Returns index of the maximum element in the specified vector
    static size_t MaxIndex( const fvector_t& vec );

    // Pads the specified 2D input (although it can be of certain depth) with the specified value
    static void AddPadding2d( const fvector_t& src, fvector_t& dst,
                              size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                              size_t depth, float_t padValue );
    static void AddPadding2d( const float_t* src, float_t* dst,
                              size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                              size_t depth, float_t padValue );

    // Removes padding from the specified 2D input (although it can be of certain depth)
    static void RemovePadding2d( const fvector_t& src, fvector_t& dst,
                                 size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                 size_t depth );
    static void RemovePadding2d( const float_t* src, float_t* dst,
                                 size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                 size_t depth );

    // Builds input to output index mapping for pooling operator - one to one mapping
    static uvector_t BuildPoolingInToOutMap( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                                             size_t poolSizeX, size_t poolSizeY,
                                             size_t horizontalStep, size_t verticalStep,
                                             BorderMode borderMode = BorderMode::Valid );

    // Builds output index to input indexes mapping for pooling operator - 1 to many mapping
    static std::vector<uvector_t> BuildPoolingOutToInMap( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                                                          size_t poolSizeX, size_t poolSizeY,
                                                          size_t horizontalStep, size_t verticalStep,
                                                          BorderMode borderMode = BorderMode::Valid );
};

// Encodes single class using one-hot encoding - a vector of all zeros except the one element set to 1,
// which index corresponds to the class value
fvector_t XDataEncodingTools::OneHotEncoding( size_t label, size_t labelsCount )
{
    fvector_t encodedClass( labelsCount, 0 );

    encodedClass[label] = 1;

    return encodedClass;
}

// Encodes a vector of classes using one-hot encoding
vector<fvector_t> XDataEncodingTools::OneHotEncoding( const uvector_t& labels, size_t labelsCount )
{
    vector<fvector_t> encodedClasses( labels.size( ) );

    for ( size_t i = 0; i < labels.size( ); i++ )
    {
        encodedClasses[i] = fvector_t( labelsCount, 0 );
        encodedClasses[i][labels[i]] = 1;
    }

    return encodedClasses;
}

// Returns index of the maximum element in the specified vector
size_t XDataEncodingTools::MaxIndex( const fvector_t& vec )
{
    size_t maxIndex = 0;
    auto   maxValue = vec[0];

    for ( size_t i = 1, n = vec.size( ); i < n; i++ )
    {
        if ( vec[i] > maxValue )
        {
            maxValue = vec[i];
            maxIndex = i;
        }
    }

    return maxIndex;
}

// Pads the specified 2D input (although it can be of certain depth) with the specified value
void XDataEncodingTools::AddPadding2d( const fvector_t& src, fvector_t& dst,
                                       size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                       size_t depth, float_t padValue )
{
    size_t dstSize = dstWidth * dstHeight * depth;

    if ( dst.size( ) != dstSize )
    {
        dst.resize( dstSize );
    }

    AddPadding2d( src.data( ), dst.data( ), srcWidth, srcHeight, dstWidth, dstHeight, depth, padValue );
}
void XDataEncodingTools::AddPadding2d( const float_t* src, float* dst,
                                       size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                       size_t depth, float_t padValue )
{
    if ( ( dstWidth >= srcWidth ) && ( dstHeight >= srcHeight ) )
    {
        size_t padWidth  = dstWidth  - srcWidth;
        size_t padHeight = dstHeight - srcHeight;
        // For even pad width/height it is distributed equally on each side.
        // However for odd value, padding goes first to right/bottom sides.
        size_t leftPad   = padWidth >> 1;
        size_t rightPad  = padWidth - leftPad;
        size_t topPad    = padHeight >> 1;
        size_t bottomPad = padHeight - topPad;

        const float* srcPtr = src;
        float*       dstPtr = dst;

        for ( size_t d = 0; d < depth; d++ )
        {
            // top padding
            for ( size_t y = 0; y < topPad; y++ )
            {
                for ( size_t x = 0; x < dstWidth; x++, dstPtr++ )
                {
                    *dstPtr = padValue;
                }
            }

            for ( size_t y = 0; y < srcHeight; y++ )
            {
                // left padding
                for ( size_t x = 0; x < leftPad; x++, dstPtr++ )
                {
                    *dstPtr = padValue;
                }

                // copying the source
                for ( size_t x = 0; x < srcWidth; x++, srcPtr++, dstPtr++ )
                {
                    *dstPtr = *srcPtr;
                }

                // right padding
                for ( size_t x = 0; x < rightPad; x++, dstPtr++ )
                {
                    *dstPtr = padValue;
                }
            }

            // bottom padding
            for ( size_t y = 0; y < bottomPad; y++ )
            {
                for ( size_t x = 0; x < dstWidth; x++, dstPtr++ )
                {
                    *dstPtr = padValue;
                }
            }
        }
    }
}

// Removes padding from the specified 2D input
void XDataEncodingTools::RemovePadding2d( const fvector_t& src, fvector_t& dst,
                                          size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                          size_t depth )
{
    size_t dstSize = dstWidth * dstHeight * depth;

    if ( dst.size( ) != dstSize )
    {
        dst.resize( dstSize );
    }

    RemovePadding2d( src.data( ), dst.data( ), srcWidth, srcHeight, dstWidth, dstHeight, depth );
}
void XDataEncodingTools::RemovePadding2d( const float_t* src, float_t* dst,
                                          size_t srcWidth, size_t srcHeight, size_t dstWidth, size_t dstHeight,
                                          size_t depth )
{
    if ( ( dstWidth <= srcWidth ) && ( dstHeight <= srcHeight ) )
    {
        size_t padWidth  = srcWidth  - dstWidth;
        size_t padHeight = srcHeight - dstHeight;
        // For even pad width/height it is distributed equally on each side.
        // However for odd value, padding goes first to right/bottom sides.
        size_t leftPad   = padWidth >> 1;
        size_t rightPad  = padWidth - leftPad;
        size_t topPad    = padHeight >> 1;
        size_t bottomPad = padHeight - topPad;

        topPad    *= srcWidth;
        bottomPad *= srcWidth;

        const float* srcPtr = src;
        float*       dstPtr = dst;

        for ( size_t d = 0; d < depth; d++ )
        {
            // skip top padding
            srcPtr += topPad;

            for ( size_t y = 0; y < dstHeight; y++ )
            {
                // skip left left padding
                srcPtr += leftPad;

                // copying the source
                for ( size_t x = 0; x < dstWidth; x++, srcPtr++, dstPtr++ )
                {
                    *dstPtr = *srcPtr;
                }

                // skip right padding
                srcPtr += rightPad;
            }

            // skip bottom padding
            srcPtr += bottomPad;
        }
    }
}

// Builds input to output index mapping for pooling operator - one to one mapping
uvector_t XDataEncodingTools::BuildPoolingInToOutMap( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                                                      size_t poolSizeX, size_t poolSizeY,
                                                      size_t horizontalStep, size_t verticalStep,
                                                      BorderMode borderMode )
{
    size_t padWidth    = 0;
    size_t padHeight   = 0;
    size_t leftPad     = 0;
    size_t topPad      = 0;

    if ( borderMode == BorderMode::Same )
    {
        padWidth     = poolSizeX - 1;
        padHeight    = poolSizeY - 1;
        leftPad      = padWidth  / 2;
        topPad       = padHeight / 2;
    }

    // calculation of output width/height as:
    //   outSize = ( inSize - kernelSize + padSize ) / step + 1
    size_t outputWidth  = ( inputWidth  - poolSizeX + padWidth )  / horizontalStep + 1;
    size_t outputHeight = ( inputHeight - poolSizeY + padHeight ) / verticalStep   + 1;

    size_t inputsCount  = inputWidth * inputHeight * inputDepth;

    // build the map providing output index for the given input index
    uvector_t inToOutMap = uvector_t( inputsCount );

    std::fill( inToOutMap.begin( ), inToOutMap.end( ), ANNT_NOT_CONNECTED );

    for ( size_t depthIndex = 0, outputIndex = 0; depthIndex < inputDepth; depthIndex++ )
    {
        for ( size_t outY = 0, inY = 0; outY < outputHeight; outY++, inY += verticalStep )
        {
            size_t inRowIndex = ( inY + depthIndex * inputHeight ) * inputWidth;

            for ( size_t outX = 0, inX = 0; outX < outputWidth; outX++, inX += horizontalStep, outputIndex++ )
            {
                size_t inStartIndex = inRowIndex + inX;

                for ( size_t poolY = 0, i = 0; poolY < poolSizeY; poolY++ )
                {
                    if ( ( inY + poolY >= topPad ) &&
                         ( inY + poolY <  topPad + inputHeight ) )
                    {
                        for ( size_t poolX = 0; poolX < poolSizeX; poolX++, i++ )
                        {
                            if ( ( inX + poolX >= leftPad ) &&
                                 ( inX + poolX <  leftPad + inputWidth ) )
                            {
                                size_t inputIndex = inStartIndex + ( poolY - topPad ) * inputWidth + poolX - leftPad;

                                inToOutMap[inputIndex] = outputIndex;
                            }
                        }
                    }
                }
            }
        }
    }

    return inToOutMap;
}

// Builds output index to input indexes mapping for pooling operator - 1 to many mapping
vector<uvector_t> XDataEncodingTools::BuildPoolingOutToInMap( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                                                              size_t poolSizeX, size_t poolSizeY,
                                                              size_t horizontalStep, size_t verticalStep,
                                                              BorderMode borderMode )
{
    size_t padWidth    = 0;
    size_t padHeight   = 0;
    size_t leftPad     = 0;
    size_t topPad      = 0;

    if ( borderMode == BorderMode::Same )
    {
        padWidth     = poolSizeX - 1;
        padHeight    = poolSizeY - 1;
        leftPad      = padWidth  / 2;
        topPad       = padHeight / 2;
    }

    // calculation of output width/height as:
    //   outSize = ( inSize - kernelSize + padSize ) / step + 1
    size_t outputWidth  = ( inputWidth  - poolSizeX + padWidth  ) / horizontalStep + 1;
    size_t outputHeight = ( inputHeight - poolSizeY + padHeight ) / verticalStep   + 1;
    size_t outputsCount = outputWidth * outputHeight * inputDepth;

    vector<uvector_t> outToInMap = vector<uvector_t>( outputsCount );

    for ( size_t depthIndex = 0, outputIndex = 0; depthIndex < inputDepth; depthIndex++ )
    {
        for ( size_t outY = 0, inY = 0; outY < outputHeight; outY++, inY += verticalStep )
        {
            size_t inRowIndex = ( inY + depthIndex * inputHeight ) * inputWidth;

            for ( size_t outX = 0, inX = 0; outX < outputWidth; outX++, inX += horizontalStep, outputIndex++ )
            {
                std::vector<size_t>& outputMap    = outToInMap[outputIndex];
                size_t               inStartIndex = inRowIndex + inX;

                for ( size_t poolY = 0, i = 0; poolY < poolSizeY; poolY++ )
                {
                    if ( ( inY + poolY >= topPad ) &&
                         ( inY + poolY <  topPad + inputHeight ) )
                    {
                        for ( size_t poolX = 0; poolX < poolSizeX; poolX++, i++ )
                        {
                            if ( ( inX + poolX >= leftPad ) &&
                                 ( inX + poolX <  leftPad + inputWidth ) )
                            {
                                size_t inputIndex = inStartIndex + ( poolY - topPad ) * inputWidth + poolX - leftPad;

                                outputMap.push_back( inputIndex );
                            }
                        }
                    }
                }
            }
        }
    }

    return outToInMap;
}


// Interface for some common operations performed on vectors
class IVectorTools
{
public:
    virtual ~IVectorTools( ) { }

    // Check if the implementation of vector tools is available on the current system
    virtual bool IsAvailable( ) const = 0;

    // Add two vectors: dst[i] += src[i]
    virtual void Add( const float*  src, float*  dst, size_t size ) const = 0;
    virtual void Add( const double* src, double* dst, size_t size ) const = 0;

    // Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
    virtual void Mul( const float*  src, float*  dst, size_t size ) const = 0;
    virtual void Mul( const double* src, double* dst, size_t size ) const = 0;

    // Dot product of two vectors: sum( vec1[i] * vec2[i] )
    virtual float  Dot( const float*  vec1, const float*  vec2, size_t size ) const = 0;
    virtual double Dot( const double* vec1, const double* vec2, size_t size ) const = 0;

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    virtual void Max( const float*  src, float  alpha, float*  dst, size_t size ) const = 0;
    virtual void Max( const double* src, double alpha, double* dst, size_t size ) const = 0;
};

class XVectorize
{
private:
    XVectorize( );

public:
    // Add two vectors: dst[i] += src[i]
    template <typename T> static inline void Add( const T* src, T* dst, size_t size )
    {
        mVectorTools->Add( src, dst, size );
    }

    // Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
    template <typename T> static inline void Mul( const T* src, T* dst, size_t size )
    {
        mVectorTools->Mul( src, dst, size );
    }

    // Dot product of two vectors: sum( vec1[i] * vec2[i] )
    template <typename T> static inline T Dot( const T* vec1, const T* vec2, size_t size )
    {
        return mVectorTools->Dot( vec1, vec2, size );
    }

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    template <typename T> static inline void Max( const T* src, T alpha, T* dst, size_t size )
    {
        mVectorTools->Max( src, alpha, dst, size );
    }

private:

    static IVectorTools* mVectorTools;
};

// Implementation of common vector routines using AVX instructions
class XAvxVectorTools : public IVectorTools
{
public:
    // Check if the implementation of vector tools is available on the current system
    bool IsAvailable( ) const override;

    // Add two vectors: dst[i] += src[i]
    void Add( const float*  src, float*  dst, size_t size ) const override;
    void Add( const double* src, double* dst, size_t size ) const override;

    // Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
    void Mul( const float*  src, float*  dst, size_t size ) const override;
    void Mul( const double* src, double* dst, size_t size ) const override;

    // Dot product of two vectors: sum( vec1[i] * vec2[i] )
    float  Dot( const float*  vec1, const float*  vec2, size_t size ) const override;
    double Dot( const double* vec1, const double* vec2, size_t size ) const override;

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    void Max( const float*  src, float  alpha, float*  dst, size_t size ) const override;
    void Max( const double* src, double alpha, double* dst, size_t size ) const override;
};

// Set of functions providing some CPU related information
class XCpu
{
private:
    XCpu( );

public:
    enum FeatureRegisters
    {
        Reg_EAX = 0,
        Reg_EBX = 1,
        Reg_ECX = 2,
        Reg_EDX = 3
    };

    // Some of the CPUID flags to check for for available instruction sets
    enum EbxFlags
    {
        Flag_AVX2   = 1 << 5,
    };

    enum EcxFlags
    {
        Flag_SSE3   = 1,
        Flag_SSSE3  = 1 << 9,
        Flag_SSE4_1 = 1 << 19,
        Flag_SSE4_2 = 1 << 20,
        Flag_AVX    = 1 << 28,
    };

    enum EdxFlags
    {
        Flag_MMX    = 1 << 24,
        Flag_SSE    = 1 << 25,
        Flag_SSE2   = 1 << 26,
    };

public:
    // Provide CPU ID - 4 32-bit registers describing CPU features
    static void CpuId( uint32_t& eax, uint32_t& ebx, uint32_t& ecx, uint32_t& edx );

    // Check if the particular feature is support by the CPU
    static bool IsFeatureSupported( uint32_t reg, uint32_t flag );

    // Get number of CPU cores provided by the system
    static uint32_t CoresCount( );
};

// Provide CPU ID - 4 32-bit registers describing CPU features
void XCpu::CpuId( uint32_t& eax, uint32_t& ebx, uint32_t& ecx, uint32_t& edx )
{
#ifdef _MSC_VER
    int cpuInfo[4];

    __cpuid( cpuInfo, 1 );

    eax = static_cast<uint32_t>( cpuInfo[0] );
    ebx = static_cast<uint32_t>( cpuInfo[1] );
    ecx = static_cast<uint32_t>( cpuInfo[2] );
    edx = static_cast<uint32_t>( cpuInfo[3] );
#elif __GNUC__
    __cpuid( 1, eax, ebx, ecx, edx );
#endif
}

// Check if the particular feature is support by the CPU
bool XCpu::IsFeatureSupported( uint32_t reg, uint32_t flag )
{
    uint32_t eax = 0, ebx = 0, ecx = 0, edx = 0;
    bool     ret = false;

    CpuId( eax, ebx, ecx, edx );

    switch ( reg )
    {
    case Reg_EAX:
        ret = ( ( eax & flag ) == flag );
        break;
    case Reg_EBX:
        ret = ( ( ebx & flag ) == flag );
        break;
    case Reg_ECX:
        ret = ( ( ecx & flag ) == flag );
        break;
    case Reg_EDX:
        ret = ( ( edx & flag ) == flag );
        break;
    }

    return ret;
}

// Get number of CPU cores provided by the system
uint32_t XCpu::CoresCount( )
{
    return std::thread::hardware_concurrency( );
}
// Helper class wrapping some AVX intrinsics
class AvxTools
{
public:
    // Add two vectors: dst[i] += src[i]
    template <typename T> static inline void Add( const T* src, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Add<T, std::true_type, std::true_type>( src, dst, size );
            }
            else
            {
                Add<T, std::true_type, std::false_type>( src, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Add<T, std::false_type, std::true_type>( src, dst, size );
            }
            else
            {
                Add<T, std::false_type, std::false_type>( src, dst, size );
            }
        }
    }

    // Multiply two vectors: dst[i] *= src[i]
    template <typename T> static inline void Mul( const T* src, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Mul<T, std::true_type, std::true_type>( src, dst, size );
            }
            else
            {
                Mul<T, std::true_type, std::false_type>( src, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Mul<T, std::false_type, std::true_type>( src, dst, size );
            }
            else
            {
                Mul<T, std::false_type, std::false_type>( src, dst, size );
            }
        }
    }

    // Dot product: sum( vec1[i] * vec2[i] )
    template <typename T> static inline T Dot( const T* vec1, const T* vec2, size_t size )
    {
        T dotProduct;

        if ( IsAligned( vec1 ) )
        {
            if ( IsAligned( vec2 ) )
            {
                dotProduct = Dot<T, std::true_type, std::true_type>( vec1, vec2, size );
            }
            else
            {
                dotProduct = Dot<T, std::true_type, std::false_type>( vec1, vec2, size );
            }
        }
        else
        {
            if ( IsAligned( vec2 ) )
            {
                dotProduct = Dot<T, std::false_type, std::true_type>( vec1, vec2, size );
            }
            else
            {
                dotProduct = Dot<T, std::false_type, std::false_type>( vec1, vec2, size );
            }
        }

        return dotProduct;
    }

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    template <typename T> static inline void Max( const T* src, T alpha, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Max<T, std::true_type, std::true_type>( src, alpha, dst, size );
            }
            else
            {
                Max<T, std::true_type, std::false_type>( src, alpha, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Max<T, std::false_type, std::true_type>( src, alpha, dst, size );
            }
            else
            {
                Max<T, std::false_type, std::false_type>( src, alpha, dst, size );
            }
        }
    }

private:
    // Unroll size for single/double precision numbers - number of those in AVX register
    template <typename T> static inline size_t UnrollSize( );

    // Load 8 single / 4 double precision numbers into AVX register
    template <typename isAligned> static inline __m256  Load( const float*  src );
    template <typename isAligned> static inline __m256d Load( const double* src );

    // Store 8 single / 4 double precision numbers from AVX register into memory
    template <typename isAligned> static inline void Store( const __m256&  value, float*  dst );
    template <typename isAligned> static inline void Store( const __m256d& value, double* dst );

    // Check if the pointer is AVX-aligned (32 byte aligned)
    template <typename T> static inline bool IsAligned( const T* ptr )
    {
        return ( ( reinterpret_cast<uintptr_t>( ptr ) % 32 ) == 0 );
    }

    // Initialize 8 single / 4 double precision numbers of AVX register with the specified value 
    static inline __m256 Set1( float value )
    {
        return _mm256_set1_ps( value );
    }
    static inline __m256d Set1( double value )
    {
        return _mm256_set1_pd( value );
    }

    // Sum 8 single / 4 double precision numbers of AVX register
    static inline float  Sum( __m256 value );
    static inline double Sum( __m256d value );

    // Add 8 single / 4 double precision numbers
    static inline __m256 Add( const __m256& value1, const __m256& value2 )
    {
        return _mm256_add_ps( value1, value2 );
    }
    static inline __m256d Add( const __m256d& value1, const __m256d& value2 )
    {
        return _mm256_add_pd( value1, value2 );
    }

    // Multiple 8 single / 4 double precision numbers
    static inline __m256 Mul( const __m256& value1, const __m256& value2 )
    {
        return _mm256_mul_ps( value1, value2 );
    }
    static inline __m256d Mul( const __m256d& value1, const __m256d& value2 )
    {
        return _mm256_mul_pd( value1, value2 );
    }

    // Multiple and Add 8 single / 4 double precision numbers: value1 * value2 + value3
    static inline __m256 MAdd( const __m256& value1, const __m256& value2, const __m256& value3 )
    {
#ifdef USE_AVX2
        return _mm256_fmadd_ps( value1, value2, value3 );
#else
        return _mm256_add_ps( _mm256_mul_ps( value1, value2 ), value3 );
#endif
    }
    static inline __m256d MAdd( const __m256d& value1, const __m256d& value2, const __m256d& value3 )
    {
#ifdef USE_AVX2
        return _mm256_fmadd_pd( value1, value2, value3 );
#else
        return _mm256_add_pd( _mm256_mul_pd( value1, value2 ), value3 );
#endif
    }

    // Maximum of 8 single / 4 double precision numbers
    static inline __m256 Max( const __m256& value1, const __m256& value2 )
    {
        return _mm256_max_ps( value1, value2 );
    }
    static inline __m256d Max( const __m256d& value1, const __m256d& value2 )
    {
        return _mm256_max_pd( value1, value2 );
    }

    // Add two vectors
    template <typename T, typename srcAligned, typename dstAligned> static void Add( const T* src, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );
            auto d0 = Load<dstAligned>(  dst );
            auto d1 = Load<dstAligned>( &dst[blockSize ] );
            auto d2 = Load<dstAligned>( &dst[blockSize2] );
            auto d3 = Load<dstAligned>( &dst[blockSize3] );

            d0 = Add( s0, d0 );
            d1 = Add( s1, d1 );
            d2 = Add( s2, d2 );
            d3 = Add( s3, d3 );

            Store<dstAligned>( d0,  dst );
            Store<dstAligned>( d1, &dst[blockSize ] );
            Store<dstAligned>( d2, &dst[blockSize2] );
            Store<dstAligned>( d3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );
            auto d = Load<dstAligned>( dst );

            d = Add( s, d );

            Store<dstAligned>( d, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst += *src;

            src++;
            dst++;
        }
    }

    // Multiply two vectors
    template <typename T, typename srcAligned, typename dstAligned> static void Mul( const T* src, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );
            auto d0 = Load<dstAligned>(  dst );
            auto d1 = Load<dstAligned>( &dst[blockSize ] );
            auto d2 = Load<dstAligned>( &dst[blockSize2] );
            auto d3 = Load<dstAligned>( &dst[blockSize3] );

            d0 = Mul( s0, d0 );
            d1 = Mul( s1, d1 );
            d2 = Mul( s2, d2 );
            d3 = Mul( s3, d3 );

            Store<dstAligned>( d0,  dst );
            Store<dstAligned>( d1, &dst[blockSize ] );
            Store<dstAligned>( d2, &dst[blockSize2] );
            Store<dstAligned>( d3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );
            auto d = Load<dstAligned>( dst );

            d = Mul( s, d );

            Store<dstAligned>( d, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst *= *src;

            src++;
            dst++;
        }
    }

    // Dot product of two vectors
    template <typename T, typename vec1Aligned, typename vec2Aligned> static T Dot( const T* vec1, const T* vec2, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;
        
        auto   sum0 = Set1( T( 0 ) );
        auto   sum1 = Set1( T( 0 ) );
        auto   sum2 = Set1( T( 0 ) );
        auto   sum3 = Set1( T( 0 ) );

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto v10 = Load<vec1Aligned>(  vec1 );
            auto v11 = Load<vec1Aligned>( &vec1[blockSize ] );
            auto v12 = Load<vec1Aligned>( &vec1[blockSize2] );
            auto v13 = Load<vec1Aligned>( &vec1[blockSize3] );
            auto v20 = Load<vec2Aligned>(  vec2 );
            auto v21 = Load<vec2Aligned>( &vec2[blockSize ] );
            auto v22 = Load<vec2Aligned>( &vec2[blockSize2] );
            auto v23 = Load<vec2Aligned>( &vec2[blockSize3] );

            sum0 = MAdd( v10, v20, sum0 );
            sum1 = MAdd( v11, v21, sum1 );
            sum2 = MAdd( v12, v22, sum2 );
            sum3 = MAdd( v13, v23, sum3 );

            vec1 += blockSize4;
            vec2 += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto v1 = Load<vec1Aligned>( vec1 );
            auto v2 = Load<vec2Aligned>( vec2 );

            sum0 = MAdd( v1, v2, sum0 );

            vec1 += blockSize;
            vec2 += blockSize;
        }

        sum0  = Add( sum0, sum1 );
        sum0  = Add( sum0, sum2 );
        sum0  = Add( sum0, sum3 );

        T sum = Sum( sum0 );

        for ( size_t i = 0; i < remainIterations; i++ )
        {
            sum += *vec1 * *vec2;

            vec1++;
            vec2++;
        }

        return sum;
    }

    // Maximum value of vector's elements and the specified alpha value
    template <typename T, typename srcAligned, typename dstAligned> static void Max( const T* src, T alpha, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        auto   alphaVec = Set1( alpha );

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );

            s0 = Max( s0, alphaVec );
            s1 = Max( s1, alphaVec );
            s2 = Max( s2, alphaVec );
            s3 = Max( s3, alphaVec );

            Store<dstAligned>( s0,  dst );
            Store<dstAligned>( s1, &dst[blockSize ] );
            Store<dstAligned>( s2, &dst[blockSize2] );
            Store<dstAligned>( s3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );

            s = Max( s, alphaVec );

            Store<dstAligned>( s, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst = ( *src > alpha ) ? *src : alpha;

            src++;
            dst++;
        }
    }
};

// Unroll size for single/double precision numbers - number of those in AVX register
template <> inline size_t AvxTools::UnrollSize<float>( )
{
    return 8;
}
template <> inline size_t AvxTools::UnrollSize<double>( )
{
    return 4;
}

// Load 8 aligned single precision numbers
template <> inline __m256 AvxTools::Load<std::true_type>( const float* src )
{
    return _mm256_load_ps( src );
}

// Load 8 unaligned single precision numbers
template <> inline __m256 AvxTools::Load<std::false_type>( const float* src )
{
    return _mm256_loadu_ps( src );
}

// Load 4 aligned double precision numbers
template <> inline __m256d AvxTools::Load<std::true_type>( const double* src )
{
    return _mm256_load_pd( src );
}

// Load 4 unaligned double precision numbers
template <> inline __m256d AvxTools::Load<std::false_type>( const double* src )
{
    return _mm256_loadu_pd( src );
}

// Store 8 signle precision numbers into aligned memory
template <> inline void AvxTools::Store<std::true_type>( const __m256& value, float* dst )
{
    _mm256_store_ps( dst, value );
}

// Store 8 signle precision numbers into unaligned memory
template <> inline void AvxTools::Store<std::false_type>( const __m256& value, float* dst )
{
    _mm256_storeu_ps( dst, value );
}

// Store 4 double precision numbers into aligned memory
template <> inline void AvxTools::Store<std::true_type>( const __m256d& value, double* dst )
{
    _mm256_store_pd( dst, value );
}

// Store 4 double precision numbers into unaligned memory
template <> inline void AvxTools::Store<std::false_type>( const __m256d& value, double* dst )
{
    _mm256_storeu_pd( dst, value );
}

// Sum 8 single / 4 double precision numbers of AVX register
inline float AvxTools::Sum( __m256 value )
{
    float mem[8];

    Store<std::false_type>( value, mem );

    return mem[0] + mem[1] + mem[2] + mem[3] + mem[4] + mem[5] + mem[6] + mem[7];
}
inline double AvxTools::Sum( __m256d value )
{
    double mem[4];

    Store<std::false_type>( value, mem );

    return mem[0] + mem[1] + mem[2] + mem[3];
}

/* ============================================================================= */

// Check if the implementation of vector tools is available on the current system
bool XAvxVectorTools::IsAvailable( ) const
{
#if defined(ANNT_USE_AVX2)
    return XCpu::IsFeatureSupported( XCpu::Reg_EBX, XCpu::Flag_AVX2 );
#elif defined(ANNT_USE_AVX)
    return XCpu::IsFeatureSupported( XCpu::Reg_ECX, XCpu::Flag_AVX );
#else
    return false;
#endif
}

// Add two vectors: dst[i] += src[i]
void XAvxVectorTools::Add( const float* src, float* dst, size_t size ) const
{
    AvxTools::Add( src, dst, size );
}
void XAvxVectorTools::Add( const double* src, double* dst, size_t size ) const
{
    AvxTools::Add( src, dst, size );
};

// Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
void XAvxVectorTools::Mul( const float*  src, float*  dst, size_t size ) const
{
    AvxTools::Mul( src, dst, size );
}
void XAvxVectorTools::Mul( const double* src, double* dst, size_t size ) const
{
    AvxTools::Mul( src, dst, size );
}

// Dot product of two vectors: sum( vec1[i] * vec2[i] )
float XAvxVectorTools::Dot( const float* vec1, const float* vec2, size_t size ) const
{
    return AvxTools::Dot( vec1, vec2, size );
}
double XAvxVectorTools::Dot( const double* vec1, const double* vec2, size_t size ) const
{
    return AvxTools::Dot( vec1, vec2, size );
}

// Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
void XAvxVectorTools::Max( const float* src, float alpha, float* dst, size_t size ) const
{
    AvxTools::Max( src, alpha, dst, size );
}
void XAvxVectorTools::Max( const double* src, double alpha, double* dst, size_t size ) const
{
    AvxTools::Max( src, alpha, dst, size );
}


// Implementation of common vector routines using SSE/SSE2 instructions
class XSseVectorTools : public IVectorTools
{
public:
    // Check if the implementation of vector tools is available on the current system
    bool IsAvailable( ) const override;

    // Add two vectors: dst[i] += src[i]
    void Add( const float*  src, float*  dst, size_t size ) const override;
    void Add( const double* src, double* dst, size_t size ) const override;

    // Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
    void Mul( const float*  src, float*  dst, size_t size ) const override;
    void Mul( const double* src, double* dst, size_t size ) const override;

    // Dot product of two vectors: sum( vec1[i] * vec2[i] )
    float  Dot( const float*  vec1, const float*  vec2, size_t size ) const override;
    double Dot( const double* vec1, const double* vec2, size_t size ) const override;

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    void Max( const float*  src, float  alpha, float*  dst, size_t size ) const override;
    void Max( const double* src, double alpha, double* dst, size_t size ) const override;
};

// Helper class wrapping some SSE intrinsics
class SseTools
{
public:
    // Add two vectors: dst[i] += src[i]
    template <typename T> static inline void Add( const T* src, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Add<T, std::true_type, std::true_type>( src, dst, size );
            }
            else
            {
                Add<T, std::true_type, std::false_type>( src, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Add<T, std::false_type, std::true_type>( src, dst, size );
            }
            else
            {
                Add<T, std::false_type, std::false_type>( src, dst, size );
            }
        }
    }

    // Multiply two vectors: dst[i] *= src[i]
    template <typename T> static inline void Mul( const T* src, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Mul<T, std::true_type, std::true_type>( src, dst, size );
            }
            else
            {
                Mul<T, std::true_type, std::false_type>( src, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Mul<T, std::false_type, std::true_type>( src, dst, size );
            }
            else
            {
                Mul<T, std::false_type, std::false_type>( src, dst, size );
            }
        }
    }

    // Dot product: sum( vec1[i] * vec2[i] )
    template <typename T> static inline T Dot( const T* vec1, const T* vec2, size_t size )
    {
        T dotProduct;

        if ( IsAligned( vec1 ) )
        {
            if ( IsAligned( vec2 ) )
            {
                dotProduct = Dot<T, std::true_type, std::true_type>( vec1, vec2, size );
            }
            else
            {
                dotProduct = Dot<T, std::true_type, std::false_type>( vec1, vec2, size );
            }
        }
        else
        {
            if ( IsAligned( vec2 ) )
            {
                dotProduct = Dot<T, std::false_type, std::true_type>( vec1, vec2, size );
            }
            else
            {
                dotProduct = Dot<T, std::false_type, std::false_type>( vec1, vec2, size );
            }
        }

        return dotProduct;
    }

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    template <typename T> static inline void Max( const T* src, T alpha, T* dst, size_t size )
    {
        if ( IsAligned( src ) )
        {
            if ( IsAligned( dst ) )
            {
                Max<T, std::true_type, std::true_type>( src, alpha, dst, size );
            }
            else
            {
                Max<T, std::true_type, std::false_type>( src, alpha, dst, size );
            }
        }
        else
        {
            if ( IsAligned( dst ) )
            {
                Max<T, std::false_type, std::true_type>( src, alpha, dst, size );
            }
            else
            {
                Max<T, std::false_type, std::false_type>( src, alpha, dst, size );
            }
        }
    }

private:

    // Unroll size for single/double precision numbers - number of those in SSE register
    template <typename T> static inline size_t UnrollSize( );

    // Load 4 single / 2 double precision numbers into SSE register
    template <typename isAligned> static inline __m128  Load( const float*  src );
    template <typename isAligned> static inline __m128d Load( const double* src );

    // Store 4 single / 2 double precision numbers from SSE register into memory
    template <typename isAligned> static inline void Store( const __m128&  value, float*  dst );
    template <typename isAligned> static inline void Store( const __m128d& value, double* dst );

    // Check if the pointer is SSE-aligned (16 byte aligned)
    template <typename T> static inline bool IsAligned( const T* ptr )
    {
        return ( ( reinterpret_cast<uintptr_t>( ptr ) % 16 ) == 0 );
    }

    // Initialize 4 single / 2 double precision numbers of SSE register with the specified value
    static inline __m128 Set1( float value )
    {
        return _mm_set1_ps( value );
    }
    static inline __m128d Set1( double value )
    {
        return _mm_set1_pd( value );
    }

    // Sum 4 single / 2 double precision numbers of SSE register
    static inline float Sum( __m128 value );
    static inline double Sum( __m128d value );

    // Add 4 single / 2 double precision numbers
    static inline __m128 Add( const __m128& value1, const __m128& value2 )
    {
        return _mm_add_ps( value1, value2 );
    }
    static inline __m128d Add( const __m128d& value1, const __m128d& value2 )
    {
        return _mm_add_pd( value1, value2 );
    }

    // Multiple 4 single / 2 double precision numbers
    static inline __m128 Mul( const __m128& value1, const __m128& value2 )
    {
        return _mm_mul_ps( value1, value2 );
    }
    static inline __m128d Mul( const __m128d& value1, const __m128d& value2 )
    {
        return _mm_mul_pd( value1, value2 );
    }

    // Multiple and Add 4 single / 2 double precision numbers: value1 * value2 + value3
    static inline __m128 MAdd( const __m128& value1, const __m128& value2, const __m128& value3 )
    {
        return _mm_add_ps( _mm_mul_ps( value1, value2 ), value3 );
    }
    static inline __m128d MAdd( const __m128d& value1, const __m128d& value2, const __m128d& value3 )
    {
        return _mm_add_pd( _mm_mul_pd( value1, value2 ), value3 );
    }

    // Maximum of 8 single / 4 double precision numbers
    static inline __m128 Max( const __m128& value1, const __m128& value2 )
    {
        return _mm_max_ps( value1, value2 );
    }
    static inline __m128d Max( const __m128d& value1, const __m128d& value2 )
    {
        return _mm_max_pd( value1, value2 );
    }

    // Add two vectors
    template <typename T, typename srcAligned, typename dstAligned> static void Add( const T* src, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );
            auto d0 = Load<dstAligned>(  dst );
            auto d1 = Load<dstAligned>( &dst[blockSize ] );
            auto d2 = Load<dstAligned>( &dst[blockSize2] );
            auto d3 = Load<dstAligned>( &dst[blockSize3] );

            d0 = Add( s0, d0 );
            d1 = Add( s1, d1 );
            d2 = Add( s2, d2 );
            d3 = Add( s3, d3 );

            Store<dstAligned>( d0,  dst );
            Store<dstAligned>( d1, &dst[blockSize ] );
            Store<dstAligned>( d2, &dst[blockSize2] );
            Store<dstAligned>( d3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );
            auto d = Load<dstAligned>( dst );

            d = Add( s, d );

            Store<dstAligned>( d, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst += *src;

            src++;
            dst++;
        }
    }

    // Multiply two vectors
    template <typename T, typename srcAligned, typename dstAligned> static void Mul( const T* src, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );
            auto d0 = Load<dstAligned>(  dst );
            auto d1 = Load<dstAligned>( &dst[blockSize ] );
            auto d2 = Load<dstAligned>( &dst[blockSize2] );
            auto d3 = Load<dstAligned>( &dst[blockSize3] );

            d0 = Mul( s0, d0 );
            d1 = Mul( s1, d1 );
            d2 = Mul( s2, d2 );
            d3 = Mul( s3, d3 );

            Store<dstAligned>( d0, dst );
            Store<dstAligned>( d1, &dst[blockSize ] );
            Store<dstAligned>( d2, &dst[blockSize2] );
            Store<dstAligned>( d3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );
            auto d = Load<dstAligned>( dst );

            d = Mul( s, d );

            Store<dstAligned>( d, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst *= *src;

            src++;
            dst++;
        }
    }

    // Dot product of two vectors
    template <typename T, typename vec1Aligned, typename vec2Aligned> static T Dot( const T* vec1, const T* vec2, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;

        auto   sum0 = Set1( T( 0 ) );
        auto   sum1 = Set1( T( 0 ) );
        auto   sum2 = Set1( T( 0 ) );
        auto   sum3 = Set1( T( 0 ) );

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto v10 = Load<vec1Aligned>(  vec1 );
            auto v11 = Load<vec1Aligned>( &vec1[blockSize ] );
            auto v12 = Load<vec1Aligned>( &vec1[blockSize2] );
            auto v13 = Load<vec1Aligned>( &vec1[blockSize3] );
            auto v20 = Load<vec2Aligned>(  vec2 );
            auto v21 = Load<vec2Aligned>( &vec2[blockSize ] );
            auto v22 = Load<vec2Aligned>( &vec2[blockSize2] );
            auto v23 = Load<vec2Aligned>( &vec2[blockSize3] );

            sum0 = MAdd( v10, v20, sum0 );
            sum1 = MAdd( v11, v21, sum1 );
            sum2 = MAdd( v12, v22, sum2 );
            sum3 = MAdd( v13, v23, sum3 );

            vec1 += blockSize4;
            vec2 += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto v1 = Load<vec1Aligned>( vec1 );
            auto v2 = Load<vec2Aligned>( vec2 );

            sum0 = MAdd( v1, v2, sum0 );

            vec1 += blockSize;
            vec2 += blockSize;
        }

        sum0  = Add( sum0, sum1 );
        sum0  = Add( sum0, sum2 );
        sum0  = Add( sum0, sum3 );

        T sum = Sum( sum0 );

        for ( size_t i = 0; i < remainIterations; i++ )
        {
            sum += *vec1 * *vec2;

            vec1++;
            vec2++;
        }

        return sum;
    }

    // Maximum value of vector's elements and the specified alpha value
    template <typename T, typename srcAligned, typename dstAligned> static void Max( const T* src, T alpha, T* dst, size_t size )
    {
        size_t blockSize        = UnrollSize<T>( );
        size_t blockSize2       = blockSize * 2;
        size_t blockSize3       = blockSize * 3;
        size_t blockSize4       = blockSize * 4;
        size_t blockIterations4 = size / blockSize4;
        size_t blockIterations  = ( size - blockIterations4 * blockSize4 ) / blockSize;
        size_t remainIterations = size - blockIterations4 * blockSize4 - blockIterations * blockSize;
        auto   alphaVec         = Set1( alpha );

        // large blocks of 4
        for ( size_t i = 0; i < blockIterations4; i++ )
        {
            auto s0 = Load<srcAligned>(  src );
            auto s1 = Load<srcAligned>( &src[blockSize ] );
            auto s2 = Load<srcAligned>( &src[blockSize2] );
            auto s3 = Load<srcAligned>( &src[blockSize3] );

            s0 = Max( s0, alphaVec );
            s1 = Max( s1, alphaVec );
            s2 = Max( s2, alphaVec );
            s3 = Max( s3, alphaVec );

            Store<dstAligned>( s0,  dst );
            Store<dstAligned>( s1, &dst[blockSize ] );
            Store<dstAligned>( s2, &dst[blockSize2] );
            Store<dstAligned>( s3, &dst[blockSize3] );

            src += blockSize4;
            dst += blockSize4;
        }

        // small blocks of 1
        for ( size_t i = 0; i < blockIterations; i++ )
        {
            auto s = Load<srcAligned>( src );

            s = Max( s, alphaVec );

            Store<dstAligned>( s, dst );

            src += blockSize;
            dst += blockSize;
        }

        // remainder for compiler to decide
        for ( size_t i = 0; i < remainIterations; i++ )
        {
            *dst = ( *src > alpha ) ? *src : alpha;

            src++;
            dst++;
        }
    }
};

// Unroll size for single/double precision numbers - number of those in SSE register
template <> inline size_t SseTools::UnrollSize<float>( )
{
    return 4;
}
template <> inline size_t SseTools::UnrollSize<double>( )
{
    return 2;
}

// Load 4 aligned single precision numbers
template <> inline __m128 SseTools::Load<std::true_type>( const float* src )
{
    return _mm_load_ps( src );
}

// Load 8 unaligned single precision numbers
template <> inline __m128 SseTools::Load<std::false_type>( const float* src )
{
    return _mm_loadu_ps( src );
}

// Load 4 aligned double precision numbers
template <> inline __m128d SseTools::Load<std::true_type>( const double* src )
{
    return _mm_load_pd( src );
}

// Load 4 unaligned double precision numbers
template <> inline __m128d SseTools::Load<std::false_type>( const double* src )
{
    return _mm_loadu_pd( src );
}

// Store 4 signle precision numbers into aligned memory
template <> inline void SseTools::Store<std::true_type>( const __m128& value, float* dst )
{
    _mm_store_ps( dst, value );
}

// Store 4 signle precision numbers into unaligned memory
template <> inline void SseTools::Store<std::false_type>( const __m128& value, float* dst )
{
    _mm_storeu_ps( dst, value );
}

// Store 2 double precision numbers into aligned memory
template <> inline void SseTools::Store<std::true_type>( const __m128d& value, double* dst )
{
    _mm_store_pd( dst, value );
}

// Store 2 double precision numbers into unaligned memory
template <> inline void SseTools::Store<std::false_type>( const __m128d& value, double* dst )
{
    _mm_storeu_pd( dst, value );
}

// Sum 4 single / 2 double precision numbers of SSE register
inline float SseTools::Sum( __m128 value )
{
    float mem[4];

    Store<std::false_type>( value, mem );

    return mem[0] + mem[1] + mem[2] + mem[3];
}
inline double SseTools::Sum( __m128d value )
{
    double mem[2];

    Store<std::false_type>( value, mem );

    return mem[0] + mem[1];
}

/* ============================================================================= */

// Check if the implementation of vector tools is available on the current system
bool XSseVectorTools::IsAvailable( ) const
{
    // the double precision part requires SSE2
#if defined(ANNT_USE_SSE)
    return XCpu::IsFeatureSupported( XCpu::Reg_EDX, XCpu::Flag_SSE2 );
#else
    return false;
#endif
}

// Add two vectors: dst[i] += src[i]
void XSseVectorTools::Add( const float* src, float* dst, size_t size ) const
{
    SseTools::Add( src, dst, size );
}
void XSseVectorTools::Add( const double* src, double* dst, size_t size ) const
{
    SseTools::Add( src, dst, size );
};

// Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
void XSseVectorTools::Mul( const float*  src, float*  dst, size_t size ) const
{
    SseTools::Mul( src, dst, size );
}
void XSseVectorTools::Mul( const double* src, double* dst, size_t size ) const
{
    SseTools::Mul( src, dst, size );
}

// Dot product of two vectors: sum( vec1[i] * vec2[i] )
float XSseVectorTools::Dot( const float* vec1, const float* vec2, size_t size ) const
{
    return SseTools::Dot( vec1, vec2, size );
}
double XSseVectorTools::Dot( const double* vec1, const double* vec2, size_t size ) const
{
    return SseTools::Dot( vec1, vec2, size );
}

// Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
void XSseVectorTools::Max( const float* src, float alpha, float* dst, size_t size ) const
{
    SseTools::Max( src, alpha, dst, size );
}
void XSseVectorTools::Max( const double* src, double alpha, double* dst, size_t size ) const
{
    SseTools::Max( src, alpha, dst, size );
}

// Helper class providing the actual implementation
class VectorToolsImpl
{
public:
    // Add two vectors: dst[i] += src[i]
    template <typename T> static inline void Add( const T* src, T* dst, size_t size )
    {
        for ( size_t i = 0; i < size; i++ )
        {
            dst[i] += src[i];
        }
    }

    // Multiply two vectors: dst[i] *= src[i]
    template <typename T> static inline void Mul( const T* src, T* dst, size_t size )
    {
        for ( size_t i = 0; i < size; i++ )
        {
            dst[i] *= src[i];
        }
    }

    // Dot product: sum( vec1[i] * vec2[i] )
    template <typename T> static inline T Dot( const T* vec1, const T* vec2, size_t size )
    {
        T dotProduct = T( 0 );

        for ( size_t i = 0; i < size; i++ )
        {
            dotProduct += vec1[i] * vec2[i];
        }

        return dotProduct;
    }

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    template <typename T> static inline void Max( const T* src, T alpha, T* dst, size_t size )
    {
        for ( size_t i = 0; i < size; i++ )
        {
            dst[i] = ( src[i] > alpha ) ? src[i] : alpha;
        }
    }
};

// Default implementation of common vector routines without using
// any extended CPU instruction sets
class XVectorTools : public IVectorTools
{
public:
    // Check if the implementation of vector tools is available on the current system
    bool IsAvailable( ) const override;

    // Add two vectors: dst[i] += src[i]
    void Add( const float*  src, float*  dst, size_t size ) const override;
    void Add( const double* src, double* dst, size_t size ) const override;

    // Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
    void Mul( const float*  src, float*  dst, size_t size ) const override;
    void Mul( const double* src, double* dst, size_t size ) const override;

    // Dot product of two vectors: sum( vec1[i] * vec2[i] )
    float  Dot( const float*  vec1, const float*  vec2, size_t size ) const override;
    double Dot( const double* vec1, const double* vec2, size_t size ) const override;

    // Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
    void Max( const float*  src, float  alpha, float*  dst, size_t size ) const override;
    void Max( const double* src, double alpha, double* dst, size_t size ) const override;
};

/* ============================================================================= */

// Get vectorization tools available on current CPU architecture
IVectorTools* GetAvailableVectorTools( )
{
    IVectorTools* vectorTools = nullptr;

#ifdef ANNT_USE_AVX
    if ( vectorTools == nullptr )
    {
        vectorTools = new XAvxVectorTools( );

        if ( !vectorTools->IsAvailable( ) )
        {
            delete vectorTools;
            vectorTools = nullptr;
        }
    }
#endif

#ifdef ANNT_USE_SSE
    if ( vectorTools == nullptr )
    {
        vectorTools = new XSseVectorTools( );

        if ( !vectorTools->IsAvailable( ) )
        {
            delete vectorTools;
            vectorTools = nullptr;
        }
    }
#endif

    if ( vectorTools == nullptr )
    {
        vectorTools = new XVectorTools( );
    }

    return vectorTools;
}

// Initialize vectorizer with what is available
IVectorTools* XVectorize::mVectorTools = GetAvailableVectorTools( );
// Check if the implementation of vector tools is available on the current system
bool XVectorTools::IsAvailable( ) const
{
    return true;
}

// Add two vectors: dst[i] += src[i]
void XVectorTools::Add( const float* src, float* dst, size_t size ) const
{
    VectorToolsImpl::Add( src, dst, size );
}
void XVectorTools::Add( const double* src, double* dst, size_t size ) const
{
    VectorToolsImpl::Add( src, dst, size );
};

// Element wise multiplication of two vectors (Hadamard product): dst[i] *= src[i]
void XVectorTools::Mul( const float*  src, float*  dst, size_t size ) const
{
    VectorToolsImpl::Mul( src, dst, size );
}
void XVectorTools::Mul( const double* src, double* dst, size_t size ) const
{
    VectorToolsImpl::Mul( src, dst, size );
}

// Dot product of two vectors: sum( vec1[i] * vec2[i] )
float XVectorTools::Dot( const float* vec1, const float* vec2, size_t size ) const
{
    return VectorToolsImpl::Dot( vec1, vec2, size );
}
double XVectorTools::Dot( const double* vec1, const double* vec2, size_t size ) const
{
    return VectorToolsImpl::Dot( vec1, vec2, size );
}

// Calculates maximum of the vector's elements and the specified value: dst[i] = max( src[i], alpha )
void XVectorTools::Max( const float* src, float alpha, float* dst, size_t size ) const
{
    VectorToolsImpl::Max( src, alpha, dst, size );
}
void XVectorTools::Max( const double* src, double alpha, double* dst, size_t size ) const
{
    VectorToolsImpl::Max( src, alpha, dst, size );
}


// Provides functions to use for paralleling for-loops
class XParallel
{
private:
    XParallel( );

public:
    // Runs the specified lambda in a parallel for loop
    template <typename Func> static inline void For( size_t size, Func func )
    {
        #ifdef ANNT_USE_OMP
        #pragma omp parallel for
        #endif
        for ( int i = 0; i < static_cast<int>( size ); i++ )
        {
            func( static_cast<size_t>( i ) );
        }
    }

    // Conditionally runs the specified lambda in a parallel for loop
    template <typename Func> static inline void For( size_t size, bool parallel, Func func )
    {
        #ifdef ANNT_USE_OMP
        if ( parallel )
        {
            #pragma omp parallel for
            for ( int i = 0; i < static_cast< int >( size ); i++ )
            {
                func( static_cast< size_t >( i ) );
            }
        }
        else
        #else
        ANNT_UNREFERENCED_PARAMETER( parallel );
        #endif
        {
            for ( size_t i = 0; i < size; i++ )
            {
                func( i );
            }
        }
    }
};


    namespace Neuro {

        class XNeuralNetwork;

        enum class LayerID
        {
            Unknown            = 0,

            FullyConnected     = 1,
            Convolution        = 2,
            RecurrentBasic     = 3,
            RecurrentLSTM      = 4,
            RecurrentGRU       = 5,

            Sigmoid            = 1000,
            Tanh               = 1001,
            Relu               = 1002,
            LeakyRelu          = 1003,
            Elu                = 1004,
            Softmax            = 1005,
            LogSoftmax         = 1006,

            MaxPooling         = 2001,
            AveragePooling     = 2002,
            DropOut            = 2003,
            BatchNormalization = 2004,
        };

        class XNeuralNetwork;
        class XNetworkInference;

        namespace Training
        {
            class XNetworkTraining;
        }

        // The class encapsulates some context passed to layers by inference/training classes
        class XNetworkContext
        {
            friend class XNetworkInference;
            friend class Training::XNetworkTraining;
            
        private:

            bool    mTrainingMode;
            size_t  mTrainingSequenceLength;   // length of sequences used to train recurrent networks
            size_t  mCurrentLayer;


            std::vector<std::vector<std::vector<void*>>> mLayersMemoryBuffers;
            std::vector<uvector_t>                       mLayersMemorySize;

        public:

            XNetworkContext( bool trainingMode ) :
                XNetworkContext( trainingMode, 1 )
            { }

            XNetworkContext( bool trainingMode, size_t sequenceLength );
            ~XNetworkContext( );

            // Checks if network is being trained
            bool IsTraining( ) const { return mTrainingMode; }

            // Get/set length of training sequences used for recurrent networks
            size_t TrainingSequenceLength( ) const
            {
                return mTrainingSequenceLength;
            }
            void SetTrainingSequenceLength( size_t sequenceLength )
            {
                mTrainingSequenceLength = sequenceLength;
            }

            // Provides specified working buffer for the sample index
            void* GetWorkingBuffer( size_t buffer, size_t sample ) const
            {
                return mLayersMemoryBuffers[mCurrentLayer][buffer][sample];
            }

        protected:

            // Allocate working buffer for layers of the network
            void AllocateWorkingBuffers( const std::shared_ptr<XNeuralNetwork>& net, size_t batchSize );

            // Clear layers' working buffers (memset zero)
            void ResetWorkingBuffers( );
            void ResetWorkingBuffers( uvector_t layersIndexes );

            // Set current layer index, so that correct working buffer could be provided
            void SetCurrentLayerIndex( size_t currentLayer )
            {
                mCurrentLayer = currentLayer;
            }

        private:

            // Free layers' working buffers
            void FreeWorkingBuffers( );
        };

        
        class ILayer
        {
            friend class XNeuralNetwork;

        protected:
            size_t mInputsCount;
            size_t mOutputsCount;

            // To be called by XNeuralNetwork to set size of layers with zero inputs/outputs,
            // which is the case for activation layers.
            virtual void Initialize( size_t inputsCount, size_t outputsCount )
            {
                mInputsCount  = inputsCount;
                mOutputsCount = outputsCount;
            }

        public:
            ILayer( size_t inputsCount, size_t outputsCount )
            {
                Initialize( inputsCount, outputsCount );
            }

            virtual ~ILayer( ) { }

            // Size of input vector expected by the layer
            size_t InputsCount( ) const
            {
                return mInputsCount;
            }

            // Size of output vector produced by the layer
            size_t OutputsCount( ) const
            {
                return mOutputsCount;
            }

            // Some of the layers may need extra memory required for processing inputs or for
            // keeping state between forward and backward pass. The method below tells
            // how many buffers are required and their size.
            //
            // For example, if a layer needs two temporary buffers of type *float_t*
            // (one vector with *inputsCount* elements and another with *outputsCount* elements),
            // it may return something like this: uvector_t( { inputsCount * sizeof( float ), outputsCount * sizeof( float ) } ).
            //
            // Each individual memory buffer is 32 byte aligned, so AVX friendly.
            //
            // Number of allocated buffers equals to number of samples in a batch.
            // 
            virtual uvector_t WorkingMemSize( bool /* trainingMode */ ) const { return uvector_t( 0 ); }

            // Reports if the layer is trainable or not (has weights/biases)
            virtual bool Trainable( ) const = 0;

            // Calculates outputs for the given inputs - forward pass
            virtual void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                        std::vector<fvector_t*>& outputs,
                                        const XNetworkContext& ctx ) = 0;

            // Propagates error to the previous layer and calculates weights/biases
            // gradients (in the case the layer is trainable)
            virtual void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                        const std::vector<fvector_t*>& outputs,
                                        const std::vector<fvector_t*>& deltas,
                                        std::vector<fvector_t*>& prevDeltas,
                                        fvector_t& gradWeights,
                                        const XNetworkContext& ctx ) = 0;

            // Saves layer's learnt parameters/weights
            virtual bool SaveLearnedParams( FILE* /* file */ ) const { return true; }
            // Loads layer's learnt parameters
            virtual bool LoadLearnedParams( FILE* /* file */ ) { return true; }

        protected:

            // Default implementation of saving layer's learned parameter, which are represented as fvector_t
            bool SaveLearnedParamsHelper( FILE* file, LayerID id, const std::vector<const fvector_t*>& params ) const
            {
                bool     ret     = false;
                uint32_t layerID = static_cast<uint32_t>( id );

                if ( fwrite( &layerID, sizeof( layerID ), 1, file ) == 1 )
                {
                    size_t i;

                    for ( i = 0; i < params.size( ); i++ )
                    {
                        uint32_t paramsCount = static_cast<uint32_t>( params[i]->size( ) );

                        if ( fwrite( &paramsCount, sizeof( paramsCount ), 1, file ) != 1 )
                        {
                            break;
                        }
                    }

                    if ( i == params.size( ) )
                    {
                        for ( i = 0; i < params.size( ); i++ )
                        {
                            if ( fwrite( params[i]->data( ), sizeof( float_t ), params[i]->size( ), file ) != params[i]->size( ) )
                            {
                                break;
                            }
                        }

                        ret = ( i == params.size( ) );
                    }
                }

                return ret;
            }

            // Default implementation of loading layer's learned parameter, which are represented as fvector_t
            bool LoadLearnedParamsHelper( FILE* file, LayerID id, std::vector<fvector_t*>& params )
            {
                bool     ret = false;
                uint32_t layerID;

                if ( ( fread( &layerID, sizeof( layerID ), 1, file ) == 1 ) &&
                    ( layerID == static_cast<uint32_t>( id ) ) )
                {
                    size_t i;

                    for ( i = 0; i < params.size( ); i++ )
                    {
                        uint32_t paramsCount;

                        if ( ( fread( &paramsCount, sizeof( paramsCount ), 1, file ) != 1 ) ||
                            ( paramsCount != static_cast<uint32_t>( params[i]->size( ) ) ) )
                        {
                            break;
                        }
                    }

                    if ( i == params.size( ) )
                    {
                        for ( i = 0; i < params.size( ); i++ )
                        {
                            if ( fread( params[i]->data( ), sizeof( float_t ), params[i]->size( ), file ) != params[i]->size( ) )
                            {
                                break;
                            }
                        }

                        ret = ( i == params.size( ) );
                    }
                }

                return ret;
            }
        };


        class IActivationLayer : public ILayer
        {
        public:
            IActivationLayer( ) : ILayer( 0, 0 )
            {
            }

            // None of the activation functions have weights/biases to train
            bool Trainable( ) const override
            {
                return false;
            }

            // Calls ForwardActivate() for individual input/output vectors passed by reference
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    ForwardActivate( inputs[i]->data( ), outputs[i]->data( ), inputs[i]->size( ) );
                } );
            }

            // Calls BackwardActivate() for individual input/output/delta vectors passed by reference
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& /* gradWeights */,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    BackwardActivate( inputs[i]->data( ), outputs[i]->data( ), deltas[i]->data( ), prevDeltas[i]->data( ), inputs[i]->size( ) );
                } );
            }

            // Applies activation function to the input vector
            virtual void ForwardActivate( const float_t* input, float_t* output, size_t len ) = 0;

            // Propagates error back to previous layer by multiplying delta with activation function's derivative
            virtual void BackwardActivate( const float_t* input, const float_t* output,
                                        const float_t* delta, float_t* prevDelta, size_t len ) = 0;
        };

        // Implementation of ELU activation function,
        //   f(x) = | x                      , x >= 0
        //          | alpha * ( exp(x) - 1 ) , x <  0
        class XEluActivation : public IActivationLayer
        {
        private:
            float   mAlpha;

        public:

            XEluActivation( float alpha = 1.0f ) : mAlpha( alpha )
            {
            }

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = ( input[i] >= float_t( 0 ) ) ? input[i] : mAlpha * ( std::exp( input[i] ) - float_t( 1 ) );
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    // derivative(ELU) = 1        , y >= 0
                    //                 = alpha + y, otherwise
                    prevDelta[i] = ( output[i] >= float_t( 0 ) ) ? delta[i] : ( mAlpha + output[i] ) * delta[i];
                }
            }
        };

        // Implementation of Leaky ReLU activation function,
        //   f(x) = | x         , x >  0
        //          | x * alpha , x <= 0
        class XLeakyReLuActivation : public IActivationLayer
        {
        private:
            float_t   mAlpha;

        public:

            XLeakyReLuActivation( float alpha = 0.01f ) : mAlpha( alpha )
            {
            }

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = ( input[i] > float_t( 0 ) ) ? input[i] : mAlpha * input[i];
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    // derivative(LeakyReLU) = 1    , y > 0
                    //                       = alpha, otherwise
                    prevDelta[i] = ( output[i] > float_t( 0 ) ) ? delta[i] : mAlpha * delta[i];
                }
            }
        };

        // Implementation of Log-SoftMax activation function - to be used with XNegativeLogLikelihoodCost,
        // which expects log-probabilities as input.
        // http://www.mlpack.org/docs/mlpack-git/doxygen/classmlpack_1_1ann_1_1LogSoftMax.html
        //
        // Using it since SoftMax+CrossEntropy can lead to NaN in gradient for unbounded activations:
        // https://github.com/Theano/Theano/issues/3162
        //
        class XLogSoftMaxActivation : public IActivationLayer
        {
        public:

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                float_t sum = 0;
                float_t max = input[0];

                for ( size_t i = 1; i < len; i++ )
                {
                    if ( input[i] > max ) max = input[i];
                }

                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = std::exp( input[i] - max );
                    sum      += output[i];
                }

                sum = std::log( sum );

                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = input[i] - max - sum;
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    prevDelta[i] = std::exp( output[i] ) +  delta[i];
                }
            }
        };

        // Implementation of ReLU activation function, f(x) = max(0, x)
        class XReLuActivation : public IActivationLayer
        {
        public:

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                XVectorize::Max( input, float_t( 0 ), output, len );
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    // derivative(ReLU) = 1, y > 0
                    //                  = 0, otherwise
                    prevDelta[i] = ( output[i] > float_t( 0 ) ) ? delta[i] : float_t( 0 );
                }
            }
        };

        // Implementation of Sigmoid activation function, f(x) = 1 / ( 1 + exp(-x) )
        class XSigmoidActivation : public IActivationLayer
        {
        public:

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = float_t( 1 ) / ( float_t( 1 ) + std::exp( -input[i] ) );
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    // derivative(Sigmoid) = y * ( 1 - y )
                    prevDelta[i] = delta[i] * output[i] * ( float_t( 1 ) - output[i] );
                }
            }
        };

        // Implementation of SoftMax activation function
        class XSoftMaxActivation : public IActivationLayer
        {
        public:

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                float_t sum = 0;
                float_t max = input[0];

                for ( size_t i = 1; i < len; i++ )
                {
                    if ( input[i] > max ) max = input[i];
                }

                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = std::exp( input[i] - max );
                    sum      += output[i];
                }

                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] /= sum;
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output,
                                const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    float_t sum = 0;

                    for ( size_t j = 0; j < len; j++ )
                    {
                        //sum += delta[j] * der[j];
                        sum += delta[j] * ( ( j == i ) ? output[i] * ( float_t( 1 ) - output[j] ) : -output[i] * output[j] );
                    }

                    prevDelta[i] = sum;
                }
            }
        };

        // Implementation of Hyperbolic Tangent activation function, f(x) = tanh(x)
        class XTanhActivation : public IActivationLayer
        {
        public:

            void ForwardActivate( const float_t* input, float_t* output, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    output[i] = std::tanh( input[i] );
                }
            }

            void BackwardActivate( const float_t* /* input */, const float_t* output, const float_t* delta, float_t* prevDelta, size_t len ) override
            {
                for ( size_t i = 0; i < len; i++ )
                {
                    // derivative(Tanh) = 1 - y^2
                    prevDelta[i] = delta[i] * ( float_t( 1 ) - output[i] * output[i] );
                }
            }
        };

        // Implementation of neural network class, which is a simple sequence of layers for now.
        // No computation/training is done here - it is moved to dedicated classes which manage
        // all the resources involved into that.
        //
        class XNeuralNetwork
        {
        private:
            std::vector<std::shared_ptr<ILayer>> mLayers;

        public:
            typedef std::vector<std::shared_ptr<ILayer>>::iterator       iterator;
            typedef std::vector<std::shared_ptr<ILayer>>::const_iterator const_iterator;

            XNeuralNetwork( )
            {
            }

            // Reports network's inputs count (inputs count of the first layer if any)
            size_t InputsCount( ) const
            {
                return ( mLayers.empty( ) ) ? 0 : mLayers.front( )->InputsCount( );
            }

            // Reports network's outputs count (outputs count of the last layer)
            size_t OutputsCount( ) const
            {
                return ( mLayers.empty( ) ) ? 0 : mLayers.back( )->OutputsCount( );
            }

            // Reports total number of layers
            size_t LayersCount( ) const
            {
                return mLayers.size( );
            }

            // Iterators to access layers
            iterator begin( )
            {
                return mLayers.begin( );
            }
            iterator end( )
            {
                return mLayers.end( );
            }
            const_iterator begin( ) const
            {
                return mLayers.begin( );
            }
            const_iterator end( ) const
            {
                return mLayers.end( );
            }

            // Provides layer at the specified index
            std::shared_ptr<ILayer> LayerAt( size_t index ) const
            {
                return ( index < mLayers.size( ) ) ? mLayers[index] : std::shared_ptr<ILayer>( nullptr );
            }

            // Adds the specified layer to the end of layers' collection
            void AddLayer( const std::shared_ptr<ILayer>& layer );

            // Saves network's learned parameters only.
            // Network structure is not saved and so same network must be constructed before loading parameters
            bool SaveLearnedParams( const std::string& fileName ) const;

            // Loads network's learned parameters.
            // A network of the same structure as saved must be created first, since this method loads only parameters/weights/biases.
            bool LoadLearnedParams( const std::string& fileName );
        };

        // Adds the specified layer to the end of layers' collection
        void XNeuralNetwork::AddLayer( const shared_ptr<ILayer>& layer )
        {
            if ( layer->InputsCount( ) == 0 )
            {
                if ( mLayers.empty( ) )
                {
                    // TODO: error
                    assert( false );
                }
                else
                {
                    size_t lastOutSize = mLayers.back( )->OutputsCount( );

                    layer->Initialize( lastOutSize, lastOutSize );
                }
            }

            if ( ( mLayers.empty( ) ) || ( layer->InputsCount( ) == mLayers.back( )->OutputsCount( ) ) )
            {
                mLayers.push_back( layer );
            }
            else
            {
                // TODO: error
                assert( false );
            }
        }

        // Saves network's learned parameters only.
        // Network structure is not saved and so same network must be constructed before loading parameters
        bool XNeuralNetwork::SaveLearnedParams( const string& fileName ) const
        {
            FILE* file = fopen( fileName.c_str( ), "wb" );
            bool  ret = false;

            if ( file != nullptr )
            {
                // float_t can be defined to other than "float", so need to save its size and stop loading
                // saves produced by incompatible builds
                uint8_t floatTypeSize = static_cast<uint8_t>( sizeof( float_t ) );

                if ( ( fwrite( "ANNT", sizeof( char ), 4, file ) == 4 ) &&
                    ( fwrite( &floatTypeSize, sizeof( floatTypeSize ), 1, file ) == 1 ) )
                {
                    ret = true;

                    for ( const_iterator layersIt = mLayers.begin( ); ( ret ) && ( layersIt != mLayers.end( ) ); layersIt++ )
                    {
                        ret = ( *layersIt )->SaveLearnedParams( file );
                    }
                }

                fclose( file );
            }

            return ret;
        }

        // Loads network's learned parameters.
        // A network of the same structure as saved must be created first, since this method loads only parameters/weights/biases.
        bool XNeuralNetwork::LoadLearnedParams( const string& fileName )
        {
            FILE* file = fopen( fileName.c_str( ), "rb" );
            bool  ret = false;

            if ( file != nullptr )
            {
                char    anntMagic[4];
                uint8_t floatTypeSize;

                if ( ( fread( anntMagic, sizeof( char ), 4, file ) == 4 ) &&
                    ( fread( &floatTypeSize, sizeof( floatTypeSize ), 1, file ) == 1 ) &&
                    ( anntMagic[0] == 'A' ) && ( anntMagic[1] == 'N' ) && ( anntMagic[2] == 'N' ) && ( anntMagic[3] == 'T' ) &&
                    ( floatTypeSize == static_cast<uint8_t>( sizeof( float_t ) ) ) )
                {
                    ret = true;

                    for ( const_iterator layersIt = mLayers.begin( ); ( ret ) && ( layersIt != mLayers.end( ) ); layersIt++ )
                    {
                        ret = ( *layersIt )->LoadLearnedParams( file );
                    }
                }

                fclose( file );
            }

            return ret;
        }

        // Implementation of artificial neural network inference - wraps
        // everything necessary to compute network's outputs for a given inputs.
        //
        class XNetworkInference
        {
        protected:
            std::shared_ptr<XNeuralNetwork>      mNetwork;

        protected:
            std::vector<std::vector<fvector_t>>  mComputeOutputsStorage;
            std::vector<std::vector<fvector_t*>> mComputeOutputs;
            std::vector<fvector_t*>              mComputeInputs;

            XNetworkContext                      mInferenceContext;

        public:
            // The passed network must be fully constructed at this point - no adding new layers
            XNetworkInference( const std::shared_ptr<XNeuralNetwork>& network );

            // Reset working buffers for all layers
            virtual void ResetState( )
            {
                mInferenceContext.ResetWorkingBuffers( );
            }

            // Reset working buffers for the specified layers
            virtual void ResetLayersState( uvector_t layersIndexes )
            {
                mInferenceContext.ResetWorkingBuffers( layersIndexes );
            }

            // Computes output vector for the given input vector
            void Compute( const fvector_t& input, fvector_t& output );

            // Runs classification for the given input - returns index of the maximum
            // element in the corresponding output vector
            size_t Classify( const fvector_t& input );

            // Tests classification for the provided inputs and target labels -
            // provides number of correctly classified samples
            size_t TestClassification( const std::vector<fvector_t>& inputs,
                                    const uvector_t& targetLabels );

        protected:

            // Helper method to compute output vectors for the given input vectors using
            // the provided storage for the intermediate outputs of all layers
            void DoCompute( const std::vector<fvector_t*>& inputs,
                            std::vector<std::vector<fvector_t*>>& outputs,
                            XNetworkContext& ctx );
        };

        XNetworkInference::XNetworkInference( const shared_ptr<XNeuralNetwork>& network ) :
            mNetwork( network ),
            mInferenceContext( false )
        {
            mComputeInputs.resize( 1 );

            // prepare output vectors for all layers and for all samples (only one sample here for now)
            for ( auto layer : *mNetwork )
            {
                mComputeOutputsStorage.push_back( vector<fvector_t>( { fvector_t( layer->OutputsCount( ) ) } ) );
                mComputeOutputs.push_back( vector<fvector_t*>( { &( mComputeOutputsStorage.back( )[0] ) } ) );
            }

            mInferenceContext.AllocateWorkingBuffers( network, 1 );
        }

        // Computes output vector for the given input vector
        void XNetworkInference::Compute( const fvector_t& input, fvector_t& output )
        {
            if ( mNetwork->LayersCount( ) != 0 )
            {
                mComputeInputs[0] = const_cast<fvector_t*>( &input );

                DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                // copy output produced by the last layer
                output = mComputeOutputsStorage.back( )[0];
            }
        }

        // Runs classification for the given input - returns index of the maximum element in the corresponding output
        size_t XNetworkInference::Classify( const fvector_t& input )
        {
            size_t classIndex = 0;

            if ( mNetwork->LayersCount( ) != 0 )
            {
                mComputeInputs[0] = const_cast<fvector_t*>( &input );

                DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                classIndex = XDataEncodingTools::MaxIndex( mComputeOutputsStorage.back( )[0] );
            }

            return classIndex;
        }

        // Tests classification for the provided inputs and target labels - provides number of correctly classified samples
        size_t XNetworkInference::TestClassification( const vector<fvector_t>& inputs, const uvector_t& targetLabels )
        {
            size_t correctLabelsCounter = 0;

            if ( mNetwork->LayersCount( ) != 0 )
            {
                for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                {
                    mComputeInputs[0] = const_cast<fvector_t*>( &( inputs[i] ) );

                    DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                    if ( XDataEncodingTools::MaxIndex( mComputeOutputsStorage.back( )[0] ) == targetLabels[i] )
                    {
                        correctLabelsCounter++;
                    }
                }
            }

            return correctLabelsCounter;
        }

        // Helper method to compute output vectors for the given input vectors
        void XNetworkInference::DoCompute( const vector<fvector_t*>& inputs,
                                        vector<vector<fvector_t*>>& outputs,
                                        XNetworkContext& ctx )
        {
            ctx.SetCurrentLayerIndex( 0 );
            mNetwork->LayerAt( 0 )->ForwardCompute( inputs, outputs[0], ctx );

            for ( size_t i = 1, layersCount = mNetwork->LayersCount( ); i < layersCount; i++ )
            {
                ctx.SetCurrentLayerIndex( i );
                mNetwork->LayerAt( i )->ForwardCompute( outputs[i - 1], outputs[i], ctx );
            }
        }

        class ITrainableLayer : public ILayer
        {
        public:
            ITrainableLayer( size_t inputsCount, size_t outputsCount ) :
                ILayer( inputsCount, outputsCount )
            {
            }

            // Reports the layer is trainable
            bool Trainable( ) const override
            {
                return true;
            }

            // Reports number of weight coefficients the layer has
            virtual size_t WeightsCount( ) const = 0;

            // Get/set layer's weights
            virtual fvector_t Weights( ) const = 0;
            virtual void SetWeights( const fvector_t& weights ) = 0;

            // Randomizes layer's weights/biases
            virtual void Randomize( ) = 0;

            // Applies updates to the layer's weights and biases
            virtual void UpdateWeights( const fvector_t& updates ) = 0;
        };


        // Implementation of fully connected layer - each neuron is connected to each input
        class XFullyConnectedLayer : public ITrainableLayer
        {
        private:
            // Weights and biases are all kept together
            fvector_t mAllWeights;

            // And here are their pointers
            float_t*  mWeights;
            float_t*  mBiases;

        public:
            XFullyConnectedLayer( size_t inputsCount, size_t outputsCount );

            // Reports number of weight coefficients the layer has
            size_t WeightsCount( ) const override
            {
                return mAllWeights.size( );
            }

            // Get/set layer's weights
            fvector_t Weights( ) const override
            {
                return mAllWeights;
            }
            void SetWeights( const fvector_t& weights ) override
            {
                mAllWeights = weights;
            }

            // Randomizes layer's weights, clears biases
            void Randomize( ) override;

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override;

            // Propagates error to the previous layer and calculates weights/biases gradients
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& gradWeights,
                                const XNetworkContext& ctx ) override;

            // Applies updates to the layer's weights and biases
            void UpdateWeights( const fvector_t& updates ) override;

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override;
            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override;
        };

        XFullyConnectedLayer::XFullyConnectedLayer( size_t inputsCount, size_t outputsCount ) :
            ITrainableLayer( inputsCount, outputsCount ),
            mAllWeights( inputsCount * outputsCount + outputsCount )
        {
            // set up weights/biases pointers
            mWeights = mAllWeights.data( );
            mBiases  = mWeights + mInputsCount * mOutputsCount;

            Randomize( );
        }

        // Randomizes layer's weights, clears biases
        void XFullyConnectedLayer::Randomize( )
        {
            float_t halfRange = sqrt( float_t( 3 ) / mInputsCount );

            for ( size_t i = 0, n = mInputsCount * mOutputsCount; i < n; i++ )
            {
                mWeights[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRange ) - halfRange;
            }
            for ( size_t i = 0; i < mOutputsCount; i++ )
            {
                mBiases[i] = 0;
            }
        }

        // Calculates outputs for the given inputs
        void XFullyConnectedLayer::ForwardCompute( const vector<fvector_t*>& inputs,
                                                vector<fvector_t*>& outputs,
                                                const XNetworkContext& ctx )
        {
            XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
            {
                const float_t* weights = mWeights;
                const float_t* input   = inputs[i]->data( );
                fvector_t&     output  = *( outputs[i] );

                for ( size_t otputIndex = 0; otputIndex < mOutputsCount; otputIndex++ )
                {
                    output[otputIndex] = XVectorize::Dot( input, weights, mInputsCount ) + mBiases[otputIndex];

                    weights += mInputsCount;
                }
            } );
        }

        // Propagates error to the previous layer and calculates weights/biases gradients
        void XFullyConnectedLayer::BackwardCompute( const vector<fvector_t*>& inputs,
                                                    const vector<fvector_t*>& /* outputs */,
                                                    const vector<fvector_t*>& deltas,
                                                    vector<fvector_t*>& prevDeltas,
                                                    fvector_t& gradWeights,
                                                    const XNetworkContext& ctx )
        {
            // set up weights/biases gradients pointers
            float_t*  gradWeightsData = gradWeights.data( );
            float_t*  gradBiasesData  = gradWeightsData + mInputsCount * mOutputsCount;

            // 1 - first propagate deltas to the previous layer
            XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
            {
                fvector_t&       prevDelta = *( prevDeltas[i] );
                const fvector_t& delta     = *( deltas[i] );

                for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                {
                    size_t  weightIndex = inputIndex;
                    float_t sum         = 0;

                    for ( size_t otputIndex = 0; otputIndex < mOutputsCount; otputIndex++, weightIndex += mInputsCount )
                    {
                        sum += delta[otputIndex] * mWeights[weightIndex];
                    }

                    prevDelta[inputIndex] = sum;
                }
            } );

            // 2 - accumulate weights' difference
            XParallel::For( mOutputsCount, ctx.IsTraining( ), [&]( size_t outputIndex )
            {
                for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                {
                    const fvector_t& input      = *( inputs[i] );
                    float_t          deltaValue = ( *( deltas[i] ) )[outputIndex];

                    for ( size_t inputIndex = 0, weightIndex = outputIndex * mInputsCount; inputIndex < mInputsCount; inputIndex++, weightIndex++ )
                    {
                        gradWeightsData[weightIndex] += deltaValue * input[inputIndex];
                    }
                }
            } );

            // 3 - accumulate baises' difference
            for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
            {
                const fvector_t& delta = *( deltas[i] );

                for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                {
                    gradBiasesData[outputIndex] += delta[outputIndex];
                }
            }
        }

        // Applies updates to the layer's weights and biases
        void XFullyConnectedLayer::UpdateWeights( const fvector_t& updates )
        {
            for ( size_t i = 0, n = mAllWeights.size( ); i < n; i++ )
            {
                mAllWeights[i] += updates[i];
            }
        }

        // Saves layer's learnt parameters/weights
        bool XFullyConnectedLayer::SaveLearnedParams( FILE* file ) const
        {
            vector<const fvector_t*> params( { &mAllWeights } );

            return SaveLearnedParamsHelper( file, LayerID::FullyConnected, params );
        }

        // Loads layer's learnt parameters
        bool XFullyConnectedLayer::LoadLearnedParams( FILE* file )
        {
            vector<fvector_t*> params( { &mAllWeights } );

            return LoadLearnedParamsHelper( file, LayerID::FullyConnected, params );
        }

        // Implementation of convolution layer - output is the result of convolving input with layer's weight (convolution kernel)
        class XConvolutionLayer : public ITrainableLayer
        {
        private:
            size_t      mInputWidth;
            size_t      mInputHeight;
            size_t      mInputDepth;
            size_t      mOutputWidth;
            size_t      mOutputHeight;
            size_t      mKernelWidth;
            size_t      mKernelHeight;
            size_t      mKernelsCount;
            size_t      mHorizontalStep;
            size_t      mVerticalStep;
            BorderMode  mBorderMode;

            std::vector<bool>   mConnectionTable;
            std::vector<size_t> mKernelOffsets;

            size_t      mPaddedWidth;
            size_t      mPaddedHeight;

            // Weights and biases are all kept together
            fvector_t   mAllWeights;
            size_t      mWeightCount; // excluding biases

            // And here are their pointers
            float_t*    mKernelsWeights;
            float_t*    mKernelsBiases;

        public:

            XConvolutionLayer( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount ) :
                XConvolutionLayer( inputWidth, inputHeight, inputDepth,
                                kernelWidth, kernelHeight, kernelsCount,
                                std::vector<bool>( ) , BorderMode::Valid, 1, 1 )
            {
            }

            XConvolutionLayer( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount,
                            BorderMode borderMode ) :
                XConvolutionLayer( inputWidth, inputHeight, inputDepth,
                                kernelWidth, kernelHeight, kernelsCount,
                                std::vector<bool>( ) , borderMode, 1, 1 )
            {
            }

            XConvolutionLayer( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount,
                            BorderMode borderMode, size_t horizontalStep, size_t verticalStep ) :
                XConvolutionLayer( inputWidth, inputHeight, inputDepth,
                                kernelWidth, kernelHeight, kernelsCount,
                                std::vector<bool>( ) , borderMode, horizontalStep, verticalStep )
            {
            }

            XConvolutionLayer( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount,
                            const std::vector<bool>& connectionTable ) :
                XConvolutionLayer( inputWidth, inputHeight, inputDepth,
                                kernelWidth, kernelHeight, kernelsCount,
                                connectionTable, BorderMode::Valid, 1, 1 )
            {
            }

            XConvolutionLayer( size_t inputWidth,  size_t inputHeight,  size_t inputDepth,
                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount,
                            const std::vector<bool>& connectionTable,
                            BorderMode borderMode, size_t horizontalStep, size_t verticalStep );

            // Reports number of weight coefficients the layer has
            size_t WeightsCount( ) const override
            {
                return mAllWeights.size( );
            }

            // Get/set layer's weights
            fvector_t Weights( ) const override
            {
                return mAllWeights;
            }
            void SetWeights( const fvector_t& weights ) override
            {
                mAllWeights = weights;
            }

            // Tells that we may need some extra memory for padding/unpadding
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize = uvector_t( 2, 0 );

                if ( mBorderMode == BorderMode::Same )
                {
                    workingMemSize[1] = workingMemSize[0] = mPaddedWidth * mPaddedHeight * mInputDepth * sizeof( float_t );
                }

                return workingMemSize;
            }

            // Randomizes layer's weights, clears biases
            void Randomize( ) override;

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override;

            // Propagates error to the previous layer and calculates weights/biases gradients
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& gradWeights,
                                const XNetworkContext& ctx ) override;

            // Applies updates to the layer's weights and biases
            void UpdateWeights( const fvector_t& updates ) override;

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override;
            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override;
        };

        XConvolutionLayer::XConvolutionLayer( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                                            size_t kernelWidth, size_t kernelHeight, size_t kernelsCount,
                                            const vector<bool>& connectionTable,
                                            BorderMode borderMode, size_t horizontalStep, size_t verticalStep ) :
                                            ITrainableLayer( 0, 0 ),
            mInputWidth( inputWidth ), mInputHeight( inputHeight ), mInputDepth( inputDepth ),
            mOutputWidth( 0 ), mOutputHeight( 0 ),
            mKernelWidth( kernelWidth ), mKernelHeight( kernelHeight ), mKernelsCount( kernelsCount ),
            mHorizontalStep( horizontalStep ), mVerticalStep( verticalStep ), mBorderMode( borderMode ),
            mConnectionTable( connectionTable ), mKernelOffsets( mInputDepth * mKernelsCount ),
            mPaddedWidth( inputWidth ), mPaddedHeight( inputHeight )
        {
            size_t padWidth = 0, padHeight = 0;

            // use input padding, if border handling mode is set to produce same size output
            // (same ouput size will be only when step size is 1; output is always smaller for larger steps)
            if ( mBorderMode == BorderMode::Same )
            {
                padWidth  = mKernelWidth  - 1;
                padHeight = mKernelHeight - 1;

                mPaddedWidth  = mInputWidth  + padWidth;
                mPaddedHeight = mInputHeight + padHeight;
            }

            // calculation of output width/height as:
            //   outSize = ( inSize - kernelSize + padSize ) / step + 1
            mOutputWidth  = ( mInputWidth  - mKernelWidth  + padWidth  ) / mHorizontalStep + 1;
            mOutputHeight = ( mInputHeight - mKernelHeight + padHeight ) / mVerticalStep   + 1;

            // total input/output size
            Initialize( mInputWidth  * mInputHeight  * mInputDepth,
                        mOutputWidth * mOutputHeight * mKernelsCount );

            // invalid or missing connections - assume all output feature maps are built using all input maps
            if ( mConnectionTable.size( ) != mInputDepth * mKernelsCount )
            {
                mConnectionTable = vector<bool>( mInputDepth * mKernelsCount, true );
            }

            // check number of kernels' weights and set offsets
            size_t totalConnectionsCount = 0;

            for ( size_t kernelIndex = 0, connectionIndex = 0; kernelIndex < mKernelsCount; kernelIndex++ )
            {
                for ( size_t inputDepthIndex = 0; inputDepthIndex < mInputDepth; inputDepthIndex++, connectionIndex++ )
                {
                    mKernelOffsets[connectionIndex] = totalConnectionsCount * mKernelWidth * mKernelHeight;

                    if ( mConnectionTable[connectionIndex] )
                    {
                        totalConnectionsCount++;
                    }
                }
            }

            // allocate vector of weights/biases
            mWeightCount = mKernelWidth * mKernelHeight * totalConnectionsCount;
            mAllWeights  = fvector_t( mWeightCount + mKernelsCount );

            // set up weights/biases pointers
            mKernelsWeights = mAllWeights.data( );
            mKernelsBiases  = mKernelsWeights + mWeightCount;

            Randomize( );
        }

        // Randomizes layer's weights, clears biases
        void XConvolutionLayer::Randomize( )
        {
            float halfRange = sqrt( 3.0f / ( mKernelWidth * mKernelHeight * mInputDepth ) );

            for ( size_t i = 0; i < mWeightCount; i++ )
            {
                mKernelsWeights[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRange ) - halfRange;
            }
            for ( size_t i = 0; i < mKernelsCount; i++ )
            {
                mKernelsBiases[i] = 0;
            }
        }

        // Calculates outputs for the given inputs
        void XConvolutionLayer::ForwardCompute( const vector<fvector_t*>& inputs,
                                                vector<fvector_t*>& outputs,
                                                const XNetworkContext& ctx )
        {
            // will be using either original input width/heigh or padded
            size_t  inputWidth   = mInputWidth;
            size_t  inputHeight  = mInputHeight;

            // decide if raw data to be used or padded
            if ( mBorderMode == BorderMode::Same )
            {
                inputWidth  = mPaddedWidth;
                inputHeight = mPaddedHeight;
            }

            size_t  inputRowInc     = inputWidth * mVerticalStep;
            // gap size after processing one input row with kernel to the next row to be processed
            size_t  inputNextRowGap = inputWidth - mKernelWidth;

            // process all samples
            XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
            {
                const float_t* inputData  = inputs[i]->data( );
                float_t*       outputData = outputs[i]->data( );

                if ( mBorderMode == BorderMode::Same )
                {
                    // get working buffer for padded inputs
                    float_t* paddedInput = static_cast<float_t*>( ctx.GetWorkingBuffer( 0, i ) );

                    XDataEncodingTools::AddPadding2d( inputData, paddedInput,
                                                    mInputWidth, mInputHeight, mPaddedWidth, mPaddedHeight,
                                                    mInputDepth, float_t( 0 ) );
                    inputData = paddedInput;
                }

                // clear the output
                fill( outputData, outputData + mOutputsCount, float_t( 0 ) );

                // go through all kernels to build output feature maps
                XParallel::For( mKernelsCount, !ctx.IsTraining( ), [&]( size_t kernelIndex )
                {
                    float_t* outputBase = outputData + kernelIndex * mOutputWidth * mOutputHeight;
                    float_t  biasValue  = mKernelsBiases[kernelIndex];

                    // go through all input layers (or feature maps produced by previous layers)
                    for ( size_t inputDepthIndex = 0; inputDepthIndex < mInputDepth; inputDepthIndex++ )
                    {
                        if ( !mConnectionTable[kernelIndex * mInputDepth + inputDepthIndex] )
                        {
                            // the input map is not used for the output feature map
                            continue;
                        }

                        const float_t* inputBase  = inputData + inputDepthIndex * inputWidth * inputHeight;
                        // get the 2D kernel for current input/output map combination
                        const float_t* kernelBase = mKernelsWeights + mKernelOffsets[kernelIndex * mInputDepth + inputDepthIndex];

                        // calculate output contributions for the current input map
                        for ( size_t oy = 0; oy < mOutputHeight; oy++ )
                        {
                            const float_t* inputRow  = inputBase + oy * inputRowInc;
                            float_t*       outputRow = outputBase + oy * mOutputWidth;

                            for ( size_t ox = 0; ox < mOutputWidth; ox++ )
                            {
                                const float_t* kernelPtr = kernelBase;
                                const float_t* inputPtr  = inputRow;
                                float_t        sum       = float_t( 0 );

                                // "convolve" input with the kernel
                                // (we actually do cross-correlation since it does not matter for CNN)
                                for ( size_t ky = 0; ky < mKernelHeight; ky++ )
                                {
                                    for ( size_t kx = 0; kx < mKernelWidth; kx++ )
                                    {
                                        sum += *inputPtr * *kernelPtr;
                                        kernelPtr++;
                                        inputPtr++;
                                    }
                                    // since input pointer was already shifted horizontally,
                                    // we need to align it back to the start of the next row
                                    inputPtr += inputNextRowGap;
                                }

                                *outputRow += sum;
                                *outputRow += biasValue;

                                // shift output/input row pointers to the next position of the sliding kernel's window
                                outputRow++;
                                inputRow += mHorizontalStep;
                            }
                        }
                    }
                } );
            } );
        }

        // Propagates error to the previous layer and calculates weights/biases gradients
        void XConvolutionLayer::BackwardCompute( const vector<fvector_t*>& inputs,
                                                const vector<fvector_t*>& /* outputs */,
                                                const vector<fvector_t*>& deltas,
                                                vector<fvector_t*>& prevDeltas,
                                                fvector_t& gradWeights,
                                                const XNetworkContext& ctx )
        {
            // set up weights/biases gradients pointers
            float_t* gradWeightsData = gradWeights.data( );
            float_t* gradBiasesData  = gradWeightsData + mWeightCount;
            size_t   outputSize      = mOutputWidth * mOutputHeight;
            // will be using either original input width/heigh or padded
            size_t   inputWidth      = mInputWidth;
            size_t   inputHeight     = mInputHeight;

            // decide if raw data to be used or padded
            if ( mBorderMode == BorderMode::Same )
            {
                inputWidth  = mPaddedWidth;
                inputHeight = mPaddedHeight;
            }

            size_t inputRowInc         = inputWidth * mVerticalStep;
            // gap size after processing one row of previous deltas to the next row to be processed
            size_t prevDeltaNextRowGap = inputWidth - mKernelWidth;
            
            // 1 - first propagate deltas to the previous layer
            XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
            {
                const float_t* deltaData     = deltas[i]->data( );
                float_t*       prevDeltaData = prevDeltas[i]->data( );

                if ( mBorderMode == BorderMode::Same )
                {
                    // get working buffer for padded previous deltas
                    prevDeltaData = static_cast<float_t*>( ctx.GetWorkingBuffer( 1, i ) );
                }

                fill( prevDeltaData, prevDeltaData + mInputsCount, float_t( 0 ) );

                // go through all input feature maps (which are the outputs of the previous layer)
                for ( size_t inputDepthIndex = 0; inputDepthIndex < mInputDepth; inputDepthIndex++ )
                {
                    float_t* prevDeltaBase = prevDeltaData + inputDepthIndex * inputWidth * inputHeight;

                    // go through all kernels, which were applied to the feature map
                    for ( size_t kernelIndex = 0; kernelIndex < mKernelsCount; kernelIndex++ )
                    {
                        if ( !mConnectionTable[kernelIndex * mInputDepth + inputDepthIndex] )
                        {
                            // the input map is not used for the output feature map
                            continue;
                        }

                        const float_t* deltaBase  = deltaData + kernelIndex * outputSize;
                        // get the 2D kernel for then current input/output map combination
                        const float_t* kernelBase = mKernelsWeights + mKernelOffsets[kernelIndex * mInputDepth + inputDepthIndex];

                        // go through the current deltas of the output produced by the current kernel
                        for ( size_t oy = 0; oy < mOutputHeight; oy++ )
                        {
                            const float_t* deltaPtr     = deltaBase     + oy * mOutputWidth;
                            float_t*       prevDeltaRow = prevDeltaBase + oy * inputRowInc;

                            for ( size_t ox = 0; ox < mOutputWidth; ox++ )
                            {
                                const float_t* kernelPtr    = kernelBase;
                                float_t*       prevDeltaPtr = prevDeltaRow;

                                // go through the kernel at current image position
                                for ( size_t ky = 0; ky < mKernelHeight; ky++ )
                                {
                                    for ( size_t kx = 0; kx < mKernelWidth; kx++ )
                                    {
                                        *prevDeltaPtr += *deltaPtr * *kernelPtr;

                                        kernelPtr++;
                                        prevDeltaPtr++;
                                    }
                                    // since previous delta pointer was already shifted horizontally,
                                    // we need to align it back to the start of the next row
                                    prevDeltaPtr += prevDeltaNextRowGap;
                                }

                                // shift current/previous delta pointers to the next position of the sliding kernel's window
                                deltaPtr++;
                                prevDeltaRow += mHorizontalStep;
                            }
                        }
                    }
                }

                if ( mBorderMode == BorderMode::Same )
                {
                    // do unpadding of previous deltas
                    XDataEncodingTools::RemovePadding2d( prevDeltaData, prevDeltas[i]->data( ), mPaddedWidth, mPaddedHeight, mInputWidth, mInputHeight, mInputDepth );
                }
            } );

            // 2 - accumulate weights' difference

            // go through all input feature maps
            XParallel::For( mInputDepth, ctx.IsTraining( ), [&]( size_t inputDepthIndex )
            {
                for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                {
                    const float_t* deltaData = deltas[i]->data( );
                    const float_t* inputData = inputs[i]->data( );

                    if ( mBorderMode == BorderMode::Same )
                    {
                        // get working buffer for padded inputs
                        inputData = static_cast<float_t*>( ctx.GetWorkingBuffer( 0, i ) );
                    }

                    const float_t* inputBase = inputData + inputDepthIndex * inputWidth * inputHeight;

                    // go through all kernels, which were applied to the feature map
                    for ( size_t kernelIndex = 0; kernelIndex < mKernelsCount; kernelIndex++ )
                    {
                        if ( !mConnectionTable[kernelIndex * mInputDepth + inputDepthIndex] )
                        {
                            // the input map is not used for the output feature map
                            continue;
                        }

                        const float_t* deltaBase = deltaData + kernelIndex * outputSize;
                        // get the 2D portion of weights' gradients for the current input/output map combination
                        float_t* gradWeightsPtr  = gradWeightsData + mKernelOffsets[kernelIndex * mInputDepth + inputDepthIndex];

                        // calculate gradients for each weight (kernel element)
                        for ( size_t ky = 0; ky < mKernelHeight; ky++ )
                        {
                            for ( size_t kx = 0; kx < mKernelWidth; kx++ )
                            {
                                float_t sum = float_t( 0 );

                                // multiply output deltas by corresponding inputs
                                for ( size_t oy = 0; oy < mOutputHeight; oy++ )
                                {
                                    const float_t* deltaPtr = deltaBase + oy * mOutputWidth;
                                    const float_t* inputPtr = inputBase + oy * inputRowInc + ky * inputWidth + kx;

                                    for ( size_t ox = 0; ox < mOutputWidth; ox++ )
                                    {
                                        sum += *deltaPtr * *inputPtr;

                                        deltaPtr++;
                                        inputPtr += mHorizontalStep;
                                    }
                                }

                                *gradWeightsPtr += sum;
                                gradWeightsPtr++;
                            }
                        }
                    }
                }
            } );

            // 3 - accumulate baises' difference
            XParallel::For( mKernelsCount, ctx.IsTraining( ), [&]( size_t kernelIndex )
            {
                float_t sum = 0;

                for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                {
                    const float_t* deltaPtr = deltas[i]->data( ) + kernelIndex * outputSize;

                    for ( size_t outputIndex = 0; outputIndex < outputSize; outputIndex++ )
                    {
                        sum += *deltaPtr;
                        deltaPtr++;
                    }
                }

                gradBiasesData[kernelIndex] += sum;
            } );
        }

        // Applies updates to the layer's weights and biases
        void XConvolutionLayer::UpdateWeights( const fvector_t& updates )
        {
            for ( size_t i = 0, n = mAllWeights.size( ); i < n; i++ )
            {
                mAllWeights[i] += updates[i];
            }
        }

        // Saves layer's learnt parameters/weights
        bool XConvolutionLayer::SaveLearnedParams( FILE* file ) const
        {
            vector<const fvector_t*> params( { &mAllWeights } );

            return SaveLearnedParamsHelper( file, LayerID::Convolution, params );
        }

        // Loads layer's learnt parameters
        bool XConvolutionLayer::LoadLearnedParams( FILE* file )
        {
            vector<fvector_t*> params( { &mAllWeights } );

            return LoadLearnedParamsHelper( file, LayerID::Convolution, params );
        }


        // Implementation of Gated Recurrent Unit (GRU) layer
        //
        // http://colah.github.io/posts/2015-08-Understanding-LSTMs/
        // https://towardsdatascience.com/understanding-gru-networks-2ef37df6c9be
        //
        // Some info on GRU backpropagation. NOTE: be careful with that. First link does not
        // provide all equations. While the second has some obvious errors.
        // Some understanding of backpropagation and derivatives math should help though.
        // https://medium.com/swlh/only-numpy-deriving-forward-feed-and-back-propagation-in-gated-recurrent-neural-networks-gru-8b6810f91bad
        // https://cran.r-project.org/web/packages/rnn/vignettes/GRU_units.html
        //
        // Update gate:
        //   Z(t) = sigmoid( Wz [X(t), H(t-1)] + Bz )
        // Reset gate:
        //   R(t) = sigmoid( Wr [X(t), H(t-1)] + Br )
        // Current memory content:
        //   H'(t) = tanh( W [X(t), R(t)*H(t-1)] + B )
        // Output/History:
        //   H(t) = (1 - Z(t)) * H(t-1) + Z(t) * H'(t)
        //
        class XGRULayer : public ITrainableLayer
        {
        private:
            XSigmoidActivation mSigmoid;
            XTanhActivation    mTanh;

            // Weights and biases are all kept together
            fvector_t mAllWeights;

            // And here are the pointers to specific weights/biases, which are used to calculate different
            // vectors from current layer's input, X(t), and its previous output/history, H(t-1):

            // 1) to calculate "update gate" vector, Z(t);
            float_t* mWeightsX2Z;
            float_t* mWeightsH2Z;
            float_t* mBiasesZ;
            // 2) to calculate "reset gate" vector, R(t);
            float_t* mWeightsX2R;
            float_t* mWeightsH2R;
            float_t* mBiasesR;
            // 3) to calculate "current memory content", H'(t);
            float_t* mWeightsX2H;
            float_t* mWeightsHR2H;
            float_t* mBiasesH;

            // --------------------------------------------------------------------------------------
            enum
            {
                // per batch
                BUFFER_INDEX_HISTORY            = 0,
                BUFFER_INDEX_HISTORY_GRAD       = 1,
                BUFFER_INDEX_DELTA              = 2,    // sum of the incoming gradient (from the next layer)
                                                        // and history gradient

                // per sample
                BUFFER_INDEX_HISTORY_PREV       = 3,    // H(t-1)
                BUFFER_INDEX_UPDATE_GATE        = 4,    // Z(t)
                BUFFER_INDEX_RESET_GATE         = 5,    // R(t)
                BUFFER_INDEX_HISTORY_PREV_RESET = 6,    // H(t-1) * R(t)
                BUFFER_INDEX_HISTORY_HAT        = 7,    // H'(t)

                BUFFER_INDEX_UPDATE_GATE_DELTA  = 8,
                BUFFER_INDEX_RESET_GATE_DELTA   = 9,
                BUFFER_INDEX_HISTORY_HAT_DELTA  = 10,
            };

        public:
            XGRULayer( size_t inputsCount, size_t outputsCount );

            // Reports number of weight coefficients the layer has
            size_t WeightsCount( ) const override
            {
                return mAllWeights.size( );
            }

            // Get/set layer's weights
            fvector_t Weights( ) const override
            {
                return mAllWeights;
            }
            void SetWeights( const fvector_t& weights ) override
            {
                mAllWeights = weights;
            }

            // Tells that we may need some extra memory for internal state/calculations
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize = uvector_t( 11, mOutputsCount * sizeof( float_t ) );

                return workingMemSize;
            }

            // Randomizes layer's weights, clears biases (forget gate biases are set to 1 though)
            void Randomize( ) override;

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override;

            // Propagates error to the previous layer and calculates weights/biases gradients
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& gradWeights,
                                const XNetworkContext& ctx ) override;

            // Applies updates to the layer's weights and biases
            void UpdateWeights( const fvector_t& updates ) override;

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override;
            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override;
        };

        XGRULayer::XGRULayer( size_t inputsCount, size_t outputsCount ) :
            ITrainableLayer( inputsCount, outputsCount ),
            mSigmoid( ), mTanh( ),
            mAllWeights( ( inputsCount * outputsCount + outputsCount * outputsCount ) * 3 + outputsCount * 3 )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights pointers
            mWeightsX2Z  = mAllWeights.data( );
            mWeightsH2Z  = mWeightsX2Z + weightsCountInputs;

            mWeightsX2R  = mWeightsH2Z + weightsCountHistory;
            mWeightsH2R  = mWeightsX2R + weightsCountInputs;

            mWeightsX2H  = mWeightsH2R + weightsCountHistory;
            mWeightsHR2H = mWeightsX2H + weightsCountInputs;

            // set up biases pointers
            mBiasesZ = mWeightsHR2H + weightsCountHistory;
            mBiasesR = mBiasesZ + mOutputsCount;
            mBiasesH = mBiasesR + mOutputsCount;

            Randomize( );
        }

        // Randomizes layer's weights, clears biases
        void XGRULayer::Randomize( )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            float halfRangeX = sqrt( 3.0f / mInputsCount );
            float halfRangeH = sqrt( 3.0f / mOutputsCount );

            for ( size_t i = 0; i < weightsCountInputs; i++ )
            {
                mWeightsX2Z[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
                mWeightsX2R[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
                mWeightsX2H[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
            }

            for ( size_t i = 0; i < weightsCountHistory; i++ )
            {
                mWeightsH2Z[i]  = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
                mWeightsH2R[i]  = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
                mWeightsHR2H[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
            }

            // See "Model Parameters" explaining why biases for Reset Gate are set to -1.0
            // https://danijar.com/tips-for-training-recurrent-neural-networks/
            for ( size_t i = 0; i < mOutputsCount; i++ )
            {
                mBiasesZ[i] =  0.0f;
                mBiasesR[i] = -1.0f;
                mBiasesH[i] =  0.0f;
            }
        }

        // Calculates outputs for the given inputs
        void XGRULayer::ForwardCompute( const vector<fvector_t*>& inputs, vector<fvector_t*>& outputs, const XNetworkContext& ctx )
        {
            size_t sequenceLen = ctx.TrainingSequenceLength( );
            size_t batchSize   = inputs.size( ) / sequenceLen;

            XParallel::For( batchSize, ctx.IsTraining( ), [&]( size_t batchIndex )
            {
                float_t* history = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY, batchIndex ) );

                for ( size_t sequenceIndex = 0; sequenceIndex < sequenceLen; sequenceIndex++ )
                {
                    size_t   sampleIndex      = batchIndex * sequenceLen + sequenceIndex;
                    float_t* input            = inputs [sampleIndex]->data( );
                    float_t* output           = outputs[sampleIndex]->data( );   // H(t)

                    float_t* historyPrev      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV, sampleIndex ) );        // H(t-1)
                    float_t* updateGate       = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_UPDATE_GATE, sampleIndex ) );         // Z(t)
                    float_t* resetGate        = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_RESET_GATE, sampleIndex ) );          // R(t)
                    float_t* historyPrevReset = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV_RESET, sampleIndex ) );  // H(t-1) * R(t)
                    float_t* historyHat       = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_HAT, sampleIndex ) );         // H'(t)

                    // remember previous history for this particular sample
                    memcpy( historyPrev, history, mOutputsCount * sizeof( float_t ) );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        // Wz [X(t), H(t-1)] + Bz
                        updateGate[outputIndex] = XVectorize::Dot( input, &( mWeightsX2Z[outputIndex * mInputsCount] ), mInputsCount ) +
                                                XVectorize::Dot( historyPrev, &( mWeightsH2Z[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                mBiasesZ[outputIndex];
                        // Wr [X(t), H(t-1)] + Br
                        resetGate[outputIndex] =  XVectorize::Dot( input, &( mWeightsX2R[outputIndex * mInputsCount] ), mInputsCount ) +
                                                XVectorize::Dot( historyPrev, &( mWeightsH2R[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                mBiasesR[outputIndex];
                        // W [X(t)] + B
                        historyHat[outputIndex] = XVectorize::Dot( input, &( mWeightsX2H[outputIndex * mInputsCount] ), mInputsCount ) +
                                                mBiasesH[outputIndex];
                    }

                    // apply activations
                    mSigmoid.ForwardActivate( updateGate, updateGate, mOutputsCount );
                    mSigmoid.ForwardActivate( resetGate, resetGate, mOutputsCount );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        // previous history multiplied by reset gate
                        historyPrevReset[outputIndex] = historyPrev[outputIndex] * resetGate[outputIndex];

                        // first part of the ouput: (1 - Z(t)) * H(t-1)
                        output[outputIndex] = historyPrev[outputIndex] * ( float_t( 1 ) - updateGate[outputIndex] );
                    }

                    // complete current memory content by adding reseted previous history ...
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        historyHat[outputIndex] += XVectorize::Dot( historyPrevReset, &( mWeightsHR2H[outputIndex * mOutputsCount] ), mOutputsCount );
                    }
                    // ... and passing through tanh() activation
                    mTanh.ForwardActivate( historyHat, historyHat, mOutputsCount );

                    // get the final output and put into history as well
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        // second part of the ouput: Z(t) * H'(t)
                        output[outputIndex] += historyHat[outputIndex] * updateGate[outputIndex];
                        history[outputIndex] = output[outputIndex];
                    }
                }
            } );
        }

        // Propagates error to the previous layer and calculates weights/biases gradients
        void XGRULayer::BackwardCompute( const vector<fvector_t*>& inputs,
                                        const vector<fvector_t*>& /* outputs */,
                                        const vector<fvector_t*>& deltas,
                                        vector<fvector_t*>& prevDeltas,
                                        fvector_t& gradWeights,
                                        const XNetworkContext& ctx )
        {
            size_t sequenceLen   = ctx.TrainingSequenceLength( );
            size_t batchSize     = inputs.size( ) / sequenceLen;

            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights gradient pointers
            float_t* gradWeightsX2Z  = gradWeights.data( );
            float_t* gradWeightsH2Z  = gradWeightsX2Z + weightsCountInputs;

            float_t* gradWeightsX2R  = gradWeightsH2Z + weightsCountHistory;
            float_t* gradWeightsH2R  = gradWeightsX2R + weightsCountInputs;

            float_t* gradWeightsX2H  = gradWeightsH2R + weightsCountHistory;
            float_t* gradWeightsHR2H = gradWeightsX2H + weightsCountInputs;

            // set up biases gradient pointers
            float_t* gradBiasesZ = gradWeightsHR2H + weightsCountHistory;
            float_t* gradBiasesR = gradBiasesZ + mOutputsCount;
            float_t* gradBiasesH = gradBiasesR + mOutputsCount;

            XParallel::For( batchSize, [&]( size_t batchIndex )
            {
                float_t* historyGrad = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_GRAD, batchIndex ) );
                float_t* delta       = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_DELTA, batchIndex ) );

                for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                {
                    size_t   sampleIndex = batchIndex * sequenceLen + sequenceIndex;
                    float_t* prevDelta   = prevDeltas[sampleIndex]->data( );

                    float_t* historyPrev = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV, sampleIndex ) );        // H(t-1)
                    float_t* updateGate  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_UPDATE_GATE, sampleIndex ) );         // Z(t)
                    float_t* resetGate   = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_RESET_GATE, sampleIndex ) );          // R(t)
                    float_t* historyHat  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_HAT, sampleIndex ) );         // H'(t)

                    float_t* dUpdateGate = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_UPDATE_GATE_DELTA, sampleIndex ) );
                    float_t* dResetGate  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_RESET_GATE_DELTA, sampleIndex ) );
                    float_t* dHistoryHat = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_HAT_DELTA, sampleIndex ) );

                    // add history gradient from the future
                    memcpy( delta, deltas[sampleIndex]->data( ), sizeof( float_t ) * mOutputsCount );
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        delta[outputIndex] += historyGrad[outputIndex];
                    }

                    // dE/dWz
                    // dH(t)/dZ(t) = H'(t) - H(t-1)
                    // delta * ( H'(t) - H(t-1) )
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dUpdateGate[outputIndex] = delta[outputIndex] * ( historyHat[outputIndex] - historyPrev[outputIndex] );
                    }
                    mSigmoid.BackwardActivate( updateGate, updateGate, dUpdateGate, dUpdateGate, mOutputsCount );

                    // dE/dWh
                    // dH(t)/dH'(t) = Z(t)
                    // delta * Z(t)
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dHistoryHat[outputIndex] = delta[outputIndex] * updateGate[outputIndex];
                    }
                    mTanh.BackwardActivate( historyHat, historyHat, dHistoryHat, dHistoryHat, mOutputsCount );

                    // progress with reset gate gradient
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        float_t weightedGradHistoryHat = float_t( 0 );

                        for ( size_t outputIndex2 = 0, weightIndex = outputIndex; outputIndex2 < mOutputsCount; outputIndex2++, weightIndex += mOutputsCount )
                        {
                            weightedGradHistoryHat += dHistoryHat[outputIndex2] * mWeightsHR2H[weightIndex];
                        }

                        // multiply with previous history to find reset gradient then
                        dResetGate[outputIndex]  = weightedGradHistoryHat * historyPrev[outputIndex];
                        // multiply with reset gate value to direct error gradient to previous history gradient
                        historyGrad[outputIndex] = weightedGradHistoryHat * resetGate[outputIndex];
                    }
                    mSigmoid.BackwardActivate( resetGate, resetGate, dResetGate, dResetGate, mOutputsCount );
                    
                    // add error gradient from output to previous history gradient
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        historyGrad[outputIndex] += ( ( 1 - updateGate[outputIndex] ) * delta[outputIndex] );
                    }

                    // input deltas for the previous layer
                    for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                    {
                        size_t  weightIndex = inputIndex;
                        float_t sum         = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mInputsCount )
                        {
                            sum += dUpdateGate[outputIndex] * mWeightsX2Z[weightIndex];
                            sum += dResetGate[outputIndex]  * mWeightsX2R[weightIndex];
                            sum += dHistoryHat[outputIndex] * mWeightsX2H[weightIndex];
                        }

                        prevDelta[inputIndex] = sum;
                    }

                    // add more to history gradient for the previous sequence of this layer
                    for ( size_t outputIndex2 = 0; outputIndex2 < mOutputsCount; outputIndex2++ )
                    {
                        size_t  weightIndex = outputIndex2;
                        float_t sum         = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mOutputsCount )
                        {
                            sum += dUpdateGate[outputIndex]  * mWeightsH2Z[weightIndex];
                            sum += dResetGate[outputIndex]   * mWeightsH2R[weightIndex];
                        }

                        historyGrad[outputIndex2] += sum;
                    }
                }
            } );

            XParallel::For( mOutputsCount, [&]( size_t outputIndex )
            {
                size_t weightIndexStartI = outputIndex * mInputsCount;
                size_t weightIndexStartH = outputIndex * mOutputsCount;

                for ( size_t batchIndex = 0; batchIndex < batchSize; batchIndex++ )
                {
                    for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                    {
                        size_t   sampleIndex      = batchIndex * sequenceLen + sequenceIndex;
                        const float_t* input      = inputs[sampleIndex]->data( );

                        float_t* historyPrev      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV, sampleIndex ) );        // H(t-1)
                        float_t* historyPrevReset = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV_RESET, sampleIndex ) );  // H(t-1) * R(t)

                        float_t* dUpdateGate      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_UPDATE_GATE_DELTA, sampleIndex ) );
                        float_t* dResetGate       = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_RESET_GATE_DELTA, sampleIndex ) );
                        float_t* dHistoryHat      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_HAT_DELTA, sampleIndex ) );

                        float_t  dUpdateGateVal   = dUpdateGate[outputIndex];
                        float_t  dResetGateVal    = dResetGate[outputIndex];
                        float_t  dHistoryHatVal   = dHistoryHat[outputIndex];

                        // accumulate gradients for inputs' weights
                        for ( size_t inputIndex = 0, weightIndex = weightIndexStartI; inputIndex < mInputsCount; inputIndex++, weightIndex++ )
                        {
                            gradWeightsX2Z[weightIndex] += dUpdateGateVal * input[inputIndex];
                            gradWeightsX2R[weightIndex] += dResetGateVal  * input[inputIndex];
                            gradWeightsX2H[weightIndex] += dHistoryHatVal * input[inputIndex];
                        }

                        // accumulate gradients for history weights
                        if ( sequenceIndex != 0 )
                        {
                            for ( size_t historyIndex = 0, weightIndex = weightIndexStartH; historyIndex < mOutputsCount; historyIndex++, weightIndex++ )
                            {
                                gradWeightsH2Z[weightIndex]  += dUpdateGateVal * historyPrev[historyIndex];
                                gradWeightsH2R[weightIndex]  += dResetGateVal  * historyPrev[historyIndex];
                                gradWeightsHR2H[weightIndex] += dHistoryHatVal * historyPrevReset[historyIndex];
                            }
                        }

                        // accumulate gradients for biases
                        gradBiasesZ[outputIndex] += dUpdateGateVal;
                        gradBiasesR[outputIndex] += dResetGateVal;
                        gradBiasesH[outputIndex] += dHistoryHatVal;
                    }
                }
            } );
        }

        // Applies updates to the layer's weights and biases
        void XGRULayer::UpdateWeights( const fvector_t& updates )
        {
            for ( size_t i = 0, n = mAllWeights.size( ); i < n; i++ )
            {
                mAllWeights[i] += updates[i];
            }
        }

        // Saves layer's learnt parameters/weights
        bool XGRULayer::SaveLearnedParams( FILE* file ) const
        {
            vector<const fvector_t*> params( { &mAllWeights } );

            return SaveLearnedParamsHelper( file, LayerID::RecurrentGRU, params );
        }

        // Loads layer's learnt parameters
        bool XGRULayer::LoadLearnedParams( FILE* file )
        {
            vector<fvector_t*> params( { &mAllWeights } );

            return LoadLearnedParamsHelper( file, LayerID::RecurrentGRU, params );
        }

        // Implementation of Long Short-Term Memory (LSTM) recurrent layer
        //
        // http://colah.github.io/posts/2015-08-Understanding-LSTMs/
        // https://medium.com/@aidangomez/let-s-do-this-f9b699de31d9
        //
        // Forget gate:
        //   F(t) = sigmoid( Wf [X(t), H(t-1)] + Bf )
        // Input gate:
        //   I(t) = sigmoid( Wi [X(t), H(t-1)] + Bi )
        // Candidate state:
        //   Z(t) = tanh( Wz [X(t), H(t-1)] + Bz )
        // State:
        //   C(t) = F(t) * C(t-1) + I(t) * Z(t)
        // Output gate:
        //   O(t) = sigmoid( Wo [X(t), H(t-1)] + Bo )
        // Output/History:
        //   H(t) = O(t) * tanh( C(t) )
        //
        class XLSTMLayer : public ITrainableLayer
        {
        private:
            XSigmoidActivation mSigmoid;
            XTanhActivation    mTanh;

            // Weights and biases are all kept together
            fvector_t mAllWeights;

            // And here are the pointers to specific weights/biases, which are used to calculate different
            // vectors from current layer's input, X(t), and its previous output/history, H(t-1):

            // 1) to calculate "forget gate" vector, F(t);
            float_t* mWeightsX2F;
            float_t* mWeightsH2F;
            float_t* mBiasesF;
            // 2) to calculate "input gate" vector, I(t);
            float_t* mWeightsX2I;
            float_t* mWeightsH2I;
            float_t* mBiasesI;
            // 3) to calculate "candidate state" vector, Z(t) (C with tilde);
            float_t* mWeightsX2Z;
            float_t* mWeightsH2Z;
            float_t* mBiasesZ;
            // 4) to calculate "output gate" vector, O(t)
            float_t* mWeightsX2O;
            float_t* mWeightsH2O;
            float_t* mBiasesO;

            // --------------------------------------------------------------------------------------
            enum
            {
                // per batch
                BUFFER_INDEX_STATE           = 0,
                BUFFER_INDEX_STATE_GRAD      = 1,
                BUFFER_INDEX_HISTORY         = 2,
                BUFFER_INDEX_HISTORY_GRAD    = 3,
                BUFFER_INDEX_DELTA           = 4,   // sum of the incoming gradient (from the next layer)
                                                    // and history gradient
                BUFFER_INDEX_STATE_DELTA     = 5,   // BUFFER_INDEX_DELTA passed backward through output gate 

                // per sample
                BUFFER_INDEX_STATE_PREV      = 6,   // C(t-1)
                BUFFER_INDEX_STATE_NEXT      = 7,   // C(t)
                BUFFER_INDEX_HISTORY_PREV    = 8,   // H(t-1)
                BUFFER_INDEX_FORGET_GATE     = 9,   // F(t)
                BUFFER_INDEX_INPUT_GATE      = 10,  // I(t)
                BUFFER_INDEX_OUTPUT_GATE     = 11,  // O(t)
                BUFFER_INDEX_CANDIDATE_STATE = 12,  // Z(t)
                BUFFER_INDEX_STATE_NEXT_TANH = 13,  // tanh(C(t))

                BUFFER_INDEX_CANDIDATE_STATE_DELTA = 14,
                BUFFER_INDEX_INPUT_GATE_DELTA      = 15,
                BUFFER_INDEX_FORGET_GATE_DELTA     = 16,
                BUFFER_INDEX_OUTPUT_GATE_DELTA     = 17,
            };

        public:

            XLSTMLayer( size_t inputsCount, size_t outputsCount );

            // Reports number of weight coefficients the layer has
            size_t WeightsCount( ) const override
            {
                return mAllWeights.size( );
            }

            // Get/set layer's weights
            fvector_t Weights( ) const override
            {
                return mAllWeights;
            }
            void SetWeights( const fvector_t& weights ) override
            {
                mAllWeights = weights;
            }

            // Tells that we may need some extra memory for internal state/calculations
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize = uvector_t( 18, mOutputsCount * sizeof( float_t ) );

                return workingMemSize;
            }

            // Randomizes layer's weights, clears biases (forget gate biases are set to 1 though)
            void Randomize( ) override;

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override;

            // Propagates error to the previous layer and calculates weights/biases gradients
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& gradWeights,
                                const XNetworkContext& ctx ) override;

            // Applies updates to the layer's weights and biases
            void UpdateWeights( const fvector_t& updates ) override;

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override;
            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override;
        };


        XLSTMLayer::XLSTMLayer( size_t inputsCount, size_t outputsCount ) :
            ITrainableLayer( inputsCount, outputsCount ),
            mSigmoid( ), mTanh( ),
            mAllWeights( ( inputsCount * outputsCount + outputsCount * outputsCount ) * 4 + outputsCount * 4 )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights pointers
            mWeightsX2F = mAllWeights.data( );
            mWeightsH2F = mWeightsX2F + weightsCountInputs;

            mWeightsX2I = mWeightsH2F + weightsCountHistory;
            mWeightsH2I = mWeightsX2I + weightsCountInputs;

            mWeightsX2Z = mWeightsH2I + weightsCountHistory;
            mWeightsH2Z = mWeightsX2Z + weightsCountInputs;

            mWeightsX2O = mWeightsH2Z + weightsCountHistory;
            mWeightsH2O = mWeightsX2O + weightsCountInputs;

            // set up biases pointers
            mBiasesF = mWeightsH2O + weightsCountHistory;
            mBiasesI = mBiasesF + mOutputsCount;
            mBiasesZ = mBiasesI + mOutputsCount;
            mBiasesO = mBiasesZ + mOutputsCount;

            Randomize( );
        }

        // Randomizes layer's weights, clears biases
        void XLSTMLayer::Randomize( )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            float halfRangeX = sqrt( 3.0f / mInputsCount );
            float halfRangeH = sqrt( 3.0f / mOutputsCount );

            for ( size_t i = 0; i < weightsCountInputs; i++ )
            {
                mWeightsX2F[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
                mWeightsX2I[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
                mWeightsX2Z[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
                mWeightsX2O[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
            }

            for ( size_t i = 0; i < weightsCountHistory; i++ )
            {
                mWeightsH2F[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
                mWeightsH2I[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
                mWeightsH2Z[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
                mWeightsH2O[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
            }

            // See "Model Parameters" explaining why biases for Forget Gate are set to 1.0
            // https://danijar.com/tips-for-training-recurrent-neural-networks/
            for ( size_t i = 0; i < mOutputsCount; i++ )
            {
                mBiasesF[i] = 1.0f;
                mBiasesI[i] = 0.0f;
                mBiasesZ[i] = 0.0f;
                mBiasesO[i] = 0.0f;
            }
        }

        // Calculates outputs for the given inputs
        void XLSTMLayer::ForwardCompute( const vector<fvector_t*>& inputs, vector<fvector_t*>& outputs, const XNetworkContext& ctx )
        {
            size_t sequenceLen = ctx.TrainingSequenceLength( );
            size_t batchSize   = inputs.size( ) / sequenceLen;

            XParallel::For( batchSize, ctx.IsTraining( ), [&]( size_t batchIndex )
            {
                float_t* state   = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE,   batchIndex ) );
                float_t* history = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY, batchIndex ) );

                for ( size_t sequenceIndex = 0; sequenceIndex < sequenceLen; sequenceIndex++ )
                {
                    size_t   sampleIndex    = batchIndex * sequenceLen + sequenceIndex;
                    float_t* input          = inputs [sampleIndex]->data( );
                    float_t* output         = outputs[sampleIndex]->data( );   // H(t)
                    float_t* statePrev      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_PREV, sampleIndex ) );      // C(t-1)
                    float_t* stateNext      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_NEXT, sampleIndex ) );      // C(t)
                    float_t* historyPrev    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV, sampleIndex ) );    // H(t-1)
                    float_t* forgetGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_FORGET_GATE, sampleIndex ) );     // F(t)
                    float_t* inputGate      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_INPUT_GATE, sampleIndex ) );      // I(t)
                    float_t* candidateState = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_CANDIDATE_STATE, sampleIndex ) ); // Z(t)
                    float_t* outputGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_OUTPUT_GATE, sampleIndex ) );     // O(t)
                    float_t* stateNextTanh  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_NEXT_TANH, sampleIndex ) ); // tanh(C(t))

                    // remember previous state/history for this particular sample
                    memcpy( statePrev, state, mOutputsCount * sizeof( float_t ) );
                    memcpy( historyPrev, history, mOutputsCount * sizeof( float_t ) );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        // Wf [X(t), H(t-1)] + Bf
                        forgetGate[outputIndex]     = XVectorize::Dot( input, &( mWeightsX2F[outputIndex * mInputsCount] ), mInputsCount ) +
                                                    XVectorize::Dot( historyPrev, &( mWeightsH2F[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                    mBiasesF[outputIndex];
                        // Wi [X(t), H(t-1)] + Bi
                        inputGate[outputIndex]      = XVectorize::Dot( input, &( mWeightsX2I[outputIndex * mInputsCount] ), mInputsCount ) +
                                                    XVectorize::Dot( historyPrev, &( mWeightsH2I[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                    mBiasesI[outputIndex];
                        // Wz [X(t), H(t-1)] + Bz
                        candidateState[outputIndex] = XVectorize::Dot( input, &( mWeightsX2Z[outputIndex * mInputsCount] ), mInputsCount ) +
                                                    XVectorize::Dot( historyPrev, &( mWeightsH2Z[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                    mBiasesZ[outputIndex];
                        // Wo [X(t), H(t-1)] + Bo
                        outputGate[outputIndex]     = XVectorize::Dot( input, &( mWeightsX2O[outputIndex * mInputsCount] ), mInputsCount ) +
                                                    XVectorize::Dot( historyPrev, &( mWeightsH2O[outputIndex * mOutputsCount] ), mOutputsCount ) +
                                                    mBiasesO[outputIndex];
                    }

                    // apply activations
                    mSigmoid.ForwardActivate( forgetGate, forgetGate, mOutputsCount );
                    mSigmoid.ForwardActivate( inputGate, inputGate, mOutputsCount );
                    mSigmoid.ForwardActivate( outputGate, outputGate, mOutputsCount );
                    mTanh.ForwardActivate( candidateState, candidateState, mOutputsCount );

                    // get the new state: C(t) = F(t) * C(t-1) + I(t) * Z(t)
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        state[outputIndex]     =
                        stateNext[outputIndex] = forgetGate[outputIndex] * statePrev[outputIndex] +
                                                inputGate[outputIndex]  * candidateState[outputIndex];
                    }

                    // get the tanh(C(t))
                    mTanh.ForwardActivate( stateNext, stateNextTanh, mOutputsCount );

                    // finally get the next output and keep it into history
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        history[outputIndex] =
                        output[outputIndex]  = outputGate[outputIndex] * stateNextTanh[outputIndex];
                    }
                }
            } );
        }

        // Propagates error to the previous layer and calculates weights/biases gradients
        void XLSTMLayer::BackwardCompute( const vector<fvector_t*>& inputs,
                                        const vector<fvector_t*>& /* outputs */,
                                        const vector<fvector_t*>& deltas,
                                        vector<fvector_t*>& prevDeltas,
                                        fvector_t& gradWeights,
                                        const XNetworkContext& ctx )
        {
            size_t sequenceLen   = ctx.TrainingSequenceLength( );
            size_t batchSize     = inputs.size( ) / sequenceLen;

            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights gradient pointers
            float_t* gradWeightsX2F = gradWeights.data( );
            float_t* gradWeightsH2F = gradWeightsX2F + weightsCountInputs;

            float_t* gradWeightsX2I = gradWeightsH2F + weightsCountHistory;
            float_t* gradWeightsH2I = gradWeightsX2I + weightsCountInputs;

            float_t* gradWeightsX2Z = gradWeightsH2I + weightsCountHistory;
            float_t* gradWeightsH2Z = gradWeightsX2Z + weightsCountInputs;

            float_t* gradWeightsX2O = gradWeightsH2Z + weightsCountHistory;
            float_t* gradWeightsH2O = gradWeightsX2O + weightsCountInputs;

            // set up biases gradient pointers
            float_t* gradBiasesF = gradWeightsH2O + weightsCountHistory;
            float_t* gradBiasesI = gradBiasesF + mOutputsCount;
            float_t* gradBiasesZ = gradBiasesI + mOutputsCount;
            float_t* gradBiasesO = gradBiasesZ + mOutputsCount;

            XParallel::For( batchSize, [&]( size_t batchIndex )
            {
                float_t* stateGrad   = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_GRAD, batchIndex ) );
                float_t* historyGrad = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_GRAD, batchIndex ) );
                float_t* delta       = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_DELTA, batchIndex ) );
                float_t* dState      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_DELTA, batchIndex ) );

                for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                {
                    size_t   sampleIndex    = batchIndex * sequenceLen + sequenceIndex;
                    float_t* prevDelta      = prevDeltas[sampleIndex]->data( );

                    float_t* statePrev      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_PREV, sampleIndex ) );      // C(t-1)
                    float_t* stateNext      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_NEXT, sampleIndex ) );      // C(t)
                    float_t* forgetGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_FORGET_GATE, sampleIndex ) );     // F(t)
                    float_t* inputGate      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_INPUT_GATE, sampleIndex ) );      // I(t)
                    float_t* candidateState = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_CANDIDATE_STATE, sampleIndex ) ); // Z(t)
                    float_t* outputGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_OUTPUT_GATE, sampleIndex ) );     // O(t)
                    float_t* stateNextTanh  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_NEXT_TANH, sampleIndex ) ); // tanh(C(t))

                    float_t* dCadidateState = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_CANDIDATE_STATE_DELTA, sampleIndex ) );
                    float_t* dInputGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_INPUT_GATE_DELTA, sampleIndex ) );
                    float_t* dForgetGate    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_FORGET_GATE_DELTA, sampleIndex ) );
                    float_t* dOutputGate    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_OUTPUT_GATE_DELTA, sampleIndex ) );

                    // add history gradient from the future
                    memcpy( delta, deltas[sampleIndex]->data( ), sizeof( float_t ) * mOutputsCount );
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        delta[outputIndex] += historyGrad[outputIndex];

                        // precalculate other history based gradients
                        dState[outputIndex]      = delta[outputIndex];
                        dOutputGate[outputIndex] = delta[outputIndex];
                    }

                    // pass deltas backward through output gate ...
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dState[outputIndex] *= outputGate[outputIndex];
                    }
                    // ... and through tanh() activation
                    mTanh.BackwardActivate( stateNext, stateNextTanh, dState, dState, mOutputsCount );

                    // add state gradient from the future
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dState[outputIndex] += stateGrad[outputIndex];

                        // precalculate other state based gradients
                        dCadidateState[outputIndex] = dState[outputIndex];
                        dInputGate[outputIndex]     = dState[outputIndex];
                        dForgetGate[outputIndex]    = dState[outputIndex];
                    }

                    // pass state gradient backward through forget gate, so it is ready for the previous sample in the time series
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        stateGrad[outputIndex] = dState[outputIndex] * forgetGate[outputIndex];
                    }
                    
                    // pass state gradient backward through input gate and tanh() activation to get candidate state gradient
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dCadidateState[outputIndex] *= inputGate[outputIndex];
                    }
                    mTanh.BackwardActivate( candidateState, candidateState, dCadidateState, dCadidateState, mOutputsCount );

                    // input gate gradients
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dInputGate[outputIndex] *= candidateState[outputIndex];
                    }
                    mSigmoid.BackwardActivate( inputGate, inputGate, dInputGate, dInputGate, mOutputsCount );

                    // forget gate gradients
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dForgetGate[outputIndex] *= statePrev[outputIndex];
                    }
                    mSigmoid.BackwardActivate( forgetGate, forgetGate, dForgetGate, dForgetGate, mOutputsCount );

                    // output gate gradients
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        dOutputGate[outputIndex] *= stateNextTanh[outputIndex];
                    }
                    mSigmoid.BackwardActivate( outputGate, outputGate, dOutputGate, dOutputGate, mOutputsCount );

                    // calculate gradients to pass to the previous layer of the network
                    for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                    {
                        size_t  weightIndex = inputIndex;
                        float_t sum         = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mInputsCount )
                        {
                            sum += dForgetGate[outputIndex]    * mWeightsX2F[weightIndex];
                            sum += dInputGate[outputIndex]     * mWeightsX2I[weightIndex];
                            sum += dOutputGate[outputIndex]    * mWeightsX2O[weightIndex];
                            sum += dCadidateState[outputIndex] * mWeightsX2Z[weightIndex];
                        }

                        prevDelta[inputIndex] = sum;
                    }

                    // calculate gradients to pass to the previous sample of the time series
                    for ( size_t outputIndex2 = 0; outputIndex2 < mOutputsCount; outputIndex2++ )
                    {
                        size_t  weightIndex = outputIndex2;
                        float_t sum         = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mOutputsCount )
                        {
                            sum += dForgetGate[outputIndex]    * mWeightsH2F[weightIndex];
                            sum += dInputGate[outputIndex]     * mWeightsH2I[weightIndex];
                            sum += dOutputGate[outputIndex]    * mWeightsH2O[weightIndex];
                            sum += dCadidateState[outputIndex] * mWeightsH2Z[weightIndex];
                        }

                        historyGrad[outputIndex2] = sum;
                    }
                }
            } );

            XParallel::For( mOutputsCount, [&]( size_t outputIndex )
            {
                size_t weightIndexStartI = outputIndex * mInputsCount;
                size_t weightIndexStartH = outputIndex * mOutputsCount;

                for ( size_t batchIndex = 0; batchIndex < batchSize; batchIndex++ )
                {
                    for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                    {
                        size_t   sampleIndex    = batchIndex * sequenceLen + sequenceIndex;
                        const float_t* input    = inputs[sampleIndex]->data( );
                        float_t* historyPrev    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_HISTORY_PREV, sampleIndex ) );    // H(t-1)

                        float_t* dCadidateState = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_CANDIDATE_STATE_DELTA, sampleIndex ) );
                        float_t* dInputGate     = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_INPUT_GATE_DELTA, sampleIndex ) );
                        float_t* dForgetGate    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_FORGET_GATE_DELTA, sampleIndex ) );
                        float_t* dOutputGate    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_OUTPUT_GATE_DELTA, sampleIndex ) );

                        float_t dInputGateVal     = dInputGate[outputIndex];
                        float_t dCadidateStateVal = dCadidateState[outputIndex];
                        float_t dForgetGateVal    = dForgetGate[outputIndex];
                        float_t dOutputGateVal    = dOutputGate[outputIndex];

                        // accumulate gradients for inputs' weights
                        for ( size_t inputIndex = 0, weightIndex = weightIndexStartI; inputIndex < mInputsCount; inputIndex++, weightIndex++ )
                        {
                            gradWeightsX2F[weightIndex] += dForgetGateVal    * input[inputIndex];
                            gradWeightsX2I[weightIndex] += dInputGateVal     * input[inputIndex];
                            gradWeightsX2Z[weightIndex] += dCadidateStateVal * input[inputIndex];
                            gradWeightsX2O[weightIndex] += dOutputGateVal    * input[inputIndex];
                        }

                        // accumulate gradients for history weights
                        if ( sequenceIndex != 0 )
                        {
                            for ( size_t historyIndex = 0, weightIndex = weightIndexStartH; historyIndex < mOutputsCount; historyIndex++, weightIndex++ )
                            {
                                gradWeightsH2F[weightIndex] += dForgetGateVal    * historyPrev[historyIndex];
                                gradWeightsH2I[weightIndex] += dInputGateVal     * historyPrev[historyIndex];
                                gradWeightsH2Z[weightIndex] += dCadidateStateVal * historyPrev[historyIndex];
                                gradWeightsH2O[weightIndex] += dOutputGateVal    * historyPrev[historyIndex];
                            }
                        }

                        // accumulate gradients for biases
                        gradBiasesF[outputIndex] += dForgetGateVal;
                        gradBiasesI[outputIndex] += dInputGateVal;
                        gradBiasesZ[outputIndex] += dCadidateStateVal;
                        gradBiasesO[outputIndex] += dOutputGateVal;
                    }
                }
            } );
        }

        // Applies updates to the layer's weights and biases
        void XLSTMLayer::UpdateWeights( const fvector_t& updates )
        {
            for ( size_t i = 0, n = mAllWeights.size( ); i < n; i++ )
            {
                mAllWeights[i] += updates[i];
            }
        }

        // Saves layer's learnt parameters/weights
        bool XLSTMLayer::SaveLearnedParams( FILE* file ) const
        {
            vector<const fvector_t*> params( { &mAllWeights } );

            return SaveLearnedParamsHelper( file, LayerID::RecurrentLSTM, params );
        }

        // Loads layer's learnt parameters
        bool XLSTMLayer::LoadLearnedParams( FILE* file )
        {
            vector<fvector_t*> params( { &mAllWeights } );

            return LoadLearnedParamsHelper( file, LayerID::RecurrentLSTM, params );
        }

        // Implementation of simple recurent layer, which perform calculations as:
        //
        //  Internal activation:  A(t) = U * X(t) + W * H(t-1) + B
        //  Output/History:       H(t) = tanh(A(t))
        //
        // See: http://www.deeplearningbook.org/contents/rnn.html
        //
        class XRecurrentLayer : public ITrainableLayer
        {
        private:
            XTanhActivation    mTanh;

            // Weights and biases are all kept together
            fvector_t   mAllWeights;

            // And here are the pointers to specific weights/biases
            float_t*    mWeightsU;
            float_t*    mWeightsW;

            float_t*    mBiasesB;
            
            // --------------------------------------------------------------------------------------
            enum
            {
                // per batch
                BUFFER_INDEX_STATE          = 0,
                BUFFER_INDEX_STATE_GRAD     = 1,

                // per sample
                BUFFER_INDEX_STATE_PREV           = 2, // H(t-1)
                BUFFER_INDEX_STATE_CURRENT        = 3, // H(t)
                BUFFER_INDEX_STATE_DELTA_CURRENT  = 4, // state delta for the current sample
            };

        public:

            XRecurrentLayer( size_t inputsCount, size_t outputsCount );

            // Reports number of weight coefficients the layer has
            size_t WeightsCount( ) const override
            {
                return mAllWeights.size( );
            }

            // Get/set layer's weights
            fvector_t Weights( ) const override
            {
                return mAllWeights;
            }
            void SetWeights( const fvector_t& weights ) override
            {
                mAllWeights = weights;
            }

            // Tells that we may need some extra memory for internal state/calculations
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize = uvector_t( 5, mOutputsCount * sizeof( float_t ) );

                return workingMemSize;
            }

            // Randomizes layer's weights, clears biases
            void Randomize( ) override;

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override;

            // Propagates error to the previous layer and calculates weights/biases gradients
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& gradWeights,
                                const XNetworkContext& ctx ) override;

            // Applies updates to the layer's weights and biases
            void UpdateWeights( const fvector_t& updates ) override;

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override;
            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override;
        };

        XRecurrentLayer::XRecurrentLayer( size_t inputsCount, size_t outputsCount ) :
            ITrainableLayer( inputsCount, outputsCount ),
            mTanh( ),
            mAllWeights( ( inputsCount + outputsCount ) * outputsCount  + outputsCount )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights pointers
            mWeightsU = mAllWeights.data( );
            mWeightsW = mWeightsU + weightsCountInputs;

            // set up biases pointers
            mBiasesB  = mWeightsW + weightsCountHistory;

            Randomize( );
        }

        // Randomizes layer's weights, clears biases
        void XRecurrentLayer::Randomize( )
        {
            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            float halfRangeX = sqrt( 3.0f / mInputsCount );
            float halfRangeH = sqrt( 3.0f / mOutputsCount );

            for ( size_t i = 0; i < weightsCountInputs; i++ )
            {
                mWeightsU[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeX ) - halfRangeX;
            }

            for ( size_t i = 0; i < weightsCountHistory; i++ )
            {
                mWeightsW[i] = ( static_cast<float_t>( rand( ) ) / RAND_MAX ) * ( float_t( 2 ) * halfRangeH ) - halfRangeH;
            }

            for ( size_t i = 0; i < mOutputsCount; i++ )
            {
                mBiasesB[i] = 0;
            }
        }

        // Calculates outputs for the given inputs
        void XRecurrentLayer::ForwardCompute( const vector<fvector_t*>& inputs,
                                            vector<fvector_t*>& outputs,
                                            const XNetworkContext& ctx )
        {
            size_t sequenceLen = ctx.TrainingSequenceLength( );
            size_t batchSize   = inputs.size( ) / sequenceLen;

            XParallel::For( batchSize, ctx.IsTraining( ), [&]( size_t batchIndex )
            {
                float_t* state = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE, batchIndex ) );

                for ( size_t sequenceIndex = 0; sequenceIndex < sequenceLen; sequenceIndex++ )
                {
                    size_t    sampleIndex  = batchIndex * sequenceLen + sequenceIndex;
                    float_t*  input        = inputs [sampleIndex]->data( );
                    float_t*  output       = outputs[sampleIndex]->data( );

                    float_t*  statePrev    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_PREV, sampleIndex ) );    // H(t-1)
                    float_t*  stateCurrent = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_CURRENT, sampleIndex ) ); // H(t)

                    // remember previous state for this particular sample
                    memcpy( statePrev, state, mOutputsCount * sizeof( float_t ) );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        state[outputIndex] =
                            // X(t) * U
                            XVectorize::Dot( input, &( mWeightsU[outputIndex * mInputsCount] ), mInputsCount ) +
                            // H(t-1) * W
                            XVectorize::Dot( statePrev, &( mWeightsW[outputIndex * mOutputsCount] ), mOutputsCount ) +
                            // B
                            mBiasesB[outputIndex];
                    }

                    // apply tanh() to get the final H(t)
                    mTanh.ForwardActivate( state, state, mOutputsCount );

                    // copy state to output
                    memcpy( output, state, mOutputsCount * sizeof( float_t ) );

                    // remember current state for this sample, to use on backward step
                    memcpy( stateCurrent, state, mOutputsCount * sizeof( float_t ) );
                }
            } );
        }

        // Propagates error to the previous layer and calculates weights/biases gradients
        void XRecurrentLayer::BackwardCompute( const vector<fvector_t*>& inputs,
                                            const vector<fvector_t*>& /* outputs */,
                                            const vector<fvector_t*>& deltas,
                                            vector<fvector_t*>& prevDeltas,
                                            fvector_t& gradWeights,
                                            const XNetworkContext& ctx )
        {
            size_t sequenceLen   = ctx.TrainingSequenceLength( );
            size_t batchSize     = inputs.size( ) / sequenceLen;

            size_t weightsCountInputs  = mInputsCount  * mOutputsCount;
            size_t weightsCountHistory = mOutputsCount * mOutputsCount;

            // set up weights gradient pointers
            float_t* gradWeightsU = gradWeights.data( );
            float_t* gradWeightsW = gradWeightsU + weightsCountInputs;
            
            // set up biases gradient pointers
            float_t* gradBiasesB = gradWeightsW + weightsCountHistory;
            
            XParallel::For( batchSize, [&]( size_t batchIndex )
            {
                // accumulated state delta
                float_t* stateGrad = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_GRAD, batchIndex ) );

                for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                {
                    size_t  sampleIndex         = batchIndex * sequenceLen + sequenceIndex;

                    const fvector_t& delta      = *( deltas[sampleIndex] );
                    fvector_t&       prevDelta  = *( prevDeltas[sampleIndex] );

                    float_t*  stateCurrent      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_CURRENT, sampleIndex ) ); // H(t)
                    float_t*  stateDeltaCurrent = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_DELTA_CURRENT, sampleIndex ) );

                    // add state gradient from the future
                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        stateDeltaCurrent[outputIndex] = delta[outputIndex] + stateGrad[outputIndex];
                    }

                    // backward pass through Tanh activation to get final state delta for current sample
                    mTanh.BackwardActivate( stateCurrent, stateCurrent, stateDeltaCurrent, stateDeltaCurrent, mOutputsCount );

                    // input deltas for the previous layer
                    for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                    {
                        size_t  weightIndex = inputIndex;
                        float_t sum = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mInputsCount )
                        {
                            sum += stateDeltaCurrent[outputIndex] * mWeightsU[weightIndex];
                        }

                        prevDelta[inputIndex] = sum;
                    }

                    // state gradients for the previous sample in the time series
                    for ( size_t outputIndex2 = 0; outputIndex2 < mOutputsCount; outputIndex2++ )
                    {
                        size_t  weightIndex = outputIndex2;
                        float_t sum = 0;

                        for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++, weightIndex += mOutputsCount )
                        {
                            sum += stateDeltaCurrent[outputIndex] * mWeightsW[weightIndex];
                        }

                        stateGrad[outputIndex2] = sum;
                    }
                }
            } );

            XParallel::For( mOutputsCount, [&]( size_t outputIndex )
            {
                size_t weightIndexStartI = outputIndex * mInputsCount;
                size_t weightIndexStartH = outputIndex * mOutputsCount;

                for ( size_t batchIndex = 0; batchIndex < batchSize; batchIndex++ )
                {
                    for ( int sequenceIndex = (int) sequenceLen - 1; sequenceIndex >= 0; sequenceIndex-- )
                    {
                        size_t  sampleIndex         = batchIndex * sequenceLen + sequenceIndex;
                        const fvector_t& input      = *( inputs[sampleIndex] );
                        float_t*  statePrev         = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_PREV, sampleIndex ) );    // H(t-1)
                        float_t*  stateDeltaCurrnet = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_STATE_DELTA_CURRENT, sampleIndex ) );

                        // accumulate weights' gradients
                        // dU
                        for ( size_t inputIndex = 0, weightIndex = weightIndexStartI; inputIndex < mInputsCount; inputIndex++, weightIndex++ )
                        {
                            gradWeightsU[weightIndex] += stateDeltaCurrnet[outputIndex] * input[inputIndex];
                        }

                        // dW
                        if ( sequenceIndex != 0 )
                        {
                            for ( size_t historyIndex = 0, weightIndex = weightIndexStartH; historyIndex < mOutputsCount; historyIndex++, weightIndex++ )
                            {
                                gradWeightsW[weightIndex] += stateDeltaCurrnet[outputIndex] * statePrev[historyIndex];
                            }
                        }

                        // accumulate biases' gradients
                        // dB
                        gradBiasesB[outputIndex] += stateDeltaCurrnet[outputIndex];
                    }
                }
            } );
        }

        // Applies updates to the layer's weights and biases
        void XRecurrentLayer::UpdateWeights( const fvector_t& updates )
        {
            for ( size_t i = 0, n = mAllWeights.size( ); i < n; i++ )
            {
                mAllWeights[i] += updates[i];
            }
        }

        // Saves layer's learnt parameters/weights
        bool XRecurrentLayer::SaveLearnedParams( FILE* file ) const
        {
            vector<const fvector_t*> params( { &mAllWeights } );

            return SaveLearnedParamsHelper( file, LayerID::RecurrentBasic, params );
        }

        // Loads layer's learnt parameters
        bool XRecurrentLayer::LoadLearnedParams( FILE* file )
        {
            vector<fvector_t*> params( { &mAllWeights } );

            return LoadLearnedParamsHelper( file, LayerID::RecurrentBasic, params );
        }

        class IProcessingLayer : public ILayer
        {
        public:
            IProcessingLayer( size_t inputsCount, size_t outputsCount ) :
                ILayer( inputsCount, outputsCount )
            {

            }

            // None of the processing layers have weights/biases to train
            bool Trainable( ) const override
            {
                return false;
            }

            // Calls BackwardProcess() to propagate error to the previous layer 
            void BackwardCompute( const std::vector<fvector_t*>& inputs,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                fvector_t& /* gradWeights */,
                                const XNetworkContext& ctx ) override
            {
                BackwardProcess( inputs, outputs, deltas, prevDeltas, ctx );
            }

            // Propagates error to the previous layer
            virtual void BackwardProcess( const std::vector<fvector_t*>& inputs,
                                        const std::vector<fvector_t*>& outputs,
                                        const std::vector<fvector_t*>& deltas,
                                        std::vector<fvector_t*>& prevDeltas,
                                        const XNetworkContext& ctx ) = 0;
        };

        // Implementation of average pooling - outputs are average values of corresponding inputs from a square window
        class XAveragePooling : public IProcessingLayer
        {
            size_t      mInputWidth;
            size_t      mInputHeight;
            size_t      mInputDepth;
            size_t      mOutputWidth;
            size_t      mOutputHeight;

            size_t      mPoolSizeX;
            size_t      mPoolSizeY;

            size_t      mHorizontalStep;
            size_t      mVerticalStep;

            BorderMode  mBorderMode;

            std::vector<uvector_t> mOutToInMap;
            uvector_t              mInToOutMap;

        public:

            XAveragePooling( size_t inputWidth, size_t inputHeight, size_t inputDepth, size_t poolSize = 2 )
                : XAveragePooling( inputWidth, inputHeight, inputDepth, poolSize, poolSize )
            { }

            XAveragePooling( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t poolSize, size_t stepSize )
                : XAveragePooling( inputWidth, inputHeight, inputDepth, poolSize, poolSize, stepSize, stepSize )
            { }

            XAveragePooling( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                            size_t poolSizeX, size_t poolSizeY,
                            size_t horizontalStep, size_t verticalStep,
                            BorderMode borderMode = BorderMode::Valid ) :
                IProcessingLayer( 0, 0 ),
                mInputWidth( inputWidth ), mInputHeight( inputHeight ), mInputDepth( inputDepth ),
                mOutputWidth( 0 ), mOutputHeight( 0 ),
                mPoolSizeX( poolSizeX ), mPoolSizeY( poolSizeY ),
                mHorizontalStep( horizontalStep ), mVerticalStep( verticalStep ),
                mBorderMode( borderMode )
            {
                size_t padWidth    = 0;
                size_t padHeight   = 0;

                if ( mBorderMode == BorderMode::Same )
                {
                    padWidth     = poolSizeX - 1;
                    padHeight    = poolSizeY - 1;
                }

                // calculation of output width/height as:
                //   outSize = ( inSize - kernelSize + padSize ) / step + 1
                mOutputWidth  = ( mInputWidth  - mPoolSizeX  + padWidth  ) / mHorizontalStep + 1;
                mOutputHeight = ( mInputHeight - mPoolSizeY  + padHeight ) / mVerticalStep   + 1;

                // total input/output size
                Initialize( mInputWidth  * mInputHeight  * mInputDepth,
                            mOutputWidth * mOutputHeight * mInputDepth );

                // build two maps:
                //   1) first tells output index for the specified input index.
                //   2) second tells indexes of inputs for the specified output;
                // An output will always have at least one input connected to it.
                // However, some inputs may not be connected at all to any of the outputs
                // (if step size is greater than pooling size).
                mInToOutMap = XDataEncodingTools::BuildPoolingInToOutMap( inputWidth, inputHeight, inputDepth, poolSizeX, poolSizeY,
                                                                        horizontalStep, verticalStep, borderMode );
                mOutToInMap = XDataEncodingTools::BuildPoolingOutToInMap( inputWidth, inputHeight, inputDepth, poolSizeX, poolSizeY,
                                                                        horizontalStep, verticalStep, borderMode );
            }

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    fvector_t& input  = *( inputs[i]  );
                    fvector_t& output = *( outputs[i] );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        const uvector_t& outputMap = mOutToInMap[outputIndex];
                        float_t          sum       = float_t( 0 );

                        for ( auto inputIndex : outputMap )
                        {
                            sum += input[inputIndex];
                        }

                        output[outputIndex] = sum / outputMap.size( );
                    }
                } );
            }
            
            // Propagates error to the previous layer
            void BackwardProcess( const std::vector<fvector_t*>& /* inputs  */,
                                const std::vector<fvector_t*>& /* outputs */,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( deltas.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    const fvector_t& delta     = *( deltas[i] );
                    fvector_t&       prevDelta = *( prevDeltas[i] );

                    for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                    {
                        if ( mInToOutMap[inputIndex] == ANNT_NOT_CONNECTED )
                        {
                            prevDelta[inputIndex] = float_t( 0 );
                        }
                        else
                        {
                            size_t outputIndex = mInToOutMap[inputIndex];

                            prevDelta[inputIndex] = delta[outputIndex] / mOutToInMap[outputIndex].size( );
                        }
                    }
                } );
            }
        };

        // Implementation of the layer performaing normalization over batch data
        class XBatchNormalization : public IProcessingLayer
        {
        private:
            bool      mFirstUpdate;

            size_t    mSpatialSize;
            size_t    mInputDepth;

            float_t   mMomentum;
            float_t   mEpsilon;

            // learn't mean and std.dev.
            fvector_t mMean;
            fvector_t mStdDev;

            enum
            {
                BUFFER_INDEX_LEARNT_VARIANCE        = 0, // must stay alive, so buffer can not be reused in backward pass
                BUFFER_INDEX_BATCH_MEAN             = 1,
                BUFFER_INDEX_BATCH_VARIANCE         = 2,
                BUFFER_INDEX_BATCH_STD_DEV          = 3,
                BUFFER_INDEX_DELTAS_DOT_OUTPUT_MEAN = 1,
                BUFFER_INDEX_DELTAS_MEAN            = 2,
            };

        public:
            XBatchNormalization( size_t inputWidth, size_t inputHeight, size_t inputDepth, float_t momentum = float_t( 0.999 ) ) :
                IProcessingLayer( inputWidth * inputHeight * inputDepth, inputWidth * inputHeight * inputDepth ),
                mFirstUpdate( true ),
                mSpatialSize( inputWidth * inputHeight ), mInputDepth( inputDepth ),
                mMomentum( momentum ), mEpsilon( float_t( 0.00001 ) )
            {
                mMean   = fvector_t( mInputDepth, float_t( 0.0f ) );
                mStdDev = fvector_t( mInputDepth, float_t( 1.0f ) );
            }

            // Tells that we may need some extra memory for keeping temporary calculations
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize( 4 );

                // these buffers we don't need per sample - only per batch, so there is a bit of memory waste

                // forward pass - learning variance
                workingMemSize[BUFFER_INDEX_LEARNT_VARIANCE] = mInputDepth *  sizeof( float_t );

                // forward pass - batch mean while training
                // backward pass - mean of deltas[i]*output[i]
                workingMemSize[BUFFER_INDEX_BATCH_MEAN] = mInputDepth *  sizeof( float_t );

                // forward pass - batch variance while training
                // backward pass - mean of deltas
                workingMemSize[BUFFER_INDEX_BATCH_VARIANCE] = mInputDepth *  sizeof( float_t );

                // forward pass - batch std.dev while training
                workingMemSize[BUFFER_INDEX_BATCH_STD_DEV] = mInputDepth *  sizeof( float_t );

                return workingMemSize;
            }

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override
            {
                float_t* batchMean      = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_BATCH_MEAN, 0 ) );
                float_t* batchVariance  = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_BATCH_VARIANCE, 0 ) );
                float_t* batchStdDEv    = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_BATCH_STD_DEV, 0 ) );

                float_t* learntMean     = mMean.data( );
                float_t* learntVariance = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_LEARNT_VARIANCE, 0 ) );
                float_t* learntStdDev   = mStdDev.data( );

                float_t* meanToUse      = ( ctx.IsTraining( ) ) ? batchMean   : learntMean;
                float_t* stdDevToUse    = ( ctx.IsTraining( ) ) ? batchStdDEv : learntStdDev;

                if ( ctx.IsTraining( ) )
                {
                    CalculateMeanAndVariance( inputs, batchMean, batchVariance );
                    CalculateStdDev( batchVariance, batchStdDEv );
                }

                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    for ( size_t depthIndex = 0; depthIndex < mInputDepth; depthIndex++ )
                    {
                        float_t meanValue   = meanToUse[depthIndex];
                        float_t stdDevValue = stdDevToUse[depthIndex];

                        const float_t* input  = inputs[i]->data( )  + depthIndex * mSpatialSize;
                        float_t*       output = outputs[i]->data( ) + depthIndex * mSpatialSize;

                        for ( size_t spatialIndex = 0; spatialIndex < mSpatialSize; spatialIndex++, input++, output++ )
                        {
                            *output = ( *input - meanValue ) / stdDevValue;
                        }
                    }
                } );

                if ( ctx.IsTraining( ) )
                {
                    if ( mFirstUpdate )
                    {
                        memcpy( learntMean, batchMean, mInputDepth * sizeof( float_t ) );
                        memcpy( learntVariance, batchVariance, mInputDepth * sizeof( float_t ) );

                        mFirstUpdate = false;
                    }
                    else
                    {
                        float_t antiMomentum = float_t( 1 ) - mMomentum;

                        // update learnt mean and variance
                        for ( size_t depthIndex = 0; depthIndex < mInputDepth; depthIndex++ )
                        {
                            learntMean[depthIndex]     = mMomentum * learntMean[depthIndex]     + antiMomentum * batchMean[depthIndex];
                            learntVariance[depthIndex] = mMomentum * learntVariance[depthIndex] + antiMomentum * batchVariance[depthIndex];
                        }
                    }

                    // calculate std dev on learnt variance
                    CalculateStdDev( learntVariance, learntStdDev );
                }
            }

            // Propagates error to the previous layer
            void BackwardProcess( const std::vector<fvector_t*>& /* inputs  */,
                                const std::vector<fvector_t*>& outputs,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                const XNetworkContext& ctx ) override
            {
                float_t* deltasDotOutputsMean = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_DELTAS_DOT_OUTPUT_MEAN, 0 ) );
                float_t* deltasMean           = static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_DELTAS_MEAN, 0 ) );

                const float_t* stdDevToUse = ( ctx.IsTraining( ) ) ? static_cast<float_t*>( ctx.GetWorkingBuffer( BUFFER_INDEX_BATCH_STD_DEV, 0 ) ) : mStdDev.data( );

                // calculate mean for delta[i]*output[i]
                XParallel::For( mInputDepth, mInputDepth != 1, [&]( size_t depthIndex )
                {
                    deltasDotOutputsMean[depthIndex] = float_t( 0 );

                    for ( size_t i = 0, n = outputs.size( ); i < n; i++ )
                    {
                        const float_t* output = outputs[i]->data( ) + depthIndex * mSpatialSize;
                        const float_t* delta  = deltas[i]->data( )  + depthIndex * mSpatialSize;

                        deltasDotOutputsMean[depthIndex] += XVectorize::Dot( delta, output, mSpatialSize ) / mSpatialSize;
                    }

                    deltasDotOutputsMean[depthIndex] /= outputs.size( );
                } );

                CalculateMean( deltas, deltasMean );

                XParallel::For( outputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    const fvector_t& output    = *( outputs[i] );
                    const fvector_t& delta     = *( deltas[i] );
                    fvector_t&       prevDelta = *( prevDeltas[i] );

                    for ( size_t depthIndex = 0; depthIndex < mInputDepth; depthIndex++ )
                    {
                        for ( size_t spatialIndex = 0, j = depthIndex * mSpatialSize; spatialIndex < mSpatialSize; spatialIndex++, j++ )
                        {
                            prevDelta[j] = ( delta[j] - deltasMean[depthIndex] - deltasDotOutputsMean[depthIndex] * output[j] ) / stdDevToUse[depthIndex];
                        }
                    }
                } );
            }

            // Saves layer's learnt parameters/weights
            bool SaveLearnedParams( FILE* file ) const override
            {
                std::vector<const fvector_t*> params( { &mMean, &mStdDev } );

                return SaveLearnedParamsHelper( file, LayerID::BatchNormalization, params );
            }

            // Loads layer's learnt parameters
            bool LoadLearnedParams( FILE* file ) override
            {
                std::vector<fvector_t*> params( { &mMean, &mStdDev } );

                return LoadLearnedParamsHelper( file, LayerID::BatchNormalization, params );
            }

        private:

            void CalculateMean( const std::vector<fvector_t*>& inputs, float_t* mean ) const
            {
                //for ( size_t depthIndex = 0; depthIndex < mInputDepth; depthIndex++ )
                XParallel::For( mInputDepth, mInputDepth != 1, [&]( size_t depthIndex )
                {
                    mean[depthIndex] = float_t( 0 );

                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        mean[depthIndex] += std::accumulate( inputs[i]->begin( ) +   depthIndex       * mSpatialSize,
                                                            inputs[i]->begin( ) + ( depthIndex + 1 ) * mSpatialSize,
                                                            float_t( 0 ) ) / mSpatialSize;
                    }

                    mean[depthIndex] /= inputs.size( );
                } );
            }

            void CalculateMeanAndVariance( const std::vector<fvector_t*>& inputs, float_t* mean, float_t* variance ) const
            {
                CalculateMean( inputs, mean );

                XParallel::For( mInputDepth, mInputDepth != 1, [&]( size_t depthIndex )
                {
                    float_t meanValue = mean[depthIndex];
                    float_t diff;

                    variance[depthIndex] = float_t( 0 );

                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        const float_t* input = inputs[i]->data( ) + depthIndex * mSpatialSize;
                        float_t        sum   = float_t( 0 );

                        for ( size_t spatialIndex = 0; spatialIndex < mSpatialSize; spatialIndex++, input++ )
                        {
                            diff = *input - meanValue;
                            sum += diff * diff;
                        }

                        variance[depthIndex] += sum / mSpatialSize;
                    }

                    variance[depthIndex] /= inputs.size( );
                } );
            }

            void CalculateStdDev( const float_t* variance, float_t* stdDev )
            {
                for ( size_t depthIndex = 0; depthIndex < mInputDepth; depthIndex++ )
                {
                    stdDev[depthIndex] = std::sqrt( variance[depthIndex] + mEpsilon );
                }
            }
        };

        class XDropOutLayer : public IProcessingLayer
        {
        private:
            float_t mDropOutRate;

            std::mt19937                            mGenerator;
            std::uniform_real_distribution<float_t> mDistribution;

        public:
            XDropOutLayer( float_t dropOutRate = float_t( 0.1f ) ) :
                IProcessingLayer( 0, 0 ),
                mDropOutRate( dropOutRate ),
                mGenerator( ), mDistribution( 0.0f, 1.0f )
            {

            }

            // Tells that we may need some extra memory for keeping drop out mask (in training mode)
            uvector_t WorkingMemSize( bool trainingMode ) const override
            {
                uvector_t workingMemSize;;

                if ( trainingMode )
                {
                    workingMemSize.push_back( mOutputsCount * sizeof( float_t ) );
                }

                return workingMemSize;
            }

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    fvector_t& input  = *( inputs[i] );
                    fvector_t& output = *( outputs[i] );

                    if ( !ctx.IsTraining( ) )
                    {
                        output = input;
                    }
                    else
                    {
                        float_t* dropOutMask = static_cast<float_t*>( ctx.GetWorkingBuffer( 0, i ) );

                        for ( size_t j = 0; j < mOutputsCount; j++ )
                        {
                            dropOutMask[j] = ( mDistribution( mGenerator ) < mDropOutRate ) ? float_t( 0.0f ) : float_t( 1.0f );
                        }

                        for ( size_t j = 0; j < mOutputsCount; j++ )
                        {
                            output[j] = input[j] * dropOutMask[j];
                        }
                    }
                } );
            }

            // Propagates error to the previous layer
            void BackwardProcess( const std::vector<fvector_t*>& /* inputs  */,
                                const std::vector<fvector_t*>& /* outputs */,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( deltas.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    const fvector_t& delta     = *( deltas[i] );
                    fvector_t&       prevDelta = *( prevDeltas[i] );

                    if ( !ctx.IsTraining( ) )
                    {
                        prevDelta = delta;
                    }
                    else
                    {
                        float_t* dropOutMask = static_cast<float_t*>( ctx.GetWorkingBuffer( 0, i ) );

                        for ( size_t j = 0; j < mOutputsCount; j++ )
                        {
                            prevDelta[j] = delta[j] * dropOutMask[j];
                        }
                    }
                } );
            }
        };


        // Implementation of maximum pooling - outputs are maximum values of corresponding inputs from a square window
        class XMaxPooling : public IProcessingLayer
        {
            size_t      mInputWidth;
            size_t      mInputHeight;
            size_t      mInputDepth;
            size_t      mOutputWidth;
            size_t      mOutputHeight;

            size_t      mPoolSizeX;
            size_t      mPoolSizeY;

            size_t      mHorizontalStep;
            size_t      mVerticalStep;

            BorderMode  mBorderMode;

            std::vector<uvector_t> mOutToInMap;
            uvector_t              mInToOutMap;

        public:

            XMaxPooling( size_t inputWidth, size_t inputHeight, size_t inputDepth, size_t poolSize = 2 )
                : XMaxPooling( inputWidth, inputHeight, inputDepth, poolSize, poolSize )
            { }

            XMaxPooling( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                        size_t poolSize, size_t stepSize )
                : XMaxPooling( inputWidth, inputHeight, inputDepth, poolSize, poolSize, stepSize, stepSize )
            { }

            XMaxPooling( size_t inputWidth, size_t inputHeight, size_t inputDepth,
                        size_t poolSizeX, size_t poolSizeY,
                        size_t horizontalStep, size_t verticalStep,
                        BorderMode borderMode = BorderMode::Valid ) :
                IProcessingLayer( 0, 0 ),
                mInputWidth( inputWidth ), mInputHeight( inputHeight ), mInputDepth( inputDepth ),
                mOutputWidth( 0 ), mOutputHeight( 0 ),
                mPoolSizeX( poolSizeX ), mPoolSizeY( poolSizeY ),
                mHorizontalStep( horizontalStep ), mVerticalStep( verticalStep ),
                mBorderMode( borderMode )
            {
                size_t padWidth    = 0;
                size_t padHeight   = 0;

                if ( mBorderMode == BorderMode::Same )
                {
                    padWidth     = poolSizeX - 1;
                    padHeight    = poolSizeY - 1;
                }

                // calculation of output width/height as:
                //   outSize = ( inSize - kernelSize + padSize ) / step + 1
                mOutputWidth  = ( mInputWidth  - mPoolSizeX  + padWidth  ) / mHorizontalStep + 1;
                mOutputHeight = ( mInputHeight - mPoolSizeY  + padHeight ) / mVerticalStep   + 1;

                // total input/output size
                Initialize( mInputWidth  * mInputHeight  * mInputDepth,
                            mOutputWidth * mOutputHeight * mInputDepth );

                // build two maps:
                //   1) first tells output index for the specified input index.
                //   2) second tells indexes of inputs for the specified output;
                // An output will always have at least one input connected to it.
                // However, some inputs may not be connected at all to any of the outputs
                // (if step size is greater than pooling size).
                mInToOutMap = XDataEncodingTools::BuildPoolingInToOutMap( inputWidth, inputHeight, inputDepth, poolSizeX, poolSizeY,
                                                                        horizontalStep, verticalStep, borderMode );
                mOutToInMap = XDataEncodingTools::BuildPoolingOutToInMap( inputWidth, inputHeight, inputDepth, poolSizeX, poolSizeY,
                                                                        horizontalStep, verticalStep, borderMode );
            }

            // Tells that we may need some extra memory for keeping indexes of maximum elements (in training mode)
            uvector_t WorkingMemSize( bool /* trainingMode */ ) const override
            {
                uvector_t workingMemSize = uvector_t( 1 );

                // we don't really need this memory when doing inference only,
                // but don't want to check that always when doing forward pass
                workingMemSize[0] = mOutputsCount * sizeof( size_t);

                return workingMemSize;
            }

            // Calculates outputs for the given inputs
            void ForwardCompute( const std::vector<fvector_t*>& inputs,
                                std::vector<fvector_t*>& outputs,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( inputs.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    fvector_t& input      = *( inputs[i] );
                    fvector_t& output     = *( outputs[i] );
                    size_t*    maxIndexes = static_cast<size_t*>( ctx.GetWorkingBuffer( 0, i ) );

                    for ( size_t outputIndex = 0; outputIndex < mOutputsCount; outputIndex++ )
                    {
                        const std::vector<size_t>& outputMap = mOutToInMap[outputIndex];
                        float_t                    maxValue  = std::numeric_limits<float_t>::lowest( );
                        size_t                     maxIndex  = 0;

                        for ( auto inputIndex : outputMap )
                        {
                            if ( input[inputIndex] > maxValue )
                            {
                                maxValue = input[inputIndex];
                                maxIndex = inputIndex;
                            }
                        }

                        output[outputIndex]     = maxValue;
                        maxIndexes[outputIndex] = maxIndex;
                    }
                } );
            }
            
            // Propagates error to the previous layer
            void BackwardProcess( const std::vector<fvector_t*>& /* inputs  */,
                                const std::vector<fvector_t*>& /* outputs */,
                                const std::vector<fvector_t*>& deltas,
                                std::vector<fvector_t*>& prevDeltas,
                                const XNetworkContext& ctx ) override
            {
                XParallel::For( deltas.size( ), ctx.IsTraining( ), [&]( size_t i )
                {
                    const fvector_t& delta      = *( deltas[i] );
                    fvector_t&       prevDelta  = *( prevDeltas[i] );
                    size_t*          maxIndexes = static_cast<size_t*>( ctx.GetWorkingBuffer( 0, i ) );

                    for ( size_t inputIndex = 0; inputIndex < mInputsCount; inputIndex++ )
                    {
                        if ( mInToOutMap[inputIndex] == ANNT_NOT_CONNECTED )
                        {
                            prevDelta[inputIndex] = float_t( 0 );
                        }
                        else
                        {
                            size_t outputIndex = mInToOutMap[inputIndex];

                            prevDelta[inputIndex] = ( maxIndexes[outputIndex] == inputIndex ) ? delta[outputIndex] : float_t( 0 );
                        }
                    }
                } );
            }
        };

        // Cost functions' interface for calculating cost and gradient
        // for a specified output/target pair
        class ICostFunction
        {
        public:
            virtual ~ICostFunction( ) { }

            // Calculates cost value of the specified output vector
            virtual float_t Cost( const fvector_t& output, const fvector_t& target ) const = 0;

            // Calculates gradient for the specified output/target pair
            virtual fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const = 0;
        };

        // Implementation of Mean Absolute Error (Manhattan) cost function
        //
        // Calculated as sum of absolute differences divided by N
        //
        class XAbsoluteCost : public ICostFunction
        {
        public:
            // Calculates cost value of the specified output vector
            float_t Cost( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t  length = output.size( );
                float_t cost   = float_t( 0 );

                for ( size_t i = 0; i < length; i++ )
                {
                    cost += std::abs( output[i] - target[i] );
                }

                cost /= length;

                return cost;
            }

            // Calculates gradient for the specified output/target pair
            fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t    length = output.size( );
                fvector_t grad( length );
                float_t   diff;

                for ( size_t i = 0; i < length; i++ )
                {
                    diff = output[i] - target[i];

                    if ( diff > float_t( 0 ) )
                    {
                        grad[i] = float_t( 1 );
                    }
                    else if ( diff < float_t( 0 ) )
                    {
                        grad[i] = float_t( -1 );
                    }
                    else
                    {
                        grad[i] = float_t( 0 );
                    }
                }

                return grad;
            }
        };

        // Implementation of binary cross entropy cost function
        class XBinaryCrossEntropyCost : public ICostFunction
        {
        public:
            // Calculates cost value of the specified output vector
            float_t Cost( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t  length = output.size( );
                float_t cost   = float_t( 0 );

                for ( size_t i = 0; i < length; i++ )
                {
                    cost += -( target[i] * std::log( output[i] ) +
                            ( float_t( 1 ) - target[i] ) * std::log( float_t( 1 ) - output[i] ) );
                }

                return cost;
            }

            // Calculates gradient for the specified output/target pair
            fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t    length = output.size( );
                fvector_t grad( length );

                for ( size_t i = 0; i < length; i++ )
                {
                    grad[i] = ( output[i] - target[i] ) / ( output[i] * ( float_t( 1 ) - output[i] ) );
                }

                return grad;
            }
        };


        // Implementation of cross entropy cost function
        class XCrossEntropyCost : public ICostFunction
        {
        public:
            // Calculates cost value of the specified output vector
            float_t Cost( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t  length = output.size( );
                float_t cost   = float_t( 0 );

                for ( size_t i = 0; i < length; i++ )
                {
                    cost += -target[i] * std::log( output[i] );
                }

                return cost;
            }

            // Calculates gradient for the specified output/target pair
            fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t    length = output.size( );
                fvector_t grad( length );

                for ( size_t i = 0; i < length; i++ )
                {
                    grad[i] = -target[i] / output[i];
                }

                return grad;
            }
        };

        // Implementation of Mean Square Error cost function
        //
        // Calculated as sum of squared differences divided by 2N
        //
        class XMSECost : public ICostFunction
        {
        public:

            // Calculates cost value of the specified output vector
            float_t Cost( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t  length = output.size( );
                float_t cost   = float_t( 0 );
                float_t diff;

                for ( size_t i = 0; i < length; i++ )
                {
                    diff  = output[i] - target[i];
                    cost += diff * diff;
                }

                cost /= ( float( 2 ) * length );

                return cost;
            }

            // Calculates gradient for the specified output/target pair
            fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t    length = output.size( );
                fvector_t grad( length );

                for ( size_t i = 0; i < length; i++ )
                {
                    grad[i] = output[i] - target[i];
                }

                return grad;
            }
        };

        // Implementation of negative log-likelihood cost function - to be used after XLogSoftMaxActivation layer,
        // which produces log-probabilities
        //
        class XNegativeLogLikelihoodCost : public ICostFunction
        {
        public:

            // Calculates cost value of the specified output vector
            float_t Cost( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t  length = output.size( );
                float_t cost   = float_t( 0 );

                for ( size_t i = 0; i < length; i++ )
                {
                    cost += -target[i] * output[i];
                }

                return cost;
            }

            // Calculates gradient for the specified output/target pair
            fvector_t Gradient( const fvector_t& output, const fvector_t& target ) const override
            {
                size_t    length = output.size( );
                fvector_t grad( length );

                for ( size_t i = 0; i < length; i++ )
                {
                    grad[i] = -target[i];
                }

                return grad;
            }
        };

        // Common interface for algorithms calculating updates for weights/biases from their gradients
        class INetworkOptimizer
        {
        protected:
            float mLearningRate;

        public:
            INetworkOptimizer( float_t learningRate ) :
                mLearningRate( learningRate )
            {
            }

            virtual ~INetworkOptimizer( ) { }

            // Learning rate to control amount of weights/biases update
            float_t LearningRate( ) const
            {
                return mLearningRate;
            }
            virtual void SetLearningRate( float_t learningRate )
            {
                mLearningRate = learningRate;
            }

            // Variables count per learning parameter (weight/bias). This can be value of the previous update
            // like in Momentum optimizer, for example, etc.
            virtual size_t ParameterVariablesCount( ) const
            {
                return 0;
            }

            // Variables count per layer. This variables replace mutable class members, so that optimization algorithms don't
            // store any state inside. For example, Adam optimizer uses them to store b1^t and b2^t values.
            virtual size_t LayerVariablesCount( ) const
            {
                return 0;
            }

            // Calculates weights/biases updates from given gradients
            virtual void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& layerVariables ) = 0;
        };

        // Implementation of Adagrad optimization algorithm
        // http://ruder.io/optimizing-gradient-descent/index.html#adagrad
        //
        class XAdagradOptimizer : public INetworkOptimizer
        {
        private:
            float_t mEpsilon;

        public:
            XAdagradOptimizer( float_t learningRate = float_t( 0.01 ) ) :
                INetworkOptimizer( learningRate ), mEpsilon( float_t( 1e-8 ) )
            {
            }

            size_t ParameterVariablesCount( ) const override
            {
                return 1;
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& /* layerVariables */ ) override
            {
                fvector_t& sqUpdatesSum = paramVariables[0];

                for ( size_t i = 0, n = updates.size( ); i < n; i++ )
                {
                    sqUpdatesSum[i] += updates[i] * updates[i];
                    updates[i]      *= -mLearningRate / std::sqrt( sqUpdatesSum[i] + mEpsilon );
                }
            }
        };

        class XAdamOptimizer : public INetworkOptimizer
        {
        private:
            float_t mEpsilon;
            float_t mB1;
            float_t mB2;

        public:
            XAdamOptimizer( float_t learningRate = float_t( 0.001 ) ) :
                INetworkOptimizer( learningRate ), mEpsilon( float_t( 1e-8 ) ),
                mB1( float_t( 0.9 ) ), mB2( float_t( 0.999 ) )
            {
            }

            // Two variables per learning parameter, m(t) and v(t)
            size_t ParameterVariablesCount( ) const override
            {
                return 2;
            }

            // Variables to keep b1^t and b2^t values
            virtual size_t LayerVariablesCount( ) const
            {
                return 3;
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& layerVariables ) override
            {
                fvector_t& mt  = paramVariables[0];
                fvector_t& vt  = paramVariables[1];
                float_t    b1t = mB1;
                float_t    b2t = mB2;

                // check if it is the first call
                if ( layerVariables[0] < float( 0.5 ) )
                {
                    layerVariables[0] = float( 1.0 );
                }
                else
                {
                    b1t = layerVariables[1];
                    b2t = layerVariables[2];
                }

                for ( size_t i = 0, n = updates.size( ); i < n; i++ )
                {
                    mt[i] = mB1 * mt[i] + ( float_t( 1 ) - mB1 ) * updates[i];
                    vt[i] = mB2 * vt[i] + ( float_t( 1 ) - mB2 ) * updates[i] * updates[i];

                    updates[i] = -mLearningRate * ( mt[i] / ( float_t( 1 ) - b1t ) ) /
                                std::sqrt( vt[i] / ( float_t( 1 ) - b2t ) + mEpsilon );
                }

                b1t *= mB1;
                b2t *= mB2;

                layerVariables[1] = b1t;
                layerVariables[2] = b2t;
            }
        };

        // Implementation of classical Stochastic Gradient Descent optimizer, which is
        //   paramUpdate(t) = -learningRate * paramGrad(t)
        // http://cs231n.github.io/neural-networks-3/#sgd
        //
        class XGradientDescentOptimizer : public INetworkOptimizer
        {
        public:
            XGradientDescentOptimizer( float_t learningRate = float_t( 0.01 ) ) :
                INetworkOptimizer( learningRate )
            {
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& /* paramVariables */, fvector_t& /* layerVariables */ ) override
            {
                for ( auto& update : updates )
                {
                    update *= -mLearningRate;
                }
            }
        };

        // Implementation of SGD with Momentum, which calculates updates as
        //   v(t) = momentum * v(t-1) + learningRate * paramGrad
        //   paramUpdate(t) = -v(t)
        // http://ruder.io/optimizing-gradient-descent/index.html#momentum
        //
        class XMomentumOptimizer : public INetworkOptimizer
        {
        private:
            float_t mMomentum;

        public:
            XMomentumOptimizer( float_t learningRate = float_t( 0.01 ), float_t momentum = float_t( 0.9 ) ) :
                INetworkOptimizer( learningRate ), mMomentum( momentum )
            {
            }

            // One variable per learning parameter - previous update value
            size_t ParameterVariablesCount( ) const override
            {
                return 1;
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& /* layerVariables */ ) override
            {
                fvector_t& vPrev = paramVariables[0];

                for ( size_t i = 0, n = updates.size( ); i < n; i++ )
                {
                    float_t vt = mMomentum * vPrev[i] + mLearningRate * updates[i];

                    updates[i] = -vt;
                    vPrev[i]   = vt;
                }
            }
        };

        // Implementation of SGD with Nesterov Momentum, which calculates updates as
        // http://cs231n.github.io/neural-networks-3/#sgd
        //
        class XNesterovMomentumOptimizer : public INetworkOptimizer
        {
        private:
            float_t mMomentum;

        public:
            XNesterovMomentumOptimizer( float_t learningRate = float_t( 0.01 ), float_t momentum = ( 0.9 ) ) :
                INetworkOptimizer( learningRate ), mMomentum( momentum )
            {
            }

            size_t ParameterVariablesCount( ) const override
            {
                return 1;
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& /* layerVariables */ ) override
            {
                fvector_t& vPrev = paramVariables[0];

                for ( size_t i = 0, n = updates.size( ); i < n; i++ )
                {
                    float_t vt = mMomentum * vPrev[i] - mLearningRate * updates[i];

                    updates[i] = -mMomentum * vPrev[i] + ( float_t( 1.0 ) + mMomentum ) * vt;
                    vPrev[i]   = vt;
                }
            }
        };



        // Implementation of RMSprop optimization algorithm
        // http://ruder.io/optimizing-gradient-descent/index.html#rmsprop
        //
        class XRMSpropOptimizer : public INetworkOptimizer
        {
        private:
            float_t mEpsilon;
            float_t mMu;

        public:
            XRMSpropOptimizer( float_t learningRate = float_t( 0.001 ), float_t mu = float_t( 0.9 ) ) :
                INetworkOptimizer( learningRate ),
                mEpsilon( float_t( 1e-8 ) ), mMu( mu )
            {
            }

            size_t ParameterVariablesCount( ) const override
            {
                return 1;
            }

            void CalculateUpdatesFromGradients( fvector_t& updates, std::vector<fvector_t>& paramVariables, fvector_t& /* layerVariables */ ) override
            {
                fvector_t& eg = paramVariables[0];

                for ( size_t i = 0, n = updates.size( ); i < n; i++ )
                {
                    eg[i]       = mMu * eg[i] + ( 1 - mMu ) * updates[i] * updates[i];
                    updates[i] *= -mLearningRate / std::sqrt( eg[i] + mEpsilon );
                }
            }
        };

        namespace Training {

            // Implementation of artificial neural network training with
            // error back propagation algorithm - wraps memory buffers and
            // infrastructure required for training a network: run forward
            // pass, calculate error, propagate it backward through the
            // network calculating gradients of weights/biases and applying
            // them.
            //
            class XNetworkTraining : public XNetworkInference
            {
            private:
                std::shared_ptr<INetworkOptimizer >  mOptimizer;
                std::shared_ptr<ICostFunction>       mCostFunction;
                bool                                 mAverageWeightGradients;

            private:
                // storage and pointers for outputs computed during training
                std::vector<std::vector<fvector_t>>  mTrainOutputsStorage;
                std::vector<std::vector<fvector_t*>> mTrainOutputs;

                // storade and pointers to compute deltas for each layer
                std::vector<std::vector<fvector_t>>  mDeltasStorage;
                std::vector<std::vector<fvector_t*>> mDeltas;

                // storage and pointers "input deltas"
                // no needed, just to make calculations consistent
                std::vector<fvector_t>               mInputDeltasStorage;
                std::vector<fvector_t*>              mInputDeltas;

                // vectors used to assemble pointers to training samples (input/outputs)
                std::vector<fvector_t*>              mTrainInputs;
                std::vector<fvector_t*>              mTargetOuputs;

                // weights/biases gradients for all layers
                std::vector<fvector_t>               mGradWeights;

                // vectors with parameter variables for optimizer
                std::vector<std::vector<fvector_t>>  mOptimizerParameterVariables;

                // vectors with layer variables for optimizer
                std::vector<fvector_t>               mOptimizerLayerVariables;

                // layers' working buffers and context for training
                XNetworkContext                      mTrainingContext;

            public:

                XNetworkTraining( const std::shared_ptr<XNeuralNetwork>& network,
                                const std::shared_ptr<INetworkOptimizer>& optimizer,
                                const std::shared_ptr<ICostFunction>& costFunction );

                // Provides access to the ANN
                std::shared_ptr<XNeuralNetwork> Network( ) const
                {
                    return mNetwork;
                }

                // Provides access to the weights/biases optimizer
                std::shared_ptr<INetworkOptimizer> Optimizer( ) const
                {
                    return mOptimizer;
                }

                // Provides access to the cost function used for error calculation
                std::shared_ptr<ICostFunction> CostFunction( ) const
                {
                    return mCostFunction;
                }

                // Average or not weights/biases gradients when running in batch mode
                bool AverageWeightGradients( ) const
                {
                    return mAverageWeightGradients;
                }
                void SetAverageWeightGradients( bool average )
                {
                    mAverageWeightGradients = average;
                }

                // Get/set length of training sequences used for recurrent networks
                size_t TrainingSequenceLength( ) const
                {
                    return mTrainingContext.TrainingSequenceLength( );
                }
                void SetTrainingSequenceLength( size_t sequenceLength )
                {
                    mTrainingContext.SetTrainingSequenceLength( sequenceLength );
                }

                // Reset working buffers for all layers
                void ResetState( ) override
                {
                    XNetworkInference::ResetState( );
                    mTrainingContext.ResetWorkingBuffers( );
                }
                // Reset working buffers for the specified layers
                void ResetLayersState( uvector_t layersIndexes ) override
                {
                    XNetworkInference::ResetLayersState( layersIndexes );
                    mTrainingContext.ResetWorkingBuffers( layersIndexes );
                }

                // Trains single input/output sample
                float_t TrainSample( const fvector_t& input, const fvector_t& targetOutput );

                // Trains single batch of prepared samples (as vectors)
                float_t TrainBatch( const std::vector<fvector_t>& inputs,
                                    const std::vector<fvector_t>& targetOutputs );

                // Trains single batch of prepared samples (as pointers to vectors)
                float_t TrainBatch( const std::vector<fvector_t*>& inputs,
                                    const std::vector<fvector_t*>& targetOutputs );

                // Trains single epoch using batches of the specified size (samples are provided as vectors)
                float_t TrainEpoch( const std::vector<fvector_t>& inputs, 
                                    const std::vector<fvector_t>& targetOutputs,
                                    size_t batchSize,
                                    bool randomPickIntoBatch = false );

                // Trains single epoch using batches of the specified size (samples are provided as pointers to vectors)
                float_t TrainEpoch( const std::vector<fvector_t*>& inputs, 
                                    const std::vector<fvector_t*>& targetOutputs,
                                    size_t batchSize,
                                    bool randomPickIntoBatch = false );

                // Tests sample - calculates real output and provides error cost
                float_t TestSample( const fvector_t& input,
                                    const fvector_t& targetOutput,
                                    fvector_t& output );

                // Tests classification for the provided inputs and target labels -
                // provides number of correctly classified samples and average cost (target outputs are used)
                size_t TestClassification( const std::vector<fvector_t>& inputs,
                                        const uvector_t& targetLabels,
                                        const std::vector<fvector_t>& targetOutputs,
                                        float_t* pAvgCost );
                size_t TestClassification( const std::vector<fvector_t*>& inputs,
                                        const uvector_t& targetLabels,
                                        const std::vector<fvector_t*>& targetOutputs,
                                        float_t* pAvgCost );

            private:

                float_t RunTraining( );
                float_t CalculateError( );
                void    DoBackwardCompute( );
                void    UpdateWeights( );
                void    AllocateTrainVectors( size_t samplesCount );
            };

            XNetworkTraining::XNetworkTraining( const shared_ptr<XNeuralNetwork>& network,
                                                const shared_ptr<INetworkOptimizer>& optimizer,
                                                const shared_ptr<ICostFunction>& costFunction ) :
                XNetworkInference( network ),
                mOptimizer( optimizer ),
                mCostFunction( costFunction ),
                mAverageWeightGradients( true ),
                mTrainingContext( true, 1 )
            {
                size_t optimizerParameterVariablesCount = mOptimizer->ParameterVariablesCount( );
                size_t optimizerLayerVariablesCount     = mOptimizer->LayerVariablesCount( );

                // allocate everything, which does not depend on batch size (number of input/output samples),
                // but only depends on layers count
                // 1) weight and bias gradients (accumulated over samples during batch);
                // 2) optimizer's variables;

                for ( auto layer : *mNetwork )
                {
                    size_t weightsCount = 0;

                    // allocate weight and bias gradients for trainable layers
                    // (for each layer, but not for each sample, since those accumulated over samples)
                    if ( layer->Trainable( ) )
                    {
                        weightsCount = static_pointer_cast<ITrainableLayer>( layer )->WeightsCount( );
                    }

                    mGradWeights.push_back( fvector_t( weightsCount ) );

                    // optimizer's variables ...
                    mOptimizerParameterVariables.push_back( vector<fvector_t>( optimizerParameterVariablesCount ) );
                    mOptimizerLayerVariables.push_back( fvector_t( optimizerLayerVariablesCount ) );

                    for ( size_t i = 0; i < optimizerParameterVariablesCount; i++ )
                    {
                        mOptimizerParameterVariables.back( )[i] = fvector_t( weightsCount );
                    }
                }
            }

            // Allocate the rest of vectors required for training - those which depend on the batch size
            void XNetworkTraining::AllocateTrainVectors( size_t samplesCount )
            {
                size_t layersCount = mNetwork->LayersCount( );

                if ( mTrainInputs.size( ) != samplesCount )
                {
                    mTrainInputs.resize( samplesCount );
                    mTargetOuputs.resize( samplesCount );

                    mTrainOutputsStorage.resize( layersCount );
                    mTrainOutputs.resize( layersCount );

                    mDeltasStorage.resize( layersCount );
                    mDeltas.resize( layersCount );

                    // prepare output vector and deltas for all samples and for all layers
                    for ( size_t layerIndex = 0; layerIndex < layersCount; layerIndex++ )
                    {
                        size_t layerOutputCount = mNetwork->LayerAt( layerIndex )->OutputsCount( );

                        mTrainOutputsStorage[layerIndex].resize( samplesCount );
                        mTrainOutputs[layerIndex].resize( samplesCount );

                        mDeltasStorage[layerIndex].resize( samplesCount );
                        mDeltas[layerIndex].resize( samplesCount );

                        for ( size_t i = 0; i < samplesCount; i++ )
                        {
                            mTrainOutputsStorage[layerIndex][i] = fvector_t( layerOutputCount );
                            mTrainOutputs[layerIndex][i]        = &( mTrainOutputsStorage[layerIndex][i] );

                            mDeltasStorage[layerIndex][i] = fvector_t( layerOutputCount );
                            mDeltas[layerIndex][i]        = &( mDeltasStorage[layerIndex][i] );
                        }
                    }

                    // to make calculations consistant, we have deltas for inputs as well ("previous" layer of the first)
                    mInputDeltasStorage.resize( samplesCount );
                    mInputDeltas.resize( samplesCount );

                    for ( size_t i = 0; i < samplesCount; i++ )
                    {
                        mInputDeltasStorage[i] = fvector_t( mNetwork->InputsCount( ) );
                        mInputDeltas[i] = &( mInputDeltasStorage[i] );
                    }

                    // allocate new buffers for layers
                    mTrainingContext.AllocateWorkingBuffers( mNetwork, samplesCount );
                }
            }

            // Calculate error of the last layer for each training sample
            float_t XNetworkTraining::CalculateError( )
            {
                vector<fvector_t>& lastOutputs = mTrainOutputsStorage.back( );
                vector<fvector_t>& lastDeltas  = mDeltasStorage.back( );
                float_t            totalCost   = 0;

                for ( size_t i = 0, n = mTrainInputs.size( ); i < n; i++ )
                {
                    fvector_t& lastDelta    = lastDeltas[i];
                    fvector_t& lastOutput   = lastOutputs[i];
                    fvector_t& targetOutput = *mTargetOuputs[i];

                    totalCost += mCostFunction->Cost( lastOutput, targetOutput );
                    lastDelta  = mCostFunction->Gradient( lastOutput, targetOutput );
                }

                totalCost /= mTrainInputs.size( );

                return totalCost;
            }

            // Propagate error through the network starting from last layer
            void XNetworkTraining::DoBackwardCompute( )
            {
                size_t  layerIndex  = mNetwork->LayersCount( ) - 1;
                
                // propagate deltas for all layers except the first one
                for ( ; layerIndex > 0; layerIndex-- )
                {
                    mTrainingContext.SetCurrentLayerIndex( layerIndex );

                    mNetwork->LayerAt( layerIndex )->
                        BackwardCompute( mTrainOutputs[layerIndex - 1], mTrainOutputs[layerIndex],
                                        mDeltas[layerIndex], mDeltas[layerIndex - 1],
                                        mGradWeights[layerIndex], mTrainingContext );
                }

                // now same for the first layer
                mTrainingContext.SetCurrentLayerIndex( 0 );

                mNetwork->LayerAt( 0 )->
                    BackwardCompute( mTrainInputs, mTrainOutputs[0],
                                    mDeltas[0], mInputDeltas,
                                    mGradWeights[0], mTrainingContext );
            }

            // Calculate weights/biases updates from gradients and apply them
            void XNetworkTraining::UpdateWeights( )
            {
                auto    itLayers          = mNetwork->begin( );
                float_t batchUpdateFactor = float_t( 1 );
                
                if ( mAverageWeightGradients )
                {
                    batchUpdateFactor /= mTrainInputs.size( );
                }

                for ( size_t i = 0, n = mNetwork->LayersCount( ); i < n; i++, ++itLayers )
                {
                    if ( (*itLayers)->Trainable( ) )
                    {
                        if ( mAverageWeightGradients )
                        {
                            std::transform( mGradWeights[i].begin( ), mGradWeights[i].end( ), mGradWeights[i].begin( ),
                                            [&]( float_t v ) -> float_t { return v * batchUpdateFactor; } );
                        }

                        mOptimizer->CalculateUpdatesFromGradients( mGradWeights[i], mOptimizerParameterVariables[i], mOptimizerLayerVariables[i] );

                        static_pointer_cast<ITrainableLayer>( *itLayers )->UpdateWeights( mGradWeights[i] );

                        // reset gradients for the next training cycle
                        fill( mGradWeights[i].begin( ), mGradWeights[i].end( ), float_t( 0 ) );
                    }
                }
            }

            // Run single training cycle
            float_t XNetworkTraining::RunTraining( )
            {
                float_t cost;

                // 1 - compute the network to get the actual output
                DoCompute( mTrainInputs, mTrainOutputs, mTrainingContext );

                // 2 - get error of the last layer
                cost = CalculateError( );

                // 3 - propagate the error backward through the network
                DoBackwardCompute( );

                // 4 - calculate weights/bias updates and apply those
                UpdateWeights( );

                return cost;
            }

            // Trains single input/output sample
            float_t XNetworkTraining::TrainSample( const fvector_t& input, const fvector_t& targetOutput )
            {
                float_t cost = 0;

                if ( mNetwork->LayersCount( ) != 0 )
                {
                    AllocateTrainVectors( 1 );

                    // get the single input/output into usable form
                    mTrainInputs[0]  = const_cast<fvector_t*>( &input );
                    mTargetOuputs[0] = const_cast<fvector_t*>( &targetOutput );

                    cost = RunTraining( );
                }

                return cost;
            }

            // Trains single batch of prepared samples (as vectors)
            float_t XNetworkTraining::TrainBatch( const vector<fvector_t>& inputs,
                                                const vector<fvector_t>& targetOutputs )
            {
                float_t cost = 0;

                if ( mNetwork->LayersCount( ) != 0 )
                {
                    AllocateTrainVectors( inputs.size( ) );

                    // prepare inputs vectors and target ouputs
                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        mTrainInputs[i]  = const_cast<fvector_t*>( &( inputs[i] ) );
                        mTargetOuputs[i] = const_cast<fvector_t*>( &( targetOutputs[i] ) );
                    }

                    cost = RunTraining( );
                }

                return cost;
            }

            // Trains single batch of prepared samples (as pointers to vectors)
            float_t XNetworkTraining::TrainBatch( const vector<fvector_t*>& inputs,
                                                const vector<fvector_t*>& targetOutputs )
            {
                float_t cost = 0;

                if ( mNetwork->LayersCount( ) != 0 )
                {
                    AllocateTrainVectors( inputs.size( ) );

                    // prepare inputs vectors and target ouputs
                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        mTrainInputs[i]  = inputs[i];
                        mTargetOuputs[i] = targetOutputs[i];
                    }

                    cost = RunTraining( );
                }

                return cost;
            }

            // Trains single epoch using batches of the specified size (samples are provided as vectors)
            float_t XNetworkTraining::TrainEpoch( const vector<fvector_t>& inputs,
                                                const vector<fvector_t>& targetOutputs,
                                                size_t batchSize, bool randomPickIntoBatch )
            {
                // It is not an average cost for all samples after completion of an epoch, since that
                // requires re-testing all samples. Instead it is an average cost over all batches.
                // However, at the end of each batch the model is updated and so may improve.
                float_t averageRunningCost = 0;
                size_t  samplesCount       = inputs.size( );

                if ( samplesCount != 0 )
                {
                    AllocateTrainVectors( batchSize );

                    if ( ( inputs.size( ) == batchSize ) && ( !randomPickIntoBatch ) )
                    {
                        averageRunningCost = TrainBatch( inputs, targetOutputs );
                    }
                    else
                    {
                        size_t iterations = ( inputs.size( ) - 1 ) / batchSize + 1;
                        
                        for ( size_t i = 0; i < iterations; i++ )
                        {
                            // prepare inputs vectors and target ouputs
                            for ( size_t j = 0; j < batchSize; j++ )
                            {
                                size_t sampleIndex;

                                if ( !randomPickIntoBatch )
                                {
                                    sampleIndex = ( i * batchSize + j ) % samplesCount;
                                }
                                else
                                {
                                    sampleIndex = rand( ) % samplesCount;
                                }

                                mTrainInputs[j]  = const_cast<fvector_t*>( &( inputs[sampleIndex] ) );
                                mTargetOuputs[j] = const_cast<fvector_t*>( &( targetOutputs[sampleIndex] ) );
                            }

                            averageRunningCost += RunTraining( );
                        }

                        averageRunningCost /= iterations;
                    }
                }

                return averageRunningCost;
            }

            // Trains single epoch using batches of the specified size (samples are provided as pointers to vectors)
            float_t XNetworkTraining::TrainEpoch( const vector<fvector_t*>& inputs,
                                                const vector<fvector_t*>& targetOutputs,
                                                size_t batchSize, bool randomPickIntoBatch )
            {
                float_t averageRunningCost = 0;
                size_t  samplesCount       = inputs.size( );

                if ( samplesCount != 0 )
                {
                    AllocateTrainVectors( batchSize );

                    if ( ( inputs.size( ) == batchSize ) && ( !randomPickIntoBatch ) )
                    {
                        averageRunningCost = TrainBatch( inputs, targetOutputs );
                    }
                    else
                    {
                        size_t iterations = ( inputs.size( ) - 1 ) / batchSize + 1;
                        
                        for ( size_t i = 0; i < iterations; i++ )
                        {
                            // prepare inputs vectors and target ouputs
                            for ( size_t j = 0; j < batchSize; j++ )
                            {
                                size_t sampleIndex;

                                if ( !randomPickIntoBatch )
                                {
                                    sampleIndex = ( i * batchSize + j ) % samplesCount;
                                }
                                else
                                {
                                    sampleIndex = rand( ) % samplesCount;
                                }

                                mTrainInputs[j]  = const_cast<fvector_t*>( inputs[sampleIndex] );
                                mTargetOuputs[j] = const_cast<fvector_t*>( targetOutputs[sampleIndex] );
                            }

                            averageRunningCost += RunTraining( );
                        }

                        averageRunningCost /= iterations;
                    }
                }

                return averageRunningCost;
            }


            // Tests sample - calculates real output and provides error cost
            float_t XNetworkTraining::TestSample( const fvector_t& input,
                                                const fvector_t& targetOutput,
                                                fvector_t& output )
            {
                float_t cost = 0;

                if ( mNetwork->LayersCount( ) != 0 )
                {
                    // compute the network to get the actual output
                    mComputeInputs[0] = const_cast<fvector_t*>( &input );
                    DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                    output = mComputeOutputsStorage.back( )[0];
                    cost   = mCostFunction->Cost( output, targetOutput );
                }

                return cost;
            }

            // Tests classification for the provided inputs and target labels - provides number of correctly classified
            // samples and average cost (target outputs are used)
            size_t XNetworkTraining::TestClassification( const std::vector<fvector_t>& inputs, const uvector_t& targetLabels,
                                                        const std::vector<fvector_t>& targetOutputs, float_t* pAvgCost )
            {
                size_t  correctLabelsCounter = 0;
                float_t cost = 0;

                if ( ( mNetwork->LayersCount( ) != 0 ) && ( inputs.size( ) != 0 ) )
                {
                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        mComputeInputs[0] = const_cast<fvector_t*>( &( inputs[i] ) );
                        DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                        cost += mCostFunction->Cost( mComputeOutputsStorage.back( )[0], targetOutputs[i] );

                        if ( XDataEncodingTools::MaxIndex( mComputeOutputsStorage.back( )[0] ) == targetLabels[i] )
                        {
                            correctLabelsCounter++;
                        }
                    }

                    cost /= inputs.size( );
                }

                if ( pAvgCost )
                {
                    *pAvgCost = cost;
                }

                return correctLabelsCounter;
            }

            size_t XNetworkTraining::TestClassification( const std::vector<fvector_t*>& inputs, const uvector_t& targetLabels,
                                                        const std::vector<fvector_t*>& targetOutputs, float_t* pAvgCost )
            {
                size_t  correctLabelsCounter = 0;
                float_t cost = 0;

                if ( ( mNetwork->LayersCount( ) != 0 ) && ( inputs.size( ) != 0 ) )
                {
                    for ( size_t i = 0, n = inputs.size( ); i < n; i++ )
                    {
                        mComputeInputs[0] = inputs[i];
                        DoCompute( mComputeInputs, mComputeOutputs, mInferenceContext );

                        cost += mCostFunction->Cost( mComputeOutputsStorage.back( )[0], *( targetOutputs[i] ) );

                        if ( XDataEncodingTools::MaxIndex( mComputeOutputsStorage.back( )[0] ) == targetLabels[i] )
                        {
                            correctLabelsCounter++;
                        }
                    }

                    cost /= inputs.size( );
                }

                if ( pAvgCost )
                {
                    *pAvgCost = cost;
                }

                return correctLabelsCounter;
            }

            // Enumeration values, which tell when to save network's parameters during training
            enum class NetworkSaveMode
            {
                NoSaving                = 0,
                OnValidationImprovement = 1,
                OnEpochEnd              = 2,
                OnTrainingEnd           = 3
            };

            /* Some helpers aimed for internal use mostly, but can be of use to custom training loop inplementations */
            namespace Helpers
            {
                // Structure to keep some common training parameters. All of those are specified by application's code,
                // but can be overridden from command line to make testing simpler without rebuild.
                typedef struct _TrainingParams
                {
                    float   LearningRate;
                    size_t  EpochsCount;
                    size_t  BatchSize;
                    bool    ShowIntermediateBatchCosts;
                    bool    RunPreTrainingTest;
                    bool    RunValidationOnly;
                    NetworkSaveMode SaveMode;
                    std::string     NetworkOutputFileName;
                    std::string     NetworkInputFileName;

                    _TrainingParams( ) :
                        LearningRate( 0.001f ), EpochsCount( 20 ), BatchSize( 48 ),
                        ShowIntermediateBatchCosts( false ), RunPreTrainingTest( true ), RunValidationOnly( false ),
                        SaveMode( NetworkSaveMode::OnValidationImprovement )
                    {
                    }
                }
                TrainingParams;

                // Parse command line extracting common training parameters
                void ParseTrainingParamsCommandLine( int argc, char** argv, TrainingParams* trainingParams );
                // Log common training parameters to stdout
                void PrintTrainingParams( const TrainingParams* trainingParams );
                // Helper function to show some progress bar on stdout
                void UpdateTrainingPogressBar( size_t lastProgress, size_t currentProgress, size_t totalSteps, size_t barLength, char barChar );
                // Prints training epoch progress (%) to stdout
                int ShowTrainingProgress( size_t currentProgress, size_t totalSteps );
                // Erases training progress from stdout (length is provided by previous ShowTrainingProgress() call)
                void EraseTrainingProgress( int stringLength );
            }

            // A helper class which encapsulates training task of a classification problem
            class XClassificationTrainingHelper
            {
            private:
                std::shared_ptr<XNetworkTraining>  mNetworkTraining;
                EpochSelectionMode                 mEpochSelectionMode;

                bool                               mRunPreTrainingTest;
                bool                               mRunValidationOnly;
                bool                               mShowIntermediateBatchCosts;

                NetworkSaveMode                    mNetworkSaveMode;
                std::string                        mNetworkOutputFileName;
                std::string                        mNetworkInputFileName;

                std::vector<fvector_t*>            mValidationInputs;
                std::vector<fvector_t*>            mValidationOutputs;
                uvector_t                          mValidationLabels;

                std::vector<fvector_t*>            mTestInputs;
                std::vector<fvector_t*>            mTestOutputs;
                uvector_t                          mTestLabels;

                int    mArgc;
                char** mArgv;

            public:
                XClassificationTrainingHelper( const std::shared_ptr<XNetworkTraining>& networkTraining,
                                            int argc = 0, char** argv = nullptr );

                // Get/set the mode of selecting data samples while running training epoch
                EpochSelectionMode SamplesSelectionMode( ) const
                {
                    return mEpochSelectionMode;
                }
                void SetSamplesSelectionMode( EpochSelectionMode selectionMode )
                {
                    mEpochSelectionMode = selectionMode;
                }

                // Run or not pre training test on training data to see the initial classification error
                bool RunPreTrainingTest( ) const
                {
                    return mRunPreTrainingTest;
                }
                void SetRunPreTrainingTest( bool runIt )
                {
                    mRunPreTrainingTest = runIt;
                }

                // Run validation only after each training epoch or classification test on training set as well
                bool RunValidationOnly( ) const
                {
                    return mRunValidationOnly;
                }
                void SetRunValidationOnly( bool validationOnly )
                {
                    mRunValidationOnly = validationOnly;
                }

                // Show cost of some training batches or progress bar
                bool ShowIntermediateBatchCosts( ) const
                {
                    return mShowIntermediateBatchCosts;
                }
                void SetShowIntermediateBatchCosts( bool showBatchCost )
                {
                    mShowIntermediateBatchCosts = showBatchCost;
                }

                // Mode of saving network's learnt parameters
                NetworkSaveMode SaveMode( ) const
                {
                    return mNetworkSaveMode;
                }
                void SetSaveMode( NetworkSaveMode saveMode )
                {
                    mNetworkSaveMode = saveMode;
                }

                // File name to save learnt paramters
                std::string OutputFileName( ) const
                {
                    return mNetworkOutputFileName;
                }
                void SetOutputFileName( const std::string outputFileName )
                {
                    mNetworkOutputFileName = outputFileName;
                }

                // File name to load learnt paramters from
                std::string InputFileName( ) const
                {
                    return mNetworkInputFileName;
                }
                void SetInputFileName( const std::string inputFileName )
                {
                    mNetworkInputFileName = inputFileName;
                }

                // Sets validation samples to use for validating classification after each training epch
                // (takes pointers of inputs/outputs, so original data must stay alive)
                void SetValidationSamples( const std::vector<fvector_t>& validationInputs,
                                        const std::vector<fvector_t>& validationOutputs,
                                        const uvector_t& validationLabels );

                // Sets test samples to use for testing classification after training is complete
                // (takes pointers of inputs/outputs, so original data must stay alive)
                void SetTestSamples( const std::vector<fvector_t>& testInputs,
                                    const std::vector<fvector_t>& testOutputs,
                                    const uvector_t& testLabels );

                // Runs training loop providing progress to stdout
                void RunTraining( size_t epochs, size_t batchSize,
                                const std::vector<fvector_t>& trainingInputs,
                                const std::vector<fvector_t>& trainingOutputs,
                                const uvector_t& trainingLabels );
            };

            namespace Helpers {

            // Parse command line extracting common training parameters
            void ParseTrainingParamsCommandLine( int argc, char** argv, TrainingParams* trainingParams )
            {
                bool showUsage = false;

                if ( argv == nullptr )
                {
                    return;
                }

                for ( int i = 1; i < argc; i++ )
                {
                    bool   parsed   = false;
                    size_t paramLen = strlen( argv[i] );

                    if ( paramLen >= 2 )
                    {
                        char* paramStart = &( argv[i][1] );

                        if ( ( argv[i][0] == '-' ) || ( argv[i][0] == '/' ) )
                        {
                            if ( ( strstr( paramStart, "bs:" ) == paramStart ) && ( paramLen > 4 ) )
                            {
                                if ( sscanf( &( argv[i][4] ), "%zu", &trainingParams->BatchSize ) == 1 )
                                {
                                    if ( trainingParams->BatchSize == 0 )
                                    {
                                        trainingParams->BatchSize = 1;
                                    }
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "ec:" ) == paramStart ) && ( paramLen > 4 ) )
                            {
                                if ( sscanf( &( argv[i][4] ), "%zu", &trainingParams->EpochsCount) == 1 )
                                {
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "lr:" ) == paramStart ) && ( paramLen > 4 ) )
                            {
                                if ( sscanf( &( argv[i][4] ), "%f", &trainingParams->LearningRate ) == 1 )
                                {
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "showBatch:" ) == paramStart ) && ( paramLen == 12 ) )
                            {
                                if ( ( argv[i][11] == '0' ) || ( argv[i][11] == '1' ) )
                                {
                                    trainingParams->ShowIntermediateBatchCosts = ( argv[i][11] == '1' );
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "runPreTrain:" ) == paramStart ) && ( paramLen == 14 ) )
                            {
                                if ( ( argv[i][13] == '0' ) || ( argv[i][13] == '1' ) )
                                {
                                    trainingParams->RunPreTrainingTest = ( argv[i][13] == '1' );
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "validateOnly:" ) == paramStart ) && ( paramLen == 15 ) )
                            {
                                if ( ( argv[i][14] == '0' ) || ( argv[i][14] == '1' ) )
                                {
                                    trainingParams->RunValidationOnly = ( argv[i][14] == '1' );
                                    parsed = true;
                                }
                            }
                            else if ( ( strstr( paramStart, "fin:" ) == paramStart ) && ( paramLen > 5 ) )
                            {
                                trainingParams->NetworkInputFileName = string( &( argv[i][5] ) );
                                parsed = true;
                            }
                            else if ( ( strstr( paramStart, "fout:" ) == paramStart ) && ( paramLen > 6 ) )
                            {
                                trainingParams->NetworkOutputFileName = string( &( argv[i][6] ) );
                                parsed = true;
                            }
                            else if ( ( strstr( paramStart, "sm:" ) == paramStart ) && ( paramLen == 5 ) )
                            {
                                int saveMode = argv[i][4] - '0';

                                if ( ( saveMode >= 1 ) && ( saveMode <= 3 ) )
                                {
                                    trainingParams->SaveMode= static_cast<NetworkSaveMode>( saveMode );
                                    parsed = true;
                                }
                            }
                        }
                    }

                    if ( !parsed )
                    {
                        showUsage = true;
                    }
                }

                if ( showUsage )
                {
                    printf( "Failed parsing some of the parameters \n\n" );

                    printf( "Available parameters are:\n" );
                    printf( "  -ec:<> - epochs count; \n" );
                    printf( "  -bs:<> - batch size; \n" );
                    printf( "  -lr:<> - learning rate; \n" );
                    printf( "  -showBatch:<0|1> - show or not intermediate batch cost; \n" );
                    printf( "  -runPreTrain:<0|1> - run or not pre training test on training data; \n" );
                    printf( "  -validateOnly:<0|1> - run test on validation data only or on test data as well after each epoch; \n" );
                    printf( "  -fin:<file name> - file to load network's parameters from; \n" );
                    printf( "  -fout:<file name> - file to save network's parameters to; \n" );
                    printf( "  -sm:<> - save mode: 1 - on validation improvement (default); \n" );
                    printf( "                      2 - at the end of each epoch; \n" );
                    printf( "                      3 - at the end of training. \n" );
                    printf( "\n" );
                }

                if ( trainingParams->NetworkOutputFileName.empty( ) )
                {
                    trainingParams->SaveMode = NetworkSaveMode::NoSaving;
                }
            }

            // Log common training parameters to stdout
            void PrintTrainingParams( const TrainingParams* trainingParams )
            {
                printf( "Learning rate: %0.4f, Epochs: %zu, Batch Size: %zu \n", trainingParams->LearningRate, trainingParams->EpochsCount, trainingParams->BatchSize );
                if ( !trainingParams->NetworkInputFileName.empty( ) )
                {
                    printf( "Network input file: %s \n", trainingParams->NetworkInputFileName.c_str( ) );
                }
                if ( ( !trainingParams->NetworkOutputFileName.empty( ) ) && ( trainingParams->SaveMode != NetworkSaveMode::NoSaving ) )
                {
                    printf( "Network output file: %s \n", trainingParams->NetworkOutputFileName.c_str( ) );
                }
                printf( "\n" );
            }

            // Helper function to show some progress bar on stdout
            void UpdateTrainingPogressBar( size_t lastProgress, size_t currentProgress, size_t totalSteps, size_t barLength, char barChar )
            {
                size_t barsDone = lastProgress    * barLength / totalSteps;
                size_t barsNeed = currentProgress * barLength / totalSteps;

                while ( barsDone++ != barsNeed )
                {
                    putchar( barChar );
                }
                fflush( stdout );
            }

            // Prints training epoch progress (%) to stdout
            int ShowTrainingProgress( size_t currentProgress, size_t totalSteps )
            {
                int printed = printf( "<%d%%>", static_cast<int>( currentProgress * 100 / totalSteps ) );
                fflush( stdout );
                return printed;
            }

            // Erases training progress from stdout (length is provided by previous ShowTrainingProgress() call)
            void EraseTrainingProgress( int stringLength )
            {
                while ( stringLength > 0 )
                {
                    printf( "\b \b" );
                    stringLength--;
                }
            }

            } // namespace Helpers

            // ========================================================================================================================

            XClassificationTrainingHelper::XClassificationTrainingHelper( const shared_ptr<XNetworkTraining>& networkTraining,
                                                                        int argc, char** argv ) :
                mNetworkTraining( networkTraining ),
                mEpochSelectionMode( EpochSelectionMode::Shuffle ),
                mRunPreTrainingTest( true ), mRunValidationOnly( false ),
                mShowIntermediateBatchCosts( false ),
                mNetworkSaveMode( NetworkSaveMode::OnValidationImprovement ), mNetworkOutputFileName( ), mNetworkInputFileName( ),
                mArgc( argc ), mArgv( argv )
            {

            }

            // Sets validation samples to use for validating classification after each training epch
            void XClassificationTrainingHelper::SetValidationSamples( const vector<fvector_t>& validationInputs,
                                                                    const vector<fvector_t>& validationOutputs,
                                                                    const uvector_t& validationLabels )
            {
                size_t samplesCount = validationInputs.size( );

                mValidationInputs.resize( validationInputs.size( ) );
                mValidationOutputs.resize( validationOutputs.size( ) );
                mValidationLabels = validationLabels;

                for ( size_t i = 0; i < samplesCount; i++ )
                {
                    mValidationInputs[i]  = const_cast<fvector_t*>( &( validationInputs[i] ) );
                    mValidationOutputs[i] = const_cast<fvector_t*>( &( validationOutputs[i] ) );
                }
            }

            // Sets test samples to use for testing classification after training is complete
            void XClassificationTrainingHelper::SetTestSamples( const std::vector<fvector_t>& testInputs,
                                                                const std::vector<fvector_t>& testOutputs,
                                                                const uvector_t& testLabels )
            {
                size_t samplesCount = testInputs.size( );

                mTestInputs.resize( testInputs.size( ) );
                mTestOutputs.resize( testOutputs.size( ) );
                mTestLabels = testLabels;

                for ( size_t i = 0; i < samplesCount; i++ )
                {
                    mTestInputs[i] = const_cast<fvector_t*>( &( testInputs[i] ) );
                    mTestOutputs[i] = const_cast<fvector_t*>( &( testOutputs[i] ) );
                }
            }

            // Runs training loop providing progress to stdout
            void XClassificationTrainingHelper::RunTraining( size_t epochs, size_t batchSize,
                                                            const vector<fvector_t>& trainingInputs,
                                                            const vector<fvector_t>& trainingOutputs,
                                                            const uvector_t& trainingLabels )
            {
                // default training parameters
                Helpers::TrainingParams     trainingParams;

                trainingParams.EpochsCount  = epochs;
                trainingParams.BatchSize    = batchSize;
                trainingParams.LearningRate = mNetworkTraining->Optimizer( )->LearningRate( );

                trainingParams.ShowIntermediateBatchCosts = mShowIntermediateBatchCosts;
                trainingParams.RunPreTrainingTest         = mRunPreTrainingTest;
                trainingParams.RunValidationOnly          = mRunValidationOnly;

                trainingParams.SaveMode              = mNetworkSaveMode;
                trainingParams.NetworkOutputFileName = mNetworkOutputFileName;
                trainingParams.NetworkInputFileName  = mNetworkInputFileName;

                // parse command line for any overrides
                Helpers::ParseTrainingParamsCommandLine( mArgc, mArgv, &trainingParams );

                // set some of the new parameters
                mNetworkTraining->Optimizer( )->SetLearningRate( trainingParams.LearningRate );

                // log current settings
                Helpers::PrintTrainingParams( &trainingParams );

                // load network parameters from the previous save file
                if ( !trainingParams.NetworkInputFileName.empty( ) )
                {
                    if ( !mNetworkTraining->Network( )->LoadLearnedParams( trainingParams.NetworkInputFileName ) )
                    {
                        printf( "Failed loading network's parameters \n\n" );
                    }
                }

                // 
                vector<fvector_t*> trainingInputsPtr( trainingInputs.size( ) );
                vector<fvector_t*> trainingOutputsPtr( trainingOutputs.size( ) );
                vector<fvector_t*> trainingInputsBatch( trainingParams.BatchSize );
                vector<fvector_t*> trainingOutputsBatch( trainingParams.BatchSize );

                size_t             samplesCount       = trainingInputs.size( );
                size_t             iterationsPerEpoch = ( samplesCount - 1 ) / trainingParams.BatchSize + 1;

                float              lastValidationAccuracy = 0.0f;
                float_t            cost;
                size_t             correct;

                steady_clock::time_point timeStartForAll = steady_clock::now( );
                steady_clock::time_point timeStart;
                long long                timeTaken;

                size_t             batchCostOutputFreq  = iterationsPerEpoch / 80;
                int                progressStringLength = 0;

                if ( batchCostOutputFreq == 0 )
                {
                    batchCostOutputFreq = 1;
                }

                // take pointers to original inputs/outputs, so those could be shuffled 
                for ( size_t i = 0; i < samplesCount; i++ )
                {
                    trainingInputsPtr[i]  = const_cast<fvector_t*>( &( trainingInputs[i]  ) );
                    trainingOutputsPtr[i] = const_cast<fvector_t*>( &( trainingOutputs[i] ) );
                }

                // check classification error before starting training
                if ( trainingParams.RunPreTrainingTest )
                {
                    timeStart = steady_clock::now( );
                    correct   = mNetworkTraining->TestClassification( trainingInputs, trainingLabels, trainingOutputs, &cost );
                    timeTaken = duration_cast<milliseconds>( steady_clock::now( ) - timeStart ).count( );

                    printf( "Before training: accuracy = %0.2f%% (%zu/%zu), cost = %0.4f, %0.3fs \n\n",
                            static_cast<float>( correct ) / trainingInputs.size( ) * 100,
                            correct, trainingInputs.size( ), static_cast<float>( cost ),
                            static_cast<float>( timeTaken ) / 1000 );
                }

                // run the specified number of epochs
                for ( size_t epoch = 0; epoch < trainingParams.EpochsCount; epoch++ )
                {
                    printf( "Epoch %3zu : ", epoch + 1 );
                    if ( !trainingParams.ShowIntermediateBatchCosts )
                    {
                        // show progress bar only
                        putchar( '[' );
                    }
                    else
                    {
                        printf( "\n" );
                    }

                    // shuffle samples if required
                    if ( mEpochSelectionMode == EpochSelectionMode::Shuffle )
                    {
                        for ( size_t i = 0; i < samplesCount / 2; i++ )
                        {
                            int swapIndex1 = rand( ) % samplesCount;
                            int swapIndex2 = rand( ) % samplesCount;

                            std::swap( trainingInputsPtr[swapIndex1],  trainingInputsPtr[swapIndex2]  );
                            std::swap( trainingOutputsPtr[swapIndex1], trainingOutputsPtr[swapIndex2] );
                        }
                    }

                    // start of epoch timing
                    timeStart = steady_clock::now( );

                    for ( size_t iteration = 0; iteration < iterationsPerEpoch; iteration++ )
                    {
                        // prepare batch inputs and ouputs
                        for ( size_t i = 0; i < trainingParams.BatchSize; i++ )
                        {
                            size_t sampleIndex;

                            if ( mEpochSelectionMode == EpochSelectionMode::RandomPick )
                            {
                                sampleIndex = rand( ) % samplesCount;
                            }
                            else
                            {
                                sampleIndex = ( iteration * trainingParams.BatchSize + i ) % samplesCount;
                            }

                            trainingInputsBatch[i]  = trainingInputsPtr[sampleIndex];
                            trainingOutputsBatch[i] = trainingOutputsPtr[sampleIndex];
                        }

                        float_t batchCost = mNetworkTraining->TrainBatch( trainingInputsBatch, trainingOutputsBatch );

                        // erase previous progress if any 
                        Helpers::EraseTrainingProgress( progressStringLength );

                        // show cost of some batches or progress bar only
                        if ( !trainingParams.ShowIntermediateBatchCosts )
                        {
                            Helpers::UpdateTrainingPogressBar( iteration, iteration + 1, iterationsPerEpoch, 50, '=' );
                        }
                        else
                        {
                            if ( ( ( iteration + 1 ) % batchCostOutputFreq ) == 0 )
                            {
                                printf( "%0.4f ", static_cast<float>( batchCost ) );

                                if ( ( ( iteration + 1 ) % ( batchCostOutputFreq * 8 ) ) == 0 )
                                {
                                    printf( "\n" );
                                }
                            }
                        }

                        // show current progress of the epoch
                        progressStringLength = Helpers::ShowTrainingProgress( iteration + 1, iterationsPerEpoch );
                    }

                    Helpers::EraseTrainingProgress( progressStringLength );
                    progressStringLength = 0;

                    // end of epoch timing
                    timeTaken = duration_cast<milliseconds>( steady_clock::now( ) - timeStart ).count( );

                    // output time spent on training
                    if ( !trainingParams.ShowIntermediateBatchCosts )
                    {
                        printf( "] " );
                    }
                    else
                    {
                        printf( "\nTime taken : " );
                    }
                    printf( "%0.3fs\n", static_cast<float>( timeTaken ) / 1000 );

                    float validationAccuracy = 0.0f;

                    // get classification error on training data after completion of an epoch
                    if ( ( !trainingParams.RunValidationOnly ) || ( mValidationInputs.size( ) == 0 ) )
                    {
                        timeStart = steady_clock::now( );
                        correct   = mNetworkTraining->TestClassification( trainingInputs, trainingLabels, trainingOutputs, &cost );
                        timeTaken = duration_cast<milliseconds>( steady_clock::now( ) - timeStart ).count( );

                        printf( "Training accuracy = %0.2f%% (%zu/%zu), cost = %0.4f, %0.3fs \n",
                                static_cast<float>( correct ) / trainingInputs.size( ) * 100,
                                correct, trainingInputs.size( ), static_cast<float>( cost ),
                                static_cast<float>( timeTaken ) / 1000 );

                        // use training accuracy, if validation data set is not provided
                        if ( mValidationInputs.size( ) == 0 )
                        {
                            validationAccuracy = static_cast<float>( correct ) / trainingInputs.size( );
                        }
                    }

                    // use validation set to check classification error on data not included into training
                    if ( mValidationInputs.size( ) != 0 )
                    {
                        timeStart = steady_clock::now( );
                        correct   = mNetworkTraining->TestClassification( mValidationInputs, mValidationLabels, mValidationOutputs, &cost );
                        timeTaken = duration_cast<milliseconds>( steady_clock::now( ) - timeStart ).count( );

                        printf( "Validation accuracy = %0.2f%% (%zu/%zu), cost = %0.4f, %0.3fs \n",
                                static_cast<float>( correct ) / mValidationInputs.size( ) * 100,
                                correct, mValidationInputs.size( ), static_cast<float>( cost ),
                                static_cast<float>( timeTaken ) / 1000 );

                        validationAccuracy = static_cast<float>( correct ) / mValidationInputs.size( );
                    }

                    // save network at the end of epoch
                    if ( trainingParams.SaveMode == NetworkSaveMode::OnEpochEnd )
                    {
                        mNetworkTraining->Network( )->SaveLearnedParams( trainingParams.NetworkOutputFileName );
                    }
                    else if ( ( trainingParams.SaveMode == NetworkSaveMode::OnValidationImprovement ) &&
                            ( validationAccuracy > lastValidationAccuracy ) )
                    {
                        mNetworkTraining->Network( )->SaveLearnedParams( trainingParams.NetworkOutputFileName );
                        lastValidationAccuracy = validationAccuracy;
                    }
                }

                // final test on test data
                if ( mTestInputs.size( ) != 0 )
                {
                    timeStart = steady_clock::now( );
                    correct   = mNetworkTraining->TestClassification( mTestInputs, mTestLabels, mTestOutputs, &cost );
                    timeTaken = duration_cast<milliseconds>( steady_clock::now( ) - timeStart ).count( );

                    printf( "\nTest accuracy = %0.2f%% (%zu/%zu), cost = %0.4f, %0.3fs \n",
                        static_cast<float>( correct ) / mTestInputs.size( ) * 100,
                        correct, mTestInputs.size( ), static_cast<float>( cost ),
                        static_cast<float>( timeTaken ) / 1000 );
                }

                // total time taken by the training
                timeTaken = duration_cast<seconds>( steady_clock::now( ) - timeStartForAll ).count( );
                printf( "\nTotal time taken : %ds (%0.2fmin) \n", static_cast<int>( timeTaken ), static_cast<float>( timeTaken ) / 60 );

                // save network when training is done
                if ( trainingParams.SaveMode == NetworkSaveMode::OnTrainingEnd )
                {
                    mNetworkTraining->Network( )->SaveLearnedParams( trainingParams.NetworkOutputFileName );
                }
            }
        }    
        XNetworkContext::XNetworkContext( bool trainingMode, size_t sequenceLength ) :
            mTrainingMode( trainingMode ), mTrainingSequenceLength( sequenceLength ), mCurrentLayer( 0 )
        {
        }

        XNetworkContext::~XNetworkContext( )
        {
            FreeWorkingBuffers( );
        }

        // Allocate working buffer for laters of the network
        void XNetworkContext::AllocateWorkingBuffers( const std::shared_ptr<XNeuralNetwork>& net, size_t batchSize )
        {
            FreeWorkingBuffers( );

            for ( auto layer : *net )
            {
                uvector_t workingMemSize = layer->WorkingMemSize( mTrainingMode );

                // filling this nice vector
                //
                //       -- for each layer
                //       |           -- for each requested buffer
                //       |           |           -- for each sample
                //       |           |           |       -- requested memory buffer
                //       |           |           |       |
                // std::vector<std::vector<std::vector<void*>>>

                mLayersMemorySize.push_back( workingMemSize );
                mLayersMemoryBuffers.push_back( vector<vector<void*>>( workingMemSize.size( ) ) );

                for ( size_t i = 0; i < workingMemSize.size( ); i++ )
                {
                    mLayersMemoryBuffers.back( )[i] = vector<void*>( );

                    for ( size_t j = 0; j < batchSize; j++ )
                    {
                        void* memBuffer = AlignedAlloc( 32, workingMemSize[i] );

                        if ( memBuffer )
                        {
                            memset( memBuffer, 0, workingMemSize[i] );
                        }

                        mLayersMemoryBuffers.back( )[i].push_back( memBuffer );
                    }
                }
            }
        }

        // Free layers' working buffers
        void XNetworkContext::FreeWorkingBuffers( )
        {
            for ( size_t i = 0; i < mLayersMemoryBuffers.size( ); i++ )
            {
                for ( size_t j = 0; j < mLayersMemoryBuffers[i].size( ); j++ )
                {
                    for ( size_t k = 0; k < mLayersMemoryBuffers[i][j].size( ); k++ )
                    {
                        AlignedFree( mLayersMemoryBuffers[i][j][k] );
                    }
                }
            }

            mLayersMemoryBuffers.clear( );
            mLayersMemorySize.clear( );
        }

        // Clear layers' working buffers (memset zero) 
        void XNetworkContext::ResetWorkingBuffers( )
        {
            for ( size_t i = 0; i < mLayersMemoryBuffers.size( ); i++ )
            {
                for ( size_t j = 0; j < mLayersMemoryBuffers[i].size( ); j++ )
                {
                    for ( size_t k = 0; k < mLayersMemoryBuffers[i][j].size( ); k++ )
                    {
                        memset( mLayersMemoryBuffers[i][j][k], 0, mLayersMemorySize[i][j] );
                    }
                }
            }
        }
        void XNetworkContext::ResetWorkingBuffers( uvector_t layersIndexes )
        {
            for ( size_t i : layersIndexes )
            {
                for ( size_t j = 0; j < mLayersMemoryBuffers[i].size( ); j++ )
                {
                    for ( size_t k = 0; k < mLayersMemoryBuffers[i][j].size( ); k++ )
                    {
                        memset( mLayersMemoryBuffers[i][j][k], 0, mLayersMemorySize[i][j] );
                    }
                }
            }
        }
    }     
}