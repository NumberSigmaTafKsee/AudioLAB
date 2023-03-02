// https://github.com/JimClay/NimbleDSP
#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "kissfft.hh"

namespace NimbleDSP {

enum FilterOperationType {STREAMING, ONE_SHOT_RETURN_ALL_RESULTS, ONE_SHOT_TRIM_TAILS};
typedef enum ParksMcClellanFilterType {PASSBAND_FILTER = 1, DIFFERENTIATOR_FILTER, HILBERT_FILTER} ParksMcClellanFilterType;

const unsigned DEFAULT_BUF_LEN = 0;

#ifndef SLICKDSP_FLOAT_TYPE
#define SLICKDSP_FLOAT_TYPE    double
#endif

#define VECTOR_TO_ARRAY(x)      (&((x)[0]))


/**
 * \brief Base class for NimbleDSP.
 *
 * Although you can instantiate objects of this type, that's not what this class is intended for.  It is the
 * base class that all of the other classes descend from which allows for a great deal of flexibility
 * through polymorphism.  It also reduces the amount of code because we don't have to replicate the same
 * functionality in each class.
 *
 * Derived classes: RealVector and ComplexVector.
 */
template <class T>
class Vector {

 public:
    /** 
     * \brief Buffer to store intermediate calculations when needed.
     */
    std::vector<T> *scratchBuf;

 protected:
    /** 
     * \brief Initializes vec to a given size and fills it with zeros.
     */
    void initSize(unsigned size) {vec = std::vector<T>(size);}
    
    /** 
     * \brief Initializes vec with the size and contents of "array".
     *
     * \param array Array to set vec equal to.
     * \param arrayLen Number of elements in array.
     */
    template <class U>
    void initArray(U *array, unsigned arrayLen);
    
 public:
    /*************************************************************************************//**
     * \brief Vector that holds the object's data.
     *
     * The class is built around this member.  Std::vector's were used because they handle the 
     * dynamic memory, have a rich set of support functions, are fast and efficient, and can
     * be accessed like a normal array when that is convenient.
     *****************************************************************************************/
	std::vector<T> vec;
    
    template <class U> friend class Vector;
    template <class U> friend class RealVector;
    template <class U> friend class ComplexVector;
    template <class U> friend class RealFirFilter;
    template <class U> friend class ComplexFirFilter;
    
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref vec and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref vec.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    Vector<T>(unsigned size = 0, std::vector<T> *scratch = NULL) {initSize(size); scratchBuf = scratch;}
    
    /**
     * \brief Vector constructor.
     *
     * Sets vec equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref vec will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    Vector<T>(std::vector<U> data, std::vector<T> *scratch = NULL) {initArray(VECTOR_TO_ARRAY(data), (unsigned) data.size()); scratchBuf = scratch;}
    
    /**
     * \brief Array constructor.
     *
     * Sets vec equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref vec will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    Vector<T>(U *data, unsigned dataLen, std::vector<T> *scratch = NULL) {initArray(data, dataLen); scratchBuf = scratch;}
    
    /**
     * \brief Copy constructor.
     */
    Vector<T>(const Vector<T>& other) {vec = other.vec; scratchBuf = other.scratchBuf;}

	/**
	 * \brief Virtual destructor.
	 */
	virtual ~Vector() {};

    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Index assignment operator.
     */
    T& operator[](unsigned index) {return vec[index];};
    
    /**
     * \brief Index operator.
     */
    const T& operator[](unsigned index) const {return vec[index];};
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Returns the size of \ref vec.
     */
    const unsigned size() const {return (const unsigned) vec.size();};
    
    /**
     * \brief Finds the first instance of "val" in \ref vec.
     *
     * \param val The value to look for in \ref vec.
     * \return Index of first instance of "val".  If there aren't any elements equal to "val"
     *      it returns -1.
     */
    const int find(const T val) const;
    
    /**
     * \brief Returns the sum of all the elements in \ref vec.
     */
	T sum() const;

};


template <class T>
template <class U>
void Vector<T>::initArray(U *array, unsigned arrayLen) {
    vec = std::vector<T>(arrayLen);
    for (unsigned i=0; i<arrayLen; i++) {
        vec[i] = (T) array[i];
    }
}

template <class T>
inline bool operator==(const Vector<T>& lhs, const Vector<T>& rhs) {
    if (lhs.size() != rhs.size())
        return false;
    
    for (unsigned i=0; i<lhs.size(); i++) {
        if (lhs[i] != rhs[i])
            return false;
    }
    return true;
}

template <class T>
inline bool operator!=(const Vector<T>& lhs, const Vector<T>& rhs) {return !(lhs == rhs);}

template <class T>
const int Vector<T>::find(const T val) const {
    for (unsigned i=0; i<size(); i++) {
        if (vec[i] == val) {
            return (int) i;
        }
    }
    return -1;
}

/**
 * \brief Finds the first instance of "val" in \ref vec.
 *
 * \param vector Buffer to operate on.
 * \param val The value to look for in \ref vec.
 * \return Index of first instance of "val".  If there aren't any elements equal to "val"
 *      it returns -1.
 */
template <class T>
const int find(Vector<T> & vector, const T val) {
    return vector.find(val);
}

template <class T>
T Vector<T>::sum() const {
	assert(vec.size() > 0);
	T vectorSum = 0;
	for (unsigned i=0; i<vec.size(); i++) {
		vectorSum += vec[i];
	}
	return vectorSum;
}

/**
 * \brief Returns the sum of all the elements in \ref vec.
 *
 * \param vector Buffer to operate on.
 */
template <class T>
T sum(const Vector<T> & vector) {
	return vector.sum();
}



Skip to content
Pull requests
Issues
Codespaces
Marketplace
Explore
@NumberSigmaTafKsee
JimClay /
NimbleDSP
Public

Fork your own copy of JimClay/NimbleDSP

Code
Issues 5
Pull requests
Actions
Projects
Wiki
Security

    Insights

NimbleDSP/src/ComplexVector.h
@JimClay
JimClay Made time/frequeny domain checks optional.
Latest commit dfddfa8 Jun 9, 2018
History
1 contributor
1796 lines (1604 sloc) 60 KB
/*
Copyright (c) 2014, James Clay
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


/**
 * @file ComplexVector.h
 *
 * Definition of the template class ComplexVector.
 */

enum DomainType {TIME_DOMAIN, FREQUENCY_DOMAIN};

/**
 * \brief Vector class for complex numbers.
 *
 * The template type should be the "plain old data" type that you want to use, not "std::complex"
 * or your own custom complex class.  The object will automatically convert the buffer type to
 * std::complex<POD_type> for you.
 */
template <class T>
class ComplexVector : public Vector< std::complex<T> > {
 public:
    template <class U> friend class ComplexFirFilter;

    /**
     * \brief Indicates whether the data in \ref buf is time domain data or frequency domain.
     */
    DomainType domain;
    
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Sets \ref domain to NimbleDSP::TIME_DOMAIN and calls Vector<T>::Vector(size, scratch).
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    ComplexVector<T>(unsigned size = DEFAULT_BUF_LEN, std::vector< std::complex<T> > *scratch = NULL) :
            Vector< std::complex<T> >(size, scratch) {domain = TIME_DOMAIN;}
            
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.  Also sets \ref domain to dataDomain.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     * \param dataDomain Indicates whether the data is time domain data or frequency domain.
     *      Valid values are NimbleDSP::TIME_DOMAIN and NimbleDSP::FREQUENCY_DOMAIN.
     */
    template <typename U>
    ComplexVector<T>(std::vector<U> data, DomainType dataDomain=TIME_DOMAIN,
                std::vector< std::complex<T> > *scratch = NULL) : Vector< std::complex<T> >(data, scratch)
                    {domain = dataDomain;}
    
    /**
     * \brief Array constructor.
     *
     * Sets buf equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.  Also sets \ref domain to dataDomain.
     * \param data Array that \ref buf will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     * \param dataDomain Indicates whether the data is time domain data or frequency domain.
     *      Valid values are NimbleDSP::TIME_DOMAIN and NimbleDSP::FREQUENCY_DOMAIN.
     */
    template <typename U>
    ComplexVector<T>(U *data, unsigned dataLen, DomainType dataDomain=TIME_DOMAIN,
                std::vector< std::complex<T> > *scratch = NULL) : Vector< std::complex<T> >(data, dataLen, scratch)
                    {domain = dataDomain;}
    
    /**
     * \brief Copy constructor.
     */
    ComplexVector<T>(const ComplexVector<T>& other)
            {this->vec = other.vec; domain = other.domain; this->scratchBuf = other.scratchBuf;}
    
    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Assignment operator from ComplexVector.
     * \return Reference to "this".
     */
    ComplexVector<T>& operator=(const ComplexVector<T>& rhs);
    
    /**
     * \brief Assignment operator from Vector.
     * \return Reference to "this".
     */
    ComplexVector<T>& operator=(const Vector<T>& rhs);
    
    /**
     * \brief Unary minus (negation) operator.
     * \return Reference to "this".
     */
    ComplexVector<T> & operator-();
    
    /**
     * \brief Add Buffer/Assignment operator.
     * \return Reference to "this".
     */
    template <class U>
    ComplexVector<T> & operator+=(const Vector<U> &rhs);
    
    /**
     * \brief Add Scalar/Assignment operator.
     * \return Reference to "this".
     */
    ComplexVector<T> & operator+=(const std::complex<T> &rhs);
    
    /**
     * \brief Subtract Buffer/Assignment operator.
     * \return Reference to "this".
     */
    template <class U>
    ComplexVector<T> & operator-=(const Vector<U> &rhs);
    
    /**
     * \brief Subtract Scalar/Assignment operator.
     * \return Reference to "this".
     */
    ComplexVector<T> & operator-=(const std::complex<T> &rhs);
    
    /**
     * \brief Multiply Buffer/Assignment operator.
     * \return Reference to "this".
     */
    template <class U>
    ComplexVector<T> & operator*=(const Vector<U> &rhs);
    
    /**
     * \brief Multiply Scalar/Assignment operator.
     * \return Reference to "this".
     */
    ComplexVector<T> & operator*=(const std::complex<T> &rhs);
    
    /**
     * \brief Divide Buffer/Assignment operator.
     * \return Reference to "this".
     */
    template <class U>
    ComplexVector<T> & operator/=(const Vector<U> &rhs);
    
    /**
     * \brief Divide Scalar/Assignment operator.
     * \return Reference to "this".
     */
    ComplexVector<T> & operator/=(const std::complex<T> &rhs);
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Sets each element of \ref vec to e^(element).
     *
     * \return Reference to "this".
     */
    //virtual Vector< std::complex<T> > & exp();
    
    /**
     * \brief Sets each element of \ref buf equal to its value to the power of "exponent".
     *
     * \param exponent Exponent to use.
     * \return Reference to "this".
     */
    ComplexVector<T> & pow(const std::complex<SLICKDSP_FLOAT_TYPE> & exponent);

    /**
     * \brief Returns the element with the maximum real component in \ref buf.
     *
     * \param maxLoc If it isn't equal to NULL the index of the maximum element
     *      will be returned via this pointer.  If more than one element is equal
     *      to the maximum value the index of the first will be returned.
     *      Defaults to NULL.
     */
    const T max(unsigned *maxLoc = NULL) const;
    
    /**
     * \brief Returns the mean (average) of the data in \ref buf.
     */
    const std::complex<SLICKDSP_FLOAT_TYPE> mean() const;
    
    /**
     * \brief Returns the variance of the data in \ref buf.
     */
    const SLICKDSP_FLOAT_TYPE var() const;
    
    /**
     * \brief Returns the standard deviation of the data in \ref buf.
     */
    const SLICKDSP_FLOAT_TYPE stdDev() const {return std::sqrt(this->var());}
    
    /**
     * \brief Sets the upper and lower limit of the values in \ref buf.
     *
     * \param val Limiting value for the data in \ref buf.  Any values that
     *      are greater than "val" are made equal to "val", and
     *      any that are less than -val are made equal to -val.  This is done
     *      independently on the real and imaginary elements of \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & saturate(const std::complex<T> & val);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & ceil(void);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & floor(void);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & round(void);
    
    /**
     * \brief Conjugates the data in \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & conj();
    
    /**
     * \brief Sets each element of \ref buf equal to its magnitude squared.
     * \return Reference to "this".
     */
    ComplexVector<T> & magSq();
    
    /**
     * \brief Sets each element of \ref buf equal to its angle.
     *
     * The angle is held in the real portion of \ref buf.
     * \return Reference to "this".
     */
    ComplexVector<T> & angle();
    
    /**
     * \brief Sets \ref buf equal to the FFT of the data in \ref buf.
     *
     * Sets \ref domain equal to NimbleDSP::FREQUENCY_DOMAIN.
     * \return Reference to "this".
     */
    ComplexVector<T> & fft();
    
    /**
     * \brief Sets \ref buf equal to the inverse FFT of the data in \ref buf.
     *
     * Sets \ref domain equal to NimbleDSP::TIME_DOMAIN.
     * \return Reference to "this".
     */
    ComplexVector<T> & ifft();
    
    /**
     * \brief Changes the elements of \ref vec to their absolute value.
     *
     * \return Reference to "this".
     */
    ComplexVector<T> & abs();
    
    /**
     * \brief Sets each element of \ref vec to e^(element).
     *
     * \return Reference to "this".
     */
    ComplexVector<T> & exp();
    
    /**
     * \brief Sets each element of \ref vec to the natural log of the element.
     *
     * \return Reference to "this".
     */
    ComplexVector<T> & log();
    
    /**
     * \brief Sets each element of \ref vec to the base 10 log of the element.
     *
     * \return Reference to "this".
     */
    ComplexVector<T> & log10();
    
    /**
     * \brief Circular rotation.
     *
     * \param numToShift Number of positions to shift in the circular rotation.  numToShift
     *      can be positive or negative.  If you visualize the 0 index value at the left and
     *      the end of the array at the right, positive numToShift values shift the array to
     *      the left, and negative values shift it to the right.
     * \return Reference to "this".
     */
    ComplexVector<T> & rotate(int numToShift);
    
    /**
     * \brief Reverses the order of the elements in \ref vec.
     *
     * \return Reference to "this".
     */
    ComplexVector<T> & reverse();
    
    /**
     * \brief Sets the length of \ref vec to "len".
     *
     * \param len The new length for \ref vec.  If len is longer than vec's current size, the
     *      new elements will be set to "val".  If len is less than vec's current size the extra
     *      elements will be cut off and the other elements will remain the same.
     * \param val The value to set any new elements to.  Defaults to 0.
     * \return Reference to "this".
     */
    ComplexVector<T> & resize(unsigned len, T val = (T) 0) {this->vec.resize(len, val); return *this;}

    /**
     * \brief Reserves "len" elements for \ref vec without actually resizing it.
     *
     * \param len The number of elements to reserve for \ref vec.
     * \return Reference to "this".
     */
    ComplexVector<T> & reserve(unsigned len) {this->vec.reserve(len); return *this;}
    
    /**
     * \brief Lengthens \ref vec by "len" elements.
     *
     * \param len The number of elements to add to \ref vec.
     * \param val The value to set the new elements to.  Defaults to 0.
     * \return Reference to "this".
     */
    ComplexVector<T> & pad(unsigned len, T val = (T) 0) {this->vec.resize(this->size()+len, val); return *this;}

    /**
     * \brief Copies the elements in the range "lower" to "upper" (inclusive) to "destination"
     *
     * \param lower The index to start copying at.
     * \param upper The last index to copy from.
     * \destination The vector to copy the vector slice to.
     */
	void slice(unsigned lower, unsigned upper, ComplexVector<T> &destination);
    
    /**
     * \brief Inserts rate-1 zeros between samples.
     *
     * \param rate Indicates how many zeros should be inserted between samples.
     * \param phase Indicates how many of the zeros should be before the samples (as opposed to
     *      after).  Valid values are 0 to "rate"-1.  Defaults to 0.
     * \return Reference to "this".
     */
    ComplexVector<T> & upsample(int rate, int phase = 0);
    
    /**
     * \brief Removes rate-1 samples out of every rate samples.
     *
     * \param rate Indicates how many samples should be removed.
     * \param phase Tells the method which sample should be the first to be kept.  Valid values
     *      are 0 to "rate"-1.  Defaults to 0.
     * \return Reference to "this".
     */
    ComplexVector<T> & downsample(int rate, int phase = 0);
    
    /**
     * \brief Replaces \ref vec with the cumulative sum of the samples in \ref vec.
     *
     * \param initialVal Initializing value for the cumulative sum.  Defaults to zero.
     * \return Reference to "this".
     */
	ComplexVector<T> & cumsum(T initialVal = 0);
    
    /**
     * \brief Replaces \ref vec with the difference between successive samples in vec.
     *
     * The resulting \ref vec is one element shorter than it was previously.
     * \return Reference to "this".
     */
	ComplexVector<T> & diff();
    
    /**
     * \brief Replaces \ref vec with the difference between successive samples in vec.
     *
     * \param previousVal The last value in the sample stream before the current contents
     *      of \ref vec.  previousVal allows the resulting vec to be the same size as the
     *      previous vec.
     * \return Reference to "this".
     */
    ComplexVector<T> & diff(std::complex<T> & previousVal);
    
    /**
     * \brief Convolution method.
     *
     * \param data The vector that will be filtered.
     * \param trimTails "False" tells the method to return the entire convolution, which is
     *      the length of "data" plus the length of "this" (the filter) - 1.  "True" tells the
     *      method to retain the size of "data" by trimming the tails at both ends of
     *      the convolution.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual ComplexVector<T> & conv(ComplexVector<T> & data, bool trimTails = false);
    
    /**
     * \brief Decimate method.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual ComplexVector<T> & decimate(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual ComplexVector<T> & interp(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Resample method.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual ComplexVector<T> & resample(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);
    
    /**
     * \brief Generates a complex tone.
     *
     * \param freq The tone frequency.
     * \param sampleFreq The sample frequency.  Defaults to 1 Hz.
     * \param phase The tone's starting phase, in radians.  Defaults to 0.
     * \param numSamples The number of samples to generate.  "0" indicates to generate
     *      this->size() samples.  Defaults to 0.
     * \return Reference to "this".
     */
    T tone(T freq, T sampleFreq = 1.0, T phase = 0.0, unsigned numSamples = 0);
    
    /**
     * \brief Modulates the data with a complex sinusoid.
     *
     * \param freq The modulating tone frequency.
     * \param sampleFreq The sample frequency of the data.  Defaults to 1 Hz.
     * \param phase The modulating tone's starting phase, in radians.  Defaults to 0.
     * \return The next phase if the tone were to continue.
     */
    T modulate(T freq, T sampleFreq = 1.0, T phase = 0.0);
};


template <class T>
ComplexVector<T>& ComplexVector<T>::operator=(const ComplexVector<T>& rhs)
{
    this->vec = rhs.vec;
    domain = rhs.domain;
    return *this;
}

template <class T>
ComplexVector<T>& ComplexVector<T>::operator=(const Vector<T> & rhs)
{
    this->vec.resize(rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = std::complex<T>(rhs[i]);
    }
    domain = TIME_DOMAIN;
    return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::operator-()
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = -(this->vec[i]);
    }
    return *this;
}

template <class T>
template <class U>
ComplexVector<T> & ComplexVector<T>::operator+=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] += rhs.vec[i];
    }
    return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::operator+=(const std::complex<T> & rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] += rhs;
    }
    return *this;
}

/**
 * \brief Buffer addition operator.
 */
template <class T, class U>
inline ComplexVector<T> operator+(ComplexVector<T> lhs, const Vector<U>& rhs)
{
    lhs += rhs;
    return lhs;
}

/**
 * \brief Scalar addition operator.
 */
template <class T>
inline ComplexVector<T> operator+(ComplexVector<T> lhs, const std::complex<T> & rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
template <class U>
ComplexVector<T> & ComplexVector<T>::operator-=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] -= rhs.vec[i];
    }
    return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::operator-=(const std::complex<T> &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] -= rhs;
    }
    return *this;
}

/**
 * \brief Buffer subtraction operator.
 */
template <class T, class U>
inline ComplexVector<T> operator-(ComplexVector<T> lhs, const Vector<U>& rhs)
{
    lhs -= rhs;
    return lhs;
}

/**
 * \brief Scalar subtraction operator.
 */
template <class T>
inline ComplexVector<T> operator-(ComplexVector<T> lhs, const std::complex<T> & rhs)
{
    lhs -= rhs;
    return lhs;
}

template <class T>
template <class U>
ComplexVector<T> & ComplexVector<T>::operator*=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= rhs.vec[i];
    }
    return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::operator*=(const std::complex<T> &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= rhs;
    }
    return *this;
}

/**
 * \brief Buffer multiplication operator.
 */
template <class T, class U>
inline ComplexVector<T> operator*(ComplexVector<T> lhs, const Vector<U>& rhs)
{
    lhs *= rhs;
    return lhs;
}

/**
 * \brief Scalar multiplication operator.
 */
template <class T>
inline ComplexVector<T> operator*(ComplexVector<T> lhs, const std::complex<T> & rhs)
{
    lhs *= rhs;
    return lhs;
}

template <class T>
template <class U>
ComplexVector<T> & ComplexVector<T>::operator/=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] /= rhs.vec[i];
    }
    return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::operator/=(const std::complex<T> &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] /= rhs;
    }
    return *this;
}

/**
 * \brief Buffer division operator.
 */
template <class T, class U>
inline ComplexVector<T> operator/(ComplexVector<T> lhs, const Vector<U> & rhs)
{
    lhs /= rhs;
    return lhs;
}

/**
 * \brief Scalar division operator.
 */
template <class T>
inline ComplexVector<T> operator/(ComplexVector<T> lhs, const std::complex<T> & rhs)
{
    lhs /= rhs;
    return lhs;
}
 /*
template <class T>
Vector< std::complex<T> > & ComplexVector<T>::exp() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = std::exp(this->vec[i]);
    }
    return *this;
}
   */
template <class T>
ComplexVector<T> & ComplexVector<T>::pow(const std::complex<SLICKDSP_FLOAT_TYPE> & exponent) {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = std::pow(this->vec[i], exponent);
    }
    return *this;
}

/**
 * \brief Sets each element of "buffer" equal to its value to the power of "exponent".
 *
 * \param buffer Buffer to operate on.
 * \param exponent Exponent to use.
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & pow(ComplexVector<T> & buffer, const std::complex<SLICKDSP_FLOAT_TYPE> exponent) {
    return buffer.pow(exponent);
}

template <class T>
const std::complex<SLICKDSP_FLOAT_TYPE> ComplexVector<T>::mean() const {
    assert(this->size() > 0);
    std::complex<SLICKDSP_FLOAT_TYPE> sum = 0;
    for (unsigned i=0; i<this->size(); i++) {
        sum += this->vec[i];
    }
    return sum / ((SLICKDSP_FLOAT_TYPE) this->size());
}

template <class T>
const T ComplexVector<T>::max(unsigned *maxLoc) const {
    assert(this->size() > 0);
    T maxVal = this->vec[0].real();
    unsigned maxIndex = 0;

    for (unsigned i=1; i<this->size(); i++) {
        if (maxVal < this->vec[i].real()) {
            maxVal = this->vec[i].real();
            maxIndex = i;
        }
    }
    if (maxLoc != NULL) {
        *maxLoc = maxIndex;
    }
    return maxVal;
}

/**
 * \brief Returns the mean (average) of the data in "buffer".
 */
template <class T>
inline const std::complex<SLICKDSP_FLOAT_TYPE> mean(ComplexVector<T> & buffer) {
    return buffer.mean();
}

template <class T>
const SLICKDSP_FLOAT_TYPE ComplexVector<T>::var() const {
    assert(this->size() > 1);
    std::complex<SLICKDSP_FLOAT_TYPE> meanVal = this->mean();
    std::complex<SLICKDSP_FLOAT_TYPE> sum = 0;
    for (unsigned i=0; i<this->size(); i++) {
        std::complex<SLICKDSP_FLOAT_TYPE> varDiff = this->vec[i] - meanVal;
        sum += varDiff * std::conj(varDiff);
    }
    return sum.real() / (this->size() - 1);
}

/**
 * \brief Returns the variance of the data in "buffer".
 */
template <class T>
inline const SLICKDSP_FLOAT_TYPE var(ComplexVector<T> & buffer) {
    return buffer.var();
}

/**
 * \brief Returns the standard deviation of the data in "buffer".
 */
template <class T>
inline const SLICKDSP_FLOAT_TYPE stdDev(ComplexVector<T> & buffer) {
    return buffer.stdDev();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::saturate(const std::complex<T> & val) {
    for (unsigned i=0; i<this->size(); i++) {
        if (this->vec[i].real() > val.real())
            this->vec[i].real(val.real());
        else if (this->vec[i].real() < -val.real())
            this->vec[i].real(-val.real());
        if (this->vec[i].imag() > val.imag())
            this->vec[i].imag(val.imag());
        else if (this->vec[i].imag() < -val.imag())
            this->vec[i].imag(-val.imag());
    }
    return *this;
}

/**
 * \brief Sets the upper and lower limit of the values in "vector".
 *
 * \param vector Data to limit.
 * \param val Limiting value for the data in "vector".  Any values that
 *      are greater than "val" are made equal to "val", and
 *      any that are less than -val are made equal to -val.  This is done
 *      independently on the real and imaginary elements of "vector".
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & saturate(ComplexVector<T> & vector, const std::complex<T> & val) {
    return vector.saturate(val);
}
    
template <class T>
ComplexVector<T> & ComplexVector<T>::fft() {
    #ifdef NIMBLEDSP_DOMAIN_CHECKS
    assert(domain == TIME_DOMAIN);
    #endif
    
    kissfft<T> fftEngine = kissfft<T>(this->size(), false);
    std::vector< std::complex<T> > fftResults(this->size());
    
    fftEngine.transform((typename kissfft_utils::traits<T>::cpx_type *) VECTOR_TO_ARRAY(this->vec),
                        (typename kissfft_utils::traits<T>::cpx_type *) VECTOR_TO_ARRAY(fftResults));
    this->vec = fftResults;
    domain = FREQUENCY_DOMAIN;
    return *this;
}

/**
 * \brief Sets "buffer" equal to the FFT of the data in buffer.
 *
 * Sets \ref domain equal to NimbleDSP::FREQUENCY_DOMAIN.
 * \param buffer Buffer to operate on.
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & fft(ComplexVector<T> &buffer) {
    return buffer.fft();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::conj() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i].imag(-this->vec[i].imag());
    }
    return *this;
}

/**
 * \brief Conjugates the data in "buffer".
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & conj(ComplexVector<T> & buffer) {
    return buffer.conj();
}

/**
 * \brief Returns the squared magnitude of "val".
 */
template <class T>
inline T magSq(const std::complex<T> &val) {
    return val.real() * val.real() + val.imag() * val.imag();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::magSq() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i].real(NimbleDSP::magSq(this->vec[i]));
        this->vec[i].imag(0);
    }
    return *this;
}

/**
 * \brief Sets each element of "buffer" equal to its magnitude squared.
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & magSq(ComplexVector<T> & buffer) {
    return buffer.magSq();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::ifft() {
    #ifdef NIMBLEDSP_DOMAIN_CHECKS
    assert(domain == FREQUENCY_DOMAIN);
    #endif
    
    kissfft<T> fftEngine = kissfft<T>(this->size(), true);
    std::vector< std::complex<T> > fftResults(this->size());
    
    fftEngine.transform((typename kissfft_utils::traits<T>::cpx_type *) VECTOR_TO_ARRAY(this->vec),
                        (typename kissfft_utils::traits<T>::cpx_type *) VECTOR_TO_ARRAY(fftResults));
    this->vec = fftResults;
    domain = TIME_DOMAIN;
    return *this;
}

/**
 * \brief Sets "buffer" equal to the inverse FFT of the data in buffer.
 *
 * Sets \ref domain equal to NimbleDSP::TIME_DOMAIN.
 * \param buffer Buffer to operate on.
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & ifft(ComplexVector<T> &buffer) {
    return buffer.ifft();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::angle() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i].real(std::arg(this->vec[i]));
        this->vec[i].imag(0);
    }
    return *this;
}

/**
 * \brief Sets each element of "buffer" equal to its angle.
 *
 * The angle is held in the real portion of "buffer".
 * \return Reference to "buffer".
 */
template <class T>
inline ComplexVector<T> & angle(ComplexVector<T> & buffer) {
    return buffer.angle();
}

/**
 * \brief Returns the angle of a single std::complex value.
 *
 * \param val Value whose angle is returned.
 */
template <class T>
inline T angle(std::complex<T> &val) {
    return std::arg(val);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::abs() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (T) std::abs(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Changes the elements of \ref vec to their absolute value.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & abs(ComplexVector<T> & vector) {
    return vector.abs();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::exp() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (std::complex<T>) std::exp(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to e^(element).
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & exp(ComplexVector<T> & vector) {
    return vector.exp();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::log() {
    for (unsigned i=0; i<this->size(); i++) {
		this->vec[i] = (std::complex<T>) std::log(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to the natural log of the element.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & log(ComplexVector<T> & vector) {
    return vector.log();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::log10() {
    for (unsigned i=0; i<this->size(); i++) {
		this->vec[i] = (std::complex<T>) std::log10(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to the base 10 log of the element.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & log10(ComplexVector<T> & vector) {
    return vector.log10();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::rotate(int numToShift) {
    while (numToShift < 0)
        numToShift += this->size();
    
    while (numToShift >= (int) this->size())
        numToShift -= this->size();
    
    if (numToShift == 0)
        return *this;

    std::rotate(this->vec.begin(), this->vec.begin()+numToShift, this->vec.end());
    return *this;
}

/**
 * \brief Circular rotation.
 *
 * \param vector Buffer to rotate.
 * \param numToShift Number of positions to shift in the circular rotation.  numToShift
 *      can be positive or negative.  If you visualize the 0 index value at the left and
 *      the end of the array at the right, positive numToShift values shift the array to
 *      the left, and negative values shift it to the right.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & rotate(ComplexVector<T> & vector, int numToShift) {
    return vector.rotate(numToShift);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::reverse() {
    std::reverse(this->vec.begin(), this->vec.end());
    return *this;
}

/**
 * \brief Reverses the order of the elements in \ref vec.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & reverse(ComplexVector<T> & vector) {
    return vector.reverse();
}

/**
 * \brief Sets the length of \ref vec to "len".
 *
 * \param vector Buffer to operate on.
 * \param len The new length for \ref vec.  If len is longer than vec's current size, the
 *      new elements will be set to "val".  If len is less than vec's current size the extra
 *      elements will be cut off and the other elements will remain the same.
 * \param val The value to set any new elements to.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & resize(ComplexVector<T> & vector, int len, T val = 0) {
    return vector.resize(len, val);
}

/**
 * \brief Lengthens \ref vec by "len" elements.
 *
 * \param vector Buffer to operate on.
 * \param len The number of elements to add to \ref vec.
 * \param val The value to set the new elements to.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & pad(ComplexVector<T> & vector, int len, T val = 0) {
    return vector.pad(len, val);
}
    
template <class T>
ComplexVector<T> & ComplexVector<T>::upsample(int rate, int phase) {
	assert(rate > 0);
	assert(phase >= 0 && phase < rate);
	if (rate == 1)
		return *this;

	int originalSize = this->vec.size();
	this->vec.resize(originalSize*rate);
	int from, to;
	for (from = originalSize - 1, to = this->size() - (rate - phase); to > 0; from--, to -= rate) {
		this->vec[to] = this->vec[from];
		this->vec[from] = 0;
	}
	return *this;
}

/**
 * \brief Inserts rate-1 zeros between samples.
 *
 * \param vector Buffer to operate on.
 * \param rate Indicates how many zeros should be inserted between samples.
 * \param phase Indicates how many of the zeros should be before the samples (as opposed to
 *      after).  Valid values are 0 to "rate"-1.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & upsample(ComplexVector<T> & vector, int rate, int phase = 0) {
    return vector.upsample(rate, phase);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::downsample(int rate, int phase) {
	assert(rate > 0);
	assert(phase >= 0 && phase < rate);
	if (rate == 1)
		return *this;

	int newSize = this->size() / rate;
	int from, to;
	for (from = phase, to = 0; to < newSize; from += rate, to++) {
		this->vec[to] = this->vec[from];
	}
	this->vec.resize(newSize);
	return *this;
}

/**
 * \brief Removes rate-1 samples out of every rate samples.
 *
 * \param vector Buffer to operate on.
 * \param rate Indicates how many samples should be removed.
 * \param phase Tells the function which sample should be the first to be kept.  Valid values
 *      are 0 to "rate"-1.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & downsample(ComplexVector<T> & vector, int rate, int phase = 0) {
    return vector.downsample(rate, phase);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::cumsum(T initialVal) {
    T sum = initialVal;
    for (unsigned i=0; i<this->size(); i++) {
        sum += this->vec[i];
        this->vec[i] = sum;
    }
    return *this;
}

/**
 * \brief Replaces "vector" with the cumulative sum of the samples in "vector".
 *
 * \param vector Data to operate on.
 * \param initialVal Initializing value for the cumulative sum.  Defaults to zero.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & cumsum(ComplexVector<T> & vector, T initialVal = 0) {
    return vector.cumsum(initialVal);
}
    
template <class T>
ComplexVector<T> & ComplexVector<T>::diff() {
	assert(this->size() > 1);
	for (unsigned i=0; i<(this->size()-1); i++) {
		this->vec[i] = this->vec[i + 1] - this->vec[i];
	}
    this->resize(this->size()-1);
    return *this;
}

/**
 * \brief Replaces \ref vec with the difference between successive samples in vec.
 *
 * The resulting \ref vec is one element shorter than it was previously.
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & diff(ComplexVector<T> & vector) {
    return vector.diff();
}

template <class T>
ComplexVector<T> & ComplexVector<T>::diff(std::complex<T> & previousVal) {
	assert(this->size() > 0);
    std::complex<T> nextPreviousVal = this->vec[this->size()-1];
	for (unsigned i=this->size()-1; i>0; i--) {
		this->vec[i] = this->vec[i] - this->vec[i - 1];
	}
    this->vec[0] = this->vec[0] - previousVal;
    previousVal = nextPreviousVal;
    return *this;
}

/**
 * \brief Replaces \ref vec with the difference between successive samples in vec.
 *
 * \param vector Buffer to operate on.
 * \param previousVal The last value in the sample stream before the current contents
 *      of \ref vec.  previousVal allows the resulting vec to be the same size as the
 *      previous vec.
 * \return Reference to "vector".
 */
template <class T>
ComplexVector<T> & diff(ComplexVector<T> & vector, std::complex<T> & previousVal) {
    return vector.diff(previousVal);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::conv(ComplexVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Convolution function.
 *
 * \param data Buffer to operate on.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution, which is
 *      the length of "data" plus the length of "filter" - 1.  "True" tells the
 *      function to retain the size of "data" by trimming the tails at both ends of
 *      the convolution.
 * \return Reference to "data", which holds the result of the convolution.
 */
template <class T>
inline ComplexVector<T> & conv(ComplexVector<T> & data, ComplexVector<T> & filter, bool trimTails = false) {
    return filter.conv(data, trimTails);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::decimate(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Decimate function.
 *
 * This function is equivalent to filtering with the \ref conv function and downsampling
 * with the \ref downsample function, but much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the decimation.
 */
template <class T>
inline ComplexVector<T> & decimate(ComplexVector<T> & data, int rate, ComplexVector<T> & filter, bool trimTails = false) {
    return filter.decimate(data, rate, trimTails);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::interp(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    else {
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Interpolation function.
 *
 * This function is equivalent to upsampling followed by filtering, but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to upsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the interpolation.
 */
template <class T>
inline ComplexVector<T> & interp(ComplexVector<T> & data, int rate, ComplexVector<T> & filter, bool trimTails = false) {
    return filter.interp(data, rate, trimTails);
}

template <class T>
ComplexVector<T> & ComplexVector<T>::resample(ComplexVector<T> & data, int interpRate, int decimateRate,  bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        int interpLen = data.size() * interpRate;
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    else {
        int interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Resample function.
 *
 * This function is equivalent to upsampling by "interpRate", filtering, and downsampling
 *      by "decimateRate", but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param interpRate Indicates how much to upsample.
 * \param decimateRate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the resampling.
 */
template <class T>
inline ComplexVector<T> & resample(ComplexVector<T> & data, int interpRate, int decimateRate,
            ComplexVector<T> & filter, bool trimTails = false) {
    return filter.resample(data, interpRate, decimateRate, trimTails);
}

template <class T>
T ComplexVector<T>::tone(T freq, T sampleFreq, T phase, unsigned numSamples) {
    assert(sampleFreq > 0.0);
    
    if (numSamples && numSamples != this->size()) {
        this->resize(numSamples);
    }
    
    T phaseInc = (freq / sampleFreq) * 2 * M_PI;
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i].real(std::cos(phase));
        this->vec[i].imag(std::sin(phase));
        phase += phaseInc;
    }
    return phase;
}

/**
 * \brief Generates a complex tone.
 *
 * \param vec The vector to put the tone in.
 * \param freq The tone frequency.
 * \param sampleFreq The sample frequency.  Defaults to 1 Hz.
 * \param phase The tone's starting phase, in radians.  Defaults to 0.
 * \param numSamples The number of samples to generate.  "0" indicates to generate
 *      this->size() samples.  Defaults to 0.
 * \return Reference to "this".
 */
template <class T>
T tone(ComplexVector<T> & vec, T freq, T sampleFreq = 1.0, T phase = 0.0, unsigned numSamples = 0) {
    return vec.tone(freq, sampleFreq, phase, numSamples);
}

template <class T>
T ComplexVector<T>::modulate(T freq, T sampleFreq, T phase) {
    assert(sampleFreq > 0.0);
    
    T phaseInc = (freq / sampleFreq) * 2 * M_PI;
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= std::complex<T>(std::cos(phase), std::sin(phase));
        phase += phaseInc;
    }
    return phase;
}

/**
 * \brief Modulates the data with a complex sinusoid.
 *
 * \param freq The modulating tone frequency.
 * \param sampleFreq The sample frequency of the data.  Defaults to 1 Hz.
 * \param phase The modulating tone's starting phase, in radians.  Defaults to 0.
 * \return The next phase if the tone were to continue.
 */
template <class T>
T modulate(ComplexVector<T> &data, T freq, T sampleFreq, T phase) {
    return data.modulate(freq, sampleFreq, phase);
}

template <class T>
void ComplexVector<T>::slice(unsigned lower, unsigned upper, ComplexVector<T> &destination) {
	assert(lower <= upper);
	assert(upper < this->size());

	destination.resize(upper - lower + 1);
	for (unsigned from=lower, to=0; from<=upper; from++, to++) {
		destination[to] = this->vec[from];
	}
}

template <class T>
ComplexVector<T> & ComplexVector<T>::ceil() {
	for (int index=0; index<this->size(); index++) {
		this->vec[index].real(std::ceil(this->vec[index].real()));
		this->vec[index].imag(std::ceil(this->vec[index].imag()));
	}
	return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::floor() {
	for (int index=0; index<this->size(); index++) {
		this->vec[index].real(std::floor(this->vec[index].real()));
		this->vec[index].imag(std::floor(this->vec[index].imag()));
	}
	return *this;
}

template <class T>
ComplexVector<T> & ComplexVector<T>::round() {
	for (int index=0; index<this->size(); index++) {
		this->vec[index].real(std::round(this->vec[index].real()));
		this->vec[index].imag(std::round(this->vec[index].imag()));
	}
	return *this;
}

/**
 * \brief Class for complex FIR filters.
 */
template <class T>
class ComplexFirFilter : public ComplexVector<T> {
 protected:
    /**
     * \brief Saved data that is used for stream filtering.
     */
    std::vector<char> savedData;
    
    /**
     * \brief Indicates how many samples are in \ref savedData.  Used for stream filtering.
     */
    int numSavedSamples;
    
    /**
     * \brief Indicates the filter phase.  Used for stream filtering.
     */
    int phase;
    
 public:
    /**
     * \brief Determines how the filter should filter.
     *
     * NimbleDSP::ONE_SHOT_RETURN_ALL_RESULTS is equivalent to "trimTails = false" of the Vector convolution methods.
     * NimbleDSP::ONE_SHOT_TRIM_TAILS is equivalent to "trimTails = true" of the Vector convolution methods.
     * NimbleDSP::STREAMING maintains the filter state from call to call so it can produce results as if it had
     *      filtered one continuous set of data.
     */
    FilterOperationType filtOperation;

    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    ComplexFirFilter<T>(unsigned size = DEFAULT_BUF_LEN, FilterOperationType operation = STREAMING, std::vector< std::complex<T> > *scratch = NULL) : ComplexVector<T>(size, scratch)
            {if (size > 0) {savedData.resize((size - 1) * sizeof(std::complex<T>)); numSavedSamples = size - 1;}
             else {savedData.resize(0); numSavedSamples = 0;} phase = 0; filtOperation = operation;}
    
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    ComplexFirFilter<T>(std::vector<U> data, FilterOperationType operation = STREAMING, std::vector<T> *scratch = NULL) : ComplexVector<T>(data, scratch)
            {savedData.resize((data.size() - 1) * sizeof(std::complex<T>)); numSavedSamples = data.size() - 1; phase = 0; filtOperation = operation;}
    
    /**
     * \brief Array constructor.
     *
     * Sets buf equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref buf will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    ComplexFirFilter<T>(U *data, unsigned dataLen, FilterOperationType operation = STREAMING, std::vector< std::complex<T> > *scratch = NULL) : ComplexVector<T>(data, dataLen, NimbleDSP::TIME_DOMAIN, scratch)
            {savedData.resize((dataLen - 1) * sizeof(std::complex<T>)); numSavedSamples = dataLen - 1; phase = 0; filtOperation = operation;}
    
    /**
     * \brief Copy constructor.
     */
    ComplexFirFilter<T>(const ComplexFirFilter<T>& other) {this->vec = other.vec; savedData = other.savedData;
            numSavedSamples = other.numSavedSamples; phase = other.phase; filtOperation = other.filtOperation;}
    
    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Assignment operator.
     */
    ComplexFirFilter<T>& operator=(const Vector<T>& rhs) {this->vec = rhs.vec; savedData.resize(this->size() - 1); phase = 0; filtOperation = STREAMING; return *this;}
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Convolution method.
     *
     * \param data The buffer that will be filtered.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual ComplexVector<T> & conv(ComplexVector<T> & data, bool trimTails = false);
    
    /**
     * \brief Decimate method.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual ComplexVector<T> & decimate(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual ComplexVector<T> & interp(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Resample method.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual ComplexVector<T> & resample(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);

    /**
     * \brief Correlation method.
     *
     * \param data The buffer that will be correlated.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual ComplexVector<T> & corr(ComplexVector<T> & data);
};


template <class T>
ComplexVector<T> & ComplexFirFilter<T>::conv(ComplexVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING:
        dataTmp->resize((this->size() - 1) + data.size());
        for (int i=0; i<this->size()-1; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + this->size() - 1] = data[i];
        }
        
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        for (int i=0; i<this->size()-1; i++) {
            savedDataArray[i] = (*dataTmp)[i + data.size()];
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & ComplexFirFilter<T>::decimate(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((data.size() + numSavedSamples - (this->size() - 1) + rate - 1)/rate);
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        int nextResultDataPoint = resultIndex * rate;
        numSavedSamples = ((int) dataTmp->size()) - nextResultDataPoint;

        for (int i=0; i<numSavedSamples; i++) {
            savedDataArray[i] = (*dataTmp)[i + nextResultDataPoint];
        }
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & ComplexFirFilter<T>::interp(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + rate - 1) / rate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * rate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((unsigned) dataTmp->size() * rate);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
                if (dataTmp->size() - dataStart == numSavedSamples) {
                    keepGoing = false;
                    phase = filterStart;
                }
            }
        }
        data.resize(resultIndex);

        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & ComplexFirFilter<T>::resample(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    int interpLen, resampLen;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + interpRate - 1) / interpRate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * interpRate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        interpLen = (unsigned) dataTmp->size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
            if (dataTmp->size() - dataStart == numSavedSamples) {
                keepGoing = false;
                phase = filterStart;
            }
        }
        data.resize(resultIndex);
        
        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & ComplexFirFilter<T>::corr(ComplexVector<T> & data) {
	this->conj();
	this->reverse();
	this->conv(data);
	this->reverse();
	this->conj();
    return data;
}

/**
 * \brief Correlation function.
 *
 * \param data Buffer to operate on.
 * \param filter The filter that will correlate with "data".
 * \param trimTails "False" tells the function to return the entire convolution, which is
 *      the length of "data" plus the length of "filter" - 1.  "True" tells the
 *      function to retain the size of "data" by trimming the tails at both ends of
 *      the convolution.
 * \return Reference to "data", which holds the result of the convolution.
 */
template <class T>
inline ComplexVector<T> & corr(ComplexVector<T> & data, ComplexFirFilter<T> & filter) {
    return filter.corr(data);
}

/**
 * \brief Class for complex IIR filters.
 */
template <class T>
class ComplexIirFilter {
 protected:
    std::vector< std::complex<T> > state;

    /** 
     * \brief Initializes numerator and denominator with the size and contents of "num" and "den".
     *
     * \param array Array to set vec equal to.
     * \param arrayLen Number of elements in array.
     */
    template <class U>
    void initArray(U *num, unsigned numLen, U *den, unsigned denLen);
    
 public:
    std::vector< std::complex<T> > numerator;
    std::vector< std::complex<T> > denominator;
    
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    ComplexIirFilter<T>(unsigned order = 2) {numerator = std::vector< std::complex<T> >(order + 1);
            denominator = std::vector< std::complex<T> >(order + 1); state = std::vector< std::complex<T> >(order + 1);}
    
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    ComplexIirFilter<T>(std::vector< std::complex<U> > & num, std::vector< std::complex<U> > & den) {
            initArray(VECTOR_TO_ARRAY(num), (unsigned) num.size(), VECTOR_TO_ARRAY(den), (unsigned) den.size());}
    
    /**
     * \brief Array constructor.
     *
     * Sets vec equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref vec will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    ComplexIirFilter<T>(U *num, unsigned numLen, U *den, unsigned denLen) {initArray(num, numLen, den, denLen);}
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Convolution method.
     *
     * \param data The vector that will be filtered.
     * \param trimTails "False" tells the method to return the entire convolution, which is
     *      the length of "data" plus the length of "this" (the filter) - 1.  "True" tells the
     *      method to retain the size of "data" by trimming the tails at both ends of
     *      the convolution.
     * \return Reference to "data", which holds the result of the convolution.
     */
    template <class U>
    ComplexVector<U> & filter(ComplexVector<U> & data);
    
};


template <class T>
template <class U>
void ComplexIirFilter<T>::initArray(U *num, unsigned numLen, U *den, unsigned denLen) {
    numerator = std::vector< std::complex<T> >(numLen);
    for (unsigned i=0; i<numLen; i++) {
        numerator[i] = (std::complex<T>) num[i];
    }
    
    denominator = std::vector< std::complex<T> >(denLen);
    for (unsigned i=0; i<denLen; i++) {
        denominator[i] = (std::complex<T>) den[i];
    }
    
    state = std::vector< std::complex<T> >(std::max(numLen, denLen));
}

template <class T>
template <class U>
ComplexVector<U> & ComplexIirFilter<T>::filter(ComplexVector<U> & data) {
    unsigned resultIndex, i;
    std::complex<U> newState0;
    
    for (resultIndex=0; resultIndex<data.size(); resultIndex++) {
        newState0 = data[resultIndex];
        
        // Update the state
        for (i=(unsigned) state.size()-1; i>=denominator.size(); i--) {
            state[i] = state[i - 1];
        }
        // Update the state and apply the feedback
        for (i=(unsigned) denominator.size()-1; i>0; i--) {
            state[i] = state[i - 1];
            newState0 -= denominator[i] * state[i];
        }
        state[0] = newState0;

        // Calculate the output
        data[resultIndex] = 0;
        for (i=0; i<numerator.size(); i++) {
            data[resultIndex] += numerator[i] * state[i];
        }
    }
    return data;
}

template <class T, class U>
ComplexVector<U> & filter(ComplexVector<U> & data, ComplexIirFilter<T> & filt) {
    return filt.filter(data);
}

#define M_2PI  6.28318530717958647692

static int NFCNS, NGRID;
static double DEV;
static double *FX, *WTX;
static bool converged;


static bool ParksMcClellan2(double *FirCoeff, int NFILT, int JTYPE, int NBANDS, double *EDGE, double *FX, double *WTX, int LGRID = 16);
static double  EFF(double FREQ, int LBAND, int JTYPE);
static double WATE(double FREQ, int LBAND, int JTYPE);
static double D(std::vector<double> &X, int K, int N, int M);
static double GEE(std::vector<double> &X, std::vector<double> &GRID, std::vector<double> &AD, std::vector<double> &Y, int K, int N);
static void Remez(std::vector<int> &IEXT, std::vector<double> &AD, std::vector<double> &ALPHA, std::vector<double> &X, std::vector<double> &Y, std::vector<double> &H, std::vector<double> &DES, std::vector<double> &GRID, std::vector<double> &WT, std::vector<double> &A, std::vector<double> &P, std::vector<double> &Q);


/* Input Values
  NFILT   FILTER LENGTH
  JTYPE   TYPE OF FILTER 1 = MULTIPLE PASSBAND/STOPBAND FILTER 2 = DIFFERENTIATOR 3 = HILBERT TRANSFORM FILTER
  NBANDS  NUMBER OF BANDS
  LGRID   GRID DENSITY, WILL BE SET TO 16 UNLESS SPECIFIED OTHERWISE BY A POSITIVE CONSTANT.
  EDGE[2*NBANDS] BANDEDGE ARRAY, LOWER AND UPPER EDGES FOR EACH BAND WITH A MAXIMUM OF 10 BANDS.
  FX[NBANDS]     DESIRED FUNCTION ARRAY (OR DESIRED SLOPE IF A DIFFERENTIATOR) FOR EACH BAND.
  WTX[NBANDS]    WEIGHT FUNCTION ARRAY IN EACH BAND.  FOR A DIFFERENTIATOR, THE WEIGHT FUNCTION IS INVERSELY PROPORTIONAL TO F.
*/
//---------------------------------------------------------------------------

bool ParksMcClellan2(double *FirCoeff, int NFILT, int JTYPE, int NBANDS, double *EDGE, double *fx, double *wtx, int LGRID)
{
 int J=0, L=0, NEG=0, NODD=0, LBAND=0;
 int NM1=0, NZ=0;
 double DELF=0.0, FUP=0.0, TEMP=0.0, CHANGE=0.0, XT=0.0;
 
 converged = true;
 
 FX = fx;
 WTX = wtx;
 
 int smallArraySize = (NFILT + 7) / 2;
 int bigArraySize = smallArraySize * 16;
    
 std::vector<int> IEXT(smallArraySize);
 std::vector<double> AD(smallArraySize), ALPHA(smallArraySize), X(smallArraySize), Y(smallArraySize), H(smallArraySize);
 std::vector<double> DES(bigArraySize), GRID(bigArraySize), WT(bigArraySize);
 std::vector<double> A(smallArraySize), P(smallArraySize), Q(smallArraySize);

 if(JTYPE == 1) NEG = 0;   // Lowpass, Bandpass, Highpass, and Notch
 else NEG = 1;             // Hilberts and Differentiators
 NODD = NFILT % 2;
 NFCNS = NFILT/2;
 if(NODD == 1 && NEG == 0) NFCNS = NFCNS + 1;

// From here on down, the code is essentially a simple Fortran to C conversion.
//  SET UP THE DENSE GRID.  THE NUMBER OF POINTS IN THE GRID IS (FILTER LENGTH + 1)*GRID DENSITY/2

	  GRID[1] = EDGE[1];
	  DELF = LGRID * NFCNS;
	  DELF = 0.5 / DELF;
	  if(NEG == 0) goto L135;
	  if(EDGE[1] < DELF) GRID[1] = DELF;
L135: J = 1;
	  L = 1;
	  LBAND = 1;
L140: FUP = EDGE[L+1];
L145: TEMP = GRID[J];

// CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT FUNCTION ON THE GRID

	  DES[J] = EFF(TEMP,LBAND,JTYPE);
	  WT[J] = WATE(TEMP,LBAND,JTYPE);
	  J=J+1;
	  GRID[J] = TEMP + DELF;
	  if(GRID[J] > FUP) goto L150;
	  goto L145;
L150: GRID[J-1] = FUP;
	  DES[J-1] = EFF(FUP,LBAND,JTYPE);
	  WT[J-1] = WATE(FUP,LBAND,JTYPE);
	  LBAND = LBAND + 1;
	  L = L + 2;
	  if(LBAND > NBANDS) goto L160;
	  GRID[J] = EDGE[L];
	  goto L140;
L160: NGRID = J-1;
	  if(NEG != NODD) goto L165;
	  if(GRID[NGRID] > (0.5-DELF)) NGRID = NGRID-1;


//SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT TO THE ORIGINAL PROBLEM

L165: if(NEG <= 0)
	   {
		if(NODD == 1) goto L200;
		for(J=1; J<=NGRID; J++)
		 {
		  CHANGE = cos(M_PI * GRID[J] );
		  DES[J] = DES[J] / CHANGE;
		  WT[J] = WT[J] * CHANGE;
		 }
		goto L200;
	   }
	  if(NODD == 1) goto L190;
	  for(J=1; J<=NGRID; J++)
	   {
		CHANGE = sin(M_PI * GRID[J]);
		DES[J] = DES[J] / CHANGE;
		WT[J] = WT[J] * CHANGE;
	   }
	  goto L200;
L190: for(J=1; J<=NGRID; J++)
	   {
		CHANGE = sin( M_2PI * GRID[J] );
		DES[J] = DES[J] / CHANGE;
		WT[J] = WT[J] * CHANGE;
	   }

// INITIAL GUESS FOR THE EXTREMAL FREQUENCIES--EQUALLY SPACED ALONG THE GRID

L200: TEMP = (double)(NGRID-1)/(double)NFCNS;
	  for(J=1; J<=NFCNS; J++)
	   {
		XT = J-1;
		IEXT[J] = XT * TEMP + 1.0;
	   }
	  IEXT[NFCNS+1] = NGRID;
	  NM1 = NFCNS-1;
	  NZ = NFCNS+1;

	  // CALL THE REMEZ EXCHANGE ALGORITHM TO DO THE APPROXIMATION PROBLEM
	  Remez(IEXT, AD, ALPHA, X, Y, H, DES, GRID, WT, A, P, Q);

	  // CALCULATE THE IMPULSE RESPONSE.
	  if(NEG > 0)goto L320;  // NEG = 1 for Hilberts and Differentiators

	  if(NODD == 0) goto L310;
	  for(J=1; J<=NM1; J++) // NEG=0 and NODD=0 for bandpass even length filter nfcns=nfilt/2
	   {
		H[J] = 0.5 * ALPHA[NZ-J];
	   }
	  H[NFCNS] = ALPHA[1];
	  goto L350;

L310: H[1] = 0.25 * ALPHA[NFCNS]; // NEG=0 and NODD=1 for bandpass odd length filter nfcns=nfilt/2 + 1
	  for(J=2; J<=NM1; J++)
	   {
		H[J] = 0.25 * (ALPHA[NZ-J] + ALPHA[NFCNS+2-J]);
	   }
	  H[NFCNS] = 0.5 * ALPHA[1] + 0.25 * ALPHA[2];
	  goto L350;


L320: if(NODD == 0) goto L330;  // NEG=1 and NODD=0 for even length Hilberts and Differentiators
	  H[1] = 0.25 * ALPHA[NFCNS];
	  H[2] = 0.25 * ALPHA[NM1];
	  for(J=3; J<=NM1; J++)
	   {
		H[J] = 0.25 * (ALPHA[NZ-J] - ALPHA[NFCNS+3-J] );
	   }
	  H[NFCNS] = 0.5 * ALPHA[1] - 0.25 * ALPHA[3];
	  H[NZ] =0.0;
	  goto L350;
L330: H[1] = 0.25 * ALPHA[NFCNS]; // NEG=1 and NODD=1 for odd length Hilberts and Differentiators
	  for(J=2; J<=NM1; J++)
	   {
		H[J] = 0.25 * (ALPHA[NZ-J] - ALPHA[NFCNS+2-J] );
	   }
	  H[NFCNS] = 0.5 * ALPHA[1] - 0.25 * ALPHA[2];




L350: // Output section. Code for the Hilberts and Differentiators was omitted.
      // Note that the Fortran arrays started at 1, but my FirCoeff array starts at 0.
	  // Write the first half of the response
	  for(J=1; J<=NFCNS; J++)
	   {
		FirCoeff[J-1] = H[J];
	   }
	  // NEG=0 and NODD=0 for bandpass even length filter nfcns=nfilt/2
	  if(NEG == 0 && NODD == 0)
	  for(J=1; J<=NFCNS; J++)
	   {
		FirCoeff[NFCNS+J-1] = H[NFCNS-J+1];
	   }
	  // NEG=0 and NODD=1 for bandpass odd length filter nfcns=nfilt/2 + 1
	  if(NEG == 0 && NODD == 1)
	  for(J=1; J<NFCNS; J++)
	   {
		FirCoeff[NFCNS+J-1] = H[NFCNS-J];
	   }
    return converged;
}

/*
FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE AS A FUNCTION OF FREQUENCY. AN ARBITRARY
FUNCTION OF FREQUENCY CAN BE APPROXIMATED if THE USER REPLACES THIS FUNCTION WITH THE
APPROPRIATE CODE TO EVALUATE THE IDEAL MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
*/


double EFF(double FREQ, int LBAND, int JTYPE)
{
 if(JTYPE == 2)  return( FX[LBAND] * FREQ );
 else return( FX[LBAND] );
}


//FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION OF FREQUENCY.  SIMILAR TO THE FUNCTION
//EFF, THIS FUNCTION CAN BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY DESIRED WEIGHTING FUNCTION.

double WATE(double FREQ, int LBAND, int JTYPE)
{
 if(JTYPE == 1 || JTYPE == 3) return(WTX[LBAND]); // JTYPE=1 Bandpass JTYPE=3 for Hilberts
 if(FX[LBAND] < 0.0001) return(WTX[LBAND]);       // JTYPE=2 for Differentiators
 return(WTX[LBAND] / FREQ);
}


/*
THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM FOR THE WEIGHTED CHEBYSHEV
APPROXIMATION OF A CONTINUOUS FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE DESIRED FUNCTION ON THIS GRID,
THE WEIGHT FUNCTION ON THE GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE EXTREMAL
FREQUENCIES. THE PROGRAM MINIMIZES THE CHEBYSHEV ERROR BY DETERMINING THE BEST LOCATION OF
THE EXTREMAL FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES THE COEFFICIENTS OF
THE BEST APPROXIMATION.
*/

void Remez(std::vector<int> &IEXT, std::vector<double> &AD, std::vector<double> &ALPHA, std::vector<double> &X, std::vector<double> &Y, std::vector<double> &H, std::vector<double> &DES, std::vector<double> &GRID, std::vector<double> &WT, std::vector<double> &A, std::vector<double> &P, std::vector<double> &Q)
{
 int J=0, ITRMAX=0, NZ=0, NZZ=0, JET=0, K=0, L=0, NU=0, JCHNGE=0, K1=0, KNZ=0, KLOW=0, NUT=0, KUP=0;
 int NUT1=0, LUCK=0, KN=0, NM1=0, KKK=0, JM1=0, JP1=0, NITER=0;
 double DNUM=0.0, DDEN=0.0, DTEMP=0.0, FT=0.0, XT=0.0, XT1=0.0, XE=0.0;
 double FSH=0.0, GTEMP=0.0, CN=0.0, DELF=0.0, AA=0.0, BB=0.0;  // SciPy declares CN as an int, which is probably inconsequential the way CN is used.
 double DEVL=0.0, COMP=0.0, YNZ=0.0, Y1 = 0.0, ERR=0.0;


	  ITRMAX = 40;      // MAXIMUM NUMBER OF ITERATIONS This was 25. It was changed by Parks in a later file.
	  DEVL = -1.0;
	  NZ = NFCNS + 1;
	  NZZ = NFCNS + 2;
	  NITER = 0;
L100: IEXT[NZZ] = NGRID + 1;
	  NITER = NITER + 1;
	  if(NITER > ITRMAX) goto L400; // L400 is where the coeff calc starts
	  for(J=1; J<=NZ; J++)
	   {
		DTEMP = GRID[ IEXT[J] ];
		X[J] = cos(DTEMP * M_2PI);
	   }

	  JET = (NFCNS-1)/15 + 1;
	  for(J=1; J<=NZ; J++)
	   {
		AD[J] = D(X,J,NZ,JET);
	   }

	  DNUM = 0.0;
	  DDEN = 0.0;
	  K = 1;
	  for(J=1; J<=NZ; J++)
	   {
		L = IEXT[J];
		DNUM += AD[J] * DES[L];
		DDEN += (double)K * AD[J]/WT[L];
		K = -K;
	   }
	  DEV = DNUM / DDEN;

//      WRITE(IOUT,131) DEV
//  131 FORMAT(1X,12HDEVIATION = ,F12.9)

	  NU = 1;
	  if(DEV > 0.0) NU = -1;
	  DEV = -(double)NU * DEV;
	  K = NU;
	  for(J=1; J<=NZ; J++)
	   {
		L = IEXT[J];
		DTEMP = (double)K * DEV/WT[L];
		Y[J] = DES[L] + DTEMP;
		K = -K;
	   }
	  if(DEV <= DEVL)
	   {
		printf("Failed to converge\n");
        converged = false;
		goto L400;
	   }
	  DEVL = DEV;
	  JCHNGE = 0;
	  K1 = IEXT[1];
	  KNZ = IEXT[NZ];
	  KLOW = 0;
	  NUT = -NU;
	  J = 1;

//SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION

L200: if(J == NZZ) YNZ = COMP;
	  if(J >= NZZ) goto L300;
	  KUP = IEXT[J+1];
	  L = IEXT[J] + 1;
	  NUT = -NUT;
	  if(J == 2) Y1 = COMP;
	  COMP = DEV;
	  if(L >= KUP) goto L220;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L220;
	  COMP = (double)NUT * ERR;
L210: L = L + 1;
	  if(L >= KUP) goto L215;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L215;
	  COMP = (double)NUT * ERR;
	  goto L210;
L215: IEXT[J] = L-1;
	  J = J + 1;
	  KLOW = L - 1;
	  JCHNGE = JCHNGE+1;
	  goto L200;
L220: L = L - 1;
L225: L = L - 1;
	  if(L <= KLOW) goto L250;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP > 0.0) goto L230;
	  if(JCHNGE <= 0) goto L225;
	  goto L260;
L230: COMP = (double)NUT * ERR;
L235: L = L - 1;
	  if(L <= KLOW) goto L240;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L240;
	  COMP = (double)NUT * ERR;
	  goto L235;
L240: KLOW = IEXT[J];
	  IEXT[J] = L + 1;
	  J = J + 1;
	  JCHNGE = JCHNGE + 1;
	  goto L200;
L250: L = IEXT[J] + 1;
	  if(JCHNGE > 0) goto L215;
L255: L = L + 1;
	  if(L >= KUP) goto L260;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = ( ERR - DES[L] ) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L255;
	  COMP = (double)NUT * ERR;
	  goto L210;
L260: KLOW = IEXT[J];
	  J = J + 1;
	  goto L200;

L300: if(J > NZZ) goto L320;
	  if(K1 > IEXT[1]) K1 = IEXT[1];
	  if(KNZ < IEXT[NZ]) KNZ = IEXT[NZ];
	  NUT1 = NUT;
	  NUT = -NU;
	  L = 0 ;
	  KUP = K1;
	  COMP = YNZ * (1.00001);
	  LUCK = 1;
L310: L = L + 1;
	  if(L >= KUP) goto L315;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L310;
	  COMP = (double)NUT * ERR;
	  J = NZZ;
	  goto L210;
L315: LUCK = 6;
	  goto L325;
L320: if(LUCK > 9) goto L350;
	  if(COMP > Y1) Y1 = COMP;
	  K1 = IEXT[NZZ];
L325: L = NGRID + 1;
	  KLOW = KNZ;
	  NUT = -NUT1;
	  COMP = Y1 * 1.00001;
L330: L = L-1;
	  if(L <= KLOW) goto L340;
	  ERR = GEE(X,GRID,AD,Y,L,NZ);
	  ERR = (ERR - DES[L]) * WT[L];
	  DTEMP = (double)NUT * ERR - COMP;
	  if(DTEMP <= 0.0) goto L330;
	  J = NZZ;
	  COMP = (double)NUT * ERR;
	  LUCK = LUCK + 10;
	  goto L235;
L340: if(LUCK == 6) goto L370;
	  for(J=1; J<=NFCNS; J++)
	   {
		IEXT[NZZ - J] = IEXT[NZ - J];
	   }
	  IEXT[1] = K1;
	  goto L100;
L350: KN = IEXT[NZZ];
	  for(J=1; J<=NFCNS; J++)
	   {
		IEXT[J] = IEXT[J+1];
	   }
	  IEXT[NZ] = KN;
	  goto L100;
L370: if(JCHNGE > 0) goto L100;

//CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION USING THE INVERSE DISCRETE FOURIER TRANSFORM

L400: NM1 = NFCNS - 1;
	  FSH = 1.0E-06;
	  GTEMP = GRID[1];
	  X[NZZ] = -2.0;
	  CN = 2 * NFCNS - 1;
	  DELF = 1.0/CN;
	  L = 1;
	  KKK = 0;

	  // Different than the SciPy code. You may find that that this condition always tests the same way.
	  if(GRID[1] < 0.01 && GRID[NGRID] > 0.49) KKK = 1;
	  // if (edge[1] == 0.0 && edge[2*nbands] == 0.5) kkk = 1;  // the SciPy code

	  if(NFCNS <= 3) KKK = 1;
	  if(KKK == 1) goto L405;
	  DTEMP = cos(M_2PI * GRID[1]);
	  DNUM = cos(M_2PI * GRID[NGRID] );
	  AA = 2.0/(DTEMP-DNUM);
	  BB = -(DTEMP+DNUM)/(DTEMP-DNUM);
L405: for(J=1; J<=NFCNS; J++)
	   {
		FT = J-1;
		FT = FT * DELF;
		XT = cos(M_2PI * FT);
		if(KKK == 1) goto L410;
		XT = (XT-BB)/AA;
		XT1 = sqrt(1.0 - XT*XT);
		FT = atan2(XT1,XT)/M_2PI;   // Different than the SciPy code
		// ft = acos(xt)/TWOPI;     // the SciPy code which was embedded in a #if .. #endif  It appears to be debug that wasn't removed.
L410:   XE = X[L];
		if(XT > XE) goto L420;
		if((XE-XT) < FSH) goto L415;
		L = L + 1;
		goto L410;
L415:   A[J] = Y[L];
		goto L425;
L420:   if((XT-XE) < FSH) goto L415;
		GRID[1] = FT;
		A[J] = GEE(X,GRID,AD,Y,1,NZ);
L425:   if(L > 1) L = L-1;
	   }
	  GRID[1] = GTEMP;
	  DDEN = M_2PI / CN;
	  for(J=1; J<=NFCNS; J++)
	   {
		DTEMP = 0.0;
		DNUM = (double)(J-1) * DDEN;
		//if(NM1 >= 1) This orig piece of code is redundant.
		for(K=1; K<=NM1; K++)
		 {
		  DTEMP += A[K+1] * cos(DNUM * (double)K);
		 }
		DTEMP = 2.0 * DTEMP + A[1];
		ALPHA[J] = DTEMP;
	   }

	  for(J=2; J<=NFCNS; J++)
	   {
		ALPHA[J] = 2.0 * ALPHA[J]/CN;
	   }
	  ALPHA[1] = ALPHA[1]/CN;
	  if(KKK == 1) goto L545;
	  P[1] = 2.0 * ALPHA[NFCNS]*BB + ALPHA[NM1];
	  P[2] = 2.0 * AA * ALPHA[NFCNS];
	  Q[1] = ALPHA[NFCNS-2] - ALPHA[NFCNS];
	  for(J=2; J<=NM1; J++)
	   {
		if(J < NM1) goto L515;
		AA = 0.5 * AA;
		BB = 0.5 * BB;
L515:   P[J+1] = 0.0;
		for(K=1; K<=J; K++)
		 {
		  A[K] = P[K];
		  P[K] = 2.0 * BB * A[K];
		 }
		P[2] += A[1] * 2.0 * AA;
		JM1 = J-1;
		for(K=1; K<=JM1; K++)
		 {
		  P[K] += Q[K] + AA * A[K+1];
		 }
		JP1 = J + 1;
		for(K=3; K<=JP1; K++)
		 {
		  P[K] += AA * A[K-1];
		 }
		if(J == NM1) goto L540;
		for(K=1; K<=J; K++)
		 {
		  Q[K] = -A[K];
		 }
		Q[1] += ALPHA[NFCNS-1-J];
	  }
L540: for(J=1; J<=NFCNS; J++)
	   {
		ALPHA[J] = P[J];
	   }
L545: if(NFCNS <= 3)
	   {
		ALPHA[NFCNS+1] = 0.0;
		ALPHA[NFCNS+2] = 0.0;
	   }
}

//-----------------------------------------------------------------------
// FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION COEFFICIENTS FOR USE IN THE FUNCTION GEE.
double D(std::vector<double> &X, int K, int N, int M)
{
 int J, L;
 double Dee, Q;
 Dee = 1.0;
 Q = X[K];
 for(L=1; L<=M; L++)
 for(J=L; J<=N; J+=M)
  {
   if(J != K)Dee = 2.0 * Dee * (Q - X[J]);
  }
 return(1.0/Dee); // Dee can and will go to zero.
}

//-----------------------------------------------------------------------
// FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE LAGRANGE INTERPOLATION FORMULA
// IN THE BARYCENTRIC FORM
double GEE(std::vector<double> &X, std::vector<double> &GRID, std::vector<double> &AD, std::vector<double> &Y, int K, int N)
{
 int j;
 double P,C,D,XF;
 P = 0.0;
 XF = GRID[K];
 XF = cos(M_2PI * XF);
 D = 0.0;
 for(j=1; j<=N; j++)
  {
   C = XF - X[j];
   C = AD[j] / C;   // C can and will go to zero.
   D = D + C;
   P = P + C*Y[j];
  }
 return(P/D);  // D can and will go to zero.
}

/*
Note:
In the 3 places noted above regarding divide by zero.
Using code similar to this helps this algorithm converge.
 #define MIN_TEST_VAL 1.0E-6
 if(fabs(C) < MIN_TEST_VAL )
  {
   if(C < 0.0)C = -MIN_TEST_VAL;
   else       C =  MIN_TEST_VAL;
  }
*/


/**
 * \brief Vector class for real numbers.
 */
template <class T>
class RealVector : public Vector<T> {
 public:
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    RealVector<T>(unsigned size = DEFAULT_BUF_LEN, std::vector<T> *scratch = NULL) : Vector<T>(size, scratch) {}
    
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealVector<T>(std::vector<U> data, std::vector<T> *scratch = NULL) : Vector<T>(data, scratch) {}
    
    /**
     * \brief Array constructor.
     *
     * Sets buf equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref buf will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealVector<T>(U *data, unsigned dataLen, std::vector<T> *scratch = NULL) : Vector<T>(data, dataLen, scratch) {}
    
    /**
     * \brief Copy constructor.
     */
    RealVector<T>(const RealVector<T>& other) {this->vec = other.vec;}
    
    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Assignment operator.
     */
    RealVector<T>& operator=(const Vector<T>& rhs) {this->vec = rhs.vec; return *this;}
    
    /**
     * \brief Unary minus (negation) operator.
     */
    RealVector<T> & operator-();
    
    /**
     * \brief Add Buffer/Assignment operator.
     */
    template <class U>
    RealVector<T> & operator+=(const Vector<U> &rhs);
    
    /**
     * \brief Add Scalar/Assignment operator.
     */
    RealVector<T> & operator+=(const T &rhs);
    
    /**
     * \brief Subtract Buffer/Assignment operator.
     */
    template <class U>
    RealVector<T> & operator-=(const Vector<U> &rhs);
    
    /**
     * \brief Subtract Scalar/Assignment operator.
     */
    RealVector<T> & operator-=(const T &rhs);
    
    /**
     * \brief Multiply Buffer/Assignment operator.
     */
    template <class U>
    RealVector<T> & operator*=(const Vector<U> &rhs);
    
    /**
     * \brief Multiply Scalar/Assignment operator.
     */
    RealVector<T> & operator*=(const T &rhs);

    /**
     * \brief Divide Buffer/Assignment operator.
     */
    template <class U>
    RealVector<T> & operator/=(const Vector<U> &rhs);
    
    /**
     * \brief Divide Scalar/Assignment operator.
     */
    RealVector<T> & operator/=(const T &rhs);
    
    /*****************************************************************************************
                                             Methods
    *****************************************************************************************/
    /**
     * \brief Sets each element of \ref buf equal to its value to the power of "exponent".
     *
     * \param exponent Exponent to use.
     * \return Reference to "this".
     */
    RealVector<T> & pow(const SLICKDSP_FLOAT_TYPE exponent);
    
    /**
     * \brief Returns the mean (average) of the data in \ref buf.
     */
    const SLICKDSP_FLOAT_TYPE mean() const;
    
    /**
     * \brief Returns the variance of the data in \ref buf.
     */
    const SLICKDSP_FLOAT_TYPE var() const;
    
    /**
     * \brief Returns the standard deviation of the data in \ref buf.
     */
    const SLICKDSP_FLOAT_TYPE stdDev() const {return std::sqrt(this->var());}
    
    /**
     * \brief Returns the median element of \ref buf.
     */
    const T median();
    
    /**
     * \brief Returns the maximum element in \ref buf.
     *
     * \param maxLoc If it isn't equal to NULL the index of the maximum element
     *      will be returned via this pointer.  If more than one element is equal
     *      to the maximum value the index of the first will be returned.
     *      Defaults to NULL.
     */
    const T max(unsigned *maxLoc = NULL) const;
    
    /**
     * \brief Returns the minimum element in \ref buf.
     *
     * \param minLoc If it isn't equal to NULL the index of the minimum element
     *      will be returned via this pointer.  If more than one element is equal
     *      to the minimum value the index of the first will be returned.
     *      Defaults to NULL.
     */
    const T min(unsigned *minLoc = NULL) const;
    
    /**
     * \brief Sets the upper and lower limit of the values in \ref buf.
     *
     * \param val Limiting value for the data in \ref buf.  Any values that
     *      are greater than "val" are made equal to "val", and
     *      any that are less than -val are made equal to -val.
     * \return Reference to "this".
     */
    RealVector<T> & saturate(T val);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    RealVector<T> & ceil(void);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    RealVector<T> & floor(void);

    /**
     * \brief Does a "ceil" operation on \ref buf.
     * \return Reference to "this".
     */
    RealVector<T> & round(void);

    /**
     * \brief Convolution method for complex data.
     *
     * \param data The buffer that will be filtered.
     * \param trimTails "False" tells the method to return the entire convolution, which is
     *      the length of "data" plus the length of "this" (the filter) - 1.  "True" tells the
     *      method to retain the size of "data" by trimming the tails at both ends of
     *      the convolution.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual ComplexVector<T> & convComplex(ComplexVector<T> & data, bool trimTails);
    
    /**
     * \brief Decimate method for complex data.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual ComplexVector<T> & decimateComplex(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method for complex data.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual ComplexVector<T> & interpComplex(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Resample method for complex data.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual ComplexVector<T> & resampleComplex(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);
    
    /**
     * \brief Changes the elements of \ref vec to their absolute value.
     *
     * \return Reference to "this".
     */
    RealVector<T> & abs();
    
    /**
     * \brief Sets each element of \ref vec to e^(element).
     *
     * \return Reference to "this".
     */
    RealVector<T> & exp();
    
    /**
     * \brief Sets each element of \ref vec to the natural log of the element.
     *
     * \return Reference to "this".
     */
    RealVector<T> & log();
    
    /**
     * \brief Sets each element of \ref vec to the base 10 log of the element.
     *
     * \return Reference to "this".
     */
    RealVector<T> & log10();

    /**
     * \brief Circular rotation.
     *
     * \param numToShift Number of positions to shift in the circular rotation.  numToShift
     *      can be positive or negative.  If you visualize the 0 index value at the left and
     *      the end of the array at the right, positive numToShift values shift the array to
     *      the left, and negative values shift it to the right.
     * \return Reference to "this".
     */
    RealVector<T> & rotate(int numToShift);
    
    /**
     * \brief Reverses the order of the elements in \ref vec.
     *
     * \return Reference to "this".
     */
    RealVector<T> & reverse();

    /**
     * \brief Sets the length of \ref vec to "len".
     *
     * \param len The new length for \ref vec.  If len is longer than vec's current size, the
     *      new elements will be set to "val".  If len is less than vec's current size the extra
     *      elements will be cut off and the other elements will remain the same.
     * \param val The value to set any new elements to.  Defaults to 0.
     * \return Reference to "this".
     */
    RealVector<T> & resize(unsigned len, T val = (T) 0) {this->vec.resize(len, val); return *this;}
    
    /**
     * \brief Lengthens \ref vec by "len" elements.
     *
     * \param len The number of elements to add to \ref vec.
     * \param val The value to set the new elements to.  Defaults to 0.
     * \return Reference to "this".
     */
    RealVector<T> & pad(unsigned len, T val = (T) 0) {this->vec.resize(this->size()+len, val); return *this;}
    
    /**
     * \brief Inserts rate-1 zeros between samples.
     *
     * \param rate Indicates how many zeros should be inserted between samples.
     * \param phase Indicates how many of the zeros should be before the samples (as opposed to
     *      after).  Valid values are 0 to "rate"-1.  Defaults to 0.
     * \return Reference to "this".
     */
    RealVector<T> & upsample(int rate, int phase = 0);
    
    /**
     * \brief Removes rate-1 samples out of every rate samples.
     *
     * \param rate Indicates how many samples should be removed.
     * \param phase Tells the method which sample should be the first to be kept.  Valid values
     *      are 0 to "rate"-1.  Defaults to 0.
     * \return Reference to "this".
     */
    RealVector<T> & downsample(int rate, int phase = 0);
    
    /**
     * \brief Replaces \ref vec with the cumulative sum of the samples in \ref vec.
     *
     * \param initialVal Initializing value for the cumulative sum.  Defaults to zero.
     * \return Reference to "this".
     */
	RealVector<T> & cumsum(T initialVal = 0);
    
    /**
     * \brief Replaces \ref vec with the difference between successive samples in vec.
     *
     * The resulting \ref vec is one element shorter than it was previously.
     * \return Reference to "this".
     */
	RealVector<T> & diff();
    
    /**
     * \brief Replaces \ref vec with the difference between successive samples in vec.
     *
     * \param previousVal The last value in the sample stream before the current contents
     *      of \ref vec.  previousVal allows the resulting vec to be the same size as the
     *      previous vec.
     * \return Reference to "this".
     */
    RealVector<T> & diff(T & previousVal);
    
    /**
     * \brief Convolution method.
     *
     * \param data The vector that will be filtered.
     * \param trimTails "False" tells the method to return the entire convolution, which is
     *      the length of "data" plus the length of "this" (the filter) - 1.  "True" tells the
     *      method to retain the size of "data" by trimming the tails at both ends of
     *      the convolution.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual RealVector<T> & conv(RealVector<T> & data, bool trimTails = false);
    
    /**
     * \brief Decimate method.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual RealVector<T> & decimate(RealVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual RealVector<T> & interp(RealVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Resample method.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The vector that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails "False" tells the method to return the entire convolution.  "True"
     *      tells the method to retain the size of "data" by trimming the tails at both
     *      ends of the convolution.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual RealVector<T> & resample(RealVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);
    
    /**
     * \brief Generates a complex tone.
     *
     * \param freq The tone frequency.
     * \param sampleFreq The sample frequency.  Defaults to 1 Hz.
     * \param phase The tone's starting phase, in radians.  Defaults to 0.
     * \param numSamples The number of samples to generate.  "0" indicates to generate
     *      this->size() samples.  Defaults to 0.
     * \return The next phase if the tone were to continue.
     */
    T tone(T freq, T sampleFreq = 1.0, T phase = 0.0, unsigned numSamples = 0);
    
    /**
     * \brief Modulates the data with a real sinusoid.
     *
     * \param freq The modulating tone frequency.
     * \param sampleFreq The sample frequency of the data.  Defaults to 1 Hz.
     * \param phase The modulating tone's starting phase, in radians.  Defaults to 0.
     * \return The next phase if the tone were to continue.
     */
    T modulate(T freq, T sampleFreq = 1.0, T phase = 0.0);
};

template <class T>
RealVector<T> & RealVector<T>::pow(const SLICKDSP_FLOAT_TYPE exponent) {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (T) std::pow(this->vec[i], exponent);
    }
    return *this;
}

/**
 * \brief Sets each element of "buffer" equal to its value to the power of "exponent".
 *
 * \param buffer The buffer to operate on.
 * \param exponent Exponent to use.
 * \return Reference to "buffer".
 */
template <class T>
RealVector<T> & pow(RealVector<T> & buffer, const SLICKDSP_FLOAT_TYPE exponent) {
    return buffer.pow(exponent);
}

template <class T>
const SLICKDSP_FLOAT_TYPE RealVector<T>::mean() const {
    assert(this->size() > 0);
    SLICKDSP_FLOAT_TYPE sum = 0;
    for (unsigned i=0; i<this->size(); i++) {
        sum += this->vec[i];
    }
    return sum / this->size();
}

/**
 * \brief Returns the mean (average) of the data in "buffer".
 * \param buffer The buffer to operate on.
 */
template <class T>
const SLICKDSP_FLOAT_TYPE mean(RealVector<T> & buffer) {
    return buffer.mean();
}

template <class T>
const SLICKDSP_FLOAT_TYPE RealVector<T>::var() const {
    assert(this->size() > 1);
    SLICKDSP_FLOAT_TYPE meanVal = mean();
    SLICKDSP_FLOAT_TYPE sum = 0;
    for (unsigned i=0; i<this->size(); i++) {
        SLICKDSP_FLOAT_TYPE varDiff = ((SLICKDSP_FLOAT_TYPE) this->vec[i]) - meanVal;
        sum += varDiff * varDiff;
    }
    return sum / (this->size() - 1);
}

/**
 * \brief Returns the variance of the data in "buffer".
 * \param buffer The buffer to operate on.
 */
template <class T>
const SLICKDSP_FLOAT_TYPE var(RealVector<T> & buffer) {
    return buffer.var();
}

/**
 * \brief Returns the standard deviation of the data in "buffer".
 * \param buffer The buffer to operate on.
 */
template <class T>
const SLICKDSP_FLOAT_TYPE stdDev(RealVector<T> & buffer) {
    return buffer.stdDev();
}

template <class T>
const T RealVector<T>::median() {
    assert(this->size() > 0);
    std::vector<T> scratchBuf = this->vec;
    std::sort(scratchBuf.begin(), scratchBuf.end());
    if (this->size() & 1) {
        // Odd number of samples
        return scratchBuf[this->size()/2];
    }
    else {
        // Even number of samples.  Average the two in the middle.
        unsigned topHalfIndex = this->size()/2;
        return (scratchBuf[topHalfIndex] + scratchBuf[topHalfIndex-1]) / ((T) 2);
    }
}

/**
 * \brief Returns the median element of "buffer".
 * \param buffer The buffer to operate on.
 */
template <class T>
const T median(RealVector<T> & buffer) {
    return buffer.median();
}

template <class T>
const T RealVector<T>::max(unsigned *maxLoc) const {
    assert(this->size() > 0);
    T maxVal = this->vec[0];
    unsigned maxIndex = 0;
    
    for (unsigned i=1; i<this->size(); i++) {
        //if (buf[i] > maxVal) {
        if (maxVal < this->vec[i]) {
            maxVal = this->vec[i];
            maxIndex = i;
        }
    }
    if (maxLoc != NULL) {
        *maxLoc = maxIndex;
    }
    return maxVal;
}

/**
 * \brief Returns the maximum element in "buffer".
 *
 * \param buffer The buffer to search.
 * \param maxLoc If it isn't equal to NULL the index of the maximum element
 *      will be returned via this pointer.  If more than one element is equal
 *      to the maximum value the index of the first will be returned.
 *      Defaults to NULL.
 */
template <class T>
const T max(RealVector<T> & buffer, unsigned *maxLoc = NULL) {
    return buffer.max(maxLoc);
}

template <class T>
const T RealVector<T>::min(unsigned *minLoc) const {
    assert(this->size() > 0);
    T minVal = this->vec[0];
    unsigned minIndex = 0;
    
    for (unsigned i=1; i<this->size(); i++) {
        if (this->vec[i] < minVal) {
            minVal = this->vec[i];
            minIndex = i;
        }
    }
    if (minLoc != NULL) {
        *minLoc = minIndex;
    }
    return minVal;
}

/**
 * \brief Returns the minimum element in "buffer".
 *
 * \param buffer The buffer to search.
 * \param minLoc If it isn't equal to NULL the index of the minimum element
 *      will be returned via this pointer.  If more than one element is equal
 *      to the minimum value the index of the first will be returned.
 *      Defaults to NULL.
 */
template <class T>
const T min(RealVector<T> & buffer, unsigned *minLoc = NULL) {
    return buffer.min(minLoc);
}

template <class T>
RealVector<T> & RealVector<T>::saturate(T val) {
    for (unsigned i=0; i<this->size(); i++) {
        if (this->vec[i] > val)
            this->vec[i] = val;
        else if (this->vec[i] < -val)
            this->vec[i] = -val;
    }
    return *this;
}

/**
 * \brief Sets the upper and lower limit of the values in "buffer".
 *
 * \param buffer The buffer to operate on.
 * \param val Limiting value for the data in \ref buf.  Any values that
 *      are greater than "val" are made equal to "val", and
 *      any that are less than -val are made equal to -val.
 * \return Reference to "buffer".
 */
template <class T>
RealVector<T> & saturate(RealVector<T> & buffer, T val) {
    return buffer.saturate(val);
}

template <class T>
ComplexVector<T> & RealVector<T>::convComplex(ComplexVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Convolution function.
 *
 * \param data Buffer to operate on.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution, which is
 *      the length of "data" plus the length of "filter" - 1.  "True" tells the
 *      function to retain the size of "data" by trimming the tails at both ends of
 *      the convolution.
 * \return Reference to "data", which holds the result of the convolution.
 */
template <class T>
inline ComplexVector<T> & conv(ComplexVector<T> & data, RealVector<T> & filter, bool trimTails = false) {
    return filter.convComplex(data, trimTails);
}

template <class T>
ComplexVector<T> & RealVector<T>::decimateComplex(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Decimate function.
 *
 * This function is equivalent to filtering with the \ref conv function and downsampling
 * with the \ref downsample function, but much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the decimation.
 */
template <class T>
inline ComplexVector<T> & decimate(ComplexVector<T> & data, int rate, RealVector<T> & filter, bool trimTails = false) {
    return filter.decimateComplex(data, rate, trimTails);
}

template <class T>
ComplexVector<T> & RealVector<T>::interpComplex(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    else {
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Interpolation function.
 *
 * This function is equivalent to upsampling followed by filtering, but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to upsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the interpolation.
 */
template <class T>
inline ComplexVector<T> & interp(ComplexVector<T> & data, int rate, RealVector<T> & filter, bool trimTails = false) {
    return filter.interpComplex(data, rate, trimTails);
}

template <class T>
ComplexVector<T> & RealVector<T>::resampleComplex(ComplexVector<T> & data, int interpRate, int decimateRate,  bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        int interpLen = data.size() * interpRate;
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    else {
        int interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Resample function.
 *
 * This function is equivalent to upsampling by "interpRate", filtering, and downsampling
 *      by "decimateRate", but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param interpRate Indicates how much to upsample.
 * \param decimateRate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the resampling.
 */
template <class T>
inline ComplexVector<T> & resample(ComplexVector<T> & data, int interpRate, int decimateRate,
            RealVector<T> & filter, bool trimTails = false) {
    return filter.resampleComplex(data, interpRate, decimateRate, trimTails);
}

template <class T>
RealVector<T> & RealVector<T>::abs() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (T) std::abs(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Changes the elements of \ref vec to their absolute value.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & abs(RealVector<T> & vector) {
    return vector.abs();
}

template <class T>
RealVector<T> & RealVector<T>::exp() {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (T) std::exp(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to e^(element).
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & exp(RealVector<T> & vector) {
    return vector.exp();
}

template <class T>
RealVector<T> & RealVector<T>::log() {
    for (unsigned i=0; i<this->size(); i++) {
		this->vec[i] = (T) std::log(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to the natural log of the element.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & log(RealVector<T> & vector) {
    return vector.log();
}

template <class T>
RealVector<T> & RealVector<T>::log10() {
    for (unsigned i=0; i<this->size(); i++) {
		this->vec[i] = (T) std::log10(this->vec[i]);
    }
    return *this;
}

/**
 * \brief Sets each element of \ref vec to the base 10 log of the element.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & log10(RealVector<T> & vector) {
    return vector.log10();
}

template <class T>
RealVector<T> & RealVector<T>::rotate(int numToShift) {
    while (numToShift < 0)
        numToShift += this->size();
    
    while (numToShift >= (int) this->size())
        numToShift -= this->size();
    
    if (numToShift == 0)
        return *this;

    std::rotate(this->vec.begin(), this->vec.begin()+numToShift, this->vec.end());
    return *this;
}

/**
 * \brief Circular rotation.
 *
 * \param vector Buffer to rotate.
 * \param numToShift Number of positions to shift in the circular rotation.  numToShift
 *      can be positive or negative.  If you visualize the 0 index value at the left and
 *      the end of the array at the right, positive numToShift values shift the array to
 *      the left, and negative values shift it to the right.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & rotate(RealVector<T> & vector, int numToShift) {
    return vector.rotate(numToShift);
}

template <class T>
RealVector<T> & RealVector<T>::reverse() {
    std::reverse(this->vec.begin(), this->vec.end());
    return *this;
}

/**
 * \brief Reverses the order of the elements in \ref vec.
 *
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & reverse(RealVector<T> & vector) {
    return vector.reverse();
}

/**
 * \brief Sets the length of \ref vec to "len".
 *
 * \param vector Buffer to operate on.
 * \param len The new length for \ref vec.  If len is longer than vec's current size, the
 *      new elements will be set to "val".  If len is less than vec's current size the extra
 *      elements will be cut off and the other elements will remain the same.
 * \param val The value to set any new elements to.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & resize(RealVector<T> & vector, int len, T val = 0) {
    return vector.resize(len, val);
}

/**
 * \brief Lengthens \ref vec by "len" elements.
 *
 * \param vector Buffer to operate on.
 * \param len The number of elements to add to \ref vec.
 * \param val The value to set the new elements to.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & pad(RealVector<T> & vector, int len, T val = 0) {
    return vector.pad(len, val);
}
    
template <class T>
RealVector<T> & RealVector<T>::upsample(int rate, int phase) {
	assert(rate > 0);
	assert(phase >= 0 && phase < rate);
	if (rate == 1)
		return *this;

	int originalSize = this->size();
	this->vec.resize(originalSize*rate);
	int from, to;
	for (from = originalSize - 1, to = this->size() - (rate - phase); to > 0; from--, to -= rate) {
		this->vec[to] = this->vec[from];
		this->vec[from] = 0;
	}
	return *this;
}

/**
 * \brief Inserts rate-1 zeros between samples.
 *
 * \param vector Buffer to operate on.
 * \param rate Indicates how many zeros should be inserted between samples.
 * \param phase Indicates how many of the zeros should be before the samples (as opposed to
 *      after).  Valid values are 0 to "rate"-1.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & upsample(RealVector<T> & vector, int rate, int phase = 0) {
    return vector.upsample(rate, phase);
}

template <class T>
RealVector<T> & RealVector<T>::downsample(int rate, int phase) {
	assert(rate > 0);
	assert(phase >= 0 && phase < rate);
	if (rate == 1)
		return *this;

	int newSize = this->size() / rate;
	int from, to;
	for (from = phase, to = 0; to < newSize; from += rate, to++) {
		this->vec[to] = this->vec[from];
	}
	this->vec.resize(newSize);
	return *this;
}

/**
 * \brief Removes rate-1 samples out of every rate samples.
 *
 * \param vector Buffer to operate on.
 * \param rate Indicates how many samples should be removed.
 * \param phase Tells the function which sample should be the first to be kept.  Valid values
 *      are 0 to "rate"-1.  Defaults to 0.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & downsample(RealVector<T> & vector, int rate, int phase = 0) {
    return vector.downsample(rate, phase);
}

template <class T>
RealVector<T> & RealVector<T>::cumsum(T initialVal) {
    T sum = initialVal;
    for (unsigned i=0; i<this->size(); i++) {
        sum += this->vec[i];
        this->vec[i] = sum;
    }
    return *this;
}

/**
 * \brief Replaces "vector" with the cumulative sum of the samples in "vector".
 *
 * \param vector Data to operate on.
 * \param initialVal Initializing value for the cumulative sum.  Defaults to zero.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & cumsum(RealVector<T> & vector, T initialVal = 0) {
    return vector.cumsum(initialVal);
}
    
template <class T>
RealVector<T> & RealVector<T>::diff() {
	assert(this->size() > 1);
	for (unsigned i=0; i<(this->size()-1); i++) {
		this->vec[i] = this->vec[i + 1] - this->vec[i];
	}
    this->resize(this->size()-1);
    return *this;
}

/**
 * \brief Replaces \ref vec with the difference between successive samples in vec.
 *
 * The resulting \ref vec is one element shorter than it was previously.
 * \param vector Buffer to operate on.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & diff(RealVector<T> & vector) {
    return vector.diff();
}

template <class T>
RealVector<T> & RealVector<T>::diff(T & previousVal) {
	assert(this->size() > 0);
    T nextPreviousVal = this->vec[this->size()-1];
	for (unsigned i=this->size()-1; i>0; i--) {
		this->vec[i] = this->vec[i] - this->vec[i - 1];
	}
    this->vec[0] = this->vec[0] - previousVal;
    previousVal = nextPreviousVal;
    return *this;
}

/**
 * \brief Replaces \ref vec with the difference between successive samples in vec.
 *
 * \param vector Buffer to operate on.
 * \param previousVal The last value in the sample stream before the current contents
 *      of \ref vec.  previousVal allows the resulting vec to be the same size as the
 *      previous vec.
 * \return Reference to "vector".
 */
template <class T>
RealVector<T> & diff(RealVector<T> & vector, T & previousVal) {
    return vector.diff(previousVal);
}

template <class T>
RealVector<T> & RealVector<T>::conv(RealVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Convolution function.
 *
 * \param data Buffer to operate on.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution, which is
 *      the length of "data" plus the length of "filter" - 1.  "True" tells the
 *      function to retain the size of "data" by trimming the tails at both ends of
 *      the convolution.
 * \return Reference to "data", which holds the result of the convolution.
 */
template <class T>
inline RealVector<T> & conv(RealVector<T> & data, RealVector<T> & filter, bool trimTails = false) {
    return filter.conv(data, trimTails);
}

template <class T>
RealVector<T> & RealVector<T>::decimate(RealVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    else {
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
    }
    return data;
}

/**
 * \brief Decimate function.
 *
 * This function is equivalent to filtering with the \ref conv function and downsampling
 * with the \ref downsample function, but much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the decimation.
 */
template <class T>
inline RealVector<T> & decimate(RealVector<T> & data, int rate, RealVector<T> & filter, bool trimTails = false) {
    return filter.decimate(data, rate, trimTails);
}

template <class T>
RealVector<T> & RealVector<T>::interp(RealVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    else {
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Interpolation function.
 *
 * This function is equivalent to upsampling followed by filtering, but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param rate Indicates how much to upsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the interpolation.
 */
template <class T>
inline RealVector<T> & interp(RealVector<T> & data, int rate, RealVector<T> & filter, bool trimTails = false) {
    return filter.interp(data, rate, trimTails);
}

template <class T>
RealVector<T> & RealVector<T>::resample(RealVector<T> & data, int interpRate, int decimateRate,  bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }
    *dataTmp = data.vec;
    
    if (trimTails) {
        int interpLen = data.size() * interpRate;
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    else {
        int interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        int resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
    }
    return data;
}

/**
 * \brief Resample function.
 *
 * This function is equivalent to upsampling by "interpRate", filtering, and downsampling
 *      by "decimateRate", but is much more efficient.
 *
 * \param data Buffer to operate on.
 * \param interpRate Indicates how much to upsample.
 * \param decimateRate Indicates how much to downsample.
 * \param filter The filter that will convolve "data".
 * \param trimTails "False" tells the function to return the entire convolution.  "True"
 *      tells the function to retain the size of "data" by trimming the tails at both
 *      ends of the convolution.
 * \return Reference to "data", which holds the result of the resampling.
 */
template <class T>
inline RealVector<T> & resample(RealVector<T> & data, int interpRate, int decimateRate,
            RealVector<T> & filter, bool trimTails = false) {
    return filter.resample(data, interpRate, decimateRate, trimTails);
}

template <class T>
RealVector<T> & RealVector<T>::operator-()
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = -this->vec[i];
    }
    return *this;
}

template <class T>
template <class U>
RealVector<T> & RealVector<T>::operator+=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] += rhs.vec[i];
    }
    return *this;
}

template <class T>
RealVector<T> & RealVector<T>::operator+=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] += rhs;
    }
    return *this;
}

template <class T, class U>
inline RealVector<T> operator+(RealVector<T> lhs, const Vector<U>& rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
inline RealVector<T> operator+(RealVector<T> lhs, const T& rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
template <class U>
RealVector<T> & RealVector<T>::operator-=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] -= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealVector<T> & RealVector<T>::operator-=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] -= rhs;
    }
    return *this;
}

template <class T, class U>
inline RealVector<T> operator-(RealVector<T> lhs, const Vector<U>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template <class T>
inline RealVector<T> operator-(RealVector<T> lhs, const T& rhs)
{
    lhs -= rhs;
    return lhs;
}

template <class T>
template <class U>
RealVector<T> & RealVector<T>::operator*=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealVector<T> & RealVector<T>::operator*=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= rhs;
    }
    return *this;
}

template <class T, class U>
inline RealVector<T> operator*(RealVector<T> lhs, const Vector<U>& rhs)
{
    lhs *= rhs;
    return lhs;
}

template <class T>
inline RealVector<T> operator*(RealVector<T> lhs, const T& rhs)
{
    lhs *= rhs;
    return lhs;
}

template <class T>
template <class U>
RealVector<T> & RealVector<T>::operator/=(const Vector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] /= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealVector<T> & RealVector<T>::operator/=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] /= rhs;
    }
    return *this;
}

template <class T, class U>
inline RealVector<T> operator/(RealVector<T> lhs, const Vector<U>& rhs)
{
    lhs /= rhs;
    return lhs;
}

template <class T>
inline RealVector<T> operator/(RealVector<T> lhs, const T& rhs)
{
    lhs /= rhs;
    return lhs;
}

template <class T>
T RealVector<T>::tone(T freq, T sampleFreq, T phase, unsigned numSamples) {
    assert(sampleFreq > 0.0);
    
    if (numSamples && numSamples != this->size()) {
        this->resize(numSamples);
    }
    
    T phaseInc = (freq / sampleFreq) * 2 * M_PI;
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = std::sin(phase);
        phase += phaseInc;
    }
    return phase;
}

/**
 * \brief Generates a complex tone.
 *
 * \param vec The vector to put the tone in.
 * \param freq The tone frequency.
 * \param sampleFreq The sample frequency.  Defaults to 1 Hz.
 * \param phase The tone's starting phase, in radians.  Defaults to 0.
 * \param numSamples The number of samples to generate.  "0" indicates to generate
 *      this->size() samples.  Defaults to 0.
 * \return Reference to "this".
 */
template <class T>
T tone(RealVector<T> & vec, T freq, T sampleFreq = 1.0, T phase = 0.0, unsigned numSamples = 0) {
    return vec.tone(freq, sampleFreq, phase, numSamples);
}

template <class T>
T RealVector<T>::modulate(T freq, T sampleFreq, T phase) {
    assert(sampleFreq > 0.0);
    
    T phaseInc = (freq / sampleFreq) * 2 * M_PI;
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] *= std::sin(phase);
        phase += phaseInc;
    }
    return phase;
}

/**
 * \brief Modulates the data with a real sinusoid.
 *
 * \param freq The modulating tone frequency.
 * \param sampleFreq The sample frequency of the data.  Defaults to 1 Hz.
 * \param phase The modulating tone's starting phase, in radians.  Defaults to 0.
 * \return The next phase if the tone were to continue.
 */
 template <class T>
T modulate(RealVector<T> &data, T freq, T sampleFreq, T phase) {
    return data.modulate(freq, sampleFreq, phase);
}

 template <class T>
 RealVector<T> & RealVector<T>::ceil() {
 	for (int index=0; index<this->size(); index++) {
 		this->vec[index] = std::ceil(this->vec[index]);
 	}
 	return *this;
 }

 template <class T>
 RealVector<T> & RealVector<T>::floor() {
 	for (int index=0; index<this->size(); index++) {
 		this->vec[index] = std::floor(this->vec[index]);
 	}
 	return *this;
 }

 template <class T>
 RealVector<T> & RealVector<T>::round() {
 	for (int index=0; index<this->size(); index++) {
 		this->vec[index] = std::round(this->vec[index]);
 	}
 	return *this;
 }


/**
 * \brief Class for real FIR filters.
 */
template <class T>
class RealFirFilter : public RealVector<T> {
 protected:
    /**
     * \brief Saved data that is used for stream filtering.
     */
    std::vector<char> savedData;
    
    /**
     * \brief Indicates how many samples are in \ref savedData.  Used for stream filtering.
     */
    int numSavedSamples;
    
    /**
     * \brief Indicates the filter phase.  Used for stream filtering.
     */
    int phase;
    
    /**
     * \brief Applies a Hamming window on the current contents of "this".
     */
    void hamming(void);
    
 public:
    /**
     * \brief Determines how the filter should filter.
     *
     * NimbleDSP::ONE_SHOT_RETURN_ALL_RESULTS is equivalent to "trimTails = false" of the Vector convolution methods.
     * NimbleDSP::ONE_SHOT_TRIM_TAILS is equivalent to "trimTails = true" of the Vector convolution methods.
     * NimbleDSP::STREAMING maintains the filter state from call to call so it can produce results as if it had
     *      filtered one continuous set of data.
     */
    FilterOperationType filtOperation;

    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    RealFirFilter<T>(unsigned size = DEFAULT_BUF_LEN, FilterOperationType operation = STREAMING, std::vector<T> *scratch = NULL) : RealVector<T>(size, scratch)
            {if (size > 0) {savedData.resize((size - 1) * sizeof(std::complex<T>)); numSavedSamples = size - 1;}
             else {savedData.resize(0); numSavedSamples = 0;} phase = 0; filtOperation = operation;}
    
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealFirFilter<T>(std::vector<U> data, FilterOperationType operation = STREAMING, std::vector<T> *scratch = NULL) : RealVector<T>(data, scratch)
            {savedData.resize((data.size() - 1) * sizeof(std::complex<T>)); numSavedSamples = data.size() - 1; phase = 0; filtOperation = operation;}
    
    /**
     * \brief Array constructor.
     *
     * Sets buf equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref buf will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealFirFilter<T>(U *data, unsigned dataLen, FilterOperationType operation = STREAMING, std::vector<T> *scratch = NULL) : RealVector<T>(data, dataLen, scratch)
            {savedData.resize((dataLen - 1) * sizeof(std::complex<T>)); numSavedSamples = dataLen - 1; phase = 0; filtOperation = operation;}
    
    /**
     * \brief Copy constructor.
     */
    RealFirFilter<T>(const RealFirFilter<T>& other) {this->vec = other.vec; savedData = other.savedData;
            numSavedSamples = other.numSavedSamples; phase = other.phase; filtOperation = other.filtOperation;}
    
    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Assignment operator.
     */
    RealFirFilter<T>& operator=(const Vector<T>& rhs) {this->vec = rhs.vec; savedData.resize(this->size() - 1); phase = 0; filtOperation = STREAMING; return *this;}
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Convolution method.
     *
     * \param data The buffer that will be filtered.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual RealVector<T> & conv(RealVector<T> & data, bool trimTails = false);
    
    /**
     * \brief Convolution method for complex data.
     *
     * \param data The buffer that will be filtered.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual ComplexVector<T> & convComplex(ComplexVector<T> & data, bool trimTails = false);
    
    /**
     * \brief Decimate method.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual RealVector<T> & decimate(RealVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Decimate method for complex data.
     *
     * This method is equivalent to filtering with the \ref conv method and downsampling
     * with the \ref downsample method, but much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the decimation.
     */
    virtual ComplexVector<T> & decimateComplex(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual RealVector<T> & interp(RealVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Interpolation method for complex data.
     *
     * This method is equivalent to upsampling followed by filtering, but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param rate Indicates how much to upsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the interpolation.
     */
    virtual ComplexVector<T> & interpComplex(ComplexVector<T> & data, int rate, bool trimTails = false);
    
    /**
     * \brief Resample method.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual RealVector<T> & resample(RealVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);
    
    /**
     * \brief Resample method for complex data.
     *
     * This method is equivalent to upsampling by "interpRate", filtering, and downsampling
     *      by "decimateRate", but is much more efficient.
     *
     * \param data The buffer that will be filtered.
     * \param interpRate Indicates how much to upsample.
     * \param decimateRate Indicates how much to downsample.
     * \param trimTails This parameter is ignored.  The operation of the filter is determined by how
     *      \ref filtOperation is set.
     * \return Reference to "data", which holds the result of the resampling.
     */
    virtual ComplexVector<T> & resampleComplex(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails = false);
    
    /**
     * \brief Parks-McClellan algorithm for generating equiripple FIR filter coefficients.
     *
     * Application of the remez algorithm to producing equiripple FIR filter coefficients.  It has issues with
     * convergence.  It can generally converge up to 128 taps- more than that and it gets iffy.
     *
     * The PM algorithm implementation is a somewhat modified version of Iowa Hills Software's port of the
     * PM algorithm from the original Fortran to C.  Much appreciation to them for their work.
     *
     * \param filterOrder Indicates that the number of taps should be filterOrder + 1.
     * \param numBands The number of pass and stop bands.  Maximum of 10 bands.
     * \param freqPoints Pairs of points specify the boundaries of the bands, thus the length of this array
     *          must be 2 * numBands.  The frequencies in between the bands are the transition bands.  "0"
     *          corresponds to 0 Hz, and "1" corresponds to the Nyquist frequency.  Example: "0.0 0.3 .5 1.0"
     *          specifies two bands, one from 0 Hz to .3 * Nyquist and the other from 0.5 Nyquist to Nyquist.
     *          .3 * Nyquist to .5 * Nyquist is the transition band.
     * \param desiredBandResponse Indicates what the desired amplitude response is in the corresponding band.
     *          This array must have "numBands" elements.
     * \param weight Indicates how much weight should be given the performance of the filter in that band.  If
     *          all of the elements are "1" then they will have equal weight which will produce a true
     *          equiripple filter (same amount of ripple in the pass and stop bands).  If the stop bands are
     *          assigned more weight than the passbands then the attenuation in the stop bands will be increased
     *          at the expense of more ripple in the pass bands.
     * \param lGrid Grid density.  Defaults to 16.  This value should generally not be set lower than 16.
     *          Setting it higher than 16 can produce a filter with a better fit to the desired response at
     *          the cost of increased computations.
     * \return Boolean that indicates whether the filter converged or not.
     */
    bool firpm(int filterOrder, int numBands, double *freqPoints, double *desiredBandResponse,
                double *weight, int lGrid = 16);
    
    /**
     * \brief Generates a filter that can delay signals by an arbitrary sub-sample time.
     *
     * \param numTaps Number of taps to use in the filter.  If the only purpose of the filter is to delay
     *          the signal then the number of taps is not crucial.  If it is doing actual stop-band filtering
     *          too then the number of taps is important.
     * \param bandwidth The approximate bandwidth of the filter, normalized to the range 0.0 to 1.0, where
     *          1.0 is the Nyquist frequency.  The bandwidth must be greater than 0 and less than 1.
     * \param delay Amount of sample time to delay.  For example, a delay value of 0.1 would indicate to delay
     *          by one tenth of a sample.  "delay" can be positive or negative.
     */
    void fractionalDelayFilter(int numTaps, double bandwidth, double delay);
    
    /**
     * \brief Correlation method.
     *
     * \param data The buffer that will be correlated.
     * \return Reference to "data", which holds the result of the convolution.
     */
    virtual RealVector<T> & corr(RealVector<T> & data);
    
    /**
     * \brief Generates a Hamming window.
     *
     * \param len The window length.
     */
     void hamming(unsigned len);
    
    /**
     * \brief Generates a Hann window.
     *
     * \param len The window length.
     */
     void hann(unsigned len);
    
    /**
     * \brief Generates a generalized Hamming window.
     *
     * \param len The window length.
     * \param alpha
     * \param beta
     */
     void generalizedHamming(unsigned len, double alpha, double beta);
    
    /**
     * \brief Generates a Blackman window.
     *
     * \param len The window length.
     */
     void blackman(unsigned len);
    
    /**
     * \brief Generates a Blackman-harris window.
     *
     * \param len The window length.
     */
     void blackmanHarris(unsigned len);
    
};


template <class T>
RealVector<T> & RealFirFilter<T>::conv(RealVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    T *savedDataArray = (T *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING:
        dataTmp->resize((this->size() - 1) + data.size());
        for (int i=0; i<this->size()-1; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + this->size() - 1] = data[i];
        }
        
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        for (int i=0; i<this->size()-1; i++) {
            savedDataArray[i] = (*dataTmp)[i + data.size()];
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & RealFirFilter<T>::convComplex(ComplexVector<T> & data, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING:
        dataTmp->resize((this->size() - 1) + data.size());
        for (int i=0; i<this->size()-1; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + this->size() - 1] = data[i];
        }
        
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        for (int i=0; i<this->size()-1; i++) {
            savedDataArray[i] = (*dataTmp)[i + data.size()];
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() + this->size() - 1);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<((int)this->size()-1) - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<(int)dataTmp->size() - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
RealVector<T> & RealFirFilter<T>::decimate(RealVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    T *savedDataArray = (T *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((data.size() + numSavedSamples - (this->size() - 1) + rate - 1)/rate);
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        int nextResultDataPoint = resultIndex * rate;
        numSavedSamples = (unsigned) dataTmp->size() - nextResultDataPoint;

        for (int i=0; i<numSavedSamples; i++) {
            savedDataArray[i] = (*dataTmp)[i + nextResultDataPoint];
        }
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & RealFirFilter<T>::decimateComplex(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((data.size() + numSavedSamples - (this->size() - 1) + rate - 1)/rate);
        for (resultIndex=0; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate, filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        int nextResultDataPoint = resultIndex * rate;
        numSavedSamples = ((int) dataTmp->size()) - nextResultDataPoint;

        for (int i=0; i<numSavedSamples; i++) {
            savedDataArray[i] = (*dataTmp)[i + nextResultDataPoint];
        }
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(((data.size() + this->size() - 1) + (rate - 1)) / rate);
        
        // Initial partial overlap
        for (resultIndex=0; resultIndex<((int)this->size()-1+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()+rate-1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - (this->size()-1), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize((data.size() + rate - 1) / rate);
        
        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0; resultIndex<(((int)this->size()-1) - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex*rate; filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size() - initialTrim + rate - 1)/rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 filterIndex>=0; dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=resultIndex*rate - ((this->size()-1) - initialTrim), filterIndex=this->size()-1;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex--) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        break;
    }
    return data;
}

template <class T>
RealVector<T> & RealFirFilter<T>::interp(RealVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    T *savedDataArray = (T *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + rate - 1) / rate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * rate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((unsigned) dataTmp->size() * rate);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
                if (dataTmp->size() - dataStart == numSavedSamples) {
                    keepGoing = false;
                    phase = filterStart;
                }
            }
        }
        data.resize(resultIndex);

        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & RealFirFilter<T>::interpComplex(ComplexVector<T> & data, int rate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + rate - 1) / rate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * rate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        data.resize((unsigned) dataTmp->size() * rate);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
                if (dataTmp->size() - dataStart == numSavedSamples) {
                    keepGoing = false;
                    phase = filterStart;
                }
            }
        }
        data.resize(resultIndex);

        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate + this->size() - 1 - (rate - 1));
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
        
        // Middle full overlap
        for (dataStart=0, filterStart=resultIndex; resultIndex<(int)dataTmp->size()*rate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int) this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        data.resize(data.size() * rate);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0; resultIndex<(int)this->size()-1 - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=initialTrim + resultIndex; filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
        }
       
        // Middle full overlap
        for (dataStart=0, filterStart=(int)this->size()-1; resultIndex<(int)dataTmp->size()*rate - initialTrim; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }

        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=rate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            ++filterStart;
            if (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= rate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}

template <class T>
RealVector<T> & RealFirFilter<T>::resample(RealVector<T> & data, int interpRate, int decimateRate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    int interpLen, resampLen;
    std::vector<T> scratch;
    std::vector<T> *dataTmp;
    T *savedDataArray = (T *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + interpRate - 1) / interpRate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * interpRate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        interpLen = (unsigned) dataTmp->size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
            if (dataTmp->size() - dataStart == numSavedSamples) {
                keepGoing = false;
                phase = filterStart;
            }
        }
        data.resize(resultIndex);
        
        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}

template <class T>
ComplexVector<T> & RealFirFilter<T>::resampleComplex(ComplexVector<T> & data, int interpRate, int decimateRate, bool trimTails) {
    int resultIndex;
    int filterIndex;
    int dataIndex;
    int dataStart, filterStart;
    int interpLen, resampLen;
    std::vector< std::complex<T> > scratch;
    std::vector< std::complex<T> > *dataTmp;
    std::complex<T> *savedDataArray = (std::complex<T> *) VECTOR_TO_ARRAY(savedData);
    
    if (data.scratchBuf == NULL) {
        dataTmp = &scratch;
    }
    else {
        dataTmp = data.scratchBuf;
    }

    switch (filtOperation) {

    case STREAMING: {
        int numTaps = (this->size() + interpRate - 1) / interpRate;
        if (numSavedSamples >= numTaps) {
            // First call to interp, have too many "saved" (really just the initial zeros) samples
            numSavedSamples = numTaps - 1;
            phase = (numTaps - 1) * interpRate;
        }
        
        dataTmp->resize(numSavedSamples + data.size());
        for (int i=0; i<numSavedSamples; i++) {
            (*dataTmp)[i] = savedDataArray[i];
        }
        for (int i=0; i<data.size(); i++) {
            (*dataTmp)[i + numSavedSamples] = data[i];
        }
        
        interpLen = (unsigned) dataTmp->size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        bool keepGoing = true;
        for (resultIndex=0, dataStart=0, filterStart=phase; keepGoing; ++resultIndex) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
            if (dataTmp->size() - dataStart == numSavedSamples) {
                keepGoing = false;
                phase = filterStart;
            }
        }
        data.resize(resultIndex);
        
        int i;
        for (i=0; dataStart<dataTmp->size(); i++, dataStart++) {
            savedDataArray[i] = (*dataTmp)[dataStart];
        }
        numSavedSamples = i;
        }
        break;

    case ONE_SHOT_RETURN_ALL_RESULTS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate + this->size() - 1 - (interpRate - 1);
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);
        
        // Initial partial overlap
        for (resultIndex=0, dataStart=0, filterStart=0; resultIndex<((int)this->size()-1+decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;

    case ONE_SHOT_TRIM_TAILS:
        *dataTmp = data.vec;
        interpLen = data.size() * interpRate;
        resampLen = (interpLen + decimateRate - 1) / decimateRate;
        data.resize(resampLen);

        // Initial partial overlap
        int initialTrim = (this->size() - 1) / 2;
        for (resultIndex=0, dataStart=0, filterStart=initialTrim;
             resultIndex<((int)this->size()-1 - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=0, filterIndex=filterStart; filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Middle full overlap
        for (; resultIndex<((int)dataTmp->size()*interpRate - initialTrim + decimateRate-1)/decimateRate; resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 filterIndex>=0; dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        
        // Final partial overlap
        for (; resultIndex<(int)data.size(); resultIndex++) {
            data[resultIndex] = 0;
            for (dataIndex=dataStart, filterIndex=filterStart;
                 dataIndex<(int)dataTmp->size(); dataIndex++, filterIndex-=interpRate) {
                data[resultIndex] += (*dataTmp)[dataIndex] * this->vec[filterIndex];
            }
            filterStart += decimateRate;
            while (filterStart >= (int)this->size()) {
                // Filter no longer overlaps with this data sample, so the first overlap sample is the next one.  We thus
                // increment the data index and decrement the filter index.
                filterStart -= interpRate;
                ++dataStart;
            }
        }
        break;
    }
    return data;
}


template <class T>
bool RealFirFilter<T>::firpm(int filterOrder, int numBands, double *freqPoints, double *desiredBandResponse,
                            double *weight, int lGrid) {
    bool converged;
    
    this->resize(filterOrder + 1);
    std::vector<double> temp(filterOrder + 1);
    
    // Need to renormalize the frequency points because for us the Nyquist frequency is 1.0, but for the
    // Iowa Hills code it is 0.5.
    for (int i=0; i<numBands*2; i++) {
        freqPoints[i] /= 2;
    }
    
    // Move the pointers back 1 (i.e. subtract one) because the ParksMcClellan code was ported from Fortran,
    // which apparently uses 1-based arrays, not 0-based arrays.
    converged = ParksMcClellan2(&(temp[0]), filterOrder + 1, PASSBAND_FILTER, numBands, freqPoints-1,
                    desiredBandResponse-1, weight-1, lGrid);
    for (int i=0; i<= filterOrder; i++) {
        (*this)[i] = (T) temp[i];
    }
    return converged;
}


template <class T>
void RealFirFilter<T>::fractionalDelayFilter(int numTaps, double bandwidth, double delay) {
    assert(bandwidth > 0 && bandwidth < 1.0);
    assert(numTaps > 0);
    
    int index;
    double tapTime;
    double timeIncrement = bandwidth * M_PI;
    
    if (numTaps % 2) {
        tapTime = (numTaps / 2 + delay) * -timeIncrement;
    }
    else {
        tapTime = (((int) numTaps / 2) - 0.5 + delay) * -timeIncrement;
    }
    
    this->resize(numTaps);
    
    // Create the delayed sinc filter
    for (index = 0; index < numTaps; index++, tapTime += timeIncrement) {
        if (tapTime != 0.0) {
            (*this)[index] = (T) (sin(tapTime) / tapTime);
        }
        else {
            (*this)[index] = (T) 1.0;
        }
    }
    
    // Window the filter
    hamming();
}

template <class T>
void RealFirFilter<T>::hamming() {
    T phase = -M_PI;
    T phaseIncrement = 2 * M_PI / (this->size() - 1);
    T alpha = 0.54;
    T beta = 0.46;
    unsigned start, end;
    
    for (start=0, end=this->size()-1; start < end; start++, end--, phase += phaseIncrement) {
        double hammingVal = alpha + beta * cos(phase);
        (*this)[start] *= hammingVal;
        (*this)[end] *= hammingVal;
    }
}

template <class T>
RealVector<T> & RealFirFilter<T>::corr(RealVector<T> & data) {
    this->reverse();
    this->conv(data);
    this->reverse();
    return data;
}

/**
 * \brief Correlation function.
 *
 * \param data Buffer to operate on.
 * \param filter The filter that will correlate with "data".
 * \return Reference to "data", which holds the result of the convolution.
 */
template <class T>
inline RealVector<T> & corr(RealVector<T> & data, RealFirFilter<T> & filter) {
    return filter.corr(data);
}

template <class T>
void RealFirFilter<T>::hamming(unsigned len)
{
    generalizedHamming(len, 0.54, 0.46);
}

template <class T>
void RealFirFilter<T>::hann(unsigned len)
{
    generalizedHamming(len, 0.5, 0.5);
}

template <class T>
void RealFirFilter<T>::generalizedHamming(unsigned len, double alpha, double beta)
{
    this->resize(len);
    double N = len - 1;
    for (unsigned index=0; index<len; index++) {
        this->vec[index] = (T) (alpha - beta * cos(2 * M_PI * ((double) index) / N));
    }
}

template <class T>
void RealFirFilter<T>::blackman(unsigned len)
{
    const double alpha[] = {0.42, 0.5, 0.08};
    
    this->resize(len);
    double N = len - 1;
    for (unsigned index=0; index<len; index++) {
        this->vec[index] = (T) (alpha[0] - alpha[1] * cos(2 * M_PI * ((double) index) / N) +
                alpha[2] * cos(4 * M_PI * ((double) index) / N));
    }
}

template <class T>
void RealFirFilter<T>::blackmanHarris(unsigned len)
{
    const double alpha[] = {0.35875, 0.48829, 0.14128, 0.01168};
    
    this->resize(len);
    double N = len - 1;
    for (unsigned index=0; index<len; index++) {
        this->vec[index] = (T) (alpha[0] - alpha[1] * cos(2 * M_PI * ((double) index) / N) +
                alpha[2] * cos(4 * M_PI * ((double) index) / N) - alpha[3] * cos(6 * M_PI * ((double) index) / N));
    }
}


/**
 * \brief Vector class for real, fixed point (i.e. short's, int's, etc.) numbers.
 */
template <class T>
class RealFixedPtVector : public RealVector<T> {
 public:
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    RealFixedPtVector<T>(unsigned size = DEFAULT_BUF_LEN, std::vector<T> *scratch = NULL) :
            RealVector<T>(size, scratch) {}
            
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealFixedPtVector<T>(std::vector<U> data, std::vector<T> *scratch = NULL) : RealVector<T>(data, scratch) {}
    
    /**
     * \brief Array constructor.
     *
     * Sets buf equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref buf will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealFixedPtVector<T>(U *data, unsigned dataLen, std::vector<T> *scratch = NULL) :
            RealVector<T>(data, dataLen, scratch) {}
    
    /**
     * \brief Copy constructor.
     */
    RealFixedPtVector<T>(const RealFixedPtVector<T>& other) {this->vec = other.vec;}
    
    /*****************************************************************************************
                                            Operators
    *****************************************************************************************/
    /**
     * \brief Assignment operator.
     */
    RealFixedPtVector<T>& operator=(const Vector<T>& rhs);
    
    /**
     * \brief Pre-increment operator.
     */
    RealFixedPtVector<T> & operator++();
    
    /**
     * \brief Post-increment operator.
     */
    RealFixedPtVector<T> operator++(int);
    
    /**
     * \brief Pre-decrement operator.
     */
    RealFixedPtVector<T> & operator--();
    
    /**
     * \brief Post-decrement operator.
     */
    RealFixedPtVector<T> operator--(int);
    
    /**
     * \brief Buffer modulo/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator%=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar modulo/assignment operator.
     */
    RealFixedPtVector<T> & operator%=(const T &rhs);
    
    /**
     * \brief Bit-wise negation operator.
     */
    RealFixedPtVector<T> & operator~();
    
    /**
     * \brief Buffer bit-wise and/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator&=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar bit-wise and/assignment operator.
     */
    RealFixedPtVector<T> & operator&=(const T &rhs);
    
    /**
     * \brief Buffer bit-wise or/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator|=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar bit-wise or/assignment operator.
     */
    RealFixedPtVector<T> & operator|=(const T &rhs);
    
    /**
     * \brief Buffer bit-wise xor/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator^=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar bit-wise xor/assignment operator.
     */
    RealFixedPtVector<T> & operator^=(const T &rhs);
    
    /**
     * \brief Buffer right shift/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator>>=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar right shift/assignment operator.
     */
    RealFixedPtVector<T> & operator>>=(const T &rhs);
    
    /**
     * \brief Buffer left shift/assignment operator.
     */
    template <class U>
    RealFixedPtVector<T> & operator<<=(const RealFixedPtVector<U> &rhs);
    
    /**
     * \brief Scalar left shift/assignment operator.
     */
    RealFixedPtVector<T> & operator<<=(const T &rhs);
    
    /*****************************************************************************************
                                             Methods
    *****************************************************************************************/
    /**
     * \brief Sets each element of \ref buf equal to its value to the power of "exponent".
     *
     * \param exponent Exponent to use.
     * \return Reference to "this".
     */
    RealVector<T> & pow(const SLICKDSP_FLOAT_TYPE exponent);
    
    /**
     * \brief Returns the mode of the data in \ref buf.
     */
    const T mode();
};


template <class T>
RealFixedPtVector<T>& RealFixedPtVector<T>::operator=(const Vector<T>& rhs)
{
    this->vec = rhs.vec;
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator++()
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = this->vec[i] + 1;
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> RealFixedPtVector<T>::operator++(int)
{
    RealFixedPtVector<T> tmp(*this);
    operator++();
    return tmp;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator--()
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = this->vec[i] - 1;
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> RealFixedPtVector<T>::operator--(int)
{
    RealFixedPtVector<T> tmp(*this);
    operator--();
    return tmp;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator%=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] %= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator%=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] %= rhs;
    }
    return *this;
}

/**
 * \brief Buffer modulo operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator%(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs %= rhs;
    return lhs;
}

/**
 * \brief Scalar modulo operator.
 */
template <class T>
inline RealFixedPtVector<T> operator%(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs %= rhs;
    return lhs;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator~()
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = ~(this->vec[i]);
    }
    return *this;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator&=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] &= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator&=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] &= rhs;
    }
    return *this;
}

/**
 * \brief Buffer bit-wise and operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator&(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs &= rhs;
    return lhs;
}

/**
 * \brief Scalar bit-wise and operator.
 */
template <class T>
inline RealFixedPtVector<T> operator&(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs &= rhs;
    return lhs;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator|=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] |= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator|=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] |= rhs;
    }
    return *this;
}

/**
 * \brief Buffer bit-wise or operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator|(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs |= rhs;
    return lhs;
}

/**
 * \brief Scalar bit-wise or operator.
 */
template <class T>
inline RealFixedPtVector<T> operator|(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs |= rhs;
    return lhs;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator^=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] ^= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator^=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] ^= rhs;
    }
    return *this;
}

/**
 * \brief Buffer bit-wise xor operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator^(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs ^= rhs;
    return lhs;
}

/**
 * \brief Scalar bit-wise xor operator.
 */
template <class T>
inline RealFixedPtVector<T> operator^(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs ^= rhs;
    return lhs;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator>>=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] >>= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator>>=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] >>= rhs;
    }
    return *this;
}

/**
 * \brief Buffer right shift operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator>>(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs >>= rhs;
    return lhs;
}

/**
 * \brief Scalar right shift operator.
 */
template <class T>
inline RealFixedPtVector<T> operator>>(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs >>= rhs;
    return lhs;
}

template <class T>
template <class U>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator<<=(const RealFixedPtVector<U> &rhs)
{
    assert(this->size() == rhs.size());
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] <<= rhs.vec[i];
    }
    return *this;
}

template <class T>
RealFixedPtVector<T> & RealFixedPtVector<T>::operator<<=(const T &rhs)
{
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] <<= rhs;
    }
    return *this;
}

/**
 * \brief Buffer left shift operator.
 */
template <class T, class U>
inline RealFixedPtVector<T> operator<<(RealFixedPtVector<T> lhs, const RealFixedPtVector<U>& rhs)
{
    lhs <<= rhs;
    return lhs;
}

/**
 * \brief Scalar left shift operator.
 */
template <class T>
inline RealFixedPtVector<T> operator<<(RealFixedPtVector<T> lhs, const T& rhs)
{
    lhs <<= rhs;
    return lhs;
}

template <class T>
RealVector<T> & RealFixedPtVector<T>::pow(const SLICKDSP_FLOAT_TYPE exponent) {
    for (unsigned i=0; i<this->size(); i++) {
        this->vec[i] = (T) std::round(std::pow(this->vec[i], exponent));
    }
    return *this;
}
    
template <class T>
const T RealFixedPtVector<T>::mode() {
    assert(this->size() > 0);
    std::vector<T> tempScratch;
    std::vector<T> *scratch;

    if (this->scratchBuf == NULL) {
        scratch = &tempScratch;
    }
    else {
        scratch = this->scratchBuf;
    }
    *scratch = this->vec;
    std::sort(scratch->begin(), scratch->end());
    
    T modeVal = 0;
    unsigned modeLen = 0;
    T currentVal = (*scratch)[0];
    unsigned currentLen = 1;
    
    for (unsigned i=1; i<this->size(); i++) {
        if ((*scratch)[i] == currentVal) {
            currentLen++;
        }
        else {
            if (currentLen > modeLen) {
                modeVal = currentVal;
                modeLen = currentLen;
            }
            currentVal = (*scratch)[i];
            currentLen = 1;
        }
    }
    if (currentLen > modeLen) {
        modeVal = currentVal;
        modeLen = currentLen;
    }
    return modeVal;
}

/**
 * \brief Returns the mode of the data in \ref buf.
 */
template <class T>
const T mode(RealFixedPtVector<T> & buffer) {
    return buffer.mode();
}

/**
 * \brief Class for real IIR filters.
 */
template <class T>
class RealIirFilter {
 protected:
    std::vector<T> state;

    /** 
     * \brief Initializes numerator and denominator with the size and contents of "num" and "den".
     *
     * \param array Array to set vec equal to.
     * \param arrayLen Number of elements in array.
     */
    template <class U>
    void initArray(U *num, unsigned numLen, U *den, unsigned denLen);
    
 public:
    std::vector<T> numerator;
    std::vector<T> denominator;
    
    /*****************************************************************************************
                                        Constructors
    *****************************************************************************************/
    /**
     * \brief Basic constructor.
     *
     * Just sets the size of \ref buf and the pointer to the scratch buffer, if one is provided.
     * \param size Size of \ref buf.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    RealIirFilter<T>(unsigned order = 2) {numerator = std::vector<T>(order + 1); denominator = std::vector<T>(order + 1);
            state = std::vector<T>(order + 1);}
    
    /**
     * \brief Vector constructor.
     *
     * Sets buf equal to the input "data" parameter and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Vector that \ref buf will be set equal to.
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealIirFilter<T>(std::vector<U> & num, std::vector<U> & den) {initArray(VECTOR_TO_ARRAY(num), (unsigned) num.size(),
            VECTOR_TO_ARRAY(den), (unsigned) den.size());}
    
    /**
     * \brief Array constructor.
     *
     * Sets vec equal to the input "data" array and sets the pointer to the scratch buffer,
     *      if one is provided.
     * \param data Array that \ref vec will be set equal to.
     * \param dataLen Length of "data".
     * \param scratch Pointer to a scratch buffer.  The scratch buffer can be shared by multiple
     *      objects (in fact, I recommend it), but if there are multiple threads then it should
     *      be shared only by objects that are accessed by a single thread.  Objects in other
     *      threads should have a separate scratch buffer.  If no scratch buffer is provided
     *      then one will be created in methods that require one and destroyed when the method
     *      returns.
     */
    template <typename U>
    RealIirFilter<T>(U *num, unsigned numLen, U *den, unsigned denLen) {initArray(num, numLen, den, denLen);}
    
    /*****************************************************************************************
                                            Methods
    *****************************************************************************************/
    /**
     * \brief Convolution method.
     *
     * \param data The vector that will be filtered.
     * \param trimTails "False" tells the method to return the entire convolution, which is
     *      the length of "data" plus the length of "this" (the filter) - 1.  "True" tells the
     *      method to retain the size of "data" by trimming the tails at both ends of
     *      the convolution.
     * \return Reference to "data", which holds the result of the convolution.
     */
    template <class U>
    Vector<U> & filter(Vector<U> & data);
    
};


template <class T>
template <class U>
void RealIirFilter<T>::initArray(U *num, unsigned numLen, U *den, unsigned denLen) {
    numerator = std::vector<T>(numLen);
    for (unsigned i=0; i<numLen; i++) {
        numerator[i] = (T) num[i];
    }
    
    denominator = std::vector<T>(denLen);
    for (unsigned i=0; i<denLen; i++) {
        denominator[i] = (T) den[i];
    }
    
    state = std::vector<T>(std::max(numLen, denLen));
}

template <class T>
template <class U>
Vector<U> & RealIirFilter<T>::filter(Vector<U> & data) {
    unsigned resultIndex, i;
    U newState0;
    
    for (resultIndex=0; resultIndex<data.size(); resultIndex++) {
        newState0 = data[resultIndex];
        
        // Update the state
        for (i=(unsigned) state.size()-1; i>=denominator.size(); i--) {
            state[i] = state[i - 1];
        }
        // Update the state and apply the feedback
        for (i=(unsigned) denominator.size()-1; i>0; i--) {
            state[i] = state[i - 1];
            newState0 -= denominator[i] * state[i];
        }
        state[0] = newState0;

        // Calculate the output
        data[resultIndex] = 0;
        for (i=0; i<numerator.size(); i++) {
            data[resultIndex] += numerator[i] * state[i];
        }
    }
    return data;
}

template <class T, class U>
Vector<U> & filter(Vector<U> & data, RealIirFilter<T> & filt) {
    return filt.filter(data);
}
};
