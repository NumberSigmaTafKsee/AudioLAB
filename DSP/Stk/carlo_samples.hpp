#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>
#include <vector>
#include <new>
#include <chrono>
#include <random>
#include <cassert>

//#include "samples/Allocator.hpp"
#include "carlo_mkl.hpp"


//template<typename T> using sample_vector = std::vector<T,Allocator::aligned_allocator<T,64>>;
//template<typename T> using sample_matrix = std::vector<sample_vector<T>,Allocator::aligned_allocator<T,64>>;

namespace Casino::MKL
{

    template<typename T> using sample_vector = MKL::Vector<T>;
    template<typename T> using sample_matrix = MKL::Matrix<T>;
    
    template<typename T>
    using complex_vector = sample_vector<std::complex<T>>;

    template<typename T>
    using complex_matrix = sample_matrix<std::complex<T>>;
}


namespace Casino
{
    using std::vector;
    template<typename T> using sample_vector  = Casino::MKL::sample_vector<T>;
    template<typename T> using sample_matrix  = Casino::MKL::sample_matrix<T>;
    template<typename T> using complex_vector = Casino::MKL::complex_vector<T>;
    template<typename T> using complex_matrix = Casino::MKL::complex_matrix<T>;
        

    template<typename T>
    struct StereoVector
    {
        sample_vector<T> samples[2];

        StereoVector(size_t n) {
            for(size_t i = 0; i < n; i++) samples[i].resize(n);
        }

        sample_vector<T>& operator[](size_t channel) { return samples[channel]; }

        sample_vector<T> get_left_channel() { return samples[0]; }
        sample_vector<T> get_right_channel() { return samples[1]; }

        void set_left_channel(sample_vector<T> & v) { samples[0] = v; }
        void set_right_channel(sample_vector<T> & v) { samples[1] = v; }
    };


    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const sample_vector<T> & left, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const sample_vector<T> & right, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const sample_vector<T> & in, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<sample_vector<T>> & in) {
        sample_vector<T> r(n*channels);
        #pragma omp parallel
        {       
        for(size_t i = 0; i < channels; i++)
            #pragma omp simd
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        }
        return r;
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<T*> & in) {
        sample_vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const sample_vector<T> & in) {
        sample_vector<sample_vector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }

    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, T * samples)
    {
        return samples[pos*num_channels + ch];
    }
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, T * samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(size_t n, const T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(size_t n, const T* & in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, size_t n, T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < n; i+=2) r[x++] = in[i];
        return r;
    }

    template<typename T>
    void set_left_channel(size_t n, const T* left, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(size_t n, const T* right, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, size_t n, const T* in, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < n; i+=2) out[i] = in[x++];
    }

    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const T** & in) {
        sample_vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        sample_vector<sample_vector<T>> r(n);
        #pragma omp parallel
        {
            for(size_t i = 0; i < channels; i++)
            {
                r[i].resize(n);
                #pragma omp simd
                for(size_t j = 0; j < n; j++)
                    r[i][j] = in[j*channels + i];
            }
        }
        return r;
    }

    template<typename T>
    bool equal_vector (sample_vector<T> & a, sample_vector<T> & b) {
        return std::equal(a.begin(),a.end(),b.end());
    }

    template<typename T>
    void copy_vector(sample_vector<T> & dst, sample_vector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy_vector(sample_vector<T> & dst, size_t n, T * src) {
        std::copy(&src[0],&src[n-1],dst.begin());
    }
    template<typename T>
    sample_vector<T> slice_vector(size_t start, size_t end, sample_vector<T> & src) {
        sample_vector<T> r(end-start);
        std::copy(src.begin()+start,src.begin()+end,r.begin());
        return r;
    }

    template<typename T>
    void copy_buffer(size_t n, T * dst, T * src) {
        memcpy(dst,src,n*sizeof(T));
    }

    template<typename T>
    sample_vector<T> slice_buffer(size_t start, size_t end, T * ptr) {
        sample_vector<T> r(end-start);    
        std::copy(&ptr[start],&ptr[end],r.begin());
        return r;
    }

    template<typename T>
    void split_stereo(size_t n, const T* input, T * left, T * right)
    {
        size_t x=0;
        for(size_t i = 0; i < n; i+=2)
        {
            left[x] = input[i];
            right[x++] = input[i+1];
        }
    }

    template<typename T>
    void split_stereo(const sample_vector<T> & input, sample_vector<T> & left, sample_vector<T> & right) {
        size_t x = input.size();
        left.resize(x/2);
        right.resize(x/2);
        split_stereo(x,input.data(),left.data(),right.data());
    }

    template<typename T>
    T insert_front(size_t n, T in, T * buffer) {
        T r = buffer[n-1];
        #pragma omp simd
        for(size_t i=0; i < n-1; i++) buffer[n+1] = buffer[n];
        buffer[0] = in;
        return r;
    }

    //============================================================
    template <class T>
    bool isEmpty(sample_vector<T> v)
    {
        return (v.size() == 0);
    }

    //============================================================
    template <class T>
    bool containsOnlyZeros(sample_vector<T> v)
    {
        bool r = true;
        if (!isEmpty(v))
        {                    
            for (int i = 0;i < v.size();i++)
            {
                if (v[i] != 0)
                {
                    return false;                    
                }
            }

            return r;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector contains only zeros" );
        }
    }

    //============================================================
    template <class T>
    bool isAllPositiveOrZero(sample_vector<T> v)
    {
        if (!isEmpty(v))
        {                 
            for (int i = 0;i < v.size();i++)
            {
                if (v[i] < 0)
                {
                    return false;
                }
            }

            return true;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector is all positive" );
        }
    }

    //============================================================
    template <class T>
    bool isAllNegativeOrZero(sample_vector<T> v)
    {
        if (!isEmpty(v))
        {                 
            for (int i = 0;i < v.size();i++)
            {
                if (v[i] > 0)
                {
                    return false;
                }
            }

            return true;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector is all negative" );
        }
    }

    //============================================================
    template <class T>
    bool contains(sample_vector<T> v, T element)
    {           
        for (int i = 0;i < v.size();i++)
        {
            if (v[i] == element)
            {
                return true;
            }
        }

        return false;
    }


    //============================================================
    template <class T>
    T max(sample_vector<T> v)
    {
        //return std::max_element(v.begin(),v.end());
        
        // if the vector isn't empty
        if (!isEmpty(v))
        {
            // set the first element as the max
            T maxVal = v[0];

            // then for each subsequent value
            #pragma omp simd
            for (int i = 1;i < v.size();i++)
            {
                // if it is larger than the max
                if (v[i] > maxVal)
                {
                    // store it as the new max value
                    maxVal = v[i];
                }
            }

            // return the maxVal
            return maxVal;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating max" );
        }    
    }

    //============================================================
    template <class T>
    int maxIndex(sample_vector<T> v)
    {
        // if the vector isn't empty
        if (!isEmpty(v))
        {
            // set the first element as the max
            T maxVal = v[0];
            int maxIndex = 0;

            // then for each subsequent value
            #pragma omp simd
            for (int i = 1;i < v.size();i++)
            {
                // if it is larger than the max
                if (v[i] > maxVal)
                {
                    // store it as the new max value
                    maxVal = v[i];

                    // store the index as the new max index
                    maxIndex = i;
                }
            }

            // return the max index
            return maxIndex;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating max index" );
        }
    }

    //============================================================
    template <class T>
    T min(sample_vector<T> v)
    {
        // if the vector isn't empty
        if (!isEmpty(v))
        {
            // set the first element as the min
            T minVal = v[0];

            // then for each subsequent value
            #pragma omp simd
            for (int i = 1;i < v.size();i++)
            {
                // if it is smaller than the min
                if (v[i] < minVal)
                {
                    // store it as the new min value
                    minVal = v[i];
                }
            }

            // return the minVal
            return minVal;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating min" );
        }
    }

    //============================================================
    template <class T>
    int minIndex(sample_vector<T> v)
    {
        // if the vector isn't empty
        if (!isEmpty(v))
        {
            // set the first element as the min
            T minVal = v[0];
            int minIndex = 0;

            // then for each subsequent value
            #pragma omp simd
            for (int i = 1;i < v.size();i++)
            {
                // if it is smaller than the min
                if (v[i] < minVal)
                {
                    // store it as the new min value
                    minVal = v[i];

                    // store the index as the new min index
                    minIndex = i;
                }
            }

            // return the min index
            return minIndex;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating minIndex" );
        }
    }

    //============================================================
    template <class T>
    void printVector(sample_vector<T> v)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            std::cout << v[i] << std::endl;
        }
    }

    //============================================================
    template <class T>
    T getLastElement(sample_vector<T> v)
    {
        // if vector is not empty
        if (v.size() > 0)
        {
            return v[v.size()-1];
        }
        else
        {
            throw std::invalid_argument( "Attempted to get last element of empty vector" );
        }
    }

    //============================================================
    template <class T>
    T getFirstElement(sample_vector<T> v)
    {
        // if vector is not empty
        if (v.size() > 0)
        {
            return v[0];
        }
        else
        {
            throw std::invalid_argument( "Attempted to get first element of empty vector" );
        }
    }


    //============================================================
    template <class T>
    sample_vector<T> getEveryNthElementStartingFromK(sample_vector<T> v,int n,int k)
    {
        if ((n >= v.size()) || (n >= v.size()))
        {
            throw std::invalid_argument( "Invalid arguments for getEveryNthElementStartingFromK()");
        }
        else
        {
            sample_vector<T> result;

            
            #pragma omp simd
            for(int i = k; i < v.size(); i += n)        
            {
                result.push_back(v[i]);            
            }

            return result;
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> getEvenElements(sample_vector<T> v)
    {
        return getEveryNthElementStartingFromK(v, 2, 0);
    }

    //============================================================
    template <class T>
    sample_vector<T> getOddElements(sample_vector<T> v)
    {
        return getEveryNthElementStartingFromK(v, 2, 1);
    }

    //============================================================
    template <class T>
    void fillVectorWith(sample_vector<T> &v,T element)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            v[i] = element;
        }
    }

    //============================================================
    template <class T>
    int countOccurrencesOf(sample_vector<T> v,T element)
    {
        int count = 0;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            if (v[i] == element)
            {
                count++;
            }
        }

        return count;
    }

    //============================================================
    template <class T>
    T sum(sample_vector<T> v)
    {
        // create a sum
        T sumVal = 0;

        // add up all elements
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            sumVal += v[i];
        }

        // return
        return sumVal;
    }

    //============================================================
    template <class T>
    T product(sample_vector<T> v)
    {
        if (!isEmpty(v))
        {
            T prod = (T) v[0];

            #pragma omp simd
            for (int i = 1;i < v.size();i++)
            {
                prod *= ((T) v[i]);
            }

            return prod;
        }
        else
        {
            throw std::invalid_argument( "Attempted to calculate the product of an empty vector" );
        }
    }

    //============================================================
    template <class T>
    T mean(sample_vector<T> v)
    {
        // if vector is not empty
        if (!isEmpty(v))
        {
            // store the length of the vector as a T
            T L = (T) v.size();

            // stor the sum of the vector as a T
            T sumVal = (T) sum(v);

            // return the mean
            return sumVal / L;
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating mean" );
        }
    }

    //============================================================
    template <class T>
    T median(sample_vector<T> v)
    {
        // if vector isn't empty
        if (!isEmpty(v))
        {
            T median;
            size_t L = v.size(); // store the size

            // sort the vector
            std::sort(v.begin(), v.end());

            // if the length is even
            if (L  % 2 == 0)
            {
                // take the average of the middle two elements
                median = ((T)(v[L / 2 - 1] + v[L / 2])) / 2.0;
            }
            else // if the length is odd
            {
                // take the middle element
                median = (T) v[(L-1) / 2];
            }

            // return the median
            return median;
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating median" );
        }
    }

    //============================================================
    template <class T>
    T variance(sample_vector<T> v)
    {
        if (!isEmpty(v))
        {
            // calculate the mean of the vector
            T mu = mean(v);

            T sumVal = 0.0;

            // sum the product of all differences from the mean
            #pragma omp simd
            for (int i = 0;i < v.size();i++)
            {
                T diff = v[i]-mu;
                sumVal += diff*diff;
            }

            // return the average of the squared differences
            return sumVal / ((T)v.size());
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating variance" );
        }
    }

    //============================================================
    template <class T>
    T standardDeviation(sample_vector<T> v)
    {
        // if vector is not empty
        if (!isEmpty(v))
        {
            // calculate the variance
            T var = variance(v);

            // if variance is non-zero
            if (var > 0)
            {
                // return the square root
                return std::sqrt(var);
            }
            else
            {
                // all differences are zero, so return 0.0
                return 0.0;
            }
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating standard deviation" );
        }
    }

    //============================================================
    template <class T>
    T norm1(sample_vector<T> v)
    {
        T sumVal = 0.0;

        // sum absolute values
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            if (v[i] > 0)
            {
                sumVal += (T) v[i];
            }
            else
            {
                sumVal += (T) (-1*v[i]);
            }
        }

        return sumVal;
    }

    //============================================================
    template <class T>
    T norm2(sample_vector<T> v)
    {
        T sumVal = 0.0;

        // sum squares
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            sumVal += (T) (v[i]*v[i]);
        }

        return std::sqrt(sumVal);
    }

    //============================================================
    template <class T>
    T magnitude(sample_vector<T> v)
    {
        // just another name for L2-norm
        return norm2(v);
    }

    //============================================================
    template <class T>
    T normP(sample_vector<T> v,T p)
    {
        T sumVal = 0.0;

        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            T val;

            if (v[i] > 0)
            {
                val = (T) v[i];
            }
            else
            {
                val = (T) (-1*v[i]);
            }

            sumVal += std::pow(val,p);
        }

        return std::pow(sumVal,1.0/p);
    }

    //============================================================
    template <class T>
    void multiplyInPlace(sample_vector<T> &v,T scalar)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            v[i] *= scalar;
        }
    }

    //============================================================
    template <class T>
    void multiplyInPlace(sample_vector<T> &v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] *= v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector multiply");
        }
    }

    //============================================================
    template <class T>
    void divideInPlace(sample_vector<T> &v,T scalar)
    {
        if (scalar != 0)
        {
            #pragma omp simd
            for (int i = 0;i < v.size();i++)
            {
                v[i] /= scalar;
            }
        }
        else
        {
            throw std::invalid_argument( "Attemted to divide a vector by a zero-valued scalar" );
        }
    }

    //============================================================
    template <class T>
    void divideInPlace(sample_vector<T> &v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            if (!contains<T>(v2, 0))
            {
                #pragma omp simd
                for (int i = 0;i < v1.size();i++)
                {
                    v1[i] /= v2[i];
                }
            }
            else
            {
                throw std::invalid_argument( "Attempted to divide by vector containing zeros");
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector divide");
        }
    }

    //============================================================
    template <class T>
    void addInPlace(sample_vector<T> &v,T value)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            v[i] += value;
        }
    }

    //============================================================
    template <class T>
    void addInPlace(sample_vector<T> &v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] += v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector add");
        }
    }

    //============================================================
    template <class T>
    void subtractInPlace(sample_vector<T> &v,T value)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] -= value;
        }
    }

    //============================================================
    template <class T>
    void subtractInPlace(sample_vector<T> &v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] -= v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector subtraction");
        }

    }

    //============================================================
    template <class T>
    void absInPlace(sample_vector<T> &v)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {        
            if ((v[i] < 0) || (v[i] == -0.0))
            {
                v[i] *= -1;
            }
        }
    }

    //============================================================
    template <class T>
    void squareInPlace(sample_vector<T> &v)
    {
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            v[i] = v[i]*v[i];
        }
    }

    //============================================================
    template <class T>
    void squareRootInPlace(sample_vector<T> &v)
    {
        if (isAllPositiveOrZero(v))
        {
            #pragma omp simd
            for (int i = 0;i < v.size();i++)
            {
                v[i] = (T) std::sqrt((T)v[i]);
            }
        }
        else
        {
            throw std::invalid_argument( "Attempted to compute square root of vector containing negative numbers");
        }
    }


    //============================================================
    template <class T>
    void sort(sample_vector<T> &v)
    {
        std::sort(v.begin(),v.end());
    }

    //============================================================
    template <class T>
    void reverse(sample_vector<T> &v)
    {
        std::reverse(v.begin(), v.end());
    }

    //============================================================
    template <class T>
    sample_vector<T> multiply(sample_vector<T> v,T scalar)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] * scalar);
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> multiply(sample_vector<T> v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            sample_vector<T> result;
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result.push_back(v1[i] * v2[i]);
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector multiply");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> divide(sample_vector<T> v,T scalar)
    {
        if (scalar != 0)
        {
            sample_vector<T> result;

            #pragma omp simd
            for (int i = 0;i < v.size();i++)
            {
                result.push_back(v[i] / scalar);
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Attemted to divide a vector by a zero-valued scalar" );
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> divide(sample_vector<T> v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            if (!contains<T>(v2, 0))
            {
                sample_vector<T> result;

                #pragma omp simd
                for (int i = 0;i < v1.size();i++)
                {
                    result.push_back(v1[i] / v2[i]);
                }

                return result;
            }
            else
            {
                throw std::invalid_argument( "Attempted to divide by vector containing zeros");
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector divide");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> add(sample_vector<T> v,T value)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] + value);
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> add(sample_vector<T> v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            sample_vector<T> result(v1.size());
            int N = v1.size();
            //#pragma omp target map(to:N) map(tofrom: result)
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result[i] = (v1[i] + v2[i]);
            }
            
            return result;
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector add");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> subtract(sample_vector<T> v,T value)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] - value);
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> subtract(sample_vector<T> v1,sample_vector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            sample_vector<T> result;
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result.push_back(v1[i] - v2[i]);
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector subtraction");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> abs(sample_vector<T> v)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            if ((v[i] < 0) || (v[i] == -0.0))
            {
                result.push_back(-1*v[i]);
            }
            else
            {
                result.push_back(v[i]);
            }
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> square(sample_vector<T> v)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i]*v[i]);
        }

        return result;
    }


    //============================================================
    template <class T>
    sample_vector<T> squareRoot(sample_vector<T> v)
    {
        if (isAllPositiveOrZero(v))
        {
            sample_vector<T> result;
            #pragma omp simd
            for (int i = 0;i < v.size();i++)
            {
                result.push_back((T) std::sqrt((T)v[i]));
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Attempted to compute square root of vector containing negative numbers");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> scale(sample_vector<T> v,T lowerLimit,T upperLimit)
    {
        sample_vector<T> result;

        T minVal = (T) min(v);
        T maxVal = (T) max(v);
        T outputRange = upperLimit - lowerLimit;
        T inputRange = maxVal - minVal;

        #pragma omp simd
        for (int i = 0;i < v.size();i++)
        {
            T value = (T) v[i];
            T scaledValue = ((value - minVal) * outputRange) / inputRange + lowerLimit;

            result.push_back(scaledValue);
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> difference(sample_vector<T> v)
    {
        sample_vector<T> result;
        #pragma omp simd
        for (int i = 1;i < v.size();i++)
        {
            result.push_back(v[i]-v[i-1]);
        }

        return result;
    }


    //============================================================
    template <class T>
    sample_vector<T> zeros(int N)
    {
        if (N >= 0)
        {
            sample_vector<T> result;
            #pragma omp simd
            for (int i = 0;i < N;i++)
            {
                result.push_back(0);
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Cannot create vector with negative length");
        }
    }

    //============================================================
    template <class T>
    sample_vector<T> ones(int N)
    {
        if (N >= 0)
        {
            sample_vector<T> result;
            #pragma omp simd
            for (int i = 0;i < N;i++)
            {
                result.push_back(1);
            }

            return result;
        }
        else
        {
            throw std::invalid_argument( "Cannot create vector with negative length");
        }
    }


    //============================================================
    template <class T>
    sample_vector<T> range(int limit1,int limit2,int step)
    {
        sample_vector<T> result;

        if (step > 0)    
        {
            #pragma omp simd
            for (T i = limit1;i < limit2;i += step)
            {
                result.push_back(i);
            }
        }
        else if (step < 0)
        {
            #pragma omp simd
            for (T i = limit1;i > limit2;i += step)
            {
                result.push_back(i);
            }
        }
        else
        {
            throw std::invalid_argument( "Cannot use a step size of 0 when creating a range of values");
        }

        return result;
    }

    //============================================================
    template <class T>
    sample_vector<T> range(int maxValue)
    {
        return range<T>(0, maxValue, 1);
    }

    //============================================================
    template <class T>
    sample_vector<T> range(int minValue,int maxValue)
    {
        return range<T>(minValue, maxValue, 1);
    }

    //============================================================
    template <class T>
    T dotProduct(sample_vector<T> v1,sample_vector<T> v2)
    {
        // if vector size is the same
        if (v1.size() == v2.size())
        {
            T sumVal = 0.0;        
            // sum the element-wise product
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {            
                sumVal += (v1[i]*v2[i]);
            }

            // return the sum as the dot product
            return sumVal;
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector dot product");
        }
    }

    //============================================================
    template <class T>
    T euclideanDistance(sample_vector<T> v1,sample_vector<T> v2)
    {
        // if vector size is the same
        if (v1.size() == v2.size())
        {
            T sumVal = 0.0;

            // sum the squared difference
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                T diff = (T) (v1[i] - v2[i]);
                sumVal += (diff*diff);
            }

            // if sum is bigger than zero
            if (sumVal > 0)
            {
                // return the square root of the sum as the Euclidean distance
                return std::sqrt(sumVal);
            }
            else // all differences were zero, so report 0.0 as Euclidean distance
            {
                return 0.0;
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in Euclidean distance calculation");
        }
    }

    //============================================================
    template <class T>
    T cosineSimilarity(sample_vector<T> v1,sample_vector<T> v2)
    {
    return dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
    }

    //============================================================
    template <class T>
    T cosineDistance(sample_vector<T> v1,sample_vector<T> v2)
    {
        return 1.0 - cosineSimilarity(v1, v2);
    }

    
    
    template<typename T>
    StereoVector<T> stereo(const sample_vector<T> & left, const sample_vector<T> & right) {
        StereoVector<T> r(left.size());
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i];
            r[1][i] = right[i];
        }
    }

    template<typename T>
    sample_vector<T> merge(const sample_vector<T> & left, const sample_vector<T> & right) {
        sample_vector<T> r(left.size()*2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[x++] = left[i];
            r[x++] = right[i];
        }
    }

    template<typename T>
    void swap(sample_vector<T> & left, sample_vector<T> & right) {
        std::swap(left,right);
    }

    template<typename T>
    bool isin(const sample_vector<T> & v, const T val) {
        return std::find(v.begin(),v.end(),val) != v.end();
    }

    template<typename T>
    StereoVector<T> pan(const sample_vector<T> & left, const sample_vector<T> & right, T amt) {
        StereoVector<T> r(left.size());
        T pan_map = ((amt+1)/2.0) * (M_PI/2.0);
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * sin(pan_map);
            r[1][i] = right[i] * cos(pan_map);
        }
        return r;
    }
    template<typename T>
    StereoVector<T> constant_power_pan(const sample_vector<T> & left, const sample_vector<T> & right, T pos) {
        StereoVector<T> r(left.size());        
        const T piover2 = 4.0*std::atan(1.0)*0.5;
        const T root2over2 = std::sqrt(2.0)*0.5;
        T thispos = pos * piover2;
        T angle   = thispos * 0.5;
        T pleft   = root2over2 * (std::cos(angle) - std::sin(angle));
        T pright  = root2over2 * (std::cos(angle) + std::sin(angle));
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * pleft;
            r[1][i] = right[i] * pright;
        }
        return r;
    }
    template<typename T>
    sample_vector<T> mix(const sample_vector<T> & a, const sample_vector<T> & b)
    {
        assert(a.size() == b.size());
        sample_vector<T> r(a.size());
        T max = -99999;
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) 
        {
            r[i] = a[i]+b[i];
            if(fabs(r[i]) > max) max = fabs(r[i]);
        }
        if(max > 0) for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }
    template<typename T>
    sample_vector<T> normalize(const sample_vector<T> & a) {
        sample_vector<T> r(a);        
        T max = std::max_element(r.begin(),r.end());
        
        if(max > 0) 
            #pragma omp simd
            for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }
    template<class A, class B>
    sample_vector<B> convert(const sample_vector<A> & v) {
        sample_vector<B> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = B(v[i]);
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, T (*f)(T value)) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = f(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, std::function<T (T)> func) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = func(v[i]);
        return r;
    }
    template<class T>
    void inplace_add(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] += func(a[i]);        
    }
    template<class T>
    void inplace_sub(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] -= func(a[i]);        
    }
    template<class T>
    void inplace_mul(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] *= func(a[i]);        
    }
    template<class T>
    void inplace_div(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] /= func(a[i]);
    }


    template<class T>
    void fill(sample_vector<T> & in, T x)
    {
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i++) in[i] = x;
    }
    template<class T>
    void zeros(sample_vector<T> & in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(sample_vector<T> & in)
    {
        fill(in,T(1));
    }

    template<typename T>
    void stdmean(T *sig_src_arr, uint32_t blockSize, T * result){
        T sum = T(0);
        uint32_t blkCnt;
        T in1,in2,in3, in4;
        assert(blockSize != 0);
        //Right shifted by 4 so divided by 4

        #pragma omp simd
        for(blkCnt = blockSize>>2U;blkCnt > 0;blkCnt--){
            in1 = *sig_src_arr++;
            in2 = *sig_src_arr++;
            in3 = *sig_src_arr++;
            in4 = *sig_src_arr++;
            sum += in1;
            sum += in2;
            sum += in3;
            sum += in4;        
        }
        

        #pragma omp simd
        for(blkCnt = blockSize% 0x4;blkCnt > 0;blkCnt--){
            sum += *sig_src_arr++;        
        }
        
        *result = sum/(T)blockSize;        
    }
    template<typename T>
    sample_vector<T> stdmean(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdMean(in.data(),in.size(),r.data());
        return r;
    }

    template<typename T>
    void stdrms(T *pSig_src_arr, uint32_t blockSize, T *pResult)
    {
        T sum = 0.0;
        uint32_t blkCnt;
        T in;
        assert(blockSize != 0);
        
        #pragma omp simd
        for(blkCnt = blockSize >>2;blkCnt > 0;blkCnt--){
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;        
        }
        
        #pragma omp simd
        for(blkCnt = blockSize%0x4;blkCnt>0;blkCnt--)
        {
            in = *pSig_src_arr++;
            sum += in*in;        
        }        
        *pResult = std::sqrt(sum/(T)blockSize);
    }

    template<typename T>
    sample_vector<T> StdRMS(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdRMS(in.data(),in.size(),r.data());
        return r;
    }

    template<typename T>
    void stddev(T * pSig_src_arr, uint32_t blockSize, T *pResult)
    {

        T sum = 0.0;
        T sumOfSquares = 0.0;
        T in;

        uint32_t blkCnt;

        T meanOfSquares, mean, squareOfMean;
        T squareOfSum = 0.0;

        T var;

        if(blockSize == 1){
            *pResult = 0;
            return;
        }

        

        #pragma omp simd
        for(blkCnt = blockSize>>2;blkCnt>0;blkCnt--){
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;        
        }

        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt>0;blkCnt--){
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        }

        meanOfSquares = sumOfSquares / ((T)blockSize-1.0);
        mean = sum/(T) blockSize;

        squareOfMean = (mean*mean) * ((T)blockSize/(T)(blockSize-1.0));

        *pResult = sqrt((meanOfSquares-squareOfMean));
    }

    template<typename T>
    sample_vector<T> stddev(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdDev(in.data(),in.size(),r.data());
        return r;
    }

    template<typename T>
    void stdvariance(T * pSig_src_arr, uint32_t blockSize, T *pResult)
    {
        T fMean, fValue;
        uint32_t blkCnt;
        T * pInput = pSig_src_arr;

        T sum = 0.0;
        T fSum = 0.0;

        T in1, in2, in3, in4;

        if(blockSize <= 1){
            *pResult = 0;
            return;
        }

        
        #pragma omp simd
        for(blkCnt = blockSize >>2U;blkCnt>0;blkCnt--){
            in1 = *pInput++;
            in2 = *pInput++;
            in3 = *pInput++;
            in4 = *pInput++;

            sum+= in1;
            sum+= in2;
            sum+= in3;
            sum+= in4;    
        }
        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt > 0;blkCnt--){
            sum += *pInput++;
        }

        fMean = sum/(T) blockSize;
        pInput = pSig_src_arr;
        
        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt > 0;blkCnt--){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
        }
        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt>0;blkCnt--){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
        }

        *pResult = fSum/(T)(blockSize-1.0);
    }

    template<typename T>
    sample_vector<T> stdvariance(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdVariation(in.data(),in.size(),r.data());
        return r;
    }

    template<typename A, typename B>
    sample_vector<A> vector_cast(sample_vector<B> & in) {
        sample_vector<A> r(in.size());
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i++)
            r[i] = (A)in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> vector_copy(T * ptr, size_t n) {
        sample_vector<T> r(n);
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = ptr[i];
        return r;
    }


    template <class T>
    void zeros(std::vector<T> & v) {
        std::fill(v.begin(),v.end(),T(0));
    }

    template<typename T>
    // r = frac
    // x = [i]
    // y = [i+1]
    T linear_interpolate(T x, T y, T r)
    {        
        return x + r*(y-x);
    }
    template<typename T>
    T cubic_interpolate(T finpos, T xm1, T x0, T x1, T x2)
    {
        //T xm1 = x [inpos - 1];
        //T x0  = x [inpos + 0];
        //T x1  = x [inpos + 1];
        //T x2  = x [inpos + 2];
        T a = (3 * (x0-x1) - xm1 + x2) / 2;
        T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
        T c = (x1 - xm1) / 2;
        return (((a * finpos) + b) * finpos + c) * finpos + x0;
    }
    // original
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite1(T x, T y0, T y1, T y2, T y3)
    {
        // 4-point, 3rd-order Hermite (x-form)
        T c0 = y1;
        T c1 = 0.5f * (y2 - y0);
        T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
        T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // james mccartney
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite2(T x, T y0, T y1, T y2, T y3)
    {
        // 4-point, 3rd-order Hermite (x-form)
        T c0 = y1;
        T c1 = 0.5f * (y2 - y0);
        T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        T c2 = y0 - y1 + c1 - c3;
        return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // james mccartney
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite3(T x, T y0, T y1, T y2, T y3)
    {
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T y0my1 = y0 - y1;
            T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
            T c2 = y0my1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // laurent de soras
    template<typename T>
    // x[i-1]
    // x[i]
    // x[i+1]
    // x[i+2]    
    inline T hermite4(T frac_pos, T xm1, T x0, T x1, T x2)
    {
        const T    c     = (x1 - xm1) * 0.5f;
        const T    v     = x0 - x1;
        const T    w     = c + v;
        const T    a     = w + v + (x2 - x0) * 0.5f;
        const T    b_neg = w + a;

        return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
    }

    template<typename T>
    sample_vector<T> mix(const sample_vector<T> & a, const sample_vector<T> & b, T f)
    {
        assert(a.size() == b.size());
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) r[i] = a[i] + f*(b[i]-a[i]);
        return r;
    }

    template<typename T>
    sample_vector<T> interp2x(const sample_vector<T> & a)
    {    
        sample_vector<T> r(a.size()*2);
        size_t n=0;
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) 
        {
            r[n++] = a[i];
            r[n++] = cubic_interpolate(T(0.5),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
        }
        return r;
    }

    template<typename T>
    sample_vector<T> interp4x(const sample_vector<T> & a)
    {    
        sample_vector<T> r(a.size()*4);
        size_t n=0;
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) 
        {
            r[n++] = a[i];
            r[n++] = cubic_interpolate(T(0.25),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
            r[n++] = cubic_interpolate(T(0.5),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
            r[n++] = cubic_interpolate(T(0.75),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
        }
        return r;
    }

    template<typename T>
    struct RingBuffer : public sample_vector<T>
    {
        size_t r=0;
        size_t w=0;

        RingBuffer(size_t n) {
            sample_vector<T>::resize(n);
        }

        void set_write_position(size_t n) {
            w = (n % sample_vector<T>::size());
        }  
        T    get() {
            return (*this)[r++];
        }
        void push(T x) {
            (*this)[w++] = x;
            w = (w % sample_vector<T>::size());
        }
        T linear() {
            T x = (*this)[r];
            T x1= (*this)[r++];
            T f = x - floor(x);
            r = r % sample_vector<T>::size();
            return linear_interpolate(x,x1,f);        
        }
        T cubic() {
            T xm1= (*this)[(r-1) % sample_vector<T>::size()];
            T x = (*this)[r];
            T x1= (*this)[(r+1) % sample_vector<T>::size()];
            T x2= (*this)[(r+2) % sample_vector<T>::size()];
            T f = x - floor(x);
            r++;
            r = r % sample_vector<T>::size();
            return cubic_interpolate(f,xm1,x,x1,x2);        
        }
    };
}

