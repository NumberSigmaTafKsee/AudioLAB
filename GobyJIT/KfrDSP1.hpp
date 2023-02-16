#pragma once

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdint>
#include <cstring>

#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>
#include <random>
#include <chrono>


#include "Undenormal.hpp"

#include "Kfr/kfrcore.hpp"
#include "Std/StdObject.h"
#include "Std/StdRandom.h"


#define KFR_NO_C_COMPLEX_TYPES
#include "kfr/capi.h"


template<typename T> using sample_vector = kfr::univector<T>;
template<typename T> using sample_matrix = kfr::univector2d<T>;
template<typename T> using complex = kfr::complex<T>;
template<typename T> using complex_vector = kfr::univector<complex<T>>;

extern double sampleRate;
extern double invSampleRate;
extern Std::RandomMersenne noise;

#define PINK_NOISE_NUM_STAGES 3

namespace DSP1
{    
    struct matrix_index {
        size_t i;
        size_t j;
    };

    template<typename T>
    struct vector_matrix : public sample_vector<T>
    {
        size_t M,N;
        vector_matrix(size_t m, size_t n) {
            M = m;
            N = n;
            resize(M,N);
        }

        size_t size() const;
        bool   empty() const;

        void   resize(size_t i, size_t j) { this->resize(i*j); M=i; N=j; }
        size_t rows() const { return M; }
        size_t cols() const { return N; }

        T& operator()(size_t i, size_t j) { return (*this)[i*N + j]; }
        T& at(size_t i, size_t j) { return (*this)[i*N + j]; }
    };

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
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const sample_vector<T> & left, sample_vector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const sample_vector<T> & right, sample_vector<T> & out) {
        size_t x = 0;
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const sample_vector<T> & in, sample_vector<T> & out) {
        size_t x = 0;
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const std::vector<sample_vector<T>> & in) {
        sample_vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const std::vector<T*> & in) {
        sample_vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    std::vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const sample_vector<T> & in) {
        std::vector<sample_vector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }

    /*
    template<typename V, typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, V & samples)
    {
        return samples[pos*num_channels + ch];
    }
    template<typename V, typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, V & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename V, typename T>
    sample_vector<T> get_left_channel(size_t n, const V& in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename V, typename T>
    sample_vector<T> get_right_channel(size_t n, const V& in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename V, typename T>
    sample_vector<T> get_channel(size_t ch, size_t n, V& in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        for(size_t i = ch; i < n; i+=2) r[x++] = in[i];
        return r;
    }

    template<typename V,typename T>
    void set_left_channel(size_t n, const T* left, V& out) {
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2) out[i] = left[x++];
    }
    template<typename V,typename T>
    void set_right_channel(size_t n, const T* right, V& out) {
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) out[i] = right[x++];
    }
    template<typename V,typename T>
    void set_channel(size_t ch, size_t n, const T* in, V& out) {
        size_t x = 0;
        for(size_t i = ch; i < n; i+=2) out[i] = in[x++];
    }

    template<typename V,typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const V& in) {
        sample_vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename V, typename T>
    std::vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const V& in) {
        std::vector<sample_vector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }
    */

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
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(size_t n, const T* & in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, size_t n, T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        for(size_t i = ch; i < n; i+=2) r[x++] = in[i];
        return r;
    }

    template<typename T>
    void set_left_channel(size_t n, const T* left, T* out) {
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(size_t n, const T* right, T* out) {
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, size_t n, const T* in, T* out) {
        size_t x = 0;
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
    std::vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        std::vector<sample_vector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
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
        if (!isEmpty(v))
        {
            for (int i = 0;i < v.size();i++)
            {
                if (v[i] != 0)
                {
                    return false;
                }
            }

            return true;
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
        // if the vector isn't empty
        if (!isEmpty(v))
        {
            // set the first element as the max
            T maxVal = v[0];

            // then for each subsequent value
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

            int i = k;

            while (i < v.size())
            {
                result.push_back(v[i]);
                i += n;
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
            sample_vector<T> result;

            for (int i = 0;i < v1.size();i++)
            {
                result.push_back(v1[i] + v2[i]);
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
            for (T i = limit1;i < limit2;i += step)
            {
                result.push_back(i);
            }
        }
        else if (step < 0)
        {
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


    template<class T>
    sample_vector<T> cos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const T value) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const sample_vector<T> value) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> exp(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log10(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> exp2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> expm1(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ilogb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log1p(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> logb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const std::vector<int> & x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const int x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const std::vector<long int> & x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const long int x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> sqrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cbrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erf(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erfc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> lgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ceil(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> floor(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> trunc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> round(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
        return r;
    }    
    template<class T>
    std::vector<long int> lround(const sample_vector<T> & v) {
        std::vector<long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
        return r;
    }    
    template<class T>
    std::vector<long long int> llround(const sample_vector<T> & v) {
        std::vector<long long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> nearbyint(const sample_vector<T> & v) {
        std::vector<long long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> remainder(const sample_vector<T> & a, const sample_vector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> copysign(const sample_vector<T> & a, const sample_vector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fdim(const sample_vector<T> & a, const sample_vector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
        return r;
    }    
    #undef fmax
    template<class T>
    sample_vector<T> fmax(const sample_vector<T> & a, const sample_vector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmin(const sample_vector<T> & a, const sample_vector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fma(const sample_vector<T> & a, const sample_vector<T> & b, const sample_vector<T> & c) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fabs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
        return r;
    }    
    
    /*
    template<class T>    
    sample_vector<T> abs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::abs((long la[i]);
        return r;
    } 
    */
       


    template<typename T>
    StereoVector<T> stereo(const sample_vector<T> & left, const sample_vector<T> & right) {
        StereoVector<T> r(left.size());
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
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * std::sin(pan_map);
            r[1][i] = right[i] * std::cos(pan_map);
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
        for(size_t i = 0; i < r.size(); i++) r[i] = a[i]+b[i];
        T max = std::max_element(r.begin(),r.end());
        if(max > 0) for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }
    template<typename T>
    sample_vector<T> normalize(const sample_vector<T> & a) {
        sample_vector<T> r(a);        
        T max = std::max_element(r.begin(),r.end());
        if(max > 0) for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }
    template<class A, class B>
    std::vector<B> convert(const std::vector<A> & v) {
        std::vector<B> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = B(v[i]);
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, T (*f)(T value)) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = f(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, std::function<T (T)> func) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = func(v[i]);
        return r;
    }
    template<class T>
    void inplace_add(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] += func(a[i]);        
    }
    template<class T>
    void inplace_sub(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] -= func(a[i]);        
    }
    template<class T>
    void inplace_mul(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] *= func(a[i]);        
    }
    template<class T>
    void inplace_div(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] /= func(a[i]);
    }

    
    template<class T>
    void fill(sample_vector<T> & in, T x)
    {
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
    void StdMean(T *sig_src_arr, uint32_t blockSize, T * result){
        T sum = T(0);
        uint32_t blkCnt;
        T in1,in2,in3, in4;
        assert(blockSize != 0);
        blkCnt = blockSize>>2U; //Right shifted by 4 so divided by 4

        while(blkCnt > 0){
            in1 = *sig_src_arr++;
            in2 = *sig_src_arr++;
            in3 = *sig_src_arr++;
            in4 = *sig_src_arr++;
            sum += in1;
            sum += in2;
            sum += in3;
            sum += in4;
            blkCnt--;
        }
        blkCnt = blockSize% 0x4;

        while(blkCnt > 0){
            sum += *sig_src_arr++;
            blkCnt--;
        }
        
        *result = sum/(T)blockSize;        
    }
    template<typename T>
    sample_vector<T> StdMean(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdMean(in.data(),in.size(),r.data());
        return r;
    }

    template<typename T>
    void StdRMS(T *pSig_src_arr, uint32_t blockSize, T *pResult)
    {
        T sum = 0.0;
        uint32_t blkCnt;
        T in;
        assert(blockSize != 0);
        blkCnt = blockSize >>2;
        while(blkCnt > 0){
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            blkCnt--;
        }

        blkCnt = blockSize%0x4;
        while(blkCnt>0)
        {
            in = *pSig_src_arr++;
            sum += in*in;
            blkCnt--;
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
    void StdDev(T * pSig_src_arr, uint32_t blockSize, T *pResult)
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

        blkCnt = blockSize>>2;

        while(blkCnt>0){
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

            blkCnt--;
        }

        blkCnt = blockSize % 0x4;

        while(blkCnt>0){
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;

            blkCnt--;
        }

        meanOfSquares = sumOfSquares / ((T)blockSize-1.0);
        mean = sum/(T) blockSize;

        squareOfMean = (mean*mean) * ((T)blockSize/(T)(blockSize-1.0));

        *pResult = sqrt((meanOfSquares-squareOfMean));
    }

    template<typename T>
    sample_vector<T> StdDev(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdDev(in.data(),in.size(),r.data());
        return r;
    }

    template<typename T>
    void StdVariance(T * pSig_src_arr, uint32_t blockSize, T *pResult)
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

        blkCnt = blockSize >>2U;
        while(blkCnt>0){
            in1 = *pInput++;
            in2 = *pInput++;
            in3 = *pInput++;
            in4 = *pInput++;

            sum+= in1;
            sum+= in2;
            sum+= in3;
            sum+= in4;

        blkCnt--;
        }
        blkCnt = blockSize % 0x4;

        while(blkCnt > 0){
            sum += *pInput++;

            blkCnt--;
        }

        fMean = sum/(T) blockSize;
        pInput = pSig_src_arr;
        blkCnt = blockSize>>2;
        while(blkCnt > 0){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;

            blkCnt--;
        }
        blkCnt = blockSize % 0x4;

        while(blkCnt>0){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;

            blkCnt--;
        }

        *pResult = fSum/(T)(blockSize-1.0);
    }

    template<typename T>
    sample_vector<T> StdVariance(sample_vector<T> & in) 
    {
        sample_vector<T> r(in.size());
        zeros(r);
        StdVariation(in.data(),in.size(),r.data());
        return r;
    }

    template<typename A, typename B>
    std::vector<A> vector_cast(std::vector<B> & in) {
        std::vector<A> r(in.size());
        for(size_t i = 0; i < in.size(); i++)
            r[i] = (A)in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> vector_copy(T * ptr, size_t n) {
        sample_vector<T> r(n);
        for(size_t i = 0; i < n; i++)
            r[i] = ptr[i];
        return r;
    }

    ///////////////////////////////////////////////////////////////
    // Utils
    ///////////////////////////////////////////////////////////////
    inline double freq_to_midi(double f) {
        return 12.0*std::log2(f/440.0) + 69;
    }

    inline double midi_to_freq(double m) {
        return std::pow(2.0, (m-69)/12)*440.0;
    }
    inline double cv_to_freq(double cv)
    {
        return std::pow(2,cv);
    }
    inline double semitone(int semi, double f)
    {
        double m = freq_to_midi(f);
        return midi_to_freq(m + semi);
    }
    inline double octave(int octave, double f) {
        double m = freq_to_midi(f);
        return midi_to_freq(m + octave*12);
    }

        ////////////////////////////////////////////////////////////////
    // KFR Complex DFT
    ////////////////////////////////////////////////////////////////
    struct CDFT
    {
        KFR_DFT_PLAN_F32 * plan;
        size_t size;

        CDFT(size_t n) {
            plan = kfr_dft_create_plan_f32(n);
            size = n;
        }
        ~CDFT() {
            if(plan) kfr_dft_delete_plan_f32(plan);
        }

        void forward(std::vector<std::complex<float>> & input, std::vector<std::complex<float>> & output)
        {
            // you cant do anything with it even if you include the c complex.h
            std::vector<kfr_c32> c_in(size*2),c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_get_temp_size_f32(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c32));
            kfr_dft_execute_f32(plan,c_out.data(),c_in.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c32));
        }
        void inverse(std::vector<std::complex<float>> & input, std::vector<std::complex<float>> & output)
        {
            std::vector<kfr_c32> c_in(size*2),c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_get_temp_size_f32(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c32));
            kfr_dft_execute_inverse_f32(plan,c_out.data(),c_in.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c32));
        }
        void dump() 
        {
            kfr_dft_dump_f32(plan);
        }
    };

    struct CDFT64
    {
        KFR_DFT_PLAN_F64 * plan;
        size_t size;

        CDFT64(size_t n) {
            plan = kfr_dft_create_plan_f64(n);
            size = n;
        }
        ~CDFT64() {
            if(plan) kfr_dft_delete_plan_f64(plan);
        }

        void forward(std::vector<std::complex<double>> & input, std::vector<std::complex<double>> & output)
        {
            std::vector<kfr_c64> c_in(size*2),c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_get_temp_size_f64(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c64));
            kfr_dft_execute_f64(plan,c_out.data(),c_in.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c64));
        }
        void inverse(std::vector<std::complex<double>> & input, std::vector<std::complex<double>> & output)
        {
            std::vector<kfr_c64> c_in(size*2),c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_get_temp_size_f64(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c64));
            kfr_dft_execute_inverse_f64(plan,c_out.data(),c_in.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c64));
        }
        void dump() 
        {
            kfr_dft_dump_f64(plan);
        }
    };

    ////////////////////////////////////////////////////////////////
    // Kfr Real DFT
    ////////////////////////////////////////////////////////////////
    struct RDFT
    {
        KFR_DFT_REAL_PLAN_F32 * plan;
        size_t size;

        RDFT(size_t n) {
            plan = kfr_dft_real_create_plan_f32(n,Perm);
            size = n;
        }
        ~RDFT() {
            if(plan) kfr_dft_real_delete_plan_f32(plan);
        }

        void forward(std::vector<float> & input, std::vector<std::complex<float>> & output)
        {
            std::vector<kfr_c32> c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_real_get_temp_size_f32(plan));            
            kfr_dft_real_execute_f32(plan,c_out.data(),input.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c32));
        }
        void inverse(std::vector<std::complex<float>> & input, std::vector<float> & output)
        {
            std::vector<kfr_c32> c_in(size*2);
            std::vector<uint8_t> temp(kfr_dft_real_get_temp_size_f32(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c32));
            kfr_dft_real_execute_inverse_f32(plan,output.data(),c_in.data(),temp.data());            
        }
        void dump() 
        {
            kfr_dft_real_dump_f32(plan);
        }
    };
    struct RDFT64
    {
        KFR_DFT_REAL_PLAN_F64 * plan;
        size_t size;

        RDFT64(size_t n) {
            plan = kfr_dft_real_create_plan_f64(n,Perm);
            size = n;
        }
        ~RDFT64() {
            if(plan) kfr_dft_real_delete_plan_f64(plan);
        }

        void forward(std::vector<double> & input, std::vector<std::complex<double>> & output)
        {
            std::vector<kfr_c64> c_out(size*2);
            std::vector<uint8_t> temp(kfr_dft_real_get_temp_size_f64(plan));            
            kfr_dft_real_execute_f64(plan,c_out.data(),input.data(),temp.data());
            output.resize(size);
            memcpy(output.data(),c_out.data(),size * sizeof(kfr_c64));
        }
        void inverse(std::vector<std::complex<double>> & input, std::vector<double> & output)
        {
            std::vector<kfr_c64> c_in(size*2);
            std::vector<uint8_t> temp(kfr_dft_real_get_temp_size_f64(plan));
            memcpy(c_in.data(),input.data(),input.size() * sizeof(kfr_c64));
            kfr_dft_real_execute_inverse_f64(plan,output.data(),c_in.data(),temp.data());            
        }
        void dump() 
        {
            kfr_dft_real_dump_f64(plan);
        }
    };

    
    ////////////////////////////////////////////////////////////////
    // Kfr FIR
    ////////////////////////////////////////////////////////////////
    struct FIR
    {
        KFR_FILTER_F32 * plan;

        FIR(const float * taps, size_t n) {
            kfr_filter_create_fir_plan_f32(taps,n);
        }
        ~FIR()
        {
            if(plan) kfr_filter_delete_plan_f32(plan);
        }
        void Process(size_t n, float * input, float * output) {
            kfr_filter_process_f32(plan,output,input,n);
        }
    };

    struct FIR64
    {
        KFR_FILTER_F64 * plan;
        FIR64(const double * taps, size_t n) {
            kfr_filter_create_fir_plan_f64(taps,n);
        }
        ~FIR64()
        {
            if(plan) kfr_filter_delete_plan_f64(plan);
        }
        void Process(size_t n, double * input, double * output) {
            kfr_filter_process_f64(plan,output,input,n);
        }
    };

    ////////////////////////////////////////////////////////////////
    // Kfr Convolution
    ////////////////////////////////////////////////////////////////
    struct CONV
    {
        KFR_FILTER_F32 * plan;
        CONV(const float * taps, size_t n, size_t block) {
            kfr_filter_create_convolution_plan_f32(taps,n,block);
        }
        ~CONV()
        {
            if(plan) kfr_filter_delete_plan_f32(plan);
        }
        void Process(size_t n, float * input, float * output) {
            kfr_filter_process_f32(plan,output,input,n);
        }
    };

    struct CONV64
    {
        KFR_FILTER_F64 * plan;
        CONV64(const double * taps, size_t n, size_t block) {
            kfr_filter_create_convolution_plan_f64(taps,n,block);
        }
        ~CONV64()
        {
            if(plan) kfr_filter_delete_plan_f64(plan);
        }
        void Process(size_t n, double * input, double * output) {
            kfr_filter_process_f64(plan,output,input,n);
        }
    };

    ////////////////////////////////////////////////////////////////
    // Kfr IIR
    ////////////////////////////////////////////////////////////////
    struct IIR
    {
        KFR_FILTER_F32 * plan;
        IIR(const float * sos, size_t n) {
            kfr_filter_create_iir_plan_f32(sos,n);
        }
        ~IIR()
        {
            if(plan) kfr_filter_delete_plan_f32(plan);
        }
        void Process(size_t n, float * input, float * output) {
            kfr_filter_process_f32(plan,output,input,n);
        }
    };

    struct IIR64
    {
        KFR_FILTER_F64 * plan;
        IIR64(const double * sos, size_t n) {
            kfr_filter_create_iir_plan_f64(sos,n);
        }
        ~IIR64()
        {
            if(plan) kfr_filter_delete_plan_f64(plan);
        }
        void Process(size_t n, double * input, double * output) {
            kfr_filter_process_f64(plan,output,input,n);
        }
    };

        ///////////////////////////////////////////////////////////////
    // Windows
    ///////////////////////////////////////////////////////////////

    enum WindowType 
    {
        rectangular     = 1,
        triangular      = 2,
        bartlett        = 3,
        cosine          = 4,
        hann            = 5,
        bartlett_hann   = 6,
        hamming         = 7,
        bohman          = 8,
        blackman        = 9,
        blackman_harris = 10,
        kaiser          = 11,
        flattop         = 12,
        gaussian        = 13,
        lanczos         = 14,
    };

    template<typename T>
    kfr::expression_pointer<T> make_window(WindowType type, size_t s, const T x = -1)
    {
        switch(type)
        {
        case hann: return DSP::make_window_hann_ptr<T>(s);
        case hamming: return DSP::make_window_hamming_ptr<T>(s);
        case blackman: return DSP::make_window_blackman_ptr<T>(s);
        case blackman_harris: return DSP::make_window_blackman_ptr<T>(s);
        case gaussian: return DSP::make_window_gaussian_ptr<T>(s);
        case triangular: return DSP::make_window_triangular_ptr<T>(s);
        case bartlett: return DSP::make_window_bartlett_ptr<T>(s);
        case cosine: return DSP::make_window_cosine_ptr<T>(s);
        case bartlett_hann: return DSP::make_window_bartlett_hann_ptr<T>(s);
        case bohman: return DSP::make_window_bohman_ptr<T>(s);
        case lanczos: return DSP::make_window_lanczos_ptr<T>(s);
        case flattop: return DSP::make_window_flattop_ptr<T>(s);
        case rectangular: return DSP::make_window_rectangular_ptr<T>(s);
        case kaiser: return DSP::make_window_kaiser_ptr<T>(s);
        }
    }
    template<typename T>
    sample_vector<T> window(WindowType type,size_t s, const T x = -1) {
        return sample_vector<T>(make_window<T>(type,s));
    }


    ///////////////////////////////////////////////////////////////
    // FIR
    ///////////////////////////////////////////////////////////////

    template<typename T>
    sample_vector<T> fir_lowpass(sample_vector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);        
        return DSP::fir_lowpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    sample_vector<T> fir_highpass(sample_vector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);        
        return DSP::fir_highpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    sample_vector<T> fir_bandpass(sample_vector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);
        return DSP::fir_bandpass(in,num_taps,window,normalize);        
    }
    template<typename T>
    sample_vector<T> fir_bandstop(sample_vector<T> in, size_t num_taps, WindowType type, bool normalize = true)
    {
        kfr::expression_pointer<T> window = make_window<T>(type,num_taps);
        return DSP::fir_bandstop(in,num_taps,window,normalize);        
    }
    ///////////////////////////////////////////////////////////////
    // DC Block
    ///////////////////////////////////////////////////////////////

    template<typename T>
    sample_vector<T> dcremove(T cutoff, sample_vector<T> in)
    {        
        return DSP::dcremove(in,cutoff);     
    }

    ///////////////////////////////////////////////////////////////
    // File I/O
    ///////////////////////////////////////////////////////////////

    template<typename T>
    sample_vector<T> load_wav(const char * filename)
    {
        return DSP::load_wav<T>(filename);        
    }
    template<typename T>
    sample_vector<T> load_mp3(const char * filename)
    {
        return DSP::load_mp3<T>(filename);        
    }
    template<typename T>
    sample_vector<T> load_flac(const char * filename)
    {
        return DSP::load_flac<T>(filename);        
    }
    template<typename T>
    void save_wav(sample_vector<T>  in, const char * filename, size_t channels, int sample_type, double sample_rate, bool use_w64=false)
    {        
        DSP::write_wav(in,filename,channels,sample_type,sample_rate,use_w64);
        
    }

    template<typename T> using WAVFileReader = DSP::WavReader<T>;
    template<typename T> using WAVFileWriter = DSP::WavWriter<T>;
    template<typename T> using MP3FileReader = DSP::MP3Reader<T>;
    template<typename T> using FLACFileReader = DSP::FlacReader<T>;

    ///////////////////////////////////////////////////////////////
    // Plot
    ///////////////////////////////////////////////////////////////
    template<typename T>
    void plot(sample_vector<T> in, std::string& name = "", std::string& options = "") {
        kfr::plot_show(name,in,options);
    }


    enum FilterType
    {
        Lowpass,
        Highpass,
        Bandpass,  
        Bandpass2, // this is used in RBJ for the cszap whatevr
        Notch,
        Bandstop,
        Allpass,
        Peak,
        Lowshelf,
        Highshelf,        
    };

    struct Biquad6DB
    {
        double a[2];
        double b[3];
        double fs,fc;
        double x1,x2,y1,y2;
        double x,y;

        FilterType filterType = Lowpass;

        Biquad6DB(FilterType type, double Fs, double Fc) {
            fs = Fs;
            fc = Fc/Fs;
            setFilter(type);
        }
        void setFilter(FilterType type) {
            filterType = type;
            switch(type) {
                case Lowpass: lowpass(fc); break;
                case Highpass: highpass(fc); break;
                case Allpass: allpass(fc); break;
            }
        }
        void setCutoff(float f) {
            fc = f/fs;
            setFilter(filterType);
        }
        void setQ(float q) {
            // does nothing right now
        }

        void lowpass(double fc)
        {
            double K = std::tan(M_PI*fc);
            b[0] = K/(K+1);
            b[1] = K/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void highpass(double fc)
        {
            double K = std::tan(M_PI*fc);
            b[0] = 1/(K+1);
            b[1] = -1/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void allpass(double fc)
        {
            double K = std::tan(M_PI*fc);
            b[0] = (K-1)/(K+1);
            b[1] = 1;
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        
        double Tick(double I)
        {
            Undenormal denormal;        
            x = I;
            y = b[0]*x + b[1]*x1 - a[0] * y1;
            x1 = x;
            y1 = y;
            return y;
        }
        
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // -12db Two Pole/Two Zero 1 section
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct Biquad12DB
    {
        double a[2];
        double b[3];
        double fs,fc,q,g;
        double x1,x2,y1,y2;
        double x,y;

        FilterType filterType = Lowpass;

        Biquad12DB(FilterType type, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            x1=x2=y1=y2=0;
            filterType = type;
            init_filter(Fc,Q);        
        }
        Biquad12DB(const kfr::biquad_params<double>& bq, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            x1=x2=y1=y2=0;        
            setCoefficients(bq);
        }
        

        void init_filter(double Fc, double Q, double gain=1)
        {
            fc = Fc/fs*0.99;        
            q = Q;
            g = gain;

            switch(filterType)
            {
                case Lowpass: lowpass(fc,q); break;
                case Highpass: highpass(fc,q); break;
                case Bandpass: bandpass(fc,q); break;
                case Notch: notch(fc,q); break;
                // fc/q dont matter q must be 0
                case Allpass: allpass(fc,0); break;
                case Peak: peak(fc,q,gain); break;
                case Lowshelf: lowshelf(fc,q); break;
                case Highshelf: highshelf(fc,q); break;
                default: assert(1==0);
            }
        }

        void setCoefficients(kfr::biquad_params<double> p)
        {
            if(p.a0 == 0) p.a0=1.0f;        
            a[0] = p.a1/p.a0;
            a[1] = p.a2/p.a0;
            b[0] = p.b0/p.a0;
            b[1] = p.b1/p.a0;
            b[2] = p.b2/p.a0;
        }
        void setCutoff(double f) {
            fc = f;
            init_filter(fc,q,g);
        }
        void setQ(double Q) {
            q  = Q;
            init_filter(fc,q,g);
        }
        void setGain(double G) {
            g = G;
            init_filter(fc,q,g);
        }


        void notch(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_notch(fc,q);  
            setCoefficients(p);
        }
        void lowpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_lowpass(fc,q);        
            setCoefficients(p);
        }
        void allpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_allpass(fc,q);        
            setCoefficients(p);
        }
        void highpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_highpass(fc,q);        
            setCoefficients(p);
        }
        void peak(double f, double Q, double gain) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_peak(fc,q, gain);        
            setCoefficients(p);
        }
        void lowshelf(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_lowshelf(fc,q);        
            setCoefficients(p);
        }
        void highshelf(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_highshelf(fc,q);        
            setCoefficients(p);
        }
        void bandpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_bandpass(fc,q);        
            setCoefficients(p);
        }

        double Tick(double I, double A = 1, double X = 0, double Y = 0)
        {
            Undenormal denormal;
            x = I;
            y = b[0]*x + b[1]*x1 + b[2]*x2 - a[0]*y1 - a[1]*y2;
            y2 = y1;
            y1 = y;
            x2 = x1;
            x1 = x;
            return y;
        }
    };




    ////////////////////////////////////////////////////////////////////////////////////////////
    // -24db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad24DB
    {
        Biquad12DB *first,*second;
        double fs,fc;
        double Q1,Q2;

        Biquad24DB(FilterType type, double Fs, double Fc, double G = 1, double Q1 = 0.54119610f, double Q2=1.3065460f)
        {
            fs = Fs;
            fc = Fc;        
            first = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            this->Q1 = Q1;
            this->Q2 = Q2;
        }
        ~Biquad24DB()
        {
            if(first) delete first;
            if(second) delete second;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
        }
        void setQ(double q1,double q2) {
            Q1 = q1;
            Q2 = q2;        
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
        }
        double Tick(double I)   {
            return second->Tick(first->Tick(I));
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////
    // -36db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad36DB
    {
        Biquad12DB *first,*second,*third;
        double fs,fc;
        double Q1,Q2,Q3;

        Biquad36DB(FilterType type, double Fs, double Fc, double G = 1, double Q2=0.70710678, double Q3=1.9318517)
        {
            fs = Fs;
            fc = Fc;
            this->Q1 = Q1;
            this->Q2 = Q2;
            this->Q3 = Q3;        
            first  = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            third  = new Biquad12DB(type,Fs,Fc,G,Q3);
        }
        ~Biquad36DB()
        {
            if(first) delete first;
            if(second) delete second;
            if(third) delete third;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
            third->lowpass(fc,Q3);
        }
        void setQ(double q1,double q2, double q3) {
            Q1 = q1;
            Q2 = q2;
            Q3 = q3;        
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
            third->lowpass(fc,Q3);
        }
        double Tick(double I)   {
            return third->Tick(second->Tick(first->Tick(I)));
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////
    // -48db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad48DB
    {
        Biquad12DB *first,*second,*third,*fourth;
        double fs,fc;
        double Q1,Q2,Q3,Q4;
        Biquad48DB(FilterType type, double Fs, double Fc, double G = 1, double Q1 = 0.50979558, double Q2=0.60134489, double Q3=0.89997622, double Q4=2.5629154)
        {
            fs = Fs;
            fc = Fc;        
            this->Q1 = Q1;
            this->Q2 = Q2;
            this->Q3 = Q3;        
            this->Q4 = Q4;        
            first = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            third = new Biquad12DB(type,Fs,Fc,G,Q3);
            fourth = new Biquad12DB(type,Fs,Fc,G,Q4);
        }
        ~Biquad48DB()
        {
            if(first) delete first;
            if(second) delete second;
            if(third) delete third;
            if(fourth) delete fourth;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
            third->lowpass(fc,Q3);
            fourth->lowpass(fc,Q4);
        }
        void setQ(double q1,double q2, double q3, double q4) {
            Q1 = q1;
            Q2 = q2;
            Q3 = q3;
            Q4 = q4;        
            first->lowpass(fc,Q1);
            second->lowpass(fc,Q2);
            third->lowpass(fc,Q3);
            fourth->lowpass(fc,Q4);
        }
        double Tick(double I)   {
            return fourth->Tick(third->Tick(second->Tick(first->Tick(I))));
        }
    };


        template<typename T>
    struct DelayLine
    {
        kfr::univector<T> delay;
        size_t read_cursor,write_cursor;
        T feedback;
        T mix;
        size_t delayLen;
        T delayTime;

        enum InterpType
        {
            None,
            NearestNeighbor,
            Linear,
            Cubic,
            Spline3,
            Spline5,
            Hermite1,
            Hermite2,
            Hermite3,
            Hermite4,
        }
        interpType;

        DelayLine() =  default;

        DelayLine(T delay_time, T delay_size = 1) {
            delay.resize(sampleRate*delay_size);
            //memset(delay.data(),0,delay.size()*sizeof(T));
            zeros(delay);
            feedback=0.5;
            mix = 0.5;
            delayTime = delay_time;        
            delayLen  = delay_time * sampleRate;        
            read_cursor  = 0;
            write_cursor = delayLen;
            interpType = Linear;
        }
        T& operator[](size_t i) {
            return delay[i % delay.size()];
        }
        
        void setDelaySize(T size) {
            delay.resize(size);
            memset(delay.data(),0,delay.size()*sizeof(T));
        }
        void reset() {
            read_cursor  = 0;
            write_cursor = delayLen;   
        }
        void setDelayTime(T f) 
        {
            delayTime = f;
            delayLen  = f;
            write_cursor = (read_cursor + delayLen) % delayLen;
            read_cursor  = 0;
            write_cursor = delayLen;
        }
        void setFeedback(T f) {
            feedback = f;
        }
        void setMix(T m) {
            mix = std::fmod(m,1);
        }
        void resize(size_t n) {
            delay.resize(n);
        }
        virtual T Tick(T I, T A = 1, T X = 1, T Y = 1) {
            
            T output;
            size_t n = read_cursor;
            delay.ringbuf_read(read_cursor,output);            
            T d1 = A*(I - Y*output*feedback);                
            T f= d1 - std::floor(d1); //(int)d1;        
            output = Interpolate(n,f);
            delay.ringbuf_write(write_cursor, I*Y + feedback*output);        
            write_cursor %= delayLen;
            read_cursor %= delayLen;
            return mix*I + (1.0-mix)*output;
        }
        size_t size() { return delay.size(); }    

        T Interpolate(size_t n, T frac)
        {
            switch(interpType)
            {
                case None: return delay[n];
                case NearestNeighbor: return NearestNeighborInterpolate(n,frac);
                case Linear: return LinearInterpolate(n,frac);
                case Cubic: return CubicInterpolate(n,frac);
                case Hermite1: return Hermite1Interpolate(n,frac);
                case Hermite2: return Hermite2Interpolate(n,frac);
                case Hermite3: return Hermite3Interpolate(n,frac);
                case Hermite4: return Hermite4Interpolate(n,frac);
                case Spline3:  return Spline3Interpolate(n,frac);
                case Spline5:  return Spline5Interpolate(n,frac);
            }
            return delay[read_cursor];
        }
        T NearestNeighborInterpolate(size_t n,T frac)
        {
            int   x  = std::round(frac);
            T x1 = delay[ (n + x) % delayLen];
            return x1;
        }
        T LinearInterpolate(size_t n,T frac)
        {            
            T x1 = delay[n];
            T x2 = delay[ (n+1) % delayLen];
            //T frac = x1 - (int)x1;
            return x1 + ((x2-x1)*frac);
        }    
        // just cubic stuff from musicdsp
        T CubicInterpolate(size_t inpos,T finpos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            T a = (3 * (x0-x1) - xm1 + x2) / 2;
            T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
            T c = (x1 - xm1) / 2;
            return (((a * finpos) + b) * finpos + c) * finpos + x0;
        }
        // just another kind of cubials it might be the same kakaloke really
        inline T Hermite1Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite2Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
            T c2 = y0 - y1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite3Interpolate(size_t inpos, T x)
        {                
                T y0 = delay[(inpos - 1) % delayLen];
                T y1 = delay[inpos + 0];
                T y2  = delay[(inpos + 1) % delayLen];
                T y3  = delay[(inpos + 2) % delayLen];
                // 4-point, 3rd-order Hermite (x-form)
                T c0 = y1;
                T c1 = 0.5f * (y2 - y0);
                T y0my1 = y0 - y1;
                T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
                T c2 = y0my1 + c1 - c3;

                return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite4Interpolate(size_t inpos, T frac_pos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            const T    c     = (x1 - xm1) * 0.5f;
            const T    v     = x0 - x1;
            const T    w     = c + v;
            const T    a     = w + v + (x2 - x0) * 0.5f;
            const T    b_neg = w + a;

            return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
        }
        T Spline3Interpolate(size_t inpos, T x)
        {            
            T L1 = delay[(inpos-1)%delayLen];
            T L0 = delay[inpos];
            T H0 = delay[(inpos + 1) % delayLen];
            T H1 = delay[(inpos + 2) % delayLen];
            return L0 + .5f* x*(H0-L1 + x*(H0 + L0*(-2) + L1 + x*( (H0 - L0)*9 + (L1 - H1)*3 + x*((L0 - H0)*15 + (H1 -  L1)*5 +  x*((H0 - L0)*6 + (L1 - H1)*2 )))));
        }
        T Spline5Interpolate(size_t inpos, T x)
        {
            /* 5-point spline*/
            int nearest_sample = delayLen;

            T p0=delay[(nearest_sample-2) % delayLen];
            T p1=delay[(nearest_sample-1) % delayLen];
            T p2=delay[nearest_sample];
            T p3=delay[(nearest_sample+1) % delayLen];
            T p4=delay[(nearest_sample+2) % delayLen];
            T p5=delay[(nearest_sample+3) % delayLen];

            return p2 + 0.04166666666*x*((p3-p1)*16.0+(p0-p4)*2.0
            + x *((p3+p1)*16.0-p0-p2*30.0- p4
            + x *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
            + x *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
            + x *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
        }
    };


    template<typename T>
    struct MultiTapDelayLine 
    {
        std::vector<size_t> tap_reads;
        size_t taps;
        size_t write_cursor, read_cursor;
        size_t delayLen;
        kfr::univector<T> delay;        
        T feedback;
        T mix;

        enum InterpType
        {
            None,
            NearestNeighbor,
            Linear,
            Cubic,
            Spline3,
            Spline5,
            Hermite1,
            Hermite2,
            Hermite3,
            Hermite4,
        }
        interpType;

        MultiTapDelayLine(T delayTime)
        {            
            delay.resize(sampleRate);
            //memset(delay.data(),0,delay.size()*sizeof(T));
            zeros(delay);
            feedback=0.5;
            mix = 0.5;             
            delayLen   = sampleRate*delayTime;            
            interpType = Linear;
        }
        void addTap(float t) {
            size_t d = t*sampleRate;            
            tap_reads.push_back(d);                    
            taps++;
        }        
        T Tick(T I, T A=1, T X=1, T Y=1) 
        {       
            T output = 0;
            T sum    = 0;            
            size_t read;
            size_t len;
            for(size_t i = 0; i < taps; i++)
            {                            
                read = tap_reads[i];                                
                delay.ringbuf_read(read,output);                        
                T f  = output - std::floor(output); //(int)output;                                        
                T x1 = output;
                read = read  % delayLen;
                T x2 = delay[read];
                output = x1 + (f*(x2-x1));
                tap_reads[i] = read;
                sum += output;                
            }         
            if( taps > 0) sum /= taps;        
                
            T x1 = 0;
            size_t r = read_cursor;
            delay.ringbuf_read(read_cursor,x1);  
            output = x1;
            T f  = output - (int)output;                                                                
            output = 0.5*(output + Interpolate(r,f));
            read_cursor = read_cursor % delayLen;
            delay.ringbuf_write(write_cursor, I + Y*feedback*output);        
            write_cursor %= delayLen;
            return mix*I + (1.0-mix)*(0.5*(sum+output));            
        }
        T Interpolate(size_t n, T frac)
        {
            switch(interpType)
            {
                case None: return delay[n];
                case NearestNeighbor: return NearestNeighborInterpolate(n,frac);
                case Linear: return LinearInterpolate(n,frac);
                case Cubic: return CubicInterpolate(n,frac);
                case Hermite1: return Hermite1Interpolate(n,frac);
                case Hermite2: return Hermite2Interpolate(n,frac);
                case Hermite3: return Hermite3Interpolate(n,frac);
                case Hermite4: return Hermite4Interpolate(n,frac);
                case Spline3:  return Spline3Interpolate(n,frac);
                case Spline5:  return Spline5Interpolate(n,frac);
            }
            return delay[read_cursor];
        }
        T NearestNeighborInterpolate(size_t n,T frac)
        {
            int   x  = std::round(frac);
            T x1 = delay[ (n + x) % delayLen];
            return x1;
        }
        T LinearInterpolate(size_t n,T frac)
        {            
            T x1 = delay[n];
            T x2 = delay[ (n+1) % delayLen];
            //T frac = x1 - (int)x1;
            return x1 + ((x2-x1)*frac);
        }    
        // just cubic stuff from musicdsp
        T CubicInterpolate(size_t inpos,T finpos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            T a = (3 * (x0-x1) - xm1 + x2) / 2;
            T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
            T c = (x1 - xm1) / 2;
            return (((a * finpos) + b) * finpos + c) * finpos + x0;
        }
        // just another kind of cubials it might be the same kakaloke really
        inline T Hermite1Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite2Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
            T c2 = y0 - y1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite3Interpolate(size_t inpos, T x)
        {                
                T y0 = delay[(inpos - 1) % delayLen];
                T y1 = delay[inpos + 0];
                T y2  = delay[(inpos + 1) % delayLen];
                T y3  = delay[(inpos + 2) % delayLen];
                // 4-point, 3rd-order Hermite (x-form)
                T c0 = y1;
                T c1 = 0.5f * (y2 - y0);
                T y0my1 = y0 - y1;
                T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
                T c2 = y0my1 + c1 - c3;

                return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite4Interpolate(size_t inpos, T frac_pos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            const T    c     = (x1 - xm1) * 0.5f;
            const T    v     = x0 - x1;
            const T    w     = c + v;
            const T    a     = w + v + (x2 - x0) * 0.5f;
            const T    b_neg = w + a;

            return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
        }
        T Spline3Interpolate(size_t inpos, T x)
        {            
            T L1 = delay[(inpos-1)%delayLen];
            T L0 = delay[inpos];
            T H0 = delay[(inpos + 1) % delayLen];
            T H1 = delay[(inpos + 2) % delayLen];
            return L0 + .5f* x*(H0-L1 + x*(H0 + L0*(-2) + L1 + x*( (H0 - L0)*9 + (L1 - H1)*3 + x*((L0 - H0)*15 + (H1 -  L1)*5 +  x*((H0 - L0)*6 + (L1 - H1)*2 )))));
        }
        T Spline5Interpolate(size_t inpos, T x)
        {
            /* 5-point spline*/
            int nearest_sample = delayLen;

            T p0=delay[(nearest_sample-2) % delayLen];
            T p1=delay[(nearest_sample-1) % delayLen];
            T p2=delay[nearest_sample];
            T p3=delay[(nearest_sample+1) % delayLen];
            T p4=delay[(nearest_sample+2) % delayLen];
            T p5=delay[(nearest_sample+3) % delayLen];

            return p2 + 0.04166666666*x*((p3-p1)*16.0+(p0-p4)*2.0
            + x *((p3+p1)*16.0-p0-p2*30.0- p4
            + x *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
            + x *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
            + x *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
        }
    };  

    template<typename T>
    struct CombFilter
    {
        DelayLine<T> delay[2];
        T x1,y,y1;
        T gain[2];
        T delayTime[2];
        
        enum {
            X_index,
            Y_index,
        };

        CombFilter(T _g1, T _g2, T _d1, T _d2)
        {
            gain[X_index] = _g1;
            gain[Y_index] = _g2;
            delayTime[0] = _d1;
            delayTime[1] = _d2;
            for(size_t i = 0; i < 2; i++)
            {
                delay[i].setDelaySize(44100);
                delay[i].setDelayTime(delayTime[i]);
            }       
            x1=y=y1=0;
        }    

        // X modulation * depth
        // Y modulation * depth
        T Tick(T I, T A=1, T X = 0, T Y=0)
        {
            T x = I;
            y = x + gain[X_index] * x1 - gain[Y_index] * y1;
            x1 = delay[X_index].Tick(x);
            y1 = delay[Y_index].Tick(y);
            return y;
        }
    };

    template<typename T>
    struct IIRCombFilter
    {
        DelayLine<T> delay;
        float g,x,y,y1;

        IIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = y1 = 0;
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x - g*y1;
            y1= delay.Tick(y);
            return y;
        }
    };

    template<typename T>
    struct FIRCombFilter
    {
        DelayLine<T> delay;
        float g,x,x1,y;

        FIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = x1 = 0;
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x + g*x1;
            x1= delay.Tick(x);
            return y;
        }
    };


    template<typename T>
    struct MultiTapCombFilter
    {
        MultiTapDelayLine<T> delay[2];
        T x1,y,y1;
        T gain[2];
        T delayTime[2];
        
        enum {
            X_index,
            Y_index,
        };

        MultiTapCombFilter(T _g1, T _g2, T _d1, T _d2)
        {
            gain[X_index] = _g1;
            gain[Y_index] = _g2;
            delayTime[0] = _d1;
            delayTime[1] = _d2;
            for(size_t i = 0; i < 2; i++)
            {
                delay[i].setDelaySize(44100);
                delay[i].setDelayTime(delayTime[i]);
            }       
            x1=y=y1=0;
        }    
        void addTap(float t) {
            delay.addTap(t);
        }
        // X modulation * depth
        // Y modulation * depth
        T Tick(T I, T A=1, T X = 0, T Y=0)
        {
            T x = I;
            y = x + gain[X_index] * x1 - gain[Y_index] * y1;
            x1 = delay[X_index].Tick(x);
            y1 = delay[Y_index].Tick(y);
            return y;
        }
    };

    template<typename T>
    struct MultiTapIIRCombFilter
    {
        MultiTapDelayLine<T> delay;
        float g,x,y,y1;

        MultiTapIIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = y1 = 0;
        }
        void addTap(float t) {
            delay.addTap(t);
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x - g*y1;
            y1= delay.Tick(y);
            return y;
        }
    };

    template<typename T>
    struct MultiTapFIRCombFilter
    {
        MultiTapDelayLine<T> delay;
        float g,x,x1,y;

        MultiTapFIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = x1 = 0;
        }
        void addTap(float t) {
            delay.addTap(t);
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x + g*x1;
            x1= delay.Tick(x);
            return y;
        }
    };

///////////////////////////////////////////////////////////////
// Convolution
///////////////////////////////////////////////////////////////

    template<typename T> using ConvolveFilter = DSP::ConvolveFilter<T>;
    
    template<typename T>
    struct ConvolutionFilter
    {
        DSP::ConvolveFilter<T> *filter;
        
        ConvolutionFilter(size_t block_size,sample_vector<T> impulse)
        {
            filter = new DSP::ConvolveFilter<T>(impulse,block_size);       
        }
        ~ConvolutionFilter() {
            if(filter) delete filter;
        }
        sample_vector<T> apply(sample_vector<T> x)
        {        
            filter->apply(x);
            return x;        
        }
        void apply(sample_vector<T> x, sample_vector<T> & out)
        {                
            filter->apply(out,x);        
        }
    };


    //////////////////////////////////////////////////////////////
    // Convolve/Correlation
    ///////////////////////////////////////////////////////////////

    template<typename T>
    sample_vector<T> convolve(sample_vector<T> a, sample_vector<T> b)
    {                
        return kfr::convolve(a,b);        
    }
    template<typename T>
    sample_vector<T> correlate(sample_vector<T> a, sample_vector<T> b)
    {                
        return kfr::correlate(a,b);        
    }
    template<typename T>
    sample_vector<T> autocorrelate(sample_vector<T> a)
    {                
        return kfr::autocorrelate(a);        
    }    


    template<typename T> using CDFTPlan = DSP::DFTPlan<T>;
    template<typename T> using RDFTPlan = DSP::DFTRealPlan<T>;
    template<typename T> using DCTPlan  = DSP::DCTPlan<T>;

    ///////////////////////////////////////////////////////////////
    // DFT/DCT
    ///////////////////////////////////////////////////////////////

    template<typename T>
    complex_vector<T> dft_forward(complex_vector<T> input)
    {        
        return DSP::run_dft(input);        
    }
    template<typename T>
    complex_vector<T> dft_inverse(complex_vector<T> input)
    {        
        return  DSP::run_idft(input);     
    }
    template<typename T>
    complex_vector<T> real_dft_forward(sample_vector<T> input)
    {        
        return run_realdft(input);     
    }
    template<typename T>
    sample_vector<T> real_dft_inverse(complex_vector<T> in)
    {        
        return run_irealdft(in);        
    }

    template<typename T>
    sample_vector<T> dct_forward(sample_vector<T> in)
    {        
        DSP::DCTPlan<T> dct(in.size());
        kfr::univector<T> out(in.size());
        dct.execute(out,in,false);
        return out;
    }
    template<typename T>
    sample_vector<T> dct_inverse(sample_vector<T> in)
    {        
        DSP::DCTPlan<T> dct(in.size());
        kfr::univector<T> out(in.size());
        dct.execute(out,in,true);
        return out;
    }

        template<typename T>
    struct FIRFilter {
    private:
        kfr::filter_fir<SampleType> * filter;
        kfr::univector<T> taps;
        
    public:
        
        // need to be able to input coefficients
        FIRFilter(size_t num_taps) { 
            taps.resize(num_taps); 
            filter = nullptr;
        }
        FIRFilter(const kfr::univector<T> & taps) {
            filter = new kfr::filter_fir<T>(taps);
        }
        ~FIRFilter() {
            if(filter != NULL) delete filter;
        }
        void set_taps(const kfr::univector<T> & taps) {
            filter = new kfr::filter_fir<T>(taps);
        }
        void bandpass(T x, T y, kfr::expression_pointer<T> & window, bool normalize=true ) {        
            kfr::fir_bandpass(taps,x,y,window,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void bandstop(T x, T y, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_bandstop(taps,x,y, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void highpass(T cutoff, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_highpass(taps,cutoff, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        void lowpass(T cutoff, kfr::expression_pointer<T> & window_type, bool normalize=true ) {        
            kfr::fir_lowpass(taps,cutoff, window_type,normalize);
            filter = new kfr::filter_fir<T>(taps);
        }
        
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }

    };

    template<typename T>
    struct FIRBandpassFilter
    {
        FIRFilter<T> * filter;

        FIRBandpassFilter(size_t num_taps, T x, T y, kfr::expression_pointer<T> & window, bool normalize = true) {
            filter = new FIRFilter<T>(num_taps);
            assert(filter != NULL);
            filter->bandpass(x,y,window,normalize);
        }
        ~FIRBandpassFilter() {
            if(filter) delete filter;
        }
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }
    };

    template<typename T>
    struct FIRBandstopFilter
    {
        FIRFilter<T> * filter;

        FIRBandstopFilter(size_t num_taps, T x, T y, kfr::expression_pointer<T> & window, bool normalize = true) {
            filter = new FIRFilter<T>(num_taps);
            assert(filter != NULL);
            filter->bandstop(x,y,window,normalize);
        }
        ~FIRBandstopFilter() {
            if(filter) delete filter;
        }
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }
    };

    template<typename T>
    struct FIRLowpassFilter
    {
        FIRFilter<T> * filter;

        FIRLowpassFilter(size_t num_taps, T x, kfr::expression_pointer<T> & window, bool normalize = true) {
            filter = new FIRFilter<T>(num_taps);
            assert(filter != NULL);
            filter->lowpass(x,window,normalize);
        }
        ~FIRLowpassFilter() {
            if(filter) delete filter;
        }
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }
    };

    template<typename T>
    struct FIRHighpassFilter
    {
        FIRFilter<T> * filter;

        FIRHighpassFilter(size_t num_taps, T x, kfr::expression_pointer<T> & window, bool normalize = true) {
            filter = new FIRFilter<T>(num_taps);
            assert(filter != NULL);
            filter->highpass(x,window,normalize);
        }
        ~FIRHighpassFilter() {
            if(filter) delete filter;
        }
        void apply(kfr::univector<T> & data) {
            filter->apply(data);
        }
        void apply(kfr::univector<T> & out, const kfr::univector<T> & in) {
            filter->apply(out,in);
        }
        void reset() { filter->reset(); }
    };


        template<typename T>
    sample_vector<T> sinewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::sinewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    sample_vector<T> squarewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::squarewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    sample_vector<T> trianglewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::trianglewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    sample_vector<T> sawtoothewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::squarewave(n,freq,sample_rate,phase);        
    }    

    template<typename T>
    struct FunctionGenerator
    {
        T phase,inc,f,sr;

        FunctionGenerator(T Fs = 44100.0f)
        {
            phase = 0;
            f     = 440.0f;
            sr    = Fs;
            inc   = f/sr;            
        }
        void Inc()
        {
            phase = fmod(phase + inc, 2*M_PI);
        }
        void rawsine()
        {
            return kfr::rawsine(phase);
        }
        void sine() {
            return kfr::sine(phase);
        }
        void sinenorm() {
            return kfr::sinenorm(phase);
        }
        void rawsquare() {
            return kfr::rawsquare(phase);
        }
        void square() {
            return kfr::square(phase);
        }
        void squarenorm() {
            return kfr::squarenorm(phase);
        }
        void rawtriangle() {
            return kfr::rawtriangle(phase);
        }
        void triangle() {
            return kfr::triangle(phase);
        }
        void trianglenorm() {
            return kfr::trianglenorm(phase);
        }
        void rawsawtooth() {
            return kfr::rawsawtooth(phase);
        }
        void sawtooth() {
            return kfr::sawtooth(phase);
        }
        void sawtoothnorm() {
            return kfr::sawtoothnorm(phase);
        }
        void isawtooth() {
            return kfr::isawtooth(phase);
        }
        void isawtoothnorm() {
            return kfr::isawtoothnorm(phase);
        }
    };


       struct BesselFilter
    {
        std::vector<Biquad12DB*> filters;
        
        double fc,fs;
        double low,high;
        int order;
        FilterType filterType;

        BesselFilter(FilterType type, int Order, double Fs, double Fc)
        {
            order = Order;
            fs    = Fs;
            fc    = Fc;
            filterType = type;
            low   = 0;
            high    = Fs/2;
            initFilter();
            
        }
        void initFilter()
        {
            switch(filterType)
            {
            case Lowpass:  lowpass(fc,fs); break;
            case Highpass: highpass(fc,fs); break;
            case Bandpass: bandpass(low,high,fs); break;
            case Bandstop: bandstop(low,high,fs); break;
            }

        }
        void setCutoff(double f) {
            fc = f;
            low = f;
            switch(filterType)
            {
                case Lowpass:  doLowpassCutoff(fc,fs); break;
                case Highpass: doHighpassCutoff(fc,fc); break;
                case Bandpass: doBandpassCutoff(low,high,fs); break;
                case Bandstop: doBandstopCutoff(low,high,fs); break;
            }
        }
        void setCutoff(double lo,double hi) {
            low = lo;
            high = hi;
            switch(filterType)
            {
                case Bandpass: doBandpassCutoff(lo,hi,fc); break;
                case Bandstop: doBandstopCutoff(lo,hi,fc); break;
            }
        }

        void doLowpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::bessel<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void lowpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::bessel<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doHighpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::bessel<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void highpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::bessel<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandpassCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::bessel<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandpass(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::bessel<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandstopCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::bessel<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandstop(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::bessel<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            double R = I;
            for(typename std::vector<Biquad12DB*>::reverse_iterator i = filters.rbegin(); i != filters.rend(); i++)
            {
                R = (*i)->Tick(R);
            }
            return R;
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////
    // Butterworth
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct ButterworthFilter
    {
        std::vector<Biquad12DB*> filters;
        
        double fc,fs;
        double low,high;
        int order;
        FilterType filterType;

        ButterworthFilter(FilterType type, int Order, double Fs, double Fc)
        {
            order = Order;
            fs    = Fs;
            fc    = Fc;
            filterType = type;
            low   = 0;
            high    = Fs/2;
            initFilter();
            
        }
        void initFilter()
        {
            switch(filterType)
            {
            case Lowpass:  lowpass(fc,fs); break;
            case Highpass: highpass(fc,fs); break;
            case Bandpass: bandpass(low,high,fs); break;
            case Bandstop: bandstop(low,high,fs); break;
            }

        }
        void setCutoff(double f) {
            fc = f;
            low = f;
            switch(filterType)
            {
                case Lowpass:  doLowpassCutoff(fc,fs); break;
                case Highpass: doHighpassCutoff(fc,fc); break;
                case Bandpass: doBandpassCutoff(low,high,fs); break;
                case Bandstop: doBandstopCutoff(low,high,fs); break;
            }
        }
        void setCutoff(double lo,double hi) {
            low = lo;
            high = hi;
            switch(filterType)
            {
                case Bandpass: doBandpassCutoff(lo,hi,fc); break;
                case Bandstop: doBandstopCutoff(lo,hi,fc); break;
            }
        }

        void doLowpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::butterworth<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void lowpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::butterworth<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doHighpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::butterworth<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void highpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::butterworth<double>(order),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandpassCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::butterworth<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandpass(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::butterworth<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandstopCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::butterworth<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandstop(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::butterworth<double>(order),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            double R = I;
            for(typename std::vector<Biquad12DB*>::reverse_iterator i = filters.rbegin(); i != filters.rend(); i++)
            {
                R = (*i)->Tick(R);
            }
            return R;
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Chebyshev1
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Chebyshev1Filter
    {
        std::vector<Biquad12DB*> filters;
        
        double fc,fs,w;
        double low,high;
        int order;
        FilterType filterType;

        Chebyshev1Filter(FilterType type, int Order, double Fs, double Fc, float W=2)
        {
            order = Order;
            fs    = Fs;
            fc    = Fc;
            filterType = type;
            low   = 0;
            high  = Fs/2;
            W     = w;
            initFilter();
            
        }
        void initFilter()
        {
            switch(filterType)
            {
            case Lowpass:  lowpass(fc,fs); break;
            case Highpass: highpass(fc,fs); break;
            case Bandpass: bandpass(low,high,fs); break;
            case Bandstop: bandstop(low,high,fs); break;
            }

        }
        void setCutoff(double f) {
            fc = f;
            low = f;
            switch(filterType)
            {
                case Lowpass:  doLowpassCutoff(fc,fs); break;
                case Highpass: doHighpassCutoff(fc,fc); break;
                case Bandpass: doBandpassCutoff(low,high,fs); break;
                case Bandstop: doBandstopCutoff(low,high,fs); break;
            }
        }
        void setCutoff(double lo,double hi) {
            low = lo;
            high = hi;
            switch(filterType)
            {
                case Bandpass: doBandpassCutoff(lo,hi,fc); break;
                case Bandstop: doBandstopCutoff(lo,hi,fc); break;
            }
        }

        void doLowpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev1<double>(order,w),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void lowpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev1<double>(order,w),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doHighpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev1<double>(order,w),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void highpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev1<double>(order,w),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandpassCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev1<double>(order,w),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandpass(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev1<double>(order,w),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandstopCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev1<double>(order,w),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandstop(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev1<double>(order,w),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            double R = I;
            for(typename std::vector<Biquad12DB*>::reverse_iterator i = filters.rbegin(); i != filters.rend(); i++)
            {
                R = (*i)->Tick(R);
            }
            return R;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Chebyshev2
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Chebyshev2Filter
    {
        std::vector<Biquad12DB*> filters;
        
        double fc,fs,W;
        double low,high;
        int order;
        
        FilterType filterType;

        Chebyshev2Filter(FilterType type, int Order, double Fs, double Fc, float w=80)
        {
            order = Order;
            fs    = Fs;
            fc    = Fc;
            filterType = type;
            low   = 0;
            high  = Fs/2;
            W     = w;
            initFilter();
            
        }
        void initFilter()
        {
            switch(filterType)
            {
            case Lowpass:  lowpass(fc,fs); break;
            case Highpass: highpass(fc,fs); break;
            case Bandpass: bandpass(low,high,fs); break;
            case Bandstop: bandstop(low,high,fs); break;
            }

        }
        void setCutoff(double f) {
            fc = f;
            low = f;
            switch(filterType)
            {
                case Lowpass:  doLowpassCutoff(fc,fs); break;
                case Highpass: doHighpassCutoff(fc,fc); break;
                case Bandpass: doBandpassCutoff(low,high,fs); break;
                case Bandstop: doBandstopCutoff(low,high,fs); break;
            }
        }
        void setCutoff(double lo,double hi) {
            low = lo;
            high = hi;
            switch(filterType)
            {
                case Bandpass: doBandpassCutoff(lo,hi,fc); break;
                case Bandstop: doBandstopCutoff(lo,hi,fc); break;
            }
        }

        void doLowpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void lowpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_lowpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doHighpassCutoff(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void highpass(double cutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_highpass(kfr::chebyshev2<double>(order,W),cutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandpassCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandpass(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandpass(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        void doBandstopCutoff(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters[i]->setCoefficients(temp);
            }            
        }
        void bandstop(double locutoff, double hicutoff, double sample_rate) {
            kfr::zpk<double> filt = kfr::iir_bandstop(kfr::chebyshev2<double>(order,W),locutoff,hicutoff,sample_rate);        
            std::vector<kfr::biquad_params<double>> bqs = kfr::to_sos<double>(filt);
            for(size_t i = 0; i < bqs.size(); i++)
            {
                kfr::biquad_params<double> temp = bqs[i];
                filters.push_back(new Biquad12DB(temp,fs,fc));
            }            
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            double R = I;
            for(typename std::vector<Biquad12DB*>::reverse_iterator i = filters.rbegin(); i != filters.rend(); i++)
            {
                R = (*i)->Tick(R);
            }
            return R;
        }
    };


     double WhiteNoise() {
      return noise.rand();
  }

  class PinkNoise {
  public:
    PinkNoise() {
    srand ( time(NULL) ); // initialize random generator
      clear();
    }

    void clear() {
      for( size_t i=0; i< PINK_NOISE_NUM_STAGES; i++ )
        state[ i ] = 0.0;
      }

    double tick() {
      static const double RMI2 = 2.0 / double(RAND_MAX); // + 1.0; // change for range [0,1)
      static const double offset = A[0] + A[1] + A[2];

    // unrolled loop
      double temp = double( noise.rand() );
      state[0] = P[0] * (state[0] - temp) + temp;
      temp = double( noise.rand() );
      state[1] = P[1] * (state[1] - temp) + temp;
      temp = double( noise.rand() );
      state[2] = P[2] * (state[2] - temp) + temp;
      return ( A[0]*state[0] + A[1]*state[1] + A[2]*state[2] )*RMI2 - offset;
    }

  protected:
    double state[ PINK_NOISE_NUM_STAGES ];
    static const double A[ PINK_NOISE_NUM_STAGES ];
    static const double P[ PINK_NOISE_NUM_STAGES ];
  };

  const double PinkNoise::A[] = { 0.02109238, 0.07113478, 0.68873558 }; // rescaled by (1+P)/(1-P)
  const double PinkNoise::P[] = { 0.3190,  0.7756,  0.9613  };

  double Pink() {
      static PinkNoise pink;
      return pink.tick();
  }
  double GaussianWhiteNoise()
  {
      double R1 = noise.rand();
      double R2 = noise.rand();

      return (double) std::sqrt( -2.0f * std::log( R1 )) * std::cos( 2.0f * M_PI * R2 );
  }
  double GaussRand (int m, double s)
  {
    static int pass = 0;
    static double y2;
    double x1, x2, w, y1;

    if (pass)
    {
        y1 = y2;
    } else  {
        do {
          x1 = 2.0f * noise.rand () - 1.0f;
          x2 = 2.0f * noise.rand () - 1.0f;
          w = x1 * x1 + x2 * x2;
        } while (w >= 1.0f);

        w = (double)std::sqrt (-2.0 * std::log (w) / w);
        y1 = x1 * w;
        y2 = x2 * w;
    }
    pass = !pass;

    return ( (y1 * s + (double) m));
  }

  // +/-0.05dB above 9.2Hz @ 44,100Hz
  class PinkingFilter
  {
    double b0, b1, b2, b3, b4, b5, b6;
  public:
    PinkingFilter() : b0(0), b1(0), b2(0), b3(0), b4(0), b5(0), b6(0) {}
    double process(const double s)
    {
      b0 = 0.99886 * b0 + s * 0.0555179;
      b1 = 0.99332 * b1 + s * 0.0750759;
      b2 = 0.96900 * b2 + s * 0.1538520;
      b3 = 0.86650 * b3 + s * 0.3104856;
      b4 = 0.55000 * b4 + s * 0.5329522;
      b5 = -0.7616 * b5 - s * 0.0168980;
      const double pink = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + (s * 0.5362)) * 0.11;
      b6 = s * 0.115926;
      return pink;
    }
  };

  class BrowningFilter
  {
  double l;
  public:
    BrowningFilter() : l(0) {}
    double process(const double s)
    {
      double brown = (l + (0.02f * s)) / 1.02f;
      l = brown;
      return brown * 3.5f; // compensate for gain
    }
  };


  double PinkNoise()
  {
      static PinkingFilter pink;
      return pink.process(noise.rand());
  }
  double BrownNoise()
  {
      static BrowningFilter brown;
      return brown.process(noise.rand());
  }


      ////////////////////////////////////////////////////////////////////////////////////////////
    // One Pole/One Zero
    // Mainly used for smoothing or blocking DC or Nyquuist
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct OnePole
    {
        double b0,a1;
        double y1;
        double x,y;
        double fc,fs;

        enum Type
        {
            Lowpass,
            Highpass
        }
        filterType = Lowpass;

        OnePole(double f, double sr)
        {
            y1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            a1 = -std::exp(-2*M_PI*fc);
            b0 = 1.0 + a1;
            
        }
        void setCutoff(double f) {
            fc = f/fs;
            if(filterType == Lowpass)
                setLowpass(f);
            else
                setHighpass(f);

        }   
        inline void setLowpass(double Fc) {
            a1 = std::exp(-2*M_PI*fc);
            b0 = 1.0 - a1;
        }

        inline void setHighpass(double Fc) {
            a1 = std::exp(-2.0 * M_PI * (0.5 - Fc));
            b0 = 1.0 - a1;
        }

        double Tick(double in) {
            x = in;
            y = b0*x - a1*y1;
            y1= y;
            return y;

        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // slightly resonant low/high pass
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct TwoPoles
    {
        double b0,a1,a2;
        double y1,y2;
        double x,y;
        double fc,fs;
        float  R;
        enum Type
        {
            Lowpass,
            Highpass,        
        }
        filterType = Lowpass;

        TwoPoles(double sr,double f)
        {
            y1 = 0.0f;   
            y2 = 0.0f;     
            fc = f/sr;
            fs = sr;        
            R  = 0.0f;
            setCutoff(f);        
        }
        // this makes the filter
        void setCutoff(double f) {
            fc = f/fs;
            if(filterType == Lowpass)
                setLowpass();
            else if(filterType == Highpass)
                setHighpass();
        }           
        void setQ(double Q) {
            if(R == 0.0f) R = 0.005f;            
            R = Q;
        }
        inline void setLowpass() {
            /* experiment
            double K = std::tan(M_PI*fc);
            double bottom = (K*K*q+K+q);
            a1 = (2*q*(K*K-1))/bottom;
            a2 = (K*K*q-K+q)/bottom;
            b0 = (K*K*q)/bottom;
            */
            b0 = (1+a1)*(1+a1);
            a1 = -2*R*std::cos(2*M_PI*fc);
            a2 = R*R;
            
        }
        inline void setHighpass() {
            /* experiment
            double K = std::tan(M_PI*fc);
            double bottom = (K*K*q+K+q);
            a1 = (2*q*(K*K-1))/bottom;
            a2 = (K*K*q-K+q)/bottom;
            b0 = (K*K*q)/bottom;
            */
            b0 = (1-a1)*(1-a1);
            a1 = -2*R*std::cos(2*M_PI*(0.5-fc));
            a2 = R*R;
            
        }

        

        double Tick(double in) {                
            x = in;
            y = b0*x - a1*y1 - a2*y2;
            y2 = y1;
            y1= y;                
            return y;
        }
    };


    	////////////////////////////////////////////////////////////////////////////////////////////
	// Rbj
	////////////////////////////////////////////////////////////////////////////////////////////
	struct RbjFilter
	{

		FilterType filter_type;

		double freq,sr,Q,dbGain;
		bool bandwidth;
		// filter coeffs
		double b0a0,b1a0,b2a0,a1a0,a2a0;

		// in/out history
		double ou1,ou2,in1,in2;

		RbjFilter(FilterType type, double sample_rate, double frequency, double q=0.707, double db_gain = 1.0f, bool q_is_bandwidth=false)
		{		
			b0a0=b1a0=b2a0=a1a0=a2a0=0.0;			
			ou1=ou2=in1=in2=0.0f;	
			filter_type = type;
			calc_filter_coeffs(frequency,sample_rate,q,db_gain,q_is_bandwidth);
		};

		double filter(double in0)
		{
			// filter
			Undenormal denormal;
			double const yn = b0a0*in0 + b1a0*in1 + b2a0*in2 - a1a0*ou1 - a2a0*ou2;

			// push in/out buffers
			in2=in1;
			in1=in0;
			ou2=ou1;
			ou1=yn;

			// return output
			return yn;
		};

		double Tick(double I, double A = 1, double X = 0, double Y = 0) {
			double cut = freq;
			double r   = Q;
			Undenormal denormal;
			calc_filter_coeffs(freq + 0.5*X*22050.0/12, sr, Q+0.5*Y, dbGain,bandwidth);
			double x   = std::tanh(A*filter(A*I));
			calc_filter_coeffs(cut, sr, r, dbGain,bandwidth);
			return x;
		}


		void calc_filter_coeffs(double frequency,double sample_rate,double q,double  db_gain,bool q_is_bandwidth)
		{
			// temp pi
			double const temp_pi=3.1415926535897932384626433832795;
			
			if(frequency < 0) frequency = 0;
			if(frequency > (sample_rate/2-1)) frequency = sample_rate/2-1;
			if(q < 0) q = 0;
			if(db_gain < 0) db_gain = 0;

			// temp coef vars
			double alpha,a0,a1,a2,b0,b1,b2;
			freq = frequency;
			sr   = sample_rate;
			Q    = q;
			dbGain = db_gain;
			bandwidth = q_is_bandwidth;
			
				
			// peaking, lowshelf and hishelf
			if(filter_type == Peak || filter_type == Lowshelf || filter_type == Highshelf)
			{
				double const A		=	std::pow(10.0,(db_gain/40.0));
				double const omega	=	2.0*temp_pi*frequency/sample_rate;
				double const tsin	=	std::sin(omega);
				double const tcos	=	std::cos(omega);
				
				if(q_is_bandwidth)
				alpha=tsin*std::sinh(std::log(2.0)/2.0*q*omega/tsin);
				else
				alpha=tsin/(2.0*q);

				double const beta	=	std::sqrt(A)/q;
				
				// peaking
				if(filter_type==Peak)
				{
					b0=double(1.0+alpha*A);
					b1=double(-2.0*tcos);
					b2=double(1.0-alpha*A);
					a0=double(1.0+alpha/A);
					a1=double(-2.0*tcos);
					a2=double(1.0-alpha/A);
				}
				
				// lowshelf
				if(filter_type== Lowshelf)
				{
					b0=double(A*((A+1.0)-(A-1.0)*tcos+beta*tsin));
					b1=double(2.0*A*((A-1.0)-(A+1.0)*tcos));
					b2=double(A*((A+1.0)-(A-1.0)*tcos-beta*tsin));
					a0=double((A+1.0)+(A-1.0)*tcos+beta*tsin);
					a1=double(-2.0*((A-1.0)+(A+1.0)*tcos));
					a2=double((A+1.0)+(A-1.0)*tcos-beta*tsin);
				}

				// hishelf
				if(filter_type==Highshelf)
				{
					b0=double(A*((A+1.0)+(A-1.0)*tcos+beta*tsin));
					b1=double(-2.0*A*((A-1.0)+(A+1.0)*tcos));
					b2=double(A*((A+1.0)+(A-1.0)*tcos-beta*tsin));
					a0=double((A+1.0)-(A-1.0)*tcos+beta*tsin);
					a1=double(2.0*((A-1.0)-(A+1.0)*tcos));
					a2=double((A+1.0)-(A-1.0)*tcos-beta*tsin);
				}
			}
			else
			{
				// other filters
				double const omega	=	2.0*temp_pi*frequency/sample_rate;
				double const tsin	=	std::sin(omega);
				double const tcos	=	std::cos(omega);

				if(q_is_bandwidth)
				alpha=tsin*std::sinh(std::log(2.0)/2.0*q*omega/tsin);
				else
				alpha=tsin/(2.0*q);

				
				// lowpass
				if(filter_type==Lowpass)
				{
					b0=(1.0-tcos)/2.0;
					b1=1.0-tcos;
					b2=(1.0-tcos)/2.0;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// hipass
				if(filter_type==Highpass)
				{
					b0=(1.0+tcos)/2.0;
					b1=-(1.0+tcos);
					b2=(1.0+tcos)/2.0;
					a0=1.0+ alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// bandpass csg
				if(filter_type==Bandpass)
				{
					b0=tsin/2.0;
					b1=0.0;
					b2=-tsin/2;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// bandpass czpg
				if(filter_type==Bandpass2)
				{
					b0=alpha;
					b1=0.0;
					b2=-alpha;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// notch
				if(filter_type==Notch || filter_type == Bandstop)
				{
					b0=1.0;
					b1=-2.0*tcos;
					b2=1.0;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}

				// allpass
				if(filter_type==Allpass)
				{
					b0=1.0-alpha;
					b1=-2.0*tcos;
					b2=1.0+alpha;
					a0=1.0+alpha;
					a1=-2.0*tcos;
					a2=1.0-alpha;
				}
			}

			// set filter coeffs
			b0a0=double(b0/a0);
			b1a0=double(b1/a0);
			b2a0=double(b2/a0);
			a1a0=double(a1/a0);
			a2a0=double(a2/a0);
		};

		
	};


    ///////////////////////////////////////////////////////////////
    // Interpolation
    ///////////////////////////////////////////////////////////////

    template<typename T>
    // r = frac
    // x = [i]
    // y = [i+1]
    T linear_interpolate(T x, T y, T r)
    {        
        return r +*x (1.0-r)*y;
        
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
    class Decimator5
    {
    private:
    float R1,R2,R3,R4,R5;
    const float h0;
    const float h1;
    const float h3;
    const float h5;
    public:
    
    Decimator5():h0(346/692.0f),h1(208/692.0f),h3(-44/692.0f),h5(9/692.0f)
    {
        R1=R2=R3=R4=R5=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R6=R5+h5x0;
        R5=R4+h3x0;
        R4=R3+h1x0;
        R3=R2+h1x0+h0*x1;
        R2=R1+h3x0;
        R1=h5x0;
        return R6;
    }
    };


    template<typename T>
    class Decimator7
    {
    private:
    float R1,R2,R3,R4,R5,R6,R7;
    const float h0,h1,h3,h5,h7;
    public:
    Decimator7():h0(802/1604.0f),h1(490/1604.0f),h3(-116/1604.0f),h5(33/1604.0f),h7(-6/1604.0f)
    {
        R1=R2=R3=R4=R5=R6=R7=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h7x0=h7*x0;
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R8=R7+h7x0;
        R7=R6+h5x0;
        R6=R5+h3x0;
        R5=R4+h1x0;
        R4=R3+h1x0+h0*x1;
        R3=R2+h3x0;
        R2=R1+h5x0;
        R1=h7x0;
        return R8;
    }
    };

    template<typename T>
    class Decimator9
    {
    private:
    float R1,R2,R3,R4,R5,R6,R7,R8,R9;
    const float h0,h1,h3,h5,h7,h9;
    public:
    Decimator9():h0(8192/16384.0f),h1(5042/16384.0f),h3(-1277/16384.0f),h5(429/16384.0f),h7(-116/16384.0f),h9(18/16384.0f)
    {
        R1=R2=R3=R4=R5=R6=R7=R8=R9=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h9x0=h9*x0;
        float h7x0=h7*x0;
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R10=R9+h9x0;
        R9=R8+h7x0;
        R8=R7+h5x0;
        R7=R6+h3x0;
        R6=R5+h1x0;
        R5=R4+h1x0+h0*x1;
        R4=R3+h3x0;
        R3=R2+h5x0;
        R2=R1+h7x0;
        R1=h9x0;
        return R10;
    }
    };


    ///////////////////////////////////////////////////////////////
    // Resampler/Upsample/Downsample
    ///////////////////////////////////////////////////////////////

    // this is interpolator/decimator
    template<typename T>
    struct Resampler
    {
        kfr::samplerate_converter<T> *resampler;

        Resampler(double insr, double outsr)
        {
            resampler = new kfr::samplerate_converter<T>(kfr::resample_quality::normal,outsr,insr);
        }
        ~Resampler() {
            if(resampler) delete resampler;
        }
        void Process(sample_vector<T> & out, sample_vector<T> & in) {
            resampler->process(out,in);        
        }
    };

    template<typename T>
    sample_vector<T> upsample2x(sample_vector<T> in)
    {
        sample_vector<T> out(in.size()*2);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*2] = in[i];
        return out;
    }
    template<typename T>
    sample_vector<T> upsample4x(sample_vector<T> in)
    {
        sample_vector<T> out(in.size()*4);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*4] = in[i];
        return out;
    }
    template<typename T>
    sample_vector<T> upsample2N(size_t n, sample_vector<T> in)
    {
        sample_vector<T> out(in.size()*2*n);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*2*n] = in[i];
        return out;
    }
    template<typename T>
    sample_vector<T> downsample2x(sample_vector<T> in)
    {
        sample_vector<T> out(in.size()/2);
        zeros(out);
        for(size_t i = 0; i < in.size()/2; i++)
            out[i] = in[i*2];
        return out;
    }
    template<typename T>
    sample_vector<T> downsample4x(sample_vector<T> in)
    {
        sample_vector<T> out(in.size()/4);
        zeros(out);
        for(size_t i = 0; i < in.size()/4; i++)
            out[i] = in[i*4];
        return out;
    }
    template<typename T>
    sample_vector<T> downsample2N(size_t n, sample_vector<T> in)
    {
        sample_vector<T> out(in.size()/(2*n));
        zeros(out);
        for(size_t i = 0; i < in.size()/(2*n); i++)
            out[i] = in[i*2*n];
        return out;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////
    // a kind of low/high pass effect
    // It can be calculaated with windowed sinc as it is a FIR filter with only 1 coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct OneZero
    {
        double b0,b1;
        double x1;
        double x,y;
        double fc,fs;
        double hp;

        enum Type
        {
            Lowpass,
            Highpass
        }
        filterType = Lowpass;

        OneZero(double f, double sr)
        {
            x1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            setLowpass(f);
        }
        void setCutoff(double f) {
            fc      = f/fs;
            double wc= 2*M_PI*fc;
            double K = std::tan(wc/2);        
            b0 = K/(K+1);
            b1 = K/(K+1);

        }
        // low pass
        inline void setLowpass(double Fc) {
            double wc      = M_PI*fc;
            double K = std::tan(wc);        
            b0 = K/(K+1);
            b1 = K/(K+1);
        }
    
        double Tick(double in) {       
            double wc= M_PI*fc;
            double K = std::tan(wc);        
            b0 =  K/(K+1);
            b1 =  K/(K+1);       
            x = in;
            y = b0*x + b1*x1;
            x1= x;        
            hp = 1 - y;
            return y;

        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // it's used to create a peculiar notch
    // It can be calculaated with windowed sinc as it is a FIR filter with only 2 coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct TwoZeros
    {
        double b0,b1,b2;
        double x1,x2;
        double x,y;
        double fc,fs,R;
        double hp;

        
        TwoZeros(double sr, double f, float r)
        {
            x1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            R  = r;
            setCutoff(f);
        }
        // this doesn't really do much
        void setCutoff(double f) {
            fc = f/fs;
            b0 = 1;
            b1 = -2*R*std::cos(2*M_PI*fc);
            b2 = R*R;
        }
        // this creates the notch
        void setQ(double q) 
        {
            R  = q;
            if(R == 0.0f) R = 0.005f;
            if(R == 1.0f) R = 0.995f;
            b0 = 1;
            b1 = -2*R*std::cos(2*M_PI*fc);
            b2 = R*R;                
        }
        double Tick(double in) {               
            x = in;
            y = b0*x + b1*x1 + b2 * x2;
            x2=x1;
            x1= x;                
            return y;

        }
    };

        ////////////////////////////////////////////////////////////////////////////////////////////
    // Zolzer
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct ZolzerBiquad
    {
        double a[2];
        double b[3];
        double fs,fc,q,g;
        double x1,x2,y1,y2;
        double x,y; 
        double res;   

        enum FilterType filter_type;

        ZolzerBiquad(FilterType type, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            res = 0;
            x1=x2=y1=y2=0;
            filter_type = type;
            init_filter(Fc,Q);        
        }

        void init_filter(double Fc, double Q, double gain=1)
        {
            fc = Fc/fs*0.5;
            
            q = Q;
            g = gain;

            switch(filter_type)
            {
                case Lowpass: lowpass(fc,q); break;
                case Highpass: highpass(fc,q); break;
                case Bandpass: bandpass(fc,q); break;
                case Notch: notch(fc,q); break;
                // fc/q dont matter q must be 0
                case Allpass: allpass(fc,0); break;
                //have to find it
                //case Peak: peak(fc,q,gain); break;
                //case Lowshelf: lowshelf(fc,q); break;
                //case Highshelf: highshelf(fc,q); break;
                default: assert(1==0);
            }
        }
        
        void setCutoff(double f) {
            fc = f;
            init_filter(fc,q,g);
        }
        void setQ(double Q) {
            q  = Q;
            init_filter(fc,q,g);
        }
        void setResonance(double R) 
        {
            res = R;
        }
        void setGain(double G) {
            g = G;
            init_filter(fc,q,g);
        }

        void notch(double f, double Q) {
            fc = f;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double Kq = Q*(1+K*K) ;
            double Kk = (K*K*Q+K+Q);        
            b[0] = Kq/Kk;
            b[1] = (2*Kq)/Kk;
            b[2] = Kq/Kk;
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        void lowpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = std::tan(M_PI*fc);
            b[0] = K/(K+1);
            b[1] = K/(K+1);
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void highpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = std::tan(M_PI*fc);
            b[0] = 1/(K+1);
            b[1] = -1/(K+1);
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void allpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = std::tan(M_PI*fc);
            b[0] = (K-1)/(K+1);
            b[1] = 1;
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void lowpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q);        
            double Kq = (K*K*Q);
            b[0] = Kq/Kk;
            b[1] = (2*Kq) /Kk;
            b[2] =  Kq / Kk;
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        void allpass(double f, double Q) {
            fc = f;                
            q  = Q;
            double K = std::tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q);        
            double Km = (K*K*Q-K+Q);
            double Kq = 2*Q*(K*K-1);
            b[0] = Km/Kk;
            b[1] = Kq/Kk;
            b[2] = 1.0f;
            a[0] = Kq/Kk;
            a[1] = Km/Kk;
        }
        void highpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q); 
            double Kq = 2*Q*(K*K-1);
            double Km = (K*K*Q-K+Q);
            b[0] = Q / Kk;
            b[1] = -(2*Q)/Kk;
            b[2] = Q / Kk;
            a[1] = Kq/Kk;
            a[2] = Km/Kk;
        }    
        void bandpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q); 
            b[0] = K / Kk;
            b[1] = 0;
            b[2] = -b[0];
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        // lowshelf
        void lfboost(double f, double G)
        {
            fc = f;
            g  = G;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);
            double Kaka1 = std::sqrt(2*V0) * K + V0*K*K;
            double Kaka2 = 1 + std::sqrt(2)*K + K*K;
            b[0] = (1+Kaka1)/Kaka2;
            b[1] = (2*(V0*K*K-1))/ Kaka2;
            b[2] = (1 - Kaka1)/Kaka2;
            a[0] = (2*(K*K-1))/Kaka2;
            a[1] = (1-std::sqrt(2)*K+K*K)/Kaka2;
        }
        // lowshelf
        void lfcut(double f, double G)
        {
            fc = f;
            g  = G;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);
            double Kaka = V0 + std::sqrt(2*V0)*K + K*K;
            b[0] = (V0*(1+std::sqrt(2)*K+K*K))/Kaka;
            b[1] = (2*V0*(K*K-1))/ Kaka;
            b[2] = (V0*(1-std::sqrt(2)*K+K*K))/Kaka;
            a[0] = (2*(K*K-V0))/Kaka;
            a[1] = (V0-std::sqrt(2*V0)*K+K*K)/Kaka;
        }
        // hishelf
        void hfboost(double f, double G)
        {
            fc = f;
            g  = G;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);            
            double Kaka = 1 + std::sqrt(2)*K + K*K;
            b[0] = (V0 + std::sqrt(2*V0)*K + K*K)/Kaka;
            b[1] = (2*(K*K-V0))/Kaka;
            b[2] = (V0 - std::sqrt(2*V0)*K + K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1-std::sqrt(2*K)+K*K)/Kaka;
        }
        // hishelf
        void hfcut(double f, double G)
        {
            fc = f;
            g  = G;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);            
            double Kaka = 1 + std::sqrt(2*V0)*K + V0*K*K;
            b[0] = (V0*(1 + std::sqrt(2)*K + K*K))/Kaka;
            b[1] = (2*V0*(K*K-1))/Kaka;
            b[2] = (V0*(1 - std::sqrt(2)*K + K*K))/Kaka;
            a[0] = (2*(V0*K*K-1))/Kaka;
            a[1] = (1-std::sqrt(2*V0)*K + V0*K*K)/Kaka;
        }
        // peak
        void boost(double f, double Q, double G)
        {
            fc = f;
            g  = G;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);            
            double Kaka = 1 + (1/Q)*K + K*K;
            b[0] = (1+(V0/Q)*K + K*K)/Kaka;
            b[1] = (2*(K*K-1))/Kaka;
            b[2] = (1- (V0/Q)*K + K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1 - (1/Q)*K + K*K)/Kaka;
        }
        //peak
        void cut(double f, double Q, double G)
        {
            fc = f;
            g  = G;
            q  = Q;
            double K = std::tan(M_PI*fc);
            double V0= std::pow(10,G/20.0);            
            double Kaka = 1 + (1/(V0*Q)*K + K*K);
            b[0] = (1 + (1/Q)*K + K*K)/Kaka;
            b[1] = (2*(K*K-1))/Kaka;
            b[2] = (1 - (1/Q)*Kaka *K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1 - (1/(V0*Q)*K +K*K))/Kaka;
        }
        
        double Tick(double I, double A = 1, double X = 0, double Y = 0)
        {
            Undenormal denormal;
            x = I; 
            y = b[0]*x + b[1]*x1 + b[2]*x2 - a[0]*y1 - a[1]*y2;        
            y2 = y1;
            y1 = y;
            x2 = x1;
            x1 = x;        
            return y;
        }
    };
}