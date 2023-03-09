#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>

#include "Kfr/kfrcore.hpp"
#include "AudioFFT/AudioFFT.h"
#include "AudioFFT/FFTConvolver.h"

namespace KfrDSP1
{    
    struct matrix_index {
        size_t i;
        size_t j;
    };

    template<typename T>
    struct vector_matrix : public kfr::univector<T>
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
        kfr::univector<T> samples[2];

        StereoVector(size_t n) {
            for(size_t i = 0; i < n; i++) samples[i].resize(n);
        }

        kfr::univector<T>& operator[](size_t channel) { return samples[channel]; }

        kfr::univector<T> get_left_channel() { return samples[0]; }
        kfr::univector<T> get_right_channel() { return samples[1]; }

        void set_left_channel(kfr::univector<T> & v) { samples[0] = v; }
        void set_right_channel(kfr::univector<T> & v) { samples[1] = v; }
    };


    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, kfr::univector<T> & samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, kfr::univector<T> & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    kfr::univector<T> get_left_channel(const kfr::univector<T> & in) {
        kfr::univector<T> r(in.size()/2);
        size_t x = 0;
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    kfr::univector<T> get_right_channel(const kfr::univector<T> & in) {
        kfr::univector<T> r(in.size()/2);
        size_t x = 0;
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    kfr::univector<T> get_channel(size_t ch, const kfr::univector<T> & in) {
        kfr::univector<T> r(in.size()/2);
        size_t x = 0;
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const kfr::univector<T> & left, kfr::univector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const kfr::univector<T> & right, kfr::univector<T> & out) {
        size_t x = 0;
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const kfr::univector<T> & in, kfr::univector<T> & out) {
        size_t x = 0;
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    kfr::univector<T> interleave(size_t n, size_t channels, const std::vector<kfr::univector<T>> & in) {
        kfr::univector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    kfr::univector<T> interleave(size_t n, size_t channels, const std::vector<T*> & in) {
        kfr::univector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    std::vector<kfr::univector<T>> deinterleave(size_t n, size_t channels, const kfr::univector<T> & in) {
        std::vector<kfr::univector<T>> r(n);
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
    kfr::univector<T> get_left_channel(size_t n, const V& in) {
        kfr::univector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename V, typename T>
    kfr::univector<T> get_right_channel(size_t n, const V& in) {
        kfr::univector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename V, typename T>
    kfr::univector<T> get_channel(size_t ch, size_t n, V& in) {
        kfr::univector<T> r(n/2);
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
    kfr::univector<T> interleave(size_t n, size_t channels, const V& in) {
        kfr::univector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename V, typename T>
    std::vector<kfr::univector<T>> deinterleave(size_t n, size_t channels, const V& in) {
        std::vector<kfr::univector<T>> r(n);
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
    kfr::univector<T> get_left_channel(size_t n, const T* in) {
        kfr::univector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    kfr::univector<T> get_right_channel(size_t n, const T* & in) {
        kfr::univector<T> r(n/2);
        size_t x = 0;
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    kfr::univector<T> get_channel(size_t ch, size_t n, T* in) {
        kfr::univector<T> r(n/2);
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
    kfr::univector<T> interleave(size_t n, size_t channels, const T** & in) {
        kfr::univector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    std::vector<kfr::univector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        std::vector<kfr::univector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }

     template<typename T>
     bool equal_vector (kfr::univector<T> & a, kfr::univector<T> & b) {
         return std::equal(a.begin(),a.end(),b.end());
     }
    template<typename T>
    void copy_vector(kfr::univector<T> & dst, kfr::univector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy_vector(kfr::univector<T> & dst, size_t n, T * src) {
        std::copy(&src[0],&src[n-1],dst.begin());
    }
    template<typename T>
    kfr::univector<T> slice_vector(size_t start, size_t end, kfr::univector<T> & src) {
        kfr::univector<T> r(end-start);
        std::copy(src.begin()+start,src.begin()+end,r.begin());
        return r;
    }

    template<typename T>
    void copy_buffer(size_t n, T * dst, T * src) {
        memcpy(dst,src,n*sizeof(T));
    }

    template<typename T>
    kfr::univector<T> slice_buffer(size_t start, size_t end, T * ptr) {
        kfr::univector<T> r(end-start);
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
    void split_stereo(const kfr::univector<T> & input, kfr::univector<T> & left, kfr::univector<T> & right) {
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
    bool isEmpty(kfr::univector<T> v)
    {
        return (v.size() == 0);
    }

    //============================================================
    template <class T>
    bool containsOnlyZeros(kfr::univector<T> v)
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
    bool isAllPositiveOrZero(kfr::univector<T> v)
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
    bool isAllNegativeOrZero(kfr::univector<T> v)
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
    bool contains(kfr::univector<T> v, T element)
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
    T max(kfr::univector<T> v)
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
    int maxIndex(kfr::univector<T> v)
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
    T min(kfr::univector<T> v)
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
    int minIndex(kfr::univector<T> v)
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
    void printVector(kfr::univector<T> v)
    {
        for (int i = 0;i < v.size();i++)
        {
            std::cout << v[i] << std::endl;
        }
    }

    //============================================================
    template <class T>
    T getLastElement(kfr::univector<T> v)
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
    T getFirstElement(kfr::univector<T> v)
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
    kfr::univector<T> getEveryNthElementStartingFromK(kfr::univector<T> v,int n,int k)
    {
        if ((n >= v.size()) || (n >= v.size()))
        {
            throw std::invalid_argument( "Invalid arguments for getEveryNthElementStartingFromK()");
        }
        else
        {
            kfr::univector<T> result;

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
    kfr::univector<T> getEvenElements(kfr::univector<T> v)
    {
        return getEveryNthElementStartingFromK(v, 2, 0);
    }

    //============================================================
    template <class T>
    kfr::univector<T> getOddElements(kfr::univector<T> v)
    {
        return getEveryNthElementStartingFromK(v, 2, 1);
    }

    //============================================================
    template <class T>
    void fillVectorWith(kfr::univector<T> &v,T element)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] = element;
        }
    }

    //============================================================
    template <class T>
    int countOccurrencesOf(kfr::univector<T> v,T element)
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
    T sum(kfr::univector<T> v)
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
    T product(kfr::univector<T> v)
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
    T mean(kfr::univector<T> v)
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
    T median(kfr::univector<T> v)
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
    T variance(kfr::univector<T> v)
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
    T standardDeviation(kfr::univector<T> v)
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
    T norm1(kfr::univector<T> v)
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
    T norm2(kfr::univector<T> v)
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
    T magnitude(kfr::univector<T> v)
    {
        // just another name for L2-norm
        return norm2(v);
    }

    //============================================================
    template <class T>
    T normP(kfr::univector<T> v,T p)
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
    void multiplyInPlace(kfr::univector<T> &v,T scalar)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] *= scalar;
        }
    }

    //============================================================
    template <class T>
    void multiplyInPlace(kfr::univector<T> &v1,kfr::univector<T> v2)
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
    void divideInPlace(kfr::univector<T> &v,T scalar)
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
    void divideInPlace(kfr::univector<T> &v1,kfr::univector<T> v2)
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
    void addInPlace(kfr::univector<T> &v,T value)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] += value;
        }
    }

    //============================================================
    template <class T>
    void addInPlace(kfr::univector<T> &v1,kfr::univector<T> v2)
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
    void subtractInPlace(kfr::univector<T> &v,T value)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] -= value;
        }
    }

    //============================================================
    template <class T>
    void subtractInPlace(kfr::univector<T> &v1,kfr::univector<T> v2)
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
    void absInPlace(kfr::univector<T> &v)
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
    void squareInPlace(kfr::univector<T> &v)
    {
        for (int i = 0;i < v.size();i++)
        {
            v[i] = v[i]*v[i];
        }
    }

    //============================================================
    template <class T>
    void squareRootInPlace(kfr::univector<T> &v)
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
    void sort(kfr::univector<T> &v)
    {
        std::sort(v.begin(),v.end());
    }

    //============================================================
    template <class T>
    void reverse(kfr::univector<T> &v)
    {
        std::reverse(v.begin(), v.end());
    }

    //============================================================
    template <class T>
    kfr::univector<T> multiply(kfr::univector<T> v,T scalar)
    {
        kfr::univector<T> result;

        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] * scalar);
        }

        return result;
    }

    //============================================================
    template <class T>
    kfr::univector<T> multiply(kfr::univector<T> v1,kfr::univector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            kfr::univector<T> result;

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
    kfr::univector<T> divide(kfr::univector<T> v,T scalar)
    {
        if (scalar != 0)
        {
            kfr::univector<T> result;

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
    kfr::univector<T> divide(kfr::univector<T> v1,kfr::univector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            if (!contains<T>(v2, 0))
            {
                kfr::univector<T> result;

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
    kfr::univector<T> add(kfr::univector<T> v,T value)
    {
        kfr::univector<T> result;

        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] + value);
        }

        return result;
    }

    //============================================================
    template <class T>
    kfr::univector<T> add(kfr::univector<T> v1,kfr::univector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            kfr::univector<T> result;

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
    kfr::univector<T> subtract(kfr::univector<T> v,T value)
    {
        kfr::univector<T> result;

        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i] - value);
        }

        return result;
    }

    //============================================================
    template <class T>
    kfr::univector<T> subtract(kfr::univector<T> v1,kfr::univector<T> v2)
    {
        if (v1.size() == v2.size())
        {
            kfr::univector<T> result;

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
    kfr::univector<T> abs(kfr::univector<T> v)
    {
        kfr::univector<T> result;

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
    kfr::univector<T> square(kfr::univector<T> v)
    {
        kfr::univector<T> result;

        for (int i = 0;i < v.size();i++)
        {
            result.push_back(v[i]*v[i]);
        }

        return result;
    }


    //============================================================
    template <class T>
    kfr::univector<T> squareRoot(kfr::univector<T> v)
    {
        if (isAllPositiveOrZero(v))
        {
            kfr::univector<T> result;

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
    kfr::univector<T> scale(kfr::univector<T> v,T lowerLimit,T upperLimit)
    {
        kfr::univector<T> result;

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
    kfr::univector<T> difference(kfr::univector<T> v)
    {
        kfr::univector<T> result;

        for (int i = 1;i < v.size();i++)
        {
            result.push_back(v[i]-v[i-1]);
        }

        return result;
    }

    //============================================================
    template <class T>
    kfr::univector<T> zeros(int N)
    {
        if (N >= 0)
        {
            kfr::univector<T> result;

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
    kfr::univector<T> ones(int N)
    {
        if (N >= 0)
        {
            kfr::univector<T> result;

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
    kfr::univector<T> range(int limit1,int limit2,int step)
    {
        kfr::univector<T> result;

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
    kfr::univector<T> range(int maxValue)
    {
        return SoundAlchemy::range<T>(0, maxValue, 1);
    }

    //============================================================
    template <class T>
    kfr::univector<T> range(int minValue,int maxValue)
    {
        return SoundAlchemy::range<T>(minValue, maxValue, 1);
    }

    //============================================================
    template <class T>
    T dotProduct(kfr::univector<T> v1,kfr::univector<T> v2)
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
    T euclideanDistance(kfr::univector<T> v1,kfr::univector<T> v2)
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
    T cosineSimilarity(kfr::univector<T> v1,kfr::univector<T> v2)
    {
    return dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
    }

    //============================================================
    template <class T>
    T cosineDistance(kfr::univector<T> v1,kfr::univector<T> v2)
    {
        return 1.0 - cosineSimilarity(v1, v2);
    }


    template<class T>
    kfr::univector<T> cos(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    kfr::univector<T> sin(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> tan(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    kfr::univector<T> acos(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    kfr::univector<T> asin(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> atan(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> atan2(const kfr::univector<T> & v, const T value) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
        return r;
    }    
    template<class T>
    kfr::univector<T> atan2(const kfr::univector<T> & v, const kfr::univector<T> value) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> cosh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    kfr::univector<T> sinh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> tanh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    kfr::univector<T> acosh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    kfr::univector<T> asinh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> atanh(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    kfr::univector<T> exp(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> log(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> log10(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> exp2(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> expm1(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> ilogb(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> log2(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> log1p(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> logb(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> scalbn(const kfr::univector<T> & v, const std::vector<int> & x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> scalbn(const kfr::univector<T> & v, const int x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
        return r;
    }    
    template<class T>
    kfr::univector<T> scalbln(const kfr::univector<T> & v, const std::vector<long int> & x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> scalbln(const kfr::univector<T> & v, const long int x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
        return r;
    }    
    template<class T>
    kfr::univector<T> pow(const kfr::univector<T> & v, const kfr::univector<T> & x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> pow(const kfr::univector<T> & v, const T x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    kfr::univector<T> pow(const T x, const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> sqrt(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> cbrt(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> hypot(const kfr::univector<T> & v, const kfr::univector<T> & x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> hypot(const kfr::univector<T> & v, const T x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
        return r;
    }    
    template<class T>
    kfr::univector<T> hypot(const T x, const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> erf(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> erfc(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> tgamma(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> lgamma(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> ceil(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> floor(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fmod(const kfr::univector<T> & v, const kfr::univector<T> & x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fmod(const kfr::univector<T> & v, const T x) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
        return r;
    }    
    template<class T>
    kfr::univector<T> fmod(const T x, const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> trunc(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> round(const kfr::univector<T> & v) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
        return r;
    }    
    template<class T>
    std::vector<long int> lround(const kfr::univector<T> & v) {
        std::vector<long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
        return r;
    }    
    template<class T>
    std::vector<long long int> llround(const kfr::univector<T> & v) {
        std::vector<long long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> nearbyint(const kfr::univector<T> & v) {
        std::vector<long long int> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> remainder(const kfr::univector<T> & a, const kfr::univector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> copysign(const kfr::univector<T> & a, const kfr::univector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fdim(const kfr::univector<T> & a, const kfr::univector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
        return r;
    }    
    #undef fmax
    template<class T>
    kfr::univector<T> fmax(const kfr::univector<T> & a, const kfr::univector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fmin(const kfr::univector<T> & a, const kfr::univector<T> & b) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fma(const kfr::univector<T> & a, const kfr::univector<T> & b, const kfr::univector<T> & c) {
        std::vector<long long int> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
        return r;
    }    
    template<class T>
    kfr::univector<T> fabs(const kfr::univector<T> & a) {
        kfr::univector<T> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
        return r;
    }    
    
    /*
    template<class T>    
    kfr::univector<T> abs(const kfr::univector<T> & a) {
        kfr::univector<T> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::abs((long la[i]);
        return r;
    } 
    */
       


    template<typename T>
    StereoVector<T> stereo(const kfr::univector<T> & left, const kfr::univector<T> & right) {
        StereoVector<T> r(left.size());
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i];
            r[1][i] = right[i];
        }
    }

    template<typename T>
    kfr::univector<T> merge(const kfr::univector<T> & left, const kfr::univector<T> & right) {
        kfr::univector<T> r(left.size()*2);
        size_t x = 0;
        for(size_t i = 0; i < left.size(); i++)
        {
            r[x++] = left[i];
            r[x++] = right[i];
        }
    }

    template<typename T>
    void swap(kfr::univector<T> & left, kfr::univector<T> & right) {
        std::swap(left,right);
    }
    
    template<typename T>
    bool isin(const kfr::univector<T> & v, const T val) {
        return std::find(v.begin(),v.end(),val) != v.end();
    }

    template<typename T>
    StereoVector<T> pan(const kfr::univector<T> & left, const kfr::univector<T> & right, T amt) {
        StereoVector<T> r(left.size());
        T pan_map = ((amt+1)/2.0) * (M_PI/2.0);
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * sin(pan_map);
            r[1][i] = right[i] * cos(pan_map);
        }
        return r;
    }
    template<typename T>
    StereoVector<T> constant_power_pan(const kfr::univector<T> & left, const kfr::univector<T> & right, T pos) {
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
    kfr::univector<T> mix(const kfr::univector<T> & a, const kfr::univector<T> & b)
    {
        assert(a.size() == b.size());
        kfr::univector<T> r(a.size());
        for(size_t i = 0; i < r.size(); i++) r[i] = a[i]+b[i];
        T max = std::max_element(r.begin(),r.end());
        if(max > 0) for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }
    template<typename T>
    kfr::univector<T> normalize(const kfr::univector<T> & a) {
        kfr::univector<T> r(a);        
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
    kfr::univector<T> kernel(const kfr::univector<T> & v, T (*f)(T value)) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = f(v[i]);
        return r;
    }
    template<class T>
    kfr::univector<T> kernel(const kfr::univector<T> & v, std::function<T (T)> func) {
        kfr::univector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = func(v[i]);
        return r;
    }
    template<class T>
    void inplace_add(const kfr::univector<T> & a, kfr::univector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] += func(a[i]);        
    }
    template<class T>
    void inplace_sub(const kfr::univector<T> & a, kfr::univector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] -= func(a[i]);        
    }
    template<class T>
    void inplace_mul(const kfr::univector<T> & a, kfr::univector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] *= func(a[i]);        
    }
    template<class T>
    void inplace_div(const kfr::univector<T> & a, kfr::univector<T> & r, std::function<T (T)> func) {        
        for(size_t i = 0; i < a.size(); i++) r[i] /= func(a[i]);
    }

    
    template<class T>
    void fill(kfr::univector<T> & in, T x)
    {
        for(size_t i = 0; i < in.size(); i++) in[i] = x;
    }
    template<class T>
    void zeros(kfr::univector<T> & in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(kfr::univector<T> & in)
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
    kfr::univector<T> StdMean(kfr::univector<T> & in) 
    {
        kfr::univector<T> r(in.size());
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
    kfr::univector<T> StdRMS(kfr::univector<T> & in) 
    {
        kfr::univector<T> r(in.size());
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
    kfr::univector<T> StdDev(kfr::univector<T> & in) 
    {
        kfr::univector<T> r(in.size());
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
    kfr::univector<T> StdVariance(kfr::univector<T> & in) 
    {
        kfr::univector<T> r(in.size());
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
    kfr::univector<T> vector_copy(T * ptr, size_t n) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
            r[i] = ptr[i];
        return r;
    }

}