#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <valarray>

namespace Casino::eigen
{   
    template<typename T>
    struct SampleVector : public Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>
    {
        SampleVector() = default;
        SampleVector(size_t n) : Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>(n) {}
        SampleVector(const SampleVector<T> & v) { *this = v; }
        SampleVector(const std::vector<T> & d) {
            resize(d.size());
            memcpy(data(),d.data(),d.size()*sizeof(T));
        }

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setLinSpaced;

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        //using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        //using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <=;


        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::fill;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cols;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setRandom;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::sum;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::dot;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cross;

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::real;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::imag;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseProduct;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseQuotient;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseSqrt;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseInverse;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseMin;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseMax; 
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs2;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::unaryExpr;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::array;
        


        

        // mean
        // min
        // max
        // stddev
        // rms
        // correlation
        // autocorrelation
        // count(n)
        // 

    };

    template<typename T>
    struct SampleMatrix : public Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
    {
        SampleMatrix() = default;
        SampleMatrix(size_t i, size_t j) : Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(i,j) {}
        SampleMatrix(const SampleMatrix<T> & m) { *this = m; }
        SampleMatrix(const std::vector<T> & d,size_t i, size_t j) {
            resize(i,j);
            for(size_t x = 0; x < i; x++)
                for(size_t y = 0; y < j; y++)
                    (*this)(x,y) = d[x*j + y];
        }

        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        //using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        //using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <=;


        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::fill;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::data;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setRandom;            

        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setIdentity;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::head;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::tail;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::segment;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::block;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::row;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::col;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rows;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cols;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::leftCols;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::middleCols;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rightCols;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::topRows;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::middleRows;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::bottomRows;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::topLeftCorner;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::topRightCorner;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::bottomLeftCorner;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::bottomRightCorner;
        /*
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::seq;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::seqN;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::A;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::v;        
        */
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::adjoint;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::transpose;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::diagonal;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::eval;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::asDiagonal;        
        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::replicate;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::reshaped;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::select;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseProduct;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseQuotient;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseSqrt;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseInverse;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseMin;        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseMax; 
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs2;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::unaryExpr;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::array;

        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::minCoeff;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::maxCoeff;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sum;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::colwise;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rowwise;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::trace;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::all;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::any;

        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::norm;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::squaredNorm;
        
        
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::real;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::imag;

        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::ldlt;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::llt;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::lu;
        //using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::qr;
        //using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::svd;
        using Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::eigenvalues;                
    };

    template<typename T>
    struct SampleArray : public Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>
    {        
        
        size_t channels;

        SampleArray() {
            channels = 1;
        }
        SampleArray(size_t size, size_t channels = 1) {
            resize(size * channels);
            this->channels = channels;
        }
        SampleArray(const std::vector<T> & s, size_t chans = 1) {
            resize(s.size());
            memcpy(data(),s.data(),s.size()*sizeof(T));
            channels = chans;
        }
        SampleArray(const T * ptr, size_t n, size_t chans = 1) {
            resize(n*chans);
            memcpy(data(),ptr,n*chans*sizeof(T));
            channels = chans;
        }
        SampleArray(const SampleArray<T> & s) {
            (*this) = s;
            channels = s.channels;
        }    
        SampleArray(const Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor> & v, size_t ch = 1) {
            (*this) = v;
            channels = ch;
        }
        SampleArray(size_t channels, const Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>& a) {
            (*this) = a;
            this->channels = channels;
        }

        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <=;

        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cols;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::fill;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setRandom;        
        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::vector;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::inverse;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::pow;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::square;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cube;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sqrt;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::rsqrt;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::exp;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log1p;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log10;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::max;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::min;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::abs;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::abs2;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sin;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cos;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::tan;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::asin;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::acos;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::atan;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sinh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cosh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::tanh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::asinh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::acosh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::atanh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::ceil;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::floor;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::round;
        
        
        SampleArray<T> clamp(T min = (T)0.0, T max = (T)1.0) {
            SampleArray<T> r(*this,channels);
            r = SampleArray<T>(r.cwiseMin(max).cwiseMax(min));
            return r;
        }
        

        T&   get_stride(size_t ch, size_t pos) { return (*this)[pos*channels + ch]; }
        void set_stride(size_t ch, size_t pos, const T val) { (*this)[pos*channels + ch] = val; } 
        
        void swap_stereo_channels() {
            assert(channels == 2);
            for(size_t i = 0; i < size(); i+= channels) {
                T temp = (*this)[i];
                (*this)[i] = (*this)[i+1];
                (*this)[i+1] = temp;
            }
        }

        void set_data(size_t n, size_t channels, T * samples)
        {
            resize(n);
            this->channels = channels;
            memcpy(data(),samples,n*sizeof(T));
        }
        void copy_data(T * samples) {
            memcpy(samples,data(),size()*sizeof(T));
        }
        void set_channel(size_t ch, T * samples) {
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = samples[x++];
        }

        
        size_t num_channels() const { return channels; } 
        size_t samples_per_channel() const { return size() / num_channels(); }
                
        
        T& operator()(size_t i, size_t ch) {
            Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor> & r = *this;
            return r[i*channels + ch];
        }
        T operator()(size_t i, size_t ch) const {
            const Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor> & r = *this;
            return r[i*channels + ch];
        }
         
        SampleArray<T> get_channel(size_t channel) {            
            SampleArray<T> r(samples_per_channel());
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];
            return r;
        }
        void get_channel(size_t channel, T * r) {                        
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];            
        }
        SampleArray<T> get_channel(size_t channel) const {            
            SampleArray<T> r(samples_per_channel());
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];
            return r;
        }
        void set_channel(const SampleArray<T> & v, size_t ch) {            
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = v[x++];
        }
        void set_channel(const T* v, size_t ch) {            
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = v[x++];
        }
        void make_stereo(const T * in) {
            set_channel(in,0);
            set_channel(in,1);
        }
        void make_stereo(const T * left, const T * right) {
            set_channel(left,0);
            set_channel(right,1);
        }

        void set_channels(size_t c) {
            channels = c;
            resize(size()*c);
        }
        void resize(size_t n) {
            Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize(n);
        }
        void resize(size_t s, size_t c) {
            channels = c;
            Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize(s * c);
        }

        SampleArray<T> get_channel_count(size_t ch, size_t pos, size_t n) {
            SampleArray<T> r(n);
            memcpy(r.data(), get_channel(ch).data()+pos, n*sizeof(T));
            return r;
        }

        bool operator == (const SampleArray<T> & s) {
            return s.channels == channels && size() == s.size();
        }

        void copy(SampleArray<T> & in, size_t block) {
            resize(in.size());
            memcpy(data(),in.data(),block*sizeof(T));
        }
        void copy(SampleArray<T> & in, size_t pos, size_t n) {
            resize((pos+n));
            memcpy(data(),in.data()+pos,n*sizeof(T));
        }
        void copy_from(const T* in, size_t block) {            
            memcpy(data(),in,block*sizeof(T));
        }
        void copy_to(T* in, size_t block) {            
            memcpy(in,data(),block*sizeof(T));
        }

        SampleArray<T> slice(size_t i, size_t n) {
            SampleArray<T> r(n);
            memcpy(r.data(), data()+i, n * sizeof(T));
            return r;
        }
        
        SampleArray<T> pan(T pan) {
            T pan_mapped = ((pan+1.0)/2.0) * (M_PI/2.0);
            SampleArray<T> r(size(),2);
            for(size_t i = 0; i < size(); i+=channels)
            {
                r[i] = (*this)[i] * std::sin(pan_mapped);
                if(num_channels() == 2)
                    r[i+1] = (*this)[i+1] * std::cos(pan_mapped);
                else
                    r[i+1] = (*this)[i] * std::cos(pan_mapped);
            }
            return r;        
        }
        SampleArray<T> stride_slice(const SampleArray<T> & in, size_t stride, size_t pos, size_t n)
        {
            SampleArray<T> r(n/stride);
            size_t x = pos;
            for(size_t i = 0; i < n; i++) {
                r[i] = (*this)[x];
                x += stride;
            } 
            return r;
        }
        
        void println() {
            std::cout << "Size=" << size() << std::endl;
            std::cout << "Channel=" << num_channels() << std::endl;
            std::cout << "Samples per channel=" << samples_per_channel() << std::endl;
            for(size_t i = 0; i < size(); i++)
                std:: cout << (*this)[i] << ",";
            std::cout << std::endl;
        }

    };

    template<typename T>
    struct SampleArray2D : public Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
    {        
        SampleArray2D() = default;
        
        SampleArray2D(size_t chan, size_t samps) {
            resize(chan,samps);
        }
        SampleArray2D(const SampleArray<T> & v) {
            resize(v.num_channels(), v.size()/v.num_channels());
            for(size_t i = 0; i < v.num_channels(); i++)             
                this->row(i) = v.get_channel(i);
        }
        SampleArray2D(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & c) {
            (*this) = c;
        }
        SampleArray2D(const SampleArray2D<T> & m) {
            (*this) = m.channels;
        }

        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <=;

        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::fill;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rows;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cols;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::row;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::col;
        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setRandom;        
        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::matrix;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::inverse;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::pow;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::square;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cube;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sqrt;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rsqrt;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::exp;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log1p;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log10;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::max;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::min;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::abs;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::abs2;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sin;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cos;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::tan;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::asin;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::acos;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::atan;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sinh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cosh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::tanh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::asinh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::acosh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::atanh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::ceil;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::floor;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::round;
        
        
        void deinterleave(const SampleArray<T> & s) {
            
            esize(s.num_channels(),s.samples_per_channel());            
            
            for(size_t i = 0; i < s.num_channels(); i++) {
                this->row(i) = s.get_channel(i);
            }
            
        }    

        SampleArray<T> interleave() { 
            SampleArray<T> r;
            r.resize(cols() * cols());
            for(size_t i = 0; i < rows(); i++)            
            {
                SampleArray<T> tmp(1,(*this)(i));            
                r.set_channel(tmp,i);
            }
            return r;
        }

        size_t num_channels() const { return Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rows(); }
                
        SampleArray<T> get_channel(size_t c) { 
            SampleArray<T> tmp(1,row(c));
            return tmp;
        }
        void set_channel(size_t c, const SampleArray<T> & s) 
        { 
            row(c) = s;
        }

        
    };
      
    template<typename T> 
    SampleArray2D<T> deinterleave(const SampleArray<T> & v) {
        SampleArray2D<T> m;
        m.deinterleave(v);
        return m;
    }

    template<typename T> 
    SampleArray<T> interleave(SampleArray2D<T> & m) {
        return m.interleave();
    }

    
    template<typename T>
    struct Window : public SampleArray<T>
    {       
        Window(size_t i) : SampleArray<T>(i) {}
        virtual ~Window() = default;     
    };

    template<typename T>
    struct Rectangle: public Window<T>
    {
        Rectangle(size_t i) : Window<T>(i) { this->fill(1); } 
    };

    template<typename T>
    struct Hamming: public Window<T>
    {
        Hamming(size_t i) : Window<T>(i) {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)(i) = 0.54 - (0.46 * std::cos(2*M_PI*i/n));
            }        
        }
    };

    template<typename T>
    struct Hanning: public Window<T>
    {
        Hanning(size_t i) : Window<T>(i) {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)(i) = 0.5*(1 - std::cos(2*M_PI*i/n));
            }        
        }
    };

    template<typename T>
    struct Blackman: public Window<T>
    {
        Blackman(size_t i) : Window<T>(i)    
        {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)                    
                (*this)(i) = 0.42 - (0.5* std::cos(2*M_PI*i/(n)) + (0.08*std::cos(4*M_PI*i/n)));        
        }
    };

    template<typename T>
    struct BlackmanHarris: public Window<T>
    {
        BlackmanHarris(size_t i) : Window<T>(i)    
        {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)            
                (*this)(i) = 0.35875 
                        - 0.48829*std::cos(2*M_PI*(i/n))
                        + 0.14128*std::cos(4.0*M_PI*(i/n)) 
                        - 0.01168*std::cos(6.0*M_PI*(i/n));
        }
    };

    template<typename T>
    struct Gaussian: public Window<T>
    {
        Gaussian(size_t i) : Window<T>(i)
        {
            T a,b,c=0.5;
            for(size_t n = 0; n < this->size(); n++)
            {
                a = (n - c*((*this).size()-1)/(std::sqrt(c)*(*this).size()-1));
                b = -c * std::sqrt(a);
                (*this)(n) = std::exp(b);
            }
        }
    };
    template<typename T>
    struct Welch: public Window<T>
    {
        Welch(size_t i) : Window<T>(i)
        {
            for(size_t i = 0; i < (*this).size(); i++)
                (*this)(i) = 1.0 - std::sqrt((2*i-(*this).size())/((*this).size()+1));        
        }
    };
    template<typename T>
    struct Parzen: public Window<T>
    {

        Parzen(size_t i) : Window<T>(i)
        {
            for(size_t i = 0; i < (*this).size(); i++)
                (*this)(i) = 1.0 - std::abs((2.0f*i-(*this).size())/((*this).size()+1.0f));        
        }    
    };
    template<typename T>
    struct Tukey: public Window<T>
    {
        Tukey(size_t i, T alpha) : Window<T>(i)
        {
            size_t num_samples = (*this).size();
            T value = (-1*(num_samples/2)) + 1;
            for(size_t i = 0; i < num_samples; i++)
            {    
                if(value >= 0 && value <= (alpha * (num_samples/2))) 
                    (*this)[i] = 1.0; 
                else if(value <= 0 && (value >= (-1*alpha*(num_samples/2)))) 
                    (*this)[i] = 1.0;
                else 
                    (*this)[i] = 0.5 * (1 + std::cos(M_PI *(((2*value)/(alpha*num_samples))-1)))        ;
                value = value + 1;
            }     
        }
    };

    template<typename T>
    void swap(SampleArray<T> & vector,size_t i, size_t j)
    {
        T x = vector(i);
        vector(i) = vector(j);
        vector(j) = x;
    }        

    template<typename T>
    void shift(SampleArray<T> & vector) 
    {
        size_t half = vector.size() / 2; 
        size_t start= half; 
        if(2*half < vector.size()) start++;
        for(size_t i = 0; i < half; i++)
        {
            swap(vector,i, i+start);
            
        }
        if(start != half)
        {
            for(size_t i = 0; i < half; i++)            
            {
                swap(vector,i+start-1,i+start);
            }
        }
    }

    template<typename T>    
    void ishift(SampleArray<T> & vector) 
    {
        size_t half = vector.size() / 2; 
        size_t start= half; 
        if(2*half < vector.size()) start++;
        for(size_t i = 0; i < half; i++)
        {
            swap(vector,i,i+start);
        }
        if(start != half)
        {
            for(size_t i = 0; i < half; i++)            
            {
                swap(vector,half,i);
            }
        }
    }

    template<typename T>    
    T quadratic_peak_pos(SampleArray<T> & vector,size_t pos)
    {
        T s0,s1,s2;
        size_t x0,x2;        
        if(pos == 0 || pos == vector.size()-1) return pos; 
        x0 = (pos < 1)? pos:pos-1;
        x2 = (pos+1 < vector.size()) ? pos+1:pos;
        if(x0 == pos) return vector(pos) <= vector(x2)? pos:x2;
        if(x2 == pos) return vector(pos) <= vector(x0)? pos:x0;
        s0 = vector(x0);
        s1 = vector(pos);
        s2 = vector(x2);
        return pos + 0.5 *(s0-s2) / (s0 - 2*s1 + s2);
    }

    template<typename T>    
    T quadratic_peak_mag(SampleArray<T> & vector,size_t pos)
    {
        T x0,x1,x2;
        size_t index = (size_t)(pos - 0.5) + 1;
        if(pos >= vector.size() || pos < (T)0) return (T)0;
        if((T)index == pos) return vector(index);
        x0 = vector(index-1);
        x1 = vector(index);
        x2 = vector(index+1);
        return x1 - 0.25 * (x0 - x2) * (pos - index);
    }

    template<typename T>    
    T median(SampleArray<T> & vector)
    {        
        size_t n = vector.size();
        size_t high,low;
        size_t median;
        size_t middle,ll,hh;

        low = 0;
        high = n-1;
        median = (low+high)/2;

        for(;;)        
        {
            if(high <= low)
                return vector(median);

            if(high == low+1)
            {
                if(vector(low) > vector(high))
                    swap(vector,low,high);
                return vector(median);
            }
            middle = (low+high)/2;
            if(vector(middle) > vector(high)) swap(vector,middle,high);
            if(vector(low) > vector(high)) swap(vector,low,high);
            if(vector(middle) > vector(low)) swap(vector,middle,low);
            swap(vector,middle,low+1);

            ll=low+1;
            hh=high;
            for(;;)
            {
                do ll++; while( vector(low) > vector(ll));
                do hh--; while( vector(hh) > vector(low));
                if(hh < ll) break;
                swap(vector,ll,hh);
            }
            swap(vector,low,hh);
            if(hh <= median) low = ll;
            if(hh >= median) high = hh-1;
        }
    }

    template<typename T>    
    T moving_threshold(SampleArray<T> & input, size_t post, size_t pre, size_t pos)
    {
        size_t length = input.size();
        size_t win_length = post+pre+1;
        SampleArray<T> vector = input;
        if(pos < post+1)
        {
            for(size_t k = 0; k < post +1 - pos; k++)
                vector(k) = 0;
            for(size_t k = post+1-pos; k < win_length; k++)
                vector(k) = vector(k+pos-post);            
        }
        else if(pos + pre < length)
        {
            for(size_t k = 0; k < win_length; k++)
                vector(k) = vector(k+pos-post);
        }
        else 
        {
            for(size_t k = 0; k < length - pos + post; k++)
                vector(k) = vector(k+pos-post);
            for(size_t k = length - pos + post; k < win_length; k++)
                vector(k) = 0;
        }        
        return median(vector);
    }

    template<typename T> 
    T zero_crossing_rate(SampleArray<T> & vector)
    {
        T zcr = 0;
        for(size_t i = 1; i < vector.size(); i++)
        {
            bool current = vector(i) > 0;
            bool prev    = vector(i-1) > 0;
            if(current != prev) zcr++;
        }    
        return zcr;   
    }
    
    template<typename T> 
    void autocorr(SampleArray<T> & vector,SampleArray<T> & output)
    {
        T tmp;
        output.resize(vector.size());
        for(size_t i = 0; i < vector.size(); i++)
        {
            tmp = (T)0;
            for(size_t j = 0; j < vector.size(); j++)
                tmp += vector(j-i) * vector(j);                
            
            output.vector(i) = tmp / (T)(vector.size()-1);
        }
    }

    
    template<typename T>
    void push(SampleArray<T> & vector,const T& new_elem)
    {
        for(size_t i = 0; i < vector.size()-1; i++) vector(i) = vector(i+1);
        vector(vector.size()-1) = new_elem;
    }

    template<typename T>
    void clamp(SampleArray<T> & vector,T absmax) { 
        for(size_t i = 0; i < vector.size(); i++)
        {
            if( vector(i) > absmax) vector(i) = absmax;
            else if(vector(i) < -absmax) vector(i) = -absmax;
        }
    }

    template<typename T> 
    void  normalize(SampleArray<T> & vector) { 
            Eigen::Matrix<T,1,Eigen::Dynamic> r = vector.matrix();
            r.normalize();
            vector = r.array();         
    }


    template<typename T> 
    bool peakpick(SampleArray<T> & vector,size_t pos)
    {
        bool tmp = false;
        tmp = (vector(pos) > vector(pos-1)) && (vector(pos) > vector(pos+1)) && vector(pos) > 0;
        return tmp;
    }
    
    template<typename T> 
    T RMS(SampleArray<T> & vector)
    {
        T sum = vector.pow(2).sum();        
        return std::sqrt(sum/vector.size());
    }

    template<typename T> 
    T peak_energy(SampleArray<T> & vector)
    {
        T peak = -10000;
        T as;         
        for(size_t i = 0; i < vector.size(); i++)
        {
            T as = std::fabs(vector(i));
            if(as > peak) peak = as;
        }
        return peak;
    }

    template<typename T> 
    T min(SampleArray<T> & vector) { 
        T min = 1e120;
        for(size_t i = 0; i < vector.size(); i++)
            if(vector(i) < min) min = vector(i);
        return min;            
        //return vector(vector.minCoeff()); 
    }

    template<typename T> 
    T max(SampleArray<T> & vector) { 
        T max = -1e120;
        for(size_t i = 0; i < vector.size(); i++)
            if(vector(i) > max) max = vector(i);
        return max;
        
    }

    template<typename T> 
    SampleArray<T>& weighted_copy(SampleArray<T> & vector,const SampleArray<T> & weight, SampleArray<T> & out)
    {
        out.vector = vector * weight.vector;
        return out;
    }    

    template<typename T> 
    T level_lin(SampleArray<T> & vector) 
    {
        T energy = vector.sqrt().sum();
        return energy/vector.size();
    }

    template<typename T> 
    T local_hfc(SampleArray<T> & vector)
    {
        T hfc = 0;
        for(size_t j = 0; j < vector.size(); j++)
            hfc += (j+1) * vector(j);
        return hfc;
    }

    template<typename T> 
    void min_removal(SampleArray<T> & vector)
    {
        T m = min(vector);
        vector -= m;
    }

    template<typename T> 
    T alpha_norm(SampleArray<T> & vector,SampleArray<T> & out, T alpha)
    {
        T tmp = vector.abs().pow(alpha).sum();
        return std::pow(tmp/vector.size(),1.0/alpha);
    }

    template<typename T> 
    void alpha_normalize(SampleArray<T> & vector,SampleArray<T> & out, T alpha)
    {
        T tmp = alpha_norm(out,alpha);
        out.vector = out.vector / tmp;         
    }    

    template<typename T> 
    T unwrap2pi(T phase)
    {
        return phase + 2*M_PI*(1.0 + std::floor(-(phase+M_PI)/(2*M_PI)));
    }

    template<typename T> 
    T quadfrac(T s0, T s1, T s2, T pf)
    {
        T tmp = s0 + (pf/2)*(pf*(s0-2.0*s1+s2)-3.0*s0+4*s1-s2);
        return tmp;
    }

    template<typename T> 
    T freqtomidi(T freq)
    {
        T midi;
        if(freq < 2.0 || freq > 100000) return 0;        
        midi = 12*std::log(midi/6.875) / 0.6931471805599453 -3;
        return midi;
    }

    template<typename T> 
    T miditofreq(T midi)
    {
        if(midi > 140) return 0;
        T freq = std::exp(((midi + 3) / 12)*0.6931471805599453)*6.875;
        return freq;
    }

    template<typename T> 
    T bintofreq(T bin, T sample_rate, T fft_size)
    {
        T freq = sample_rate / fft_size;
        return freq * std::max(bin,(T)0);
    }


    template<typename T> 
    T bintomidi(T bin, T sample_rate, T fft_size)
    {
        T midi = bintofreq(bin,sample_rate,fft_size);
        return freqtomidi(midi);
    }

    template<typename T> 
    T freqtobin(T freq, T sample_rate, T fft_size)
    {
        T bin = fft_size / sample_rate;
        return std::max(freq,(T)0)*bin;
    }

    template<typename T> 
    T miditobin(T midi, T sample_rate, T fft_size)
    {
        T freq = miditofeq(midi);
        return freqtobin(freq,sample_rate, fft_size);
    }
    template<typename T> 
    bool is_power_of_two(uint64_t a)
    {
        if((a & (a-1))==0)
            return true;
        return false;
    }
    template<typename T> 
    uint64_t next_power_of_two(uint64_t a)
    {
        uint64_t i = 1;
        while(i < a) i <<= 1;
        return i;
    }

    template<typename T> 
    T hztomel(T freq)
    {
        T lin_space = (T)3/(T)200;
        T split_hz = (T)1000;
        T split_mel = split_hz * lin_space;
        T log_space = (T)27 / std::log(6.4);
        assert(freq >= (T)0);
        if(freq < split_hz)        
            return freq * lin_space;        
        else
            return split_mel + log_space * std::log(freq/split_hz);
    }
    template<typename T> 
    T hztomel_htk(T freq)
    {        
        T split_hz = (T)700;        
        T log_space = (T)1127;
        assert(freq >= (T)0);
        return log_space * std::log(1 + freq/split_hz);
    }
        
    template<typename T> 
    T meltohz(T mel)
    {
        T lin_space = (T)200/(T)3;
        T split_hz  = (T)1000;
        T split_mel = split_hz / lin_space;
        T log_spacing = std::pow(6.4,1.0/27.0);
        assert(mel >= 0);
        if( mel < split_mel) return lin_space * mel;
        return split_hz * std::pow(log_spacing, mel-split_mel);
    }

    template<typename T> 
    T meltohz_htk(T mel)
    {
        T split_hz = (T)700;
        T log_space = (T)1/(T)1127;
        assert(mel >= 0);
        return split_hz * (std::exp(mel *log_space) -(T)1);
    }

    template<typename T> 
    T db_spl(SampleArray<T> & vector) { return 10 * std::log10(level_lin(vector)); }

    template<typename T> 
    T level_detect(SampleArray<T> & vector,T threshold) { 
        T db = db_spl(vector);
        if(db < threshold) return 1;
        return db;
    }

    template<typename T> 
    size_t size(SampleArray<T> & v) { return v.size(); }
    
    template<typename T> 
    SampleArray<T> random(size_t n, size_t channels=1) { 
        SampleArray<T> r(n,channels);
        r.random();
        return r;
    }

    template<typename T> 
    void fill(SampleArray<T> & v, T value) { v.fill(value); }

    template<typename T> 
    void resize(SampleArray<T> & v, size_t n, size_t channels=1) { v.resize(n,channels); }

    template<typename T> 
    bool    silence_detect(SampleArray<T> & vector,T threshold) { return db_spl(vector) < threshold; }

    template<typename T>
    void set_left_channel(const std::vector<T> & left, SampleArray<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_left_channel(const SampleArray<T> & left, SampleArray<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_left_channel(const SampleArray<T> & left, std::vector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const std::vector<T> & right, SampleArray<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_right_channel(const SampleArray<T> & right, SampleArray<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_right_channel(const SampleArray<T> & right, std::vector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const std::vector<T> & in, SampleArray<T> & out) {
        size_t x = 0;
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    SampleArray<T> interleave(size_t n, size_t channels, const std::vector<SampleArray<T>> & in) {
        std::vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    std::vector<SampleArray<T>> deinterleave(size_t n, size_t channels, const SampleArray<T> & in) {
        std::vector<SampleArray<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }
};