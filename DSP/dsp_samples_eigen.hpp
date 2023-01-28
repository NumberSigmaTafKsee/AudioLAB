#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <valarray>

typedef float SampleType; 

//#include "TinyEigen.hpp"

template<typename T>
using EigenArray   = Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>;
template<typename T>
using EigenArray2D = Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

namespace Casino::TinyEigen
{        
    template<typename T>
    struct SampleVector : public EigenArray<T>
    {        
        
        size_t channels;

        SampleVector() {
            channels = 1;
        }
        SampleVector(size_t size, size_t channels = 1) {
            resize(size * channels);
            this->channels = channels;
        }
        SampleVector(const std::vector<T> & s, size_t chans = 1) {
            resize(s.size());
            memcpy(data(),s.data(),s.size()*sizeof(T));
            channels = chans;
        }
        SampleVector(const T * ptr, size_t n, size_t chans = 1) {
            resize(n*chans);
            memcpy(data(),ptr,n*chans*sizeof(T));
            channels = chans;
        }
        SampleVector(const SampleVector<T> & s) {
            (*this) = s;
            channels = s.channels;
        }    
        SampleVector(const EigenArray<T> & v, size_t ch = 1) {
            (*this) = v;
            channels = ch;
        }
        SampleVector(size_t channels, const EigenArray<T>& a) {
            (*this) = a;
            this->channels = channels;
        }

        SampleVector<T> clamp(T min = (T)0.0, T max = (T)1.0) {
            SampleVector<T> r(*this,channels);
            r = SampleVector<T>(r.cwiseMin(max).cwiseMax(min));
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

        size_t size() const { return EigenArray<T>::size(); }
        size_t num_channels() const { return channels; } 
        void zero() { EigenArray<T>::setZero(); }
        void ones() { EigenArray<T>::setOnes(); }
        void random() { EigenArray<T>::setRandom(); }
        void fill(T s) { EigenArray<T>::fill(s); }
        size_t samples_per_channel() const { return size() / num_channels(); }
        T sum() { return EigenArray<T>::sum(); }
        T min() { return EigenArray<T>::minCoeff(); }
        T max() { return EigenArray<T>::maxCoeff(); }
        size_t min_index() { size_t i; EigenArray<T>::minCoeff(&i); return i; }
        size_t max_index() { size_t i; EigenArray<T>::maxCoeff(&i); return i; }
        T* data() { return EigenArray<T>::data(); }
        
        void normalize() { EigenArray<T>::matrix().normalize(); }
        SampleVector<T> normalized() { return SampleVector<T>(EigenArray<T>::matrix().normalized().array()); }

        SampleVector<T>& operator = (const SampleVector<T> & v) {
            (*this) = v;
            channels  = v.channels;
            return *this;
        }    
        T& operator()(size_t i, size_t ch) {
            EigenArray<T> & r = *this;
            return r[i*channels + ch];
        }
        T operator()(size_t i, size_t ch) const {
            const EigenArray<T> & r = *this;
            return r[i*channels + ch];
        }
        T& operator[](size_t i) {
            EigenArray<T> & r = *this;
            return r[i];
        }
        T operator[](size_t i) const {
            const EigenArray<T> & r = *this;
            return r[i];
        }
         
        SampleVector<T> get_channel(size_t channel) {            
            SampleVector<T> r(samples_per_channel());
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
        SampleVector<T> get_channel(size_t channel) const {            
            SampleVector<T> r(samples_per_channel());
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];
            return r;
        }
        void set_channel(const SampleVector<T> & v, size_t ch) {            
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
            EigenArray<T>::resize(n);
        }
        void resize(size_t s, size_t c) {
            channels = c;
            EigenArray<T>::resize(s * c);
        }

        SampleVector<T> get_channel_count(size_t ch, size_t pos, size_t n) {
            SampleVector<T> r(n);
            memcpy(r.data(), get_channel(ch).data()+pos, n*sizeof(T));
            return r;
        }

        bool operator == (const SampleVector<T> & s) {
            return s.channels == channels && size() == s.size();
        }

        void copy(SampleVector<T> & in, size_t block) {
            resize(in.size());
            memcpy(data(),in.data(),block*sizeof(T));
        }
        void copy(SampleVector<T> & in, size_t pos, size_t n) {
            resize((pos+n));
            memcpy(data(),in.data()+pos,n*sizeof(T));
        }
        void copy_from(const T* in, size_t block) {            
            memcpy(data(),in,block*sizeof(T));
        }
        void copy_to(T* in, size_t block) {            
            memcpy(in,data(),block*sizeof(T));
        }

        SampleVector<T> slice(size_t i, size_t n) {
            SampleVector<T> r(n);
            memcpy(r.data(), data()+i, n * sizeof(T));
            return r;
        }
        
        SampleVector<T> pan(SampleType pan) {
            SampleType pan_mapped = ((pan+1.0)/2.0) * (M_PI/2.0);
            SampleVector<T> r(size(),2);
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
        SampleVector<T> stride_slice(const SampleVector<T> & in, size_t stride, size_t pos, size_t n)
        {
            SampleVector<T> r(n/stride);
            size_t x = pos;
            for(size_t i = 0; i < n; i++) {
                r[i] = (*this)[x];
                x += stride;
            } 
            return r;
        }
        SampleVector<T> operator + (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r += s;
            return r;
        }
        SampleVector<T> operator - (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r -= s;
            return r;
        }
        SampleVector<T> operator * (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r *= s;
            return r;
        }
        SampleVector<T> operator / (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r /= s;
            return r;
        }
        
        SampleVector<T> operator + (const T s) {    
            SampleVector<T> r(*this);
            r += s;
            return r;
        }
        SampleVector<T> operator - (const T s) {        
            SampleVector<T> r(*this);
            r -= s;
            return r;
        }
        SampleVector<T> operator / (const T s) {        
            SampleVector<T> r(*this);
            r /= s;
            return r;
        }
        SampleVector<T> operator * (const T s) {        
            SampleVector<T> r(*this);
            r *= s;
            return r;
        }
        

        SampleVector<T>& operator += (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) += s;
            return *this;
        }
        SampleVector<T>& operator -= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) -= s;
            return *this;
        }
        SampleVector<T>& operator *= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) *= s;
            return *this;
        }
        SampleVector<T>& operator /= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) /= s;
            return *this;
        }
        
        SampleVector<T> operator += (const T s) {
            (*this) += s;
            return *this;
        }
        SampleVector<T> operator -= (const T s) {
            (*this) -= s;
            return *this;
        }
        SampleVector<T> operator /= (const T s) {
            (*this) /= s;
            return *this;
        }
        SampleVector<T> operator *= (const T s) {
            (*this) *= s;
            return *this;
        }
        
        
        SampleVector<T> abs() { 
            SampleVector<T> r(*this);
            r = this->abs();
            return r;
        }
        SampleVector<T> abs2() { 
            SampleVector<T> r(*this);
            r = this->abs2();
            return r;
        }
        SampleVector<T> inverse() { 
            SampleVector<T> r(*this);
            r = this->inverse();
            return r;
        }
        SampleVector<T> exp() { 
            SampleVector<T> r(*this);
            r = this->exp();
            return r;
        }
        SampleVector<T> log() { 
            SampleVector<T> r(*this);
            r = this->log();
            return r;
        }
        SampleVector<T> log1p() { 
            SampleVector<T> r(*this);
            r = this->log1p();
            return r;
        }
        SampleVector<T> log10() { 
            SampleVector<T> r(*this);
            r = this->log10();
            return r;
        }
        SampleVector<T> pow(const SampleVector<T> & s) { 
            SampleVector<T> r(*this);
            r = this->pow(s);
            return r;
        }
        SampleVector<T> pow(const T s) { 
            SampleVector<T> r(*this);
            r = this->pow(s);        
            return r;        
        }
        SampleVector<T> sqrt() {
            SampleVector<T> r(*this);
            r = this->sqrt();
            return r;
        }
        SampleVector<T> rsqrt() {
            SampleVector<T> r(*this);
            r = this->rsqrt();
            return r;
        }
        SampleVector<T> square() {
            SampleVector<T> r(*this);
            r = this->square();
            return r;
        }
        SampleVector<T> sin() {
            SampleVector<T> r(*this);
            r = this->sin();
            return r;
        }
        SampleVector<T> cos() {
            SampleVector<T> r(*this);
            r = this->cos();
            return r;
        }
        SampleVector<T> tan() {
            SampleVector<T> r(*this);
            r = this->tan();
            return r;
        }
        SampleVector<T> asin() {
            SampleVector<T> r(*this);
            r = this->asin();
            return r;
        }
        SampleVector<T> acos() {
            SampleVector<T> r(*this);
            r = this->acos();
            return r;
        }
        SampleVector<T> atan() {
            SampleVector<T> r(*this);
            r = this->atan();
            return r;
        }        
        SampleVector<T> sinh() {
            SampleVector<T> r(*this);
            r = this->sinh();
            return r;
        }
        SampleVector<T> cosh() {
            SampleVector<T> r(*this);
            r = this->cosh();
            return r;
        }
        SampleVector<T> tanh() {
            SampleVector<T> r(*this);
            r = this->tanh();
            return r;
        }
        SampleVector<T> ceil() { 
            SampleVector<T> r(*this);
            r = this->ceil();
            return r;        
        }
        SampleVector<T> floor() { 
            SampleVector<T> r(*this);
            r = this->floor();
            return r;
        }
        SampleVector<T> round() { 
            SampleVector<T> r(*this);
            r = this->round();
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
    struct SampleMatrix : public EigenArray2D<T>
    {        
        SampleMatrix() = default;
        
        SampleMatrix(size_t chan, size_t samps) {
            resize(chan,samps);
        }
        SampleMatrix(const SampleVector<T> & v) {
            resize(v.num_channels(), v.size()/v.num_channels());
            for(size_t i = 0; i < v.num_channels(); i++)             
                this->row(i) = v.get_channel(i);
        }
        SampleMatrix(const EigenArray2D<T> & c) {
            (*this) = c;
        }
        SampleMatrix(const SampleMatrix<T> & m) {
            (*this) = m.channels;
        }

        T * data() { return EigenArray2D<T>::data(); }
        
        void deinterleave(const SampleVector<T> & s) {
            
            esize(s.num_channels(),s.samples_per_channel());            
            
            for(size_t i = 0; i < s.num_channels(); i++) {
                this->row(i) = s.get_channel(i);
            }
            
        }    

        SampleVector<T> interleave() { 
            SampleVector<T> r;
            r.resize(cols() * cols());
            for(size_t i = 0; i < rows(); i++)            
            {
                SampleVector<T> tmp(1,(*this)(i));            
                r.set_channel(tmp,i);
            }
            return r;
        }

        size_t num_channels() const { return EigenArray2D<T>::rows(); }
        size_t rows() const { return EigenArray2D<T>::rows(); }
        size_t cols() const { return EigenArray2D<T>::cols(); }
        size_t size() const { return rows()*cols(); }
        void   resize(size_t r, size_t c) { EigenArray2D<T>::resize(r,c); }
        
        SampleVector<T>& row(size_t i) {
            return EigenArray2D<T>::row(i);
        }
        SampleVector<T>& col(size_t i) {
            return EigenArray2D<T>::col(i);
        }
        SampleVector<T> operator[](size_t ch) { 
            SampleVector<T> tmp(1,row(ch));
            return tmp;
        }
        SampleVector<T> get_channel(size_t c) { 
            SampleVector<T> tmp(1,row(c));
            return tmp;
        }
        void set_channel(size_t c, const SampleVector<T> & s) 
        { 
            row(c) = s;
        }

        SampleMatrix<T>& operator = (const SampleMatrix<T> & v) {        
            (*this)  = v;
            return *this;
        }
        T& operator()(size_t c, size_t n) {
            return (*this)(c,n);
        }
        bool operator == (const SampleMatrix<T> & b) {
            return rows() == b.rows() && cols() == b.cols();
        }
        SampleMatrix<T> operator + (const SampleMatrix<T> & m) {        
            SampleMatrix<T> r(*this);
            r = r + m;
            return r;
        }
        SampleMatrix<T> operator - (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r = r - m;
            return r;
        }
        SampleMatrix<T> operator * (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r = (r.matrix() * m.matrix()).array();
            return r;
        }
        SampleMatrix<T> operator / (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r =  r / m;
            return r;
        }
        
        SampleMatrix<T> operator + (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r + s;        
            return r;
        }
        SampleMatrix<T> operator - (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r - s;
            return r;
        }
        SampleMatrix<T> operator * (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r * s;
            return r;
        }
        SampleMatrix<T> operator / (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r / s;
            return r;
        }
        

        SampleMatrix<T> abs() { 
            SampleMatrix<T> r(*this);        
            r = abs();
            return r;
        }
        SampleMatrix<T> abs2() { 
            SampleMatrix<T> r(*this);        
            r = abs2();
            return r;
        }
        SampleMatrix<T> inverse() { 
            SampleMatrix<T> r(*this);        
            r = inverse();
            return r;
        }
        SampleMatrix<T> exp() { 
            SampleMatrix<T> r(*this);        
            r = exp();
            return r;
        }
        SampleMatrix<T> log() { 
            SampleMatrix<T> r(*this);        
            r = log();
            return r;
        }
        SampleMatrix<T> log1p() { 
            SampleMatrix<T> r(*this);        
            r = log1p();
            return r;
        }
        SampleMatrix<T> log10() { 
            SampleMatrix<T> r(*this);        
            r = log10();
            return r;
        }
        SampleMatrix<T> pow(const T s) { 
            SampleMatrix<T> r(*this);        
            r = pow(s);
            return r;
        }
        SampleMatrix<T> sqrt() {
            SampleMatrix<T> r(*this);        
            r = sqrt();
            return r;
        }
        SampleMatrix<T> rsqrt() {
            SampleMatrix<T> r(*this);        
            r = rsqrt();
            return r;
        }
        SampleMatrix<T> square() {
            SampleMatrix<T> r(*this);        
            r = square();
            return r;
        }
        SampleMatrix<T> sin() {
            SampleMatrix<T> r(*this);        
            r = sin();
            return r;
        }
        SampleMatrix<T> cos() {
            SampleMatrix<T> r(*this);        
            r = cos();
            return r;
        }
        SampleMatrix<T> tan() {
            SampleMatrix<T> r(*this);        
            r = tan();
            return r;
        }
        SampleMatrix<T> asin() {
            SampleMatrix<T> r(*this);        
            r = asin();
            return r;
        }
        SampleMatrix<T> acos() {
            SampleMatrix<T> r(*this);
            r = acos();
            return r;
        }
        SampleMatrix<T> atan() {
            SampleMatrix<T> r(*this);
            r = atan();
            return r;
        }
        SampleMatrix<T> sinh() {
            SampleMatrix<T> r(*this);        
            r = sinh();
            return r;
        }
        SampleMatrix<T> cosh() {
            SampleMatrix<T> r(*this);        
            r = cosh();
            return r;
        }
        SampleMatrix<T> tanh() {
            SampleMatrix<T> r(*this);        
            r = tanh();
            return r;
        }
        SampleMatrix<T> ceil() {
            SampleMatrix<T> r(*this);        
            r = ceil();
            return r;
        }
        SampleMatrix<T> floor() {
            SampleMatrix<T> r(*this);        
            r = floor();
            return r;
        }
        SampleMatrix<T> round() {
            SampleMatrix<T> r(*this);        
            r = round();
            return r;
        }
    };

    template<typename T>  SampleVector<T> abs( SampleVector<T> & m) { return m.abs(); }
    template<typename T>  SampleVector<T> abs2( SampleVector<T> & m) { return m.abs2(); }
    template<typename T>  SampleVector<T> inverse( SampleVector<T> & m) { return m.inverse(); }
    template<typename T>  SampleVector<T> exp( SampleVector<T> & m) { return m.exp(); }
    template<typename T>  SampleVector<T> log( SampleVector<T> & m) { return m.log(); }
    template<typename T>  SampleVector<T> log1p( SampleVector<T> & m) { return m.log1p(); }
    template<typename T>  SampleVector<T> log10( SampleVector<T> & m) { return m.log10(); }
    template<typename T>  SampleVector<T> pow( SampleVector<T> & m,  SampleVector<T> & p) { return m.pow(p); }
    template<typename T>  SampleVector<T> pow( SampleVector<T> & m,  T p) { return m.pow(p); }
    template<typename T>  SampleVector<T> sqrt( SampleVector<T> & m) { return m.sqrt(); }
    template<typename T>  SampleVector<T> rsqrt( SampleVector<T> & m) { return m.rsqrt(); }
    template<typename T>  SampleVector<T> square( SampleVector<T> & m) { return m.square(); }
    template<typename T>  SampleVector<T> sin( SampleVector<T> & m) { return m.sin(); }
    template<typename T>  SampleVector<T> cos( SampleVector<T> & m) { return m.cos(); }
    template<typename T>  SampleVector<T> tan( SampleVector<T> & m) { return m.tan(); }
    template<typename T>  SampleVector<T> asin( SampleVector<T> & m) { return m.asin(); }
    template<typename T>  SampleVector<T> acos( SampleVector<T> & m) { return m.acos(); }
    template<typename T>  SampleVector<T> atan( SampleVector<T> & m) { return m.atan(); }
    template<typename T>  SampleVector<T> sinh( SampleVector<T> & m) { return m.sinh(); }
    template<typename T>  SampleVector<T> cosh( SampleVector<T> & m) { return m.cosh(); }
    template<typename T>  SampleVector<T> tanh( SampleVector<T> & m) { return m.tanh(); }
    template<typename T>  SampleVector<T> ceil( SampleVector<T> & m) { return m.ceil(); }
    template<typename T>  SampleVector<T> floor( SampleVector<T> & m) { return m.floor(); }
    template<typename T>  SampleVector<T> round( SampleVector<T> & m) { return m.round(); }

    template<typename T>  SampleMatrix<T> abs( SampleMatrix<T> & m) { return m.abs(); } 
    template<typename T>  SampleMatrix<T> abs2( SampleMatrix<T> & m) { return m.abs2(); }
    template<typename T>  SampleMatrix<T> inverse( SampleMatrix<T> & m) { return m.inverse(); }
    template<typename T>  SampleMatrix<T> exp( SampleMatrix<T> & m) { return m.exp(); }
    template<typename T>  SampleMatrix<T> log( SampleMatrix<T> & m) { return m.log(); }
    template<typename T>  SampleMatrix<T> log1p( SampleMatrix<T> & m) { return m.log1p(); }
    template<typename T>  SampleMatrix<T> log10( SampleMatrix<T> & m) { return m.log10(); }
    template<typename T>  SampleMatrix<T> pow( SampleMatrix<T> & m,  T &p) { return m.pow(p); }
    template<typename T>  SampleMatrix<T> sqrt( SampleMatrix<T> & m) { return m.sqrt(); }
    template<typename T>  SampleMatrix<T> rsqrt( SampleMatrix<T> & m) { return m.rsqrt(); }
    template<typename T>  SampleMatrix<T> square( SampleMatrix<T> & m) { return m.square(); }
    template<typename T>  SampleMatrix<T> sin( SampleMatrix<T> & m) { return m.sin(); }
    template<typename T>  SampleMatrix<T> cos( SampleMatrix<T> & m) { return m.cos(); }
    template<typename T>  SampleMatrix<T> tan( SampleMatrix<T> & m) { return m.tan(); }
    template<typename T>  SampleMatrix<T> asin( SampleMatrix<T> & m) { return m.asin(); }
    template<typename T>  SampleMatrix<T> acos( SampleMatrix<T> & m) { return m.acos(); }
    template<typename T>  SampleMatrix<T> atan( SampleMatrix<T> & m) { return m.atan(); }
    template<typename T>  SampleMatrix<T> sinh( SampleMatrix<T> & m) { return m.sinh(); }
    template<typename T>  SampleMatrix<T> cosh( SampleMatrix<T> & m) { return m.cosh(); }
    template<typename T>  SampleMatrix<T> tanh( SampleMatrix<T> & m) { return m.tanh(); }
    template<typename T>  SampleMatrix<T> ceil( SampleMatrix<T> & m) { return m.ceil(); }
    template<typename T>  SampleMatrix<T> floor( SampleMatrix<T> & m) { return m.floor(); }
    template<typename T>  SampleMatrix<T> round( SampleMatrix<T> & m) { return m.round(); }    
};