#pragma once
#include <valarray>

template<typename T>
struct ArrayVector
{
    // this will be replaced with eigen (or arrayfire, or another backend) later.
    std::valarray<T> vector;
    size_t channels;

    ArrayVector() {
        channels = 1;
    }
    ArrayVector(size_t size, size_t channels = 1) {
        vector.resize(size * channels);
        this->channels = channels;
    }
    ArrayVector<T>(const std::vector<T> & s, size_t chans = 1) {
        vector = s;
        channels = chans;
    }
    ArrayVector<T>(const ArrayVector<T><T> & s) {
        vector = s.vector;
        channels = s.channels;
    }

    size_t size() const { return vector.size(); }
    size_t num_channels() const { return channels; } 

    ArrayVector<T>& operator = (const ArrayVector<T> & v) {
        vector = v.vector;
        channels  = v.channels;
        return *this;
    }
    ArrayType& operator()(size_t i, size_t ch) {
        return vector[i*channels + ch];
    }

    ArrayVector<T> get_channel(size_t channel) {
        //ArrayVector<T> r(*this);
        //if(channels == 1) return r;
        ArrayVector<T> r(*this);
        r.vector = vector[std::slice(0,vector.size(),channel-1)];
        return r;
    }
    void set_channel(const ArrayVector<T> & v, size_t ch) {
        //size_t x = 0;
        //for(size_t i = ch; i < vector.size(); i+=channels) {
        //    vector[i] = v[x++];
        //}
        std::valarray<std::size_t> idx;
        idx.resize(size()/num_channels);
        size_t x = 0;
        for(size_t i = ch-1; i < size(); i+=channels) idx[x++] = i;
        vector[idx] = v.vector;
    }

    bool operator == (const ArrayVector<T> & s) {
        return s.channels == channels && vector.size() == s.vector.size();
    }

    ArrayVector<T> operator + (const ArrayVector<T> &s) {
        assert(*this == s);
        ArrayVector<T> r(*this);
        r.vector += s.vector;
        return r;
    }
    ArrayVector<T> operator - (const ArrayVector<T> &s) {
        assert(*this == s);
        ArrayVector<T> r(*this);
        r.vector -= s.vector;
        return r;
    }
    ArrayVector<T> operator * (const ArrayVector<T> &s) {
        assert(*this == s);
        ArrayVector<T> r(*this);
        r.vector *= s.vector;
        return r;
    }
    ArrayVector<T> operator / (const ArrayVector<T> &s) {
        assert(*this == s);
        ArrayVector<T> r(*this);
        r.vector /= s.vector;
        return r;
    }
    ArrayVector<T> operator % (const ArrayVector<T> &s) {
        assert(*this == s);
        ArrayVector<T> r(*this);
        r.vector %= s.vector;
        return r;
    }

    ArrayVector<T> operator + (const ArrayType s) {    
        ArrayVector<T> r(*this);
        r.vector += s;
        return r;
    }
    ArrayVector<T> operator - (const ArrayType s) {        
        ArrayVector<T> r(*this);
        r.vector -= s;
        return r;
    }
    ArrayVector<T> operator / (const ArrayType s) {        
        ArrayVector<T> r(*this);
        r.vector /= s;
        return r;
    }
    ArrayVector<T> operator * (const ArrayType s) {        
        ArrayVector<T> r(*this);
        r.vector *= s;
        return r;
    }
    ArrayVector<T> operator % (const ArrayType s) {        
        ArrayVector<T> r(*this);
        r.vector %= s;
        return r;
    }


    ArrayVector<T>& operator += (const ArrayVector<T>& s) {        
        assert(*this == s);
        vector += s.vector;
        return *this;
    }
    ArrayVector<T>& operator -= (const ArrayVector<T>& s) {        
        assert(*this == s);
        vector -= s.vector;
        return *this;
    }
    ArrayVector<T>& operator *= (const ArrayVector<T>& s) {        
        assert(*this == s);
        vector *= s.vector;
        return *this;
    }
    ArrayVector<T>& operator /= (const ArrayVector<T>& s) {        
        assert(*this == s);
        vector /= s.vector;
        return *this;
    }
    ArrayVector<T>& operator %= (const ArrayVector<T>& s) {        
        assert(*this == s);
        vector %= s.vector;
        return *this;
    }

    ArrayVector<T> operator += (const ArrayType s) {
        vector += s;
        return *this;
    }
    ArrayVector<T> operator -= (const ArrayType s) {
        vector -= s;
        return *this;
    }
    ArrayVector<T> operator /= (const ArrayType s) {
        vector /= s;
        return *this;
    }
    ArrayVector<T> operator *= (const ArrayType s) {
        vector *= s;
        return *this;
    }
    ArrayVector<T> operator %= (const ArrayType s) {        
        vector %= s;
        return *this;
    }

    ArrayType sum() { return vector.sum(); }
    ArrayType min() { return vector.min(); }
    ArrayType max() { return vector.max(); }

    void resize(size_t s, size_t c) {
        channels = c;
        vector.resize(s * c);
    }

    ArrayVector<T> abs() { 
        ArrayVector<T> r(*this);
        r.vector = std::abs(vector);
        return r;
    }
    ArrayVector<T> exp() { 
        ArrayVector<T> r(*this);
        r.vector = std::exp(vector);
        return r;
    }
    ArrayVector<T> log() { 
        ArrayVector<T> r(*this);
        r.vector = std::log(vector);
        return r;
    }
    ArrayVector<T> log10() { 
        ArrayVector<T> r(*this);
        r.vector = std::log10(vector);
        return r;
    }
    ArrayVector<T> pow(const ArrayVector<T> & s) { 
        ArrayVector<T> r(*this);
        r.vector = std::pow(vector,s.vector);
        return r;
    }
    ArrayVector<T> pow(const ArrayType s) { 
        ArrayVector<T> r(*this);
        r.vector = std::pow(vector,s);        
        return r;        
    }
    ArrayVector<T> sqrt() {
        ArrayVector<T> r(*this);
        r.vector = std::sqrt(vector);
        return r;
    }
    ArrayVector<T> sin() {
        ArrayVector<T> r(*this);
        r.vector = std::sin(vector);
        return r;
    }
    ArrayVector<T> cos() {
        ArrayVector<T> r(*this);
        r.vector = std::cos(vector);
        return r;
    }
    ArrayVector<T> tan() {
        ArrayVector<T> r(*this);
        r.vector = std::tan(vector);
        return r;
    }
    ArrayVector<T> asin() {
        ArrayVector<T> r(*this);
        r.vector = std::asin(vector);
        return r;
    }
    ArrayVector<T> acos() {
        ArrayVector<T> r(*this);
        r.vector = std::acos(vector);
        return r;
    }
    ArrayVector<T> atan() {
        ArrayVector<T> r(*this);
        r.vector = std::atan(vector);
        return r;
    }
    ArrayVector<T> atan2(const ArrayType s) {
        ArrayVector<T> r(*this);
        r.vector = std::atan2(vector,s);
        return r;
    }
    ArrayVector<T> sinh() {
        ArrayVector<T> r(*this);
        r.vector = std::sinh(vector);
        return r;
    }
    ArrayVector<T> cosh() {
        ArrayVector<T> r(*this);
        r.vector = std::cosh(vector);
        return r;
    }
    ArrayVector<T> tanh() {
        ArrayVector<T> r(*this);
        r.vector = std::tanh(vector);
        return r;
    }
};

template<typename T>
struct ArrayMatrix {
    std::vector<ArrayVector<T>> channels;
    ArrayMatrix<T>(size_t chan, size_t samps) {
        channels.resize(chan);
        for(size_t i = 0; i < chan; i++)
            channels[i].resize(samps);
    }
    ArrayMatrix<T>(const ArrayVector<T> & v) {
        channels.resize(v.num_channels());
        for(size_t i = 0; i < v.num_channels(); i++) 
            channels[i] = v.get_channel(i);
    }

    ArrayVector<T> deinterleave() { 
        ArrayVector<T> r;
        r.resize(channels[0].size() * channels.size());
        for(size_t i = 0; i < channels.size(); i++)
            r.set_channel(channels[i]);
        return r;
    }

    size_t num_channels() const { return channels.size(); }
    ArrayVector<T>& operator[](size_t ch) { return channels[ch]; }
    ArrayVector<T>& get_channel(size_t c) { return channels[c]; }
    void set_channel(size_t c, const ArrayVector<T> & s) { channels[c] = s; }

    ArrayMatrix<T>& operator = (const ArrayMatrix<T> & v) {        
        channels  = v.channels;
        return *this;
    }
    ArrayType& operator()(size_t c, size_t n) {
        return channels[c][n];
    }

    bool check(const ArrayMatrix<T> & b) {
        if(channels.size() != b.channels.size()) return false;
        for(size_t i = 0; i < channels.size(); i++)
            if(channels[i].size() != b.channels[i].size()) return false;
        return true;
    }
    bool operator == (const ArrayMatrix<T> & b) {
        return check(b);
    }
    ArrayMatrix<T> operator + (const ArrayMatrix<T> & m) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] +  m.channels[i];
        return r;
    }
    ArrayMatrix<T> operator - (const ArrayMatrix<T> & m) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] -  m.channels[i];
        return r;
    }
    ArrayMatrix<T> operator * (const ArrayMatrix<T> & m) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] *  m.channels[i];
        return r;
    }
    ArrayMatrix<T> operator / (const ArrayMatrix<T> & m) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] / m.channels[i];
        return r;
    }
    ArrayMatrix<T> operator % (const ArrayMatrix<T> & m) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] % m.channels[i];
        return r;
    }

    ArrayMatrix<T> operator + (const ArrayType s) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] + s;
        return r;
    }
    ArrayMatrix<T> operator - (const ArrayType s) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] - s;
        return r;
    }
    ArrayMatrix<T> operator * (const ArrayType s) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] * s;
        return r;
    }
    ArrayMatrix<T> operator / (const ArrayType s) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] / s;
        return r;
    }
    ArrayMatrix<T> operator % (const ArrayType s) {
        assert(*this == m);
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < m.channels.size(); i++)
            r.channels[i] = channels[i] % s;
        return r;
    }


     ArrayMatrix<T> abs() { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].abs();
        return r;
    }
    ArrayMatrix<T> exp() { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].exp();
        return r;
    }
    ArrayMatrix<T> log() { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].log();
        return r;
    }
    ArrayMatrix<T> log10() { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].log10();
        return r;
    }
    ArrayMatrix<T> pow(const ArrayVector<T> & s) { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].pow(s);
        return r;
    }
    ArrayMatrix<T> pow(const ArrayType s) { 
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].pow(s);
        return r;
    }
    ArrayMatrix<T> sqrt() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].sqrt();
        return r;
    }
    ArrayMatrix<T> sin() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].sin();
        return r;
    }
    ArrayMatrix<T> cos() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].cos();
        return r;
    }
    ArrayMatrix<T> tan() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].tan();
        return r;
    }
    ArrayMatrix<T> asin() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].asin();
        return r;
    }
    ArrayMatrix<T> acos() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].acos();
        return r;
    }
    ArrayMatrix<T> atan() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].atan();
        return r;
    }
    ArrayMatrix<T> atan2(const ArrayType s) {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].atan2(s);
        return r;
    }
    ArrayMatrix<T> sinh() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].sinh();
        return r;
    }
    ArrayMatrix<T> cosh() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].cosh();
        return r;
    }
    ArrayMatrix<T> tanh() {
        ArrayMatrix<T> r(*this);
        for(size_t i = 0; i < channels.size(); i++)
            r.channels[i] = channels[i].tanh();
        return r;
    }
};

inline ArrayMatrix<T> interleave(const std::vector<ArrayVector<T>> & c) {
    ArrayMatrix<T> m(c.size(),c[0].size);
    for(size_t i = 0; i < c.size(); i++)
        m.set_channel(i,c[i]);
    return m;
}

