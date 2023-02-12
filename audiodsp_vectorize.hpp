#pragma once

    void vectorize(size_t n, DspFloatType g,DspFloatType * in, DspFloatType * out, DspFloatType (*func)(DspFloatType,DspFloatType))
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) out[i] = func(in[i],g);
    }
    void add(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a[i] + b[i];
    }
    void sub(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a[i] - b[i];
    }
    void mul(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a[i] * b[i];
    }
    void div(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a[i] / b[i];
    }
    void mod(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::fmod(a[i],b[i]);
    }

    void add(size_t n, DspFloatType a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a + b[i];
    }
    void sub(size_t n, DspFloatType a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a - b[i];
    }
    void mul(size_t n, DspFloatType a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a * b[i];
    }
    void div(size_t n, DspFloatType a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a / b[i];
    }
    void div(size_t n, DspFloatType *a, DspFloatType b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = a[i] / b;
    }
    void mod(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::fmod(a[i],b[i]);
    }
    
    void pow(size_t n, DspFloatType * a, DspFloatType b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::pow(a[i],b);
    }
    void pow(size_t n, DspFloatType * a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::pow(a[i],b[i]);
    }
    void pow(size_t n, DspFloatType a, DspFloatType *b, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::pow(a,b[i]);
    }
    void sqrt(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::sqrt(a[i]);
    }
    void sin(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::sin(a[i]);
    }
    void cos(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::cos(a[i]);
    }
    void tan(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::tan(a[i]);
    }
    void asin(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::asin(a[i]);
    }
    void acos(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::acos(a[i]);
    }
    void atan(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::atan(a[i]);
    }
    void sinh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::sinh(a[i]);
    }
    void cosh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::cosh(a[i]);
    }
    void tanh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::tanh(a[i]);
    }
    void asinh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::asinh(a[i]);
    }
    void acosh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::acosh(a[i]);
    }
    void atanh(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::atanh(a[i]);
    }
    void exp(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::exp(a[i]);
    }
    void log(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::log(a[i]);
    }
    void log10(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::log10(a[i]);
    }
    void abs(size_t n, DspFloatType * a, DspFloatType * r) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = std::abs(a[i]);
    }


    void matmul(int m, int n, int k, DspFloatType * a, DspFloatType * b, DspFloatType * r) {
        assert(n == k);    
        #pragma omp simd
        for(size_t i = 0; i < m; i++)
        {        
            for(size_t j = 0; j < n; j++)
            {
                DspFloatType sum=0;
                for(size_t h = 0; h < k; h++ )
                    sum += a[i*n+j] * b[k*n+j];
                r[i*n+j] = sum;
            }
        }
    }

    void transpose(int m, int n, DspFloatType * a, DspFloatType * r)
    {
        #pragma omp simd
        for(size_t i = 0; i < m; i++)
            for(size_t j = 0; j < n; j++)
                r[i*n+j] = a[j*m + i];
    }