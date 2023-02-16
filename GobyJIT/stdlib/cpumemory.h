#pragma once
#include <complex.h>
#include <cstdint>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <assert.h>

typedef struct
{
    void * ptr;
    int    n_elements;
    int    e_size;
}
CPUMemory;

//typedef float _Complex scomplex;
//typedef double _Complex dcomplex;

#ifdef __cplusplus
extern "C" {
#endif 


CPUMemory* cm_new(int num, int element_size);
void       cm_free(CPUMemory * mem);

float  cm_get_float(CPUMemory * p, int index);
double cm_get_double(CPUMemory * p, int index);
//scomplex cm_get_scomplex(CPUMemory * p, int index);
//dcomplex cm_get_domplex(CPUMemory * p, int index);
uint8_t cm_get_byte(CPUMemory * p, int index);
uint16_t cm_get_word(CPUMemory * p, int index);
uint32_t cm_get_dword(CPUMemory * p, int index);
uint64_t cm_get_qword(CPUMemory * p, int index);
int8_t  cm_get_int8(CPUMemory * p, int index);
int16_t cm_get_int16(CPUMemory * p, int index);
int32_t cm_get_int32(CPUMemory * p, int index);
int64_t cm_get_int64(CPUMemory * p, int index);


void cm_set_byte(CPUMemory * p, int index, uint8_t value);
void cm_set_word(CPUMemory * p, int index, uint16_t value);
void cm_set_dword(CPUMemory * p, int index, uint32_t value);
void cm_set_qword(CPUMemory * p, int index, uint64_t value);
void cm_set_uint8(CPUMemory * p, int index, uint8_t value);
void cm_set_uint16(CPUMemory * p, int index, uint16_t value);
void cm_set_uint32(CPUMemory * p, int index, uint32_t value);
void cm_set_uint64(CPUMemory * p, int index, uint64_t value);
void cm_set_int8(CPUMemory * p, int index,  int8_t value);
void cm_set_int16(CPUMemory * p, int index, int16_t value);
void cm_set_int32(CPUMemory * p, int index, int32_t value);
void cm_set_int64(CPUMemory * p, int index, int64_t value);

void cm_set_float(CPUMemory * p, int index, float value);
void cm_set_double(CPUMemory * p, int index, double value);
//void cm_set_scomplex(CPUMemory * p, int index, scomplex value);
//void cm_set_dcomplex(CPUMemory * p, int index, dcomplex value);

#ifdef __cplusplus
} 
#endif 


template<typename T>
struct CpuMemory 
{
    CPUMemory * mem;

    CpuMemory(int size) {
        mem = cm_new(size, sizeof(T));
    }
    CpuMemory(int size, int element_size) {
        mem = cm_new(size,element_size);
    }
    ~CpuMemory() {
        if(mem) cm_free(mem);
    }

    T& operator[](size_t index) { return ((T*)mem->ptr)[index]; }
    T     __getitem(size_t index) { return ((T*)mem->ptr)[index]; }
    void  __setitem(size_t index, const T value) { ((T*)mem->ptr)[index] = value; }

    double   get_double(size_t i) { return cm_get_double(mem,i); }
    float    get_float(size_t i) { return cm_get_float(mem,i); }
    //scomplex get_scomplex(size_t i) { return cm_get_scomplex(mem,i); }
    //dcomplex get_dcomplex(size_t i) { return cm_get_dcomplex(mem,i); }
    uint8_t  get_byte(size_t i) { return cm_get_byte(mem,i); }
    uint16_t get_word(size_t i) { return cm_get_word(mem,i); }
    uint32_t get_dword(size_t i) { return cm_get_dword(mem,i); }
    uint64_t get_qword(size_t i) { return cm_get_qword(mem,i); }
    int8_t   get_int8(size_t i) { return cm_get_int8(mem,i); }
    int16_t  get_int16(size_t i) { return cm_get_int16(mem,i); }
    int32_t  get_int32(size_t i) { return cm_get_int32(mem,i); }
    int64_t  get_int64(size_t i) { return cm_get_int64(mem,i); }

    void    set_double(size_t i, double v) { cm_set_double(mem,i,v); }
    void    set_float(size_t i, float v) { cm_set_float(mem,i,v); }
    //void    set_scomplex(size_t i, scomplex v) { cm_set_scomplex(mem,i,v); }
    //void    set_dcomplex(size_t i, dcomplex v) { cm_set_dcomplex(mem,i,v); }
    void    set_byte(size_t i, uint8_t v) { cm_set_byte(mem,i,v); }
    void    set_word(size_t i, uint16_t v) { cm_set_word(mem,i,v); }
    void    set_dword(size_t i, uint32_t v) { cm_set_dword(mem,i,v); }
    void    set_qword(size_t i, uint64_t v) { cm_set_qword(mem,i,v); }
    void    set_int8(size_t i, int8_t v) { cm_set_int8(mem,i,v); }
    void    set_int16(size_t i, int16_t v) { cm_set_int16(mem,i,v); }
    void    set_int32(size_t i, int32_t v) { cm_set_int32(mem,i,v); }
    void    set_int64(size_t i, int64_t v) { cm_set_int64(mem,i,v); }

};



inline int sizeof_float() { return sizeof(float); }
inline int sizeof_double() { return sizeof(double); }
//inline int sizeof_scomplex() { return sizeof(scomplex); }
//inline int sizeof_dcomplex() { return sizeof(dcomplex); }
inline int sizeof_byte() { return sizeof(uint8_t); }
inline int sizeof_word() { return sizeof(uint16_t); }
inline int sizeof_dword() { return sizeof(uint32_t); }
inline int sizeof_qword() { return sizeof(uint64_t); }
inline int sizeof_int8() { return sizeof(int8_t); }
inline int sizeof_int16() { return sizeof(int16_t); }
inline int sizeof_int32() { return sizeof(uint32_t); }
inline int sizeof_int64() { return sizeof(uint64_t); }
