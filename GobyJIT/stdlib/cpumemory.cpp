#include "cpumemory.h"

CPUMemory* cm_new(int num, int element_size) {
    CPUMemory * p = (CPUMemory*)calloc(1,sizeof(CPUMemory));
    assert(p!=NULL);
    p->ptr = calloc(num,element_size);
    assert(p->ptr != NULL);
    p->n_elements = num;
    p->e_size     = element_size;
    return p;
}
void cm_free(CPUMemory * mem) {
    assert(mem != NULL);
    if(mem->ptr) free(mem->ptr);
    free(mem);
}
float  cm_get_float(CPUMemory * p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((float*)p->ptr)[index];
}
double  cm_get_double(CPUMemory * p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((double*)p->ptr)[index];
}
/*
scomplex  cm_get_scomplex(CPUMemory * p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((scomplex*)p->ptr)[index];
}
dcomplex  cm_get_dcomplex(CPUMemory * p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((dcomplex*)p->ptr)[index];
}
*/
uint8_t cm_get_byte(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((uint8_t*)p->ptr)[index];
}
uint16_t cm_get_word(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((uint16_t*)p->ptr)[index];
}
uint32_t cm_get_dword(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((uint32_t*)p->ptr)[index];
}
uint64_t cm_get_qword(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((uint64_t*)p->ptr)[index];
}

int8_t cm_get_int8(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((int8_t*)p->ptr)[index];
}
int16_t cm_get_int16(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((int16_t*)p->ptr)[index];
}
int32_t cm_get_int32(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((int32_t*)p->ptr)[index];
}
int64_t cm_get_int64(CPUMemory *p, int index) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    return ((int64_t*)p->ptr)[index];
}

void cm_set_float(CPUMemory * p, int index, float value) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((float*)p->ptr)[index]=value;
}
void cm_set_double(CPUMemory * p, int index, double value) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((double*)p->ptr)[index]=value;
}
/*
void cm_set_scomplex(CPUMemory * p, int index, scomplex value) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((scomplex*)p->ptr)[index]=value;
}
void cm_set_dcomplex(CPUMemory * p, int index, dcomplex value) {
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((dcomplex*)p->ptr)[index]=value;
}
*/


void cm_set_byte(CPUMemory * p, int index, uint8_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint8_t*)p->ptr)[index]=value;
}
void cm_set_word(CPUMemory * p, int index, uint16_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint16_t*)p->ptr)[index]=value;
}
void cm_set_dword(CPUMemory * p, int index, uint32_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint32_t*)p->ptr)[index]=value;
}
void cm_set_qword(CPUMemory * p, int index, uint64_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint64_t*)p->ptr)[index]=value;
}
void cm_set_uint8(CPUMemory * p, int index, uint8_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint8_t*)p->ptr)[index]=value;
}
void cm_set_uint16(CPUMemory * p, int index, uint16_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint16_t*)p->ptr)[index]=value;
}
void cm_set_uint32(CPUMemory * p, int index, uint32_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint32_t*)p->ptr)[index]=value;
}
void cm_set_uint64(CPUMemory * p, int index, uint64_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((uint64_t*)p->ptr)[index]=value;
}
void cm_set_int8(CPUMemory * p, int index,  int8_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((int8_t*)p->ptr)[index]=value;
}
void cm_set_int16(CPUMemory * p, int index, int16_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((int16_t*)p->ptr)[index]=value;
}
void cm_set_int32(CPUMemory * p, int index, int32_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((int32_t*)p->ptr)[index]=value;
}
void cm_set_int64(CPUMemory * p, int index, int64_t value)
{
    assert(p != NULL);
    assert(p->ptr != NULL);
    while(index < 0) index += p->n_elements;
    assert(index < p->n_elements);
    ((int64_t*)p->ptr)[index]=value;
}

