#include "bitvector.h"
#include "cpumemory.h"
#include <stdio.h>
#include <stlib.h>
#include <string.h>
#include <assert.h>

BitVector* bv_encode_memory(CPUMemory * p) {
    BitVector* ptr = bv_ptr(p->ptr, p->n_elements * p->e_size * 8);
    assert(ptr != NULL);
    return ptr;
}
CPUMemory* bv_decode_memory(BitVector * bits) {
    CPUMemory * p = cm_new(bits->num_bits/sizeof(double), sizeof(double));
    assert(p != NULL);
    return p;
}

BitVector * bv_new(int num_bits) {
    BitVector * p = (BitVector*)calloc(1,sizeof(BitVector));
    assert(p != NULL);
    p->bits = (uint8_t*)calloc(num_bits/8+1,sizeof(uint8_t));
    assert(p->bits != NULL);
    p->size = num_bits/8+1;
    p->num_bits = num_bits;
    return p;
}

BitVector * bv_ptr(void * ptr, int size) {
    BitVector * p = (BitVector*)calloc(1,sizeof(BitVector));
    assert(p != NULL);
    int num_bits = size*sizeof(uint8_t);
    p->bits = (uint8_t*)calloc(num_bits/8+1,sizeof(uint8_t));
    memcpy(p->bits,ptr,size*sizeof(uint8_t));
    assert(p->bits != NULL);
    p->size = num_bits/8+1;
    p->num_bits = num_bits;
    return p;
}
int bv_get(BitVector * p, int index)
{
    assert(index < p->num_bits);    
    assert(p != NULL);
    assert(p->bits != NULL);
    assert(index < p->size);
    int i = index/sizeof(uint8_t);
    int l = index % 8;
    int bit = p->bits[i] & (1 << l);
    return bit;    
}

void bv_set(BitVector * p, int index) {
    assert(index < p->num_bits);    
    assert(p != NULL);
    assert(p->bits != NULL);
    assert(index < p->size);
    int i = index/sizeof(uint8_t);
    int l = index % 8;
    int mask = (1 << l);
    p->bits[i] |= mask;
}

void bv_flip(BitVector * p, int index) {
    assert(index < p->num_bits);    
    assert(p != NULL);
    assert(p->bits != NULL);
    assert(index < p->size);
    int i = index/sizeof(uint8_t);
    int l = index % 8;
    int mask = (1 << l);
    int bit = p->bits[i] & (1 << l);
    p->bits[i] &= ~bit;    
}

#define BV_ASSERT(x) { assert((x) != NULL); assert((x)->bits != NULL); }



BitVector* bv_copy(BitVector * p) {
    assert(p != NULL);
    assert(p->bits != NULL);
    BitVector * r = bv_new(p->num_bits);
    memcpy(r->bits,p->bits,p->size);
    r->num_bits = p->num_bits;
    r->size = p->size;
    return r;
}

BitVector* bv_op(BitVector * a, BitVector * b, BV_OP op) 
{
    uint8_t x,y,r;
    BV_ASSERT(a);
    BV_ASSERT(b);
    assert(a->num_bits == b->num_bits);
    BitVector * p = bv_copy(a);
    for(int i = 0; i < p->size; i++) {
        x = a->bits[i];
        y = b->bits[i];
        if(op == BV_NOR || op == BV_NAND || op == BV_NXOR) {
            x = ~x;
            y = ~y;
        }
        switch(op) {
            case BV_AND:
            case BV_NAND: r = x & y; break;
            case BV_OR:
            case BV_NOR: r = x | y; break;
            case BV_XOR:
            case BV_NXOR: r = x ^ y; break;
        }
        p->bits[i] = r;
    }
    return p;
}
BitVector* bv_not(BitVector * a) {
    BitVector * p = bv_copy(a);
    for(int i = 0; i < p->size; i++)
        p->bits[i] = ~a->bits[i];
    return p;
}
void  bv_resize(BitVector * p, int bits) {
    BV_ASSERT(p);
    uint8_t * n = (uint8_t*)calloc(bits/8+1,sizeof(uint8_t));
    assert(n != NULL);
    memcpy(n,p->bits,bits/8);
    free(p->bits);
    p->bits = n;
    p->num_bits = bits;
    p->size = bits/8+1;
}
