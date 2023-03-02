#ifndef __BITVECTOR_H
#define __BITVECTOR_H

typedef struct CPUMemory;

typedef enum {
    BV_AND,
    BV_OR,
    BV_XOR,
    BV_NAND,
    BV_NOR,
    BV_NXOR,
}
BV_OP;

typedef struct {
    uint8_t * bits;
    int       size;
    int       num_bits;
}
BitVector;



#ifdef __cplusplus
extern "C" { 
#endif 

BitVector * bv_new(int num_bits);
BitVector * bv_ptr(void * ptr, int size);
int         bv_get(BitVector * p, int index);
void        bv_set(BitVector * p, int index);
void        bv_flip(BitVector * p, int index);
BitVector* bv_encode_memory(CPUMemory * p);
CPUMemory* bv_decode_memory(BitVector * bits);
BitVector* bv_copy(BitVector * p);
BitVector* bv_not(BitVector * a);
BitVector* bv_op(BitVector * a, BitVector * b, BV_OP op);
void        bv_resize(BitVector * p, int bits);

#ifdef __cplusplus
}
#endif

#endif