#ifndef __CTYPE_H
#define __CTYPE_H

#include <cstdint>

union UCType {
    int8_t i8;
    uint8_t u8;
    int16_t i16;
    uint16_t u16;
    int32_t i32;
    uint32_t u32;
    int64_t i64;
    uint64_t u64;
    float f32;
    double f64;
    const char * str;
};

enum CType_t {
    CTYPE_I8,
    CTYPE_U8,
    CTYPE_I16,
    CTYPE_U16,
    CTYPE_I32,
    CTYPE_U32,
    CTYPE_I64,
    CTYPE_U64,
    CTYPE_F32,
    CTYPE_F64,
    CTYPE_STR,    
};

struct CType {
    UCType  t;
    CType_t type;
};

#endif
