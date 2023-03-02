%module bitcrush
%{
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

typedef float SPFLOAT;
#define SP_OK 0

typedef struct {
    SPFLOAT bitdepth;
    SPFLOAT srate;
    uint32_t sr;
    SPFLOAT incr;
    SPFLOAT index;
    int32_t sample_index;
    SPFLOAT value;
} sp_bitcrush;

typedef struct bitcrush_t
{
    sp_bitcrush * crush;
}
BitCrush;

%}

%include "stdint.i"

typedef struct {
    SPFLOAT bitdepth;
    SPFLOAT srate;
    uint32_t sr;
    SPFLOAT incr;
    SPFLOAT index;
    int32_t sample_index;
    SPFLOAT value;
} sp_bitcrush;

typedef struct bitcrush_t
{
    sp_bitcrush * crush;
}
BitCrush;

%inline 
%{
int sp_bitcrush_create(sp_bitcrush **p)
{
    *p = malloc(sizeof(sp_bitcrush));
    return SP_OK;
}

int sp_bitcrush_destroy(sp_bitcrush **p)
{
    free(*p);
    return SP_OK;
}

int sp_bitcrush_init(uint32_t sr, sp_bitcrush *p)
{
    p->bitdepth = 8;
    p->srate = 10000;
    p->sr = 44100;
    p->incr = 1000;
    p->sample_index = 0;
    p->index = 0.0;
    p->value = 0.0;
    return SP_OK;
}

int sp_bitcrush_compute(sp_bitcrush *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT bits = pow(2, floor(p->bitdepth));
    SPFLOAT foldamt = p->sr / p->srate;
    SPFLOAT sig;
    *out = *in * 65536.0;
    *out += 32768;
    *out *= (bits / 65536.0);
    *out = floor(*out);
    *out = *out * (65536.0 / bits) - 32768;
    sig = *out;
    p->incr = foldamt;

    /* apply downsampling */
    if (p->index < (SPFLOAT)p->sample_index) {
        p->index += p->incr;
        p->value = sig;
    }

    *out = p->value;

    p->sample_index++;

    *out /= 65536.0;
    return SP_OK;
}

BitCrush* BitCrush_New(uint32_t sr)
{
    BitCrush * p = (BitCrush*)calloc(1,sizeof(BitCrush));
    assert( p != NULL);
    sp_bitcrush_create(&p->crush);
    assert( p->crush != NULL );
    sp_bitcrush_init(sr, p->crush);
    return p;
}
float BitCrush_Tick(BitCrush * c, float in)
{
    float out = 0.0f;
    sp_bitcrush_compute(c->crush, &in, &out);
    return out;
}
%}