%module brown 
%{
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>

#define SP_OK 0
typedef float SPFLOAT;


typedef struct {
    SPFLOAT brown;
    uint32_t sr;
} sp_brown;

typedef struct 
{
    sp_brown * brown;
}
Brown;

%}

%include "stdint.i"


typedef struct {
    SPFLOAT brown;
    uint32_t sr;
} sp_brown;


typedef struct 
{
    sp_brown * brown;
}
Brown;

/*
 * Brown
 *
 * Brownian noise algorithm based on implementation found here:
 * http://vellocet.com/dsp/noise/VRand.h
 *
 *
 */

%inline 
%{

int sp_brown_create(sp_brown **p)
{
    *p = malloc(sizeof(sp_brown));
    return SP_OK;
}

int sp_brown_destroy(sp_brown **p)
{
    free(*p);
    return SP_OK;
}

int sp_brown_init(uint32_t sr, sp_brown *p)
{
    p->brown = 0.0;
    p->sr = sr;
    return SP_OK;
}

int sp_brown_compute(sp_brown *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT r;
    while (1) {
        r = (rand() % RAND_MAX) / (SPFLOAT)(RAND_MAX);
        r = ((r * 2) - 1) * 0.5;
        p->brown += r;
        if (p->brown < -8.0f || p->brown > 8.0f) {
            p->brown -= r;
        } else {
            break;
        }
    }
    *out = p->brown * 0.0625;
    return SP_OK;
}

Brown* Brown_New(uint32_t sr)
{
    Brown* p = (Brown*)calloc(1,sizeof(Brown));
    assert(p != NULL);
    srand(time(NULL));
    sp_brown_create(&p->brown);
    assert(p->brown != NULL);
    sp_brown_init(sr, p->brown);
    return p;
}

float Brown_Tick(Brown * p)
{
    float in = 0.0;
    float out = 0.0;
    sp_brown_compute(p->brown,&in,&out);
    return out;
}

%}