#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

typedef float SPFLOAT;

const int SP_OK=0;

enum { CLEAR, ATTACK, DECAY, SUSTAIN, RELEASE };

typedef struct {
    SPFLOAT atk;
    SPFLOAT dec;
    SPFLOAT sus;
    SPFLOAT rel;
    uint32_t timer;
    uint32_t atk_time;
    SPFLOAT a;
    SPFLOAT b;
    SPFLOAT y;
    SPFLOAT x;
    SPFLOAT prev;
    uint32_t sr;
    int mode;
} sp_adsr;

typedef struct tauadsr_t
{
    sp_adsr *adsr;
}
TauADSR;

TauADSR* TauADSR_New(uint32_t sr);
void     TauADSR_Free(TauADSR * p);
SPFLOAT  TauADSR_Tick(TauADSR * p,SPFLOAT gate);


static int sp_adsr_create(sp_adsr **p)
{
    *p = malloc(sizeof(sp_adsr));
    return SP_OK;
}

static int sp_adsr_destroy(sp_adsr **p)
{
    free(*p);
    return SP_OK;
}

static int sp_adsr_init(uint32_t sr, sp_adsr *p)
{
    p->atk = 0.1;
    p->dec = 0.1;
    p->sus = 0.5;
    p->rel = 0.3;
    p->timer = 0;
    p->a = 0;
    p->b = 0;
    p->y = 0;
    p->x = 0;
    p->prev = 0;
    p->atk_time = p->atk * sr;
    p->mode = CLEAR;
    p->sr = sr;
    return SP_OK;
}

static SPFLOAT tau2pole(sp_adsr *p, SPFLOAT tau)
{
    return exp(-1.0 / (tau * p->sr));
}

static SPFLOAT adsr_filter(sp_adsr *p)
{
    p->y = p->b * p->x  + p->a * p->y;
    return p->y;
}

static int sp_adsr_compute(sp_adsr *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT pole;
    if (p->prev < *in && p->mode != DECAY) {
        p->mode = ATTACK;
        p->timer = 0;
        /* quick fix: uncomment if broken */
        /* pole = tau2pole(sp, p, p->atk * 0.75); */
        /* p->atk_time = p->atk * sp->sr * 1.5; */
        pole = tau2pole( p, p->atk * 0.6);
        p->atk_time = p->atk * p->sr;
        p->a = pole;
        p->b = 1 - pole;
    } else if (p->prev > *in) {
        p->mode = RELEASE;
        pole = tau2pole(p, p->rel);
        p->a = pole;
        p->b = 1 - pole;
    }

    p->x = *in;
    p->prev = *in;

    switch (p->mode) {
        case CLEAR:
            *out = 0;
            break;
        case ATTACK:
            p->timer++;
            *out = adsr_filter(p);
            /* quick fix: uncomment if broken */
            /* if(p->timer > p->atk_time) { */
            if (*out > 0.99) {
                p->mode = DECAY;
                pole = tau2pole(p, p->dec);
                p->a = pole;
                p->b = 1 - pole;
            }
            break;
        case DECAY:
        case RELEASE:
            p->x *= p->sus;
            *out = adsr_filter(p);
        default:
            break;
    }

    return SP_OK;
}


void     TauADSR_Free(TauADSR * p)
{
    if(p->adsr) sp_adsr_destroy(&p->adsr);
    free(p);
}

TauADSR* TauADSR_New(uint32_t sr)
{
    TauADSR * p = (TauADSR*)calloc(1,sizeof(TauADSR));
    assert(p != NULL);
    sp_adsr_create(&p->adsr);
    assert(p->adsr != NULL);
    sp_adsr_init(sr, p->adsr);
    return p;
}

SPFLOAT TauADSR_Tick(TauADSR * p,SPFLOAT gate)
{
    SPFLOAT out = 0.0;
    sp_adsr_compute(p->adsr,&gate,&out);
    return out;
}