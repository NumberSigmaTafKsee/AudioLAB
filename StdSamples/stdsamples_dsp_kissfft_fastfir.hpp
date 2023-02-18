#pragma once

#include "_kiss_fft_guts.h"


/*
 Some definitions that allow real or complex filtering
*/
#ifdef REAL_FASTFIR
#define MIN_FFT_LEN 2048
#include "kiss_fftr.h"
typedef kiss_fft_scalar kffsamp_t;
typedef kiss_fftr_cfg kfcfg_t;
#define FFT_ALLOC kiss_fftr_alloc
#define FFTFWD kiss_fftr
#define FFTINV kiss_fftri
#else
#define MIN_FFT_LEN 1024
typedef kiss_fft_cpx kffsamp_t;
typedef kiss_fft_cfg kfcfg_t;
#define FFT_ALLOC kiss_fft_alloc
#define FFTFWD kiss_fft
#define FFTINV kiss_fft
#endif

typedef struct kiss_fastfir_state *kiss_fastfir_cfg;



kiss_fastfir_cfg kiss_fastfir_alloc(const kffsamp_t * imp_resp,size_t n_imp_resp,
        size_t * nfft,void * mem,size_t*lenmem);

/* see do_file_filter for usage */
size_t kiss_fastfir( kiss_fastfir_cfg cfg, kffsamp_t * inbuf, kffsamp_t * outbuf, size_t n, size_t *offset);



static int verbose=0;


struct kiss_fastfir_state{
    size_t nfft;
    size_t ngood;
    kfcfg_t fftcfg;
    kfcfg_t ifftcfg;
    kiss_fft_cpx * fir_freq_resp;
    kiss_fft_cpx * freqbuf;
    size_t n_freq_bins;
    kffsamp_t * tmpbuf;
};


kiss_fastfir_cfg kiss_fastfir_alloc(
        const kffsamp_t * imp_resp,size_t n_imp_resp,
        size_t *pnfft, /* if <= 0, an appropriate size will be chosen */
        void * mem,size_t*lenmem)
{
    kiss_fastfir_cfg st = NULL;
    size_t len_fftcfg,len_ifftcfg;
    size_t memneeded = sizeof(struct kiss_fastfir_state);
    char * ptr;
    size_t i;
    size_t nfft=0;
    float scale;
    int n_freq_bins;
    if (pnfft)
        nfft=*pnfft;

    if (nfft<=0) {
        /* determine fft size as next power of two at least 2x 
         the impulse response length*/
        i=n_imp_resp-1;
        nfft=2;
        do{
             nfft<<=1;
        }while (i>>=1);
#ifdef MIN_FFT_LEN
        if ( nfft < MIN_FFT_LEN )
            nfft=MIN_FFT_LEN;
#endif        
    }
    if (pnfft)
        *pnfft = nfft;

#ifdef REAL_FASTFIR
    n_freq_bins = nfft/2 + 1;
#else
    n_freq_bins = nfft;
#endif
    /*fftcfg*/
    FFT_ALLOC (nfft, 0, NULL, &len_fftcfg);
    memneeded += len_fftcfg;  
    /*ifftcfg*/
    FFT_ALLOC (nfft, 1, NULL, &len_ifftcfg);
    memneeded += len_ifftcfg;  
    /* tmpbuf */
    memneeded += sizeof(kffsamp_t) * nfft;
    /* fir_freq_resp */
    memneeded += sizeof(kiss_fft_cpx) * n_freq_bins;
    /* freqbuf */
    memneeded += sizeof(kiss_fft_cpx) * n_freq_bins;
    
    if (lenmem == NULL) {
        st = (kiss_fastfir_cfg) malloc (memneeded);
    } else {
        if (*lenmem >= memneeded)
            st = (kiss_fastfir_cfg) mem;
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;

    st->nfft = nfft;
    st->ngood = nfft - n_imp_resp + 1;
    st->n_freq_bins = n_freq_bins;
    ptr=(char*)(st+1);

    st->fftcfg = (kfcfg_t)ptr;
    ptr += len_fftcfg;

    st->ifftcfg = (kfcfg_t)ptr;
    ptr += len_ifftcfg;

    st->tmpbuf = (kffsamp_t*)ptr;
    ptr += sizeof(kffsamp_t) * nfft;

    st->freqbuf = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;
    
    st->fir_freq_resp = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;

    FFT_ALLOC (nfft,0,st->fftcfg , &len_fftcfg);
    FFT_ALLOC (nfft,1,st->ifftcfg , &len_ifftcfg);

    memset(st->tmpbuf,0,sizeof(kffsamp_t)*nfft);
    /*zero pad in the middle to left-rotate the impulse response 
      This puts the scrap samples at the end of the inverse fft'd buffer */
    st->tmpbuf[0] = imp_resp[ n_imp_resp - 1 ];
    for (i=0;i<n_imp_resp - 1; ++i) {
        st->tmpbuf[ nfft - n_imp_resp + 1 + i ] = imp_resp[ i ];
    }

    FFTFWD(st->fftcfg,st->tmpbuf,st->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    scale = 1.0 / st->nfft;

    for ( i=0; i < st->n_freq_bins; ++i ) {
#ifdef USE_SIMD
        st->fir_freq_resp[i].r *= _mm_set1_ps(scale);
        st->fir_freq_resp[i].i *= _mm_set1_ps(scale);
#else
        st->fir_freq_resp[i].r *= scale;
        st->fir_freq_resp[i].i *= scale;
#endif
    }
    return st;
}

static void fastconv1buf(const kiss_fastfir_cfg st,const kffsamp_t * in,kffsamp_t * out)
{
    size_t i;
    /* multiply the frequency response of the input signal by
     that of the fir filter*/
    FFTFWD( st->fftcfg, in , st->freqbuf );
    for ( i=0; i<st->n_freq_bins; ++i ) {
        kiss_fft_cpx tmpsamp; 
        C_MUL(tmpsamp,st->freqbuf[i],st->fir_freq_resp[i]);
        st->freqbuf[i] = tmpsamp;
    }

    /* perform the inverse fft*/
    FFTINV(st->ifftcfg,st->freqbuf,out);
}

/* n : the size of inbuf and outbuf in samples
   return value: the number of samples completely processed
   n-retval samples should be copied to the front of the next input buffer */
static size_t kff_nocopy(
        kiss_fastfir_cfg st,
        const kffsamp_t * inbuf, 
        kffsamp_t * outbuf,
        size_t n)
{
    size_t norig=n;
    while (n >= st->nfft ) {
        fastconv1buf(st,inbuf,outbuf);
        inbuf += st->ngood;
        outbuf += st->ngood;
        n -= st->ngood;
    }
    return norig - n;
}

static
size_t kff_flush(kiss_fastfir_cfg st,const kffsamp_t * inbuf,kffsamp_t * outbuf,size_t n)
{
    size_t zpad=0,ntmp;

    ntmp = kff_nocopy(st,inbuf,outbuf,n);
    n -= ntmp;
    inbuf += ntmp;
    outbuf += ntmp;

    zpad = st->nfft - n;
    memset(st->tmpbuf,0,sizeof(kffsamp_t)*st->nfft );
    memcpy(st->tmpbuf,inbuf,sizeof(kffsamp_t)*n );
    
    fastconv1buf(st,st->tmpbuf,st->tmpbuf);
    
    memcpy(outbuf,st->tmpbuf,sizeof(kffsamp_t)*( st->ngood - zpad ));
    return ntmp + st->ngood - zpad;
}

size_t kiss_fastfir(
        kiss_fastfir_cfg vst,
        kffsamp_t * inbuf,
        kffsamp_t * outbuf,
        size_t n_new,
        size_t *offset)
{
    size_t ntot = n_new + *offset;
    if (n_new==0) {
        return kff_flush(vst,inbuf,outbuf,ntot);
    }else{
        size_t nwritten = kff_nocopy(vst,inbuf,outbuf,ntot);
        *offset = ntot - nwritten;
        /*save the unused or underused samples at the front of the input buffer */
        memcpy( inbuf , inbuf+nwritten , *offset * sizeof(kffsamp_t) );
        return nwritten;
    }
}