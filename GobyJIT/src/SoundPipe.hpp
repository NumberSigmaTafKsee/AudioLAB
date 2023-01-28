#pragma once

extern "C"
{
#include "SoundPipe/soundpipe.h"
/*
#include "SoundPipe/adsr.h"
#include "SoundPipe/allpass.h"
#include "SoundPipe/atone.h"
#include "SoundPipe/autowah.h"
#include "SoundPipe/bal.h"
#include "SoundPipe/bar.h"
#include "SoundPipe/base.h"
#include "SoundPipe/biquad.h"
#include "SoundPipe/biscale.h"
#include "SoundPipe/bitcrush.h"
#include "SoundPipe/blsaw.h"
#include "SoundPipe/blsquare.h"
#include "SoundPipe/bltriangle.h"
#include "SoundPipe/brown.h"
#include "SoundPipe/butbp.h"
#include "SoundPipe/butbr.h"
#include "SoundPipe/buthp.h"
#include "SoundPipe/butlp.h"
#include "SoundPipe/clip.h"
#include "SoundPipe/clock.h"
#include "SoundPipe/comb.h"
#include "SoundPipe/compressor.h"
#include "SoundPipe/conv.h"
#include "SoundPipe/count.h"
#include "SoundPipe/crossfade.h"
#include "SoundPipe/dcblock.h"
#include "SoundPipe/delay.h"
#include "SoundPipe/diode.h"
#include "SoundPipe/diskin.h"
#include "SoundPipe/dist.h"
#include "SoundPipe/dmetro.h"
#include "SoundPipe/drip.h"
#include "SoundPipe/dtrig.h"
#include "SoundPipe/dust.h"
#include "SoundPipe/eqfil.h"
#include "SoundPipe/expon.h"
#include "SoundPipe/fftwrapper.h"
#include "SoundPipe/fof.h"
#include "SoundPipe/fofilt.h"
#include "SoundPipe/fog.h"
#include "SoundPipe/fold.h"
#include "SoundPipe/foo.h"
#include "SoundPipe/fosc.h"
#include "SoundPipe/ftbl.h"
#include "SoundPipe/gbuzz.h"
#include "SoundPipe/hilbert.h"
#include "SoundPipe/in.h"
#include "SoundPipe/incr.h"
#include "SoundPipe/jack.h"
#include "SoundPipe/jcrev.h"
#include "SoundPipe/jitter.h"
#include "SoundPipe/line.h"
#include "SoundPipe/lpc.h"
#include "SoundPipe/lpf18.h"
#include "SoundPipe/maygate.h"
#include "SoundPipe/metro.h"
#include "SoundPipe/mincer.h"
#include "SoundPipe/mode.h"
#include "SoundPipe/moogladder.h"
#include "SoundPipe/noise.h"
#include "SoundPipe/nsmp.h"
#include "SoundPipe/osc.h"
#include "SoundPipe/oscmorph.h"
#include "SoundPipe/padsynth.h"
#include "SoundPipe/pan2.h"
#include "SoundPipe/panst.h"
#include "SoundPipe/pareq.h"
#include "SoundPipe/paulstretch.h"
#include "SoundPipe/pdhalf.h"
#include "SoundPipe/peaklim.h"
#include "SoundPipe/phaser.h"
#include "SoundPipe/phasor.h"
#include "SoundPipe/pinknoise.h"
#include "SoundPipe/pitchamdf.h"
#include "SoundPipe/pluck.h"
#include "SoundPipe/port.h"
#include "SoundPipe/posc3.h"
#include "SoundPipe/progress.h"
#include "SoundPipe/prop.h"
#include "SoundPipe/pshift.h"
#include "SoundPipe/ptrack.h"
#include "SoundPipe/randh.h"
#include "SoundPipe/randi.h"
#include "SoundPipe/randmt.h"
#include "SoundPipe/random.h"
#include "SoundPipe/reson.h"
#include "SoundPipe/reverse.h"
#include "SoundPipe/revsc.h"
#include "SoundPipe/rms.h"
#include "SoundPipe/rpi.h"
#include "SoundPipe/rpt.h"
#include "SoundPipe/rspline.h"
#include "SoundPipe/samphold.h"
#include "SoundPipe/saturator.h"
#include "SoundPipe/scale.h"
#include "SoundPipe/scrambler.h"
#include "SoundPipe/sdelay.h"
#include "SoundPipe/slice.h"
#include "SoundPipe/smoothdelay.h"
#include "SoundPipe/spa.h"
#include "SoundPipe/sparec.h"
#include "SoundPipe/sp_base.h"
#include "SoundPipe/streson.h"
#include "SoundPipe/switch.h"
#include "SoundPipe/tabread.h"
#include "SoundPipe/tadsr.h"
#include "SoundPipe/talkbox.h"
#include "SoundPipe/tblrec.h"
#include "SoundPipe/tbvcf.h"
#include "SoundPipe/tdiv.h"
#include "SoundPipe/tenv.h"
#include "SoundPipe/tenv2.h"
#include "SoundPipe/tenvx.h"
#include "SoundPipe/tevent.h"
#include "SoundPipe/tgate.h"
#include "SoundPipe/thresh.h"
#include "SoundPipe/timer.h"
#include "SoundPipe/tin.h"
#include "SoundPipe/tone.h"
#include "SoundPipe/trand.h"
#include "SoundPipe/tseg.h"
#include "SoundPipe/tseq.h"
#include "SoundPipe/vdelay.h"
#include "SoundPipe/voc.h"
#include "SoundPipe/vocoder.h"
#include "SoundPipe/waveset.h"
#include "SoundPipe/wavin.h"
#include "SoundPipe/wavout.h"
#include "SoundPipe/wpkorg35.h"
#include "SoundPipe/zitarev.h"
*/
}

namespace SoundPipe
{
    inline sp_data* create_soundpipe(int sample_rate)
    {
        sp_data *sp;
        sp_create(&sp);
        sp->sr = sample_rate;
        return sp;
    }

    inline void destroy_soundpipe(sp_data * sp) 
    {
        sp_destroy(&sp);
    }

    struct SoundPipe
    {
        sp_data * sp;
        SoundPipe(sp_data * d) : sp(d) {
            
        }
        virtual float Tick(float I, float A=1, float X=0, float Y=0) = 0;
    };

    struct ADSR : public SoundPipe
    {
        
        sp_adsr * adsr;
        float     gate;
        enum { CLEAR, ATTACK, DECAY, SUSTAIN, RELEASE };
        ADSR(sp_data * data, float a, float d, float s, float r) : SoundPipe(data) 
        {
            sp_adsr_create(&adsr);            
            sp_adsr_init(sp, adsr);    
            adsr->atk = a;
            adsr->dec = d;
            adsr->sus = s;
            adsr->rel = r;
            gate = 0;
        }
        ~ADSR() {
            if(adsr) sp_adsr_destroy(&adsr);
        }
        void noteOn() { gate = 1.0f; }
        void noteOff() { gate = 0.0f; }
        void setAttack(float a) { adsr->atk = a; }
        void setDecay(float d) { adsr->dec = d; }
        void setSustain(float s) { adsr->sus = s; }
        void setRelease(float r) { adsr->rel = r; }

        float Tick(float I=0, float A = 1, float X = -1, float Y = 1) {            
            float out= 0;
            sp_adsr_compute(sp,adsr,&gate,&out);            
            return out;
        }
    };

    struct AllPass : public SoundPipe
    {
        sp_allpass * obj;

        AllPass(sp_data * data, float looptime) : SoundPipe(data) {
            sp_allpass_create(&obj);            
            sp_allpass_init(sp,obj,looptime);
        }
        ~AllPass() {
            if(obj) sp_allpass_destroy(&obj);
        }
        float Tick(float I, float A = 1, float X = -1, float Y = 1) {
            float out = 0.0f;
            float in  = I;
            if(in < X) in = X;
            if(in > Y) in = Y;
            sp_allpass_compute(sp,obj,&in,&out);
            return A*out;
        }
    };
}        