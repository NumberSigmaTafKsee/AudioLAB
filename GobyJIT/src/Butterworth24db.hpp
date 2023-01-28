#pragma once
#include <cmath>
#include "Undenormal.hpp"
#include "SoundObject.hpp"
#undef saturate

namespace Filters
{
    class CFilterButterworth24db : public FilterProcessor
    {
    public:
        CFilterButterworth24db(void);
        ~CFilterButterworth24db(void);
        void SetSampleRate(float fs);
        void Set(float cutoff, float q);
        float Run(float input);

        double Tick(double I, double A=1, double X=0, double Y=0) {
            return Run(I);
        }

    private:
        float t0, t1, t2, t3;
        float coef0, coef1, coef2, coef3;
        float history1, history2, history3, history4;
        float gain;
        float min_cutoff, max_cutoff;
    };


    #define BUDDA_Q_SCALE 6.f

    inline CFilterButterworth24db::CFilterButterworth24db(void) : FilterProcessor()
    {
        this->history1 = 0.f;
        this->history2 = 0.f;
        this->history3 = 0.f;
        this->history4 = 0.f;

        this->SetSampleRate(44100.f);
        this->Set(22050.f, 0.0);
    }

    inline CFilterButterworth24db::~CFilterButterworth24db(void)
    {
    }

    inline void CFilterButterworth24db::SetSampleRate(float fs)
    {
        float pi = 4.f * atanf(1.f);

        this->t0 = 4.f * fs * fs;
        this->t1 = 8.f * fs * fs;
        this->t2 = 2.f * fs;
        this->t3 = pi / fs;

        this->min_cutoff = fs * 0.01f;
        this->max_cutoff = fs * 0.45f;
    }

    inline void CFilterButterworth24db::Set(float cutoff, float q)
    {
        if (cutoff < this->min_cutoff)
                cutoff = this->min_cutoff;
        else if(cutoff > this->max_cutoff)
                cutoff = this->max_cutoff;


        float wp = this->t2 * tanf(this->t3 * cutoff);
        float bd, bd_tmp, b1, b2;

        b1 = (0.765367f / q) / wp;
        b2 = 1.f / (wp * wp);

        bd_tmp = this->t0 * b2 + 1.f;

        bd = 1.f / (bd_tmp + this->t2 * b1);

        this->gain = bd * 0.5f;

        this->coef2 = (2.f - this->t1 * b2);

        this->coef0 = this->coef2 * bd;
        this->coef1 = (bd_tmp - this->t2 * b1) * bd;

        b1 = (1.847759f / q) / wp;

        bd = 1.f / (bd_tmp + this->t2 * b1);

        this->gain *= bd;
        this->coef2 *= bd;
        this->coef3 = (bd_tmp - this->t2 * b1) * bd;
    }

    inline float CFilterButterworth24db::Run(float input)
    {
        Undenormal denormal;
        float output = input * this->gain;
        float new_hist;


        output -= this->history1 * this->coef0;
        new_hist = output - this->history2 * this->coef1;

        output = new_hist + this->history1 * 2.f;
        output += this->history2;

        this->history2 = this->history1;
        this->history1 = new_hist;

        output -= this->history3 * this->coef2;
        new_hist = output - this->history4 * this->coef3;

        output = new_hist + this->history3 * 2.f;
        output += this->history4;

        this->history4 = this->history3;
        this->history3 = new_hist;

        return output;
    }   
}
