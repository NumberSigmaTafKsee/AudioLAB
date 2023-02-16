#pragma once

#include "AudioProcessor.hpp"

namespace Filters
{
    struct ButterWorthLPFSection
    {
        double a[3];
        double b[2];
        double gain;
        double radius = 0.5;
        double x,y,x1,x2,y1,y2;
        double Fc,Fs;
        double k,n;
        bool   odd=false;
        
        ButterWorthLPFSection(const ButterWorthLPFSection & s, bool odd = false) {
            memcpy(a,s.a,sizeof(double)*3);
            memcpy(b,s.b,sizeof(double)*2);
            gain = s.gain;
            radius = s.radius;
            Fc = s.Fc;
            Fs = s.Fs;  
            k  = s.k;
            n  = s.n;      
            this->odd = odd;
        }
        ButterWorthLPFSection(double cutoff, double R, double k, double n, double sr)
        {            
            
            x = y = x1 = x2 = y1 = x2 = 0;
            Fc = cutoff;
            Fs = sr;
            radius = R;
            this->k = k;
            this->n = n;
            setCutoff(Fc);
        }
        ButterWorthLPFSection& operator = (const ButterWorthLPFSection & s)
        {
            memcpy(a,s.a,sizeof(double)*3);
            memcpy(b,s.b,sizeof(double)*2);
            gain = s.gain;
            radius = s.radius;
            Fc = s.Fc;
            Fs = s.Fs;
        }
        void setCutoff(float f)
        {

            double omegac = 2.0 * Fs * tan(M_PI*f/Fs);
            double zeta   = -radius*cos(M_PI*(2.0*k+n-1)/(2.0*n));                
            Fc = f;

            if(odd)
            {
                a[0] = (omegac/2.0)/(1+(omegac/2.0));
                a[1] = (omegac/2.0)/(1+(omegac/2.0));
                a[2] = 0;
                b[0] = -(1-(omegac/2.0))/(1+(omegac/2.0));
                b[1] = 0;
                b[0] /= (b[0] + b[1]);
                b[1] /= (b[0] + b[1]);
            }
            else {
                a[0] = omegac*omegac;
                a[1] = 2.0 * omegac * omegac;
                a[2] = a[1];

                double b0 = (4.0 * Fs * Fs) + (4.0 * Fs * zeta * omegac) + (omegac * omegac);
                b[0] = ((2.0 * omegac * omegac) - (8.0 * Fs * Fs)) / -b0;
                b[1] = ((4.0 * Fs * Fs) - (4.0 * Fs * zeta * omegac) + (omegac*omegac))/-b0;
                gain = 1.0/b0;
            }
        }
        void setResonance(float r)
        {
            radius = r;
        }
        double Tick(double input, double A=1, double X=0, double Y=0)
        {
            x = input * gain;
            y = a[0]* x + a[1]*x1 + a[2]*x2 + b[0]*y1 + b[1]*y2;
            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;
            return A*y;
        }
    };

    struct ButterWorthLPF
    {
        std::vector<ButterWorthLPFSection*> sections;

        ButterWorthLPF(double cutoff, double radius, int numSections, double sr)
        {
            sections.resize(numSections);
            for(size_t i = 0; i < numSections; i++)
            {
                sections[i] = new ButterWorthLPFSection(cutoff,radius,i+1,numSections*2,sr);            
            }
        }
        ~ButterWorthLPF() {
            for(size_t i = 0; i < sections.size(); i++)
            {
                delete sections[i];
            }
        }
        void setCutoff(float f)
        {
            printf("%f\n",f);
            for(size_t i = 0; i < sections.size(); i++)
                sections[i]->setCutoff(f);
        }
        void setResonance(float r) {
            if( r < 0 || r >= 1.0) return;
            for(size_t i = 0; i < sections.size(); i++)
                sections[i]->setResonance(1-r);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            double output = I;
            for(size_t i = 0; i < sections.size(); i++)
                output = sections[i]->Tick(output,A,X,Y);
            return output;
        }
    };


    struct ButterworthAudioProcessor : public MonoAudioProcessor
    {

    };
}