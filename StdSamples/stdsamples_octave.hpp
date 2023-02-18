// Octave is experimental
// It is too slow to real-time audio

#pragma once

#include "Octopus.hpp"


Octopus::Octopus interp;

kfr3::Filters::BiquadSOS oct_butterlp(int order, DspFloatType Q=1.0)
{
    ValueList values;
    values(0) = order;
    values(1) = 1.0;
    values(2) = 's';
    values = Octopus::Functions::butter(values,2);
    values = Octopus::Functions::tf2sos(values,2);
    MatrixXd m = values(0).matrix_value();
    kfr3::Filters::BiquadSOS sos;
    int n = order/2;
    if(n == 0) n = 1;
    for(size_t i = 0; i < n; i++)
    {
        kfr3::Filters::BiquadSection c;
        c.z[0] = m(i,0);
        c.z[1] = m(i,1);
        c.z[2] = m(i,2);
        c.p[0] = m(i,3);
        c.p[1] = (1.0/Q)*m(i,4);
        c.p[2] = m(i,5);
        sos.push_back(c);
    }
    return sos;
}
kfr3::Filters::BiquadSOS oct_butterhp(int order, DspFloatType Q=1.0)
{
    ValueList values;
    values(0) = order;
    values(1) = 1.0;
    values(2) = "high";
    values(3) = 's';
    values = Octopus::Functions::butter(values,2);
    values = Octopus::Functions::tf2sos(values,2);
    MatrixXd m = values(0).matrix_value();
    kfr3::Filters::BiquadSOS sos;
    int n = order/2;
    if(n == 0) n = 1;
    for(size_t i = 0; i < n; i++)
    {
        kfr3::Filters::BiquadSection c;
        c.z[0] = m(i,0);
        c.z[1] = m(i,1);
        c.z[2] = m(i,2);
        c.p[0] = m(i,3);
        c.p[1] = (1.0/Q)*m(i,4);
        c.p[2] = m(i,5);
        sos.push_back(c);
    }
    return sos;
}    

kfr3::Filters::BiquadSOS oct_butterbp(int order, DspFloatType Q=1.0)
{
    ValueList values;
    VectorXd temp(2);
    temp(0) = 0.0;
    temp(1) = 1.0;
    values(0) = order;
    values(1) = temp;
    values(2) = "bandpass";
    values(3) = 's';
    values = Octopus::Functions::butter(values,2);
    values = Octopus::Functions::tf2sos(values,2);
    MatrixXd m = values(0).matrix_value();
    kfr3::Filters::BiquadSOS sos;
    int n = order/2;
    if(n == 0) n = 1;
    for(size_t i = 0; i < n; i++)
    {
        kfr3::Filters::BiquadSection c;
        c.z[0] = m(i,0);
        c.z[1] = m(i,1);
        c.z[2] = m(i,2);
        c.p[0] = m(i,3);
        c.p[1] = (1.0/Q)*m(i,4);
        c.p[2] = m(i,5);
        sos.push_back(c);
    }
    return sos;
} 


kfr3::Filters::BiquadSOS oct_butterbs(int order, DspFloatType Q=1.0)
{
    ValueList values;
    VectorXd temp(2);
    temp(0) = 0.0;
    temp(1) = 1.0;
    values(0) = order;
    values(1) = temp;
    values(2) = "stop";
    values(3) = 's';
    values = Octopus::Functions::butter(values,2);
    values = Octopus::Functions::tf2sos(values,2);
    MatrixXd m = values(0).matrix_value();
    kfr3::Filters::BiquadSOS sos;
    int n = order/2;
    if(n == 0) n = 1;
    for(size_t i = 0; i < n; i++)
    {
        kfr3::Filters::BiquadSection c;
        c.z[0] = m(i,0);
        c.z[1] = m(i,1);
        c.z[2] = m(i,2);
        c.p[0] = m(i,3);
        c.p[1] = (1.0/Q)*m(i,4);
        c.p[2] = m(i,5);
        sos.push_back(c);
    }
    return sos;
} 

struct OctaveButterworthFilter
{
    enum {
        LOW,
        HIGH,
        BAND,
        REJECT
    };
    int type = LOW;
    IIRFilters::BiquadFilterCascade cascade;
    DspFloatType sampleRate,cutLow,cutHigh;
    OctaveButterworth(int Type, float sr) {
        type = Type;
        sampleRate = sr;
    }
    IIRFilters::BiquadSOS Filterize(OctopusComplexMatrix & zeros, OctopusComplexMatrix & poles, DspFloatType K)
    {
        IIRFilters::BiquadSOS sos;    
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;

            std::complex<DspFloatType> z1(0,0);
            if(zeros.cols() > 0)
                z1 = zeros(order/2,0);
            std::complex<DspFloatType> p1  = poles(order/2,0);        
            
            DspFloatType y1 = abs(z1);
            DspFloatType y2 = 0;
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = K*1.0/x1;
            c.z[1] = K*y1/x1;
            c.z[2] = K*y2/x1;
            c.p[0] = 1.0;
            c.p[1] = 1.0/x1;
            c.p[2] = x2/x1;
            sos.push_back(c);        
            n++;
        }            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> z1(0,0),z2(0,0);
            if(zeros.cols() > 0) {
                z1 = zeros(order,n);
                z2 = zeros(order/2,n);
            std::complex<DspFloatType> p1  = poles(i-n,0);
            std::complex<DspFloatType> p2  = poles(order-(i-n)-1,0);
            
            DspFloatType y1 = abs(z1*z2);
            DspFloatType y2 = abs(-z1-z2);
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = K*y1/x1;
            c.z[1] = K*y2/x1;
            c.z[2] = type == LOW? 0 : K/x1;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }
    void setCutoff(int order, DspFloatType fc1, DspFloatType sr) {
        OctopusValueList l;
        OctopusComplexMatrix    poles,zeros;
        IIRFilters::BiquadSOS sos;
        DspFloatType k;
        
        switch(type) {
            case LOW:
                l(0) = order;
                l(1) = 1.0;
                l(2) = 's';
                l = Octopus::butter(l,3);
                zeros = l(0).complex_matrix_value();
                poles = l(1).complex_matrix_value();
                k     = std::abs(l(2).complex_scalar_value());
                sos = Filterize(zeros,poles,k);
                break;                
        }
        auto c = IIRFilters::AnalogBiquadCascade(sos,fc1,sr);
        setCoefficients(c);
    }
    void setCoefficients(IIRFilters::BiquadSOS & sos) {
        cascade.setCoefficients(sos);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return cascade.Tick(I,A,X,Y);
    }    
};
