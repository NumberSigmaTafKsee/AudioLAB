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
