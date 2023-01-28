#pragma once

class TSClipper
{
public:
    
    TSClipper();
    
    ~TSClipper();
    
    void prepare(float newFs);

    float processSample(float Vi);
    
    void setKnob(float newDrive);
    
private:
    
    const float eta = 1.f; // Change for non-ideal diode
    const float Is = 1e-15;
    const float Vt = 26e-3;
    
    float Fs = 48000.f;
    float Ts = 1.f/Fs;
    
    float C1 = 47e-9;
    float R1 = Ts/(2.f*C1);
    float C2 = 51e-12;
    float R2 = Ts/(2.*C2);
    float drivePot = 1.f;
    float P1 = drivePot*500e3;
    float R3 = 51000.f + P1;
    float R4 = 4700.f;

    // Combined Resistances
    float G1 = (1.f + R4/R1);
    float G4 = (1.f + R1/R4);
    
    // States
    float x1 = 0.f;
    float x2 = 0.f;
    float Vd = 0.f;
    
    float thr = 0.00000000001f;
    
    void updateCoefficients();
};

TSClipper::TSClipper(){}

TSClipper::~TSClipper(){}

void TSClipper::prepare(float newFs)
{
    if (Fs != newFs)
    {
        Fs = newFs;
        updateCoefficients();
    }
}

float TSClipper::processSample(float Vi)
{
    float p = -Vi/(G4*R4) + R1/(G4*R4)*x1 - x2;
    int iter = 1;
    float b = 1.f;
    float fd = p + Vd/R2 + Vd/R3 + 2.f*Is * sinh(Vd/(eta*Vt));
    
    while (iter < 50 && abs(fd) > thr)
    {
        float den = 2.f*Is/(eta*Vt) * cosh(Vd/(eta*Vt)) + 1.f/R2 + 1.f/R3;
        float Vnew = Vd - b*fd/den;
        float fn = p + Vnew/R2 + Vnew/R3 + 2.f*Is * sinh(Vnew/(eta*Vt));
        
        if (abs(fn) < abs(fd))
        {
            Vd = Vnew;
            b = 1.f;
        }
        else
            b *= 0.5f;
    
        fd = p + Vd/R2 + Vd/R3 + 2.f*Is * sinh(Vd/(eta*Vt));
        iter++;
    }
    
    float Vo = Vd + Vi;
    x1 = (2.f/R1)*(Vi/G1 + x1*R4/G1) - x1;
    x2 = (2.f/R2)*(Vd) - x2;
    
    return Vo;
}

void TSClipper::setKnob(float newDrive)
{
    if (drivePot != newDrive)
    {
        drivePot = newDrive;
        updateCoefficients();
    }
}

void TSClipper::updateCoefficients()
{
    Ts = 1.f/Fs;
    
    R1 = Ts/(2.f*C1);
    R2 = Ts/(2.f*C2);
    P1 = drivePot * 500e3;
    R3 = 51000.f + P1;
    
    // Combined Resistances
    G1 = (1.f + R4/R1);
    G4 = (1.f + R1/R4);
}