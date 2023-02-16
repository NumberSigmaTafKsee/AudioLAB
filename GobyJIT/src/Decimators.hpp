#pragma once

#include <cmath>
#include "SoundObject.hpp"

namespace DSP
{
  class Decimator5 : public FunctionProcessor
  {
    private:
    double R1,R2,R3,R4,R5;
    const double h0;
    const double h1;
    const double h3;
    const double h5;
    double prev=0;
    public:
    
    Decimator5(): FunctionProcessor(),h0(346/692.0f),h1(208/692.0f),h3(-44/692.0f),h5(9/692.0f)
    {
      R1=R2=R3=R4=R5=0.0f;
    }
    double Calc(const double x0,const double x1)
    {
      double h5x0=h5*x0;
      double h3x0=h3*x0;
      double h1x0=h1*x0;
      double R6=R5+h5x0;
      R5=R4+h3x0;
      R4=R3+h1x0;
      R3=R2+h1x0+h0*x1;
      R2=R1+h3x0;
      R1=h5x0;
      return R6;
    }
    double Tick(double I, double A=1, double X=0, double Y=0) {
      double out = Calc(prev,I);
      prev = I;
      return out;
    }
  };

  class Decimator7 : public FunctionProcessor
  {
    private:
    double R1,R2,R3,R4,R5,R6,R7;
    const double h0,h1,h3,h5,h7;
    double prev=0;
    public:
    Decimator7():FunctionProcessor(),h0(802/1604.0f),h1(490/1604.0f),h3(-116/1604.0f),h5(33/1604.0f),h7(-6/1604.0f)
    {
      R1=R2=R3=R4=R5=R6=R7=0.0f;
    }
    double Calc(const double x0,const double x1)
    {
      double h7x0=h7*x0;
      double h5x0=h5*x0;
      double h3x0=h3*x0;
      double h1x0=h1*x0;
      double R8=R7+h7x0;
      R7=R6+h5x0;
      R6=R5+h3x0;
      R5=R4+h1x0;
      R4=R3+h1x0+h0*x1;
      R3=R2+h3x0;
      R2=R1+h5x0;
      R1=h7x0;
      return R8;
    }
    double Tick(double I, double A=1, double X=0, double Y=0) {
      double out = Calc(prev,I);
      prev = I;
      return out;
    }
  };
  
  class Decimator9 : public FunctionProcessor
  {
    private:
    double R1,R2,R3,R4,R5,R6,R7,R8,R9;
    const double h0,h1,h3,h5,h7,h9;
    double prev=0;
    public:
    Decimator9(): FunctionProcessor(), h0(8192/16384.0f),h1(5042/16384.0f),h3(-1277/16384.0f),h5(429/16384.0f),h7(-116/16384.0f),h9(18/16384.0f)
    {
      R1=R2=R3=R4=R5=R6=R7=R8=R9=0.0f;
    }
    double Calc(const double x0,const double x1)
    {
      double h9x0=h9*x0;
      double h7x0=h7*x0;
      double h5x0=h5*x0;
      double h3x0=h3*x0;
      double h1x0=h1*x0;
      double R10=R9+h9x0;
      R9=R8+h7x0;
      R8=R7+h5x0;
      R7=R6+h3x0;
      R6=R5+h1x0;
      R5=R4+h1x0+h0*x1;
      R4=R3+h3x0;
      R3=R2+h5x0;
      R2=R1+h7x0;
      R1=h9x0;
      return R10;
    }
    double Tick(double I, double A=1, double X=0, double Y=0) {
      double out = Calc(prev,I);
      prev = I;
      return out;
    }
  };
}