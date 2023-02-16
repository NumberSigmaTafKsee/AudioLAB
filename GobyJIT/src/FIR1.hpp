#include <Eigen/Dense>
#include <Eigen/Core>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

#define BANDPASS 0
#define LOW 1
#define HIGH 2
#define DC0 3

int FIRTypeCheck(ArrayXd Wn){
    int Type;
    if (Wn.size() == 1)
        Type = LOW;
    else if (Wn.size() == 2)
        Type = HIGH;
    else
        Type = DC0;
    return Type;
}

ArrayXd sinc(ArrayXd t){
    double pi = 3.1415926535;
    ArrayXd res(t.size());
    for (int ii = 0; ii < t.size(); ii++){
        res(ii) = t(ii) == 0? 1:sin(pi * t(ii)) / pi / t(ii);
    }
    return res;
}

// Linear-phase FIR filter design using least-squares error minimization
ArrayXd firls(int N, ArrayXd F, ArrayXd M){
    double pi = 3.1415926535;
    ArrayXd W = ArrayXd::Ones(2);
    N++; // 滤波器长度
    F = F / 2; W = sqrt(W);
    ArrayXd dF(F.size() - 1);
    bool alldF2e0 = true;
    for (int ii = 0; ii < dF.size(); ii++){
        dF(ii) = F(ii + 1) - F(ii);
        if (ii % 1 == 0 && dF(ii) != 0)
            alldF2e0 = false;
    }

    int fullband;
    int constant_weights = 1;
    if (dF.size() > 0 && alldF2e0)
        fullband = 1;
    else
        fullband = 0;
    for (int ii = 0; ii < W.size(); ii++)
        if (W(0) != W(ii)){
            constant_weights = 0;
            break;
        }

    int L = (N - 1) / 2;
    int Nodd = N % 2;

    ArrayXd m = ArrayXd::LinSpaced(L + 1,0,L);
    if (!Nodd){
        m = m + 5;
    }

    ArrayXd k = m;
    bool need_matrix = !fullband || !constant_weights;
    if (need_matrix){
        // TODO
    }
    double b0;
    if (Nodd){
       k = m.segment(1, k.size() - 1);
       b0 = 0;
    }
    ArrayXd b = ArrayXd::Zero(k.size());

    for (int s = 0; s < F.size(); s+=2){
        double _m_ = (M(s+1)-M(s))/(F(s+1)-F(s));
        double b1 = M(s) - _m_ * F(s);
        if (Nodd){
            b0 += (b1 * (F(s+1) - F(s)) + _m_ / 2.0 * (F(s+1)*F(s+1)-F(s)*F(s))) * abs(pow(W((s+1)/2), 2));
        }

        ArrayXd temp = 1.0 / pow(k, 2);

        b = b + (_m_/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s))) * temp) * abs(pow(W((s+1)/2), 2));
        b = b + (F(s+1)*(_m_*F(s+1)+b1)*sinc(2*k*F(s+1)) - F(s)*(_m_*F(s)+b1)*sinc(2*k*F(s))) * abs(pow(W((s+1)/2), 2));
    }

    ArrayXd _b_;
    if (Nodd){
        _b_ = ArrayXd(1 + b.size());
        _b_(0) = b0;
        _b_.segment(1, b.size()) = b;
        b = _b_;
    }

    ArrayXd a = (W(0) * W(0)) * 4 * b;
    if (Nodd)
        a(0) = a(0) / 2;
    ArrayXd h;
    if (Nodd){
        h = ArrayXd(2 * L + 1);
        for (int ii = 0; ii < 2 * L + 1; ii++){
            if (ii < L)
                h(ii) = a(L - ii) / 2;
            else if (ii == L)
                h(ii) = a(0);
            else{
                h(ii) = a(ii - L - 1) / 2;
            }
        }
    }else{
        h = ArrayXd(2 * L + 2);
        for (int ii = 0; ii < 2 * L + 1; ii++){
            if (ii < L + 1){
                h(ii) = a(a.size() - ii - 1) / 2;
            }else{
                h(ii) = a(ii - L - 1) / 2;
            }
        }
    }
    return h;
}

ArrayXd hammingWin(int N){
    ArrayXd ret(N);
    double pi = 3.1415926535;
    for(int n = 0; n < N; n++)
    {
        ret(n) = 0.54 - 0.46 * cos (2 * pi *  ( double )n / ( N - 1 ));
    }
    return ret;
}

ArrayXd fir1(int N, ArrayXd Wn){
    ArrayXd _res_;
    for (int ii = 0; ii < Wn.size(); ii++) if (Wn(ii) >= 1.0) return _res_;

    int Type = FIRTypeCheck(Wn);
    // 查看这个是带通还是高通、带阻 xxxx 滤波器 ¯\_(ツ)_/¯

    // 找到对应的频率 ¯\_(ツ)_/¯
    int nbands = Wn.size() + 1;
    if (nbands > 2 && Type == BANDPASS){
        Type = DC0;
    }

    ArrayXXd ff(2, Wn.size() + 1);
    ff(0, 0) = 0;
    ff.row(0).segment(1, Wn.size()) = Wn;
    ff.row(1).segment(0, Wn.size()) = Wn;
    ff(1, Wn.size()) = 1;

    ArrayXd fff(2 * (Wn.size() + 1));
    for (int ii = 0; ii < ff.cols(); ii++)
        fff.segment(ii * ff.rows(), ff.rows()) = ff.col(ii);

    // 找到想要的幅值 ¯\_(ツ)_/¯
    int First_Band = Type != DC0 && Type != HIGH;
    ArrayXd bands = ArrayXd::LinSpaced(nbands, 0, nbands - 1);
    ArrayXi mags(bands.size());
    for (int ii = 0; ii < mags.size(); ii++)
        mags(ii) = (First_Band + int(bands(ii))) % 2;
    ArrayXd aa(2 * bands.size());
    for (int ii = 0; ii < mags.size(); ii++){
        int offset = ii * 2;
        aa(offset) = mags(ii);
        aa(offset + 1) = mags(ii);
    }

    int L = N + 1; // FIR滤波器窗口长度
    ArrayXd Wind = hammingWin(L); // 构建汉明窗口
    ArrayXd hh = firls(L-1, fff, aa);
    return Wind * hh;
}