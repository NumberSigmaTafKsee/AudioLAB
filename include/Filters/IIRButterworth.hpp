// this is all old 
// if you really want to do it the hard way

#pragma once

#include "IIRAnalog.hpp"
#include "IIRFilters.hpp"
#include "IIRAnalogFilter.hpp"


namespace Filters::IIR::ButterworthFilters
{
    // k = pole position
    // n = order or number of poles
    std::complex<DspFloatType> butterworthpole(DspFloatType k, DspFloatType n)
    {
        DspFloatType p = M_PI * ((2 * k + n - 1) / (2 * n));
        return std::complex<DspFloatType>(-std::cos(p), std::sin(p));
    }
    std::complex<DspFloatType> cbutterworthpole(DspFloatType k, DspFloatType n)
    {
        // wc * exp(j PI * (2*k+n-1)/2n)
        std::complex<DspFloatType> p(0, M_PI * (2 * k + n - 1) / 2 * n);
        return exp(p);
    }
    std::complex<DspFloatType> ButterworthPoles(DspFloatType K, DspFloatType N)
    {
        DspFloatType theta = ((2*K+N-1)/(2*N))*M_PI;
        DspFloatType sk = cos(theta);
        DspFloatType omegak = sin(theta);
        std::complex<DspFloatType> p(sk,omegak);
        return p;
    }    
    
   

    BiquadSOS ButterworthLowPassSolver(int n)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {
            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = 1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 1;
            sos[index].z[1] = 0;
            sos[index].z[2] = 0;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthLowPassSolverRadius(int n, DspFloatType R)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = 1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -R * 2 * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 1;
            sos[index].z[1] = 0;
            sos[index].z[2] = 0;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthLowPassSolverQ(int n, DspFloatType Q)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = 1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (1 / (2 * Q)) * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 1;
            sos[index].z[1] = 0;
            sos[index].z[2] = 0;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthLowPassSolverRQ(int n, DspFloatType R, DspFloatType Q)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = 1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (R / (2 * Q)) * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 1;
            sos[index].z[1] = 0;
            sos[index].z[2] = 0;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }

    BiquadSOS ButterworthHighPassSolver(int n)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = -1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = 1;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthHighPassSolveQ(int n,DspFloatType Q)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = -1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (1.0/(2.0*Q)) * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = 1;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthHighPassSolveR(int n,DspFloatType R)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = -1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * R * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = 1;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthHighPassSolveRQ(int n,DspFloatType R, DspFloatType Q)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            sos[0].z[0] = 1;
            sos[0].z[1] = 0;
            sos[0].z[2] = 0;
            sos[0].p[0] = 1;
            sos[0].p[1] = -1;
            sos[0].p[2] = 0;
            total = (n - 1) / 2;
            index = 1;
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (R/(2.0*Q)) * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = 1;
            sos[index].p[0] = 1;
            sos[index].p[1] = X;
            sos[index].p[2] = 1;
            index++;
        }
        return sos;
    }   
    BiquadSOS ButterworthBandPassSolver(int n, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandpass must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType wu = Wu*Wl;
            DspFloatType wl = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandPassSolverQ(int n, DspFloatType Q, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandpass must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType wu = Wu*Wl;
            DspFloatType wl = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (1/2*Q)*std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandPassSolverR(int n, DspFloatType R, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandpass must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType wu = Wu*Wl;
            DspFloatType wl = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * R*std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandPassSolverRQ(int n, DspFloatType R, DspFloatType Q, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        int total = n / 2;
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandpass must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {
            DspFloatType wu = Wu*Wl;
            DspFloatType wl = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (R/(2*Q))*std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandStopSolver(int n, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        
        int total = n / 2;

        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandstop must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {            
            DspFloatType wl = Wu*Wl;
            DspFloatType wu = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandStopSolverQ(int n, DspFloatType Q, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        
        int total = n / 2;

        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandstop must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {            
            DspFloatType wl = Wu*Wl;
            DspFloatType wu = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (1.0/(2.0*Q))*std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandStopSolverR(int n, DspFloatType R, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        
        int total = n / 2;

        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandstop must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {            
            DspFloatType wl = Wu*Wl;
            DspFloatType wu = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * R* std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    BiquadSOS ButterworthBandStopSolverRQ(int n, DspFloatType R, DspFloatType Q, DspFloatType Wu, DspFloatType Wl)
    {
        BiquadSOS sos;
        
        int total = n / 2;

        if (n % 2 != 0)
            sos.resize(n / 2 + 1);
        else
            sos.resize(n / 2);

        
        int index = 0;
        if (n % 2 != 0)
        {

            printf("Bandstop must be even order\n");
            assert(1==0);
        }

        for (int k = 0; k < total; k++)
        {            
            DspFloatType wl = Wu*Wl;
            DspFloatType wu = Wu-Wl;
            DspFloatType wu2= wu*wu;
            DspFloatType wl2= wl*wl;

            DspFloatType K = (DspFloatType)(k + 1);
            DspFloatType N = (DspFloatType)n;
            DspFloatType X1 = (2 * K + N - 1);
            DspFloatType X2 = 2 * N;
            // DspFloatType X  = -2*(R/Q)*std::cos(M_PI*X1/X2);
            DspFloatType X = -2 * (R/(2.0*Q))*std::cos(M_PI * X1 / X2);
            sos[index].z[0] = 0;
            sos[index].z[1] = 0;
            sos[index].z[2] = wl2/wu2;
            sos[index].p[0] = 1;
            sos[index].p[1] = X*(-sqrt(2)*wl*wu-2*wu)/wu2;
            sos[index].p[2] = (wl2 + sqrt(2)*wl + 1)/wu2;            
            index++;
        }
        return sos;
    }
    
    struct ButterworthLowPassFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, Q;
        

        ButterworthLowPassFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q  = 0.5;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = q;
        }        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;

            BiquadSOS s           = ButterworthLowPassSolver(2);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthResonantLowPassFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, R;
        

        ButterworthResonantLowPassFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            R  = 0.5;
            setCutoff(1000.0f);
        }
        void setResonance(DspFloatType r)
        {
            R = r;
            setCutoff(Fc);
        }        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            BiquadSOS s           = ButterworthLowPassSolverRadius(2,R);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
                default: printf("No port %d\n",port);
            }
        }

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthDampedLowPassFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, Q;
        

        ButterworthDampedLowPassFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = q;
            setCutoff(Fc);
        }        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            BiquadSOS s           = ButterworthLowPassSolverQ(2,Q);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };


    struct ButterworthLowPassCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, Q;
        size_t order;
        BiquadSOS sos;

        ButterworthLowPassCascadeFilter(size_t n,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q  = 0.5;
            order = n;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = 1;
            for(size_t i = 0; i < order; i++)
                Q *= q;
        }
        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if( Fc >= Sr/2.0)
                return;

            BiquadSOS bw = ButterworthLowPassSolver(order);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr, order % 2 != 0);            
            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthResonantLowPassCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, R;
        size_t order;
        BiquadSOS sos;

        ButterworthResonantLowPassCascadeFilter(size_t n,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            R = 0.5;
            order = n;
            setCutoff(1000.0f);
        }
        void setResonance(DspFloatType r)
        {
            R = r;
            setCutoff(Fc);
        }
        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if( Fc >= Sr/2.0)
                return;

            BiquadSOS bw =  ButterworthLowPassSolverRadius(order,R);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr);

            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
            }
        }

        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    
    struct ButterworthDampedLowPassCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, Q;
        size_t order;
        BiquadSOS sos;

        ButterworthDampedLowPassCascadeFilter(size_t n,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            order = n;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = 1;
            for(size_t i = 0; i < order; i++)
                Q *= q;
        }        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            BiquadSOS bw = ButterworthLowPassSolverQ(order,Q);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr);

            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }

        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };


    struct ButterworthHighPassFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, Q;
        

        ButterworthHighPassFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = q;
        }        
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            BiquadSOS s = ButterworthHighPassSolver(2);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }

        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthHighPassCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, Q;
        size_t order;
        BiquadSOS sos;

        ButterworthHighPassCascadeFilter(size_t n,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            order = n;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = q;            
        }
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            BiquadSOS bw = ButterworthHighPassSolver(order);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr);

            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }

        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthBandPassFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, Q;
        DspFloatType bw = 0.2;

        ButterworthBandPassFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;            
            Q = 0.5;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q) {
            Q = q;
        }
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            DspFloatType Wu = Fc+bw;
            DspFloatType Wl = Fc-bw;
            if(Wu > 1) Wu = 0.995;
            if(Wl < 0.005) Wl = 0.005;

            
            BiquadSOS s = ButterworthBandPassSolver(2,Wu,Wl);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);         
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };
    
    struct ButterworthBandPassCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, Q,bw;
        size_t order;
        BiquadSOS sos;

        ButterworthBandPassCascadeFilter(size_t n, DspFloatType BW = 0.2,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            bw = BW;
            order = n;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = 1;
            for(size_t i = 0; i < order; i++)
                Q *= q;
        }
        void setBandWidth(DspFloatType BW) {
            bw = BW;
        }
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;
                            DspFloatType Wu = Fc+bw;
            DspFloatType Wl = Fc-bw;
            if(Wu > 1) Wu = 0.995;
            if(Wl < 0.005) Wl = 0.005;

            BiquadSOS bw = ButterworthBandPassSolver(order,Wu,Wl);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr);

            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_BW,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_BW: setBandWidth(v); break;
                default: printf("No port %d\n",port);
            }
        }

        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };

    struct ButterworthBandStopFilter12db
    {

        BiquadTypeI bc;
        DspFloatType Fc, Sr, Q;
        DspFloatType bw = 0.2;

        ButterworthBandStopFilter12db(DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q = 0.5;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q) {
            Q = q;
        }
        void setBandWidth(DspFloatType BW) {
            bw = BW;
        }
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;
                
            DspFloatType Wu = Fc+bw;
            DspFloatType Wl = Fc-bw;
            if(Wu > 1) Wu = 0.995;
            if(Wl < 0.005) Wl = 0.005;

            
            BiquadSOS s = ButterworthBandStopSolver(2,Wu,Wl);
            BiquadSection section = AnalogBiquadSection(s[0], fc, Sr);

            bc.setCoefficients(section);
        }
        enum {
            PORT_CUTOFF,
            PORT_BW,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_BW: setBandWidth(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };
    
    struct ButterworthBandStopCascadeFilter
    {

        BiquadTypeICascade bc;
        DspFloatType Fc, Sr, Q,bw;
        size_t order;
        BiquadSOS sos;

        

        ButterworthBandStopCascadeFilter(size_t n, DspFloatType BW = 0.2,DspFloatType sampleRate=44100.0)
        {            
            Sr = sampleRate;
            Q  = 0.5;
            bw = BW;
            order = n;
            setCutoff(1000.0f);
        }
        void setQ(DspFloatType q)
        {
            Q = 1;
            for(size_t i = 0; i < order; i++)
                Q *= q;
        }
        void setBandWidth(DspFloatType BW) {
            bw = BW;
        }
        void setCutoff(DspFloatType fc)
        {
            Fc = fc;

            if (Fc == 0)
                return;
            if (Fc >= Sr/2.0)
                return;

            DspFloatType Wu = Fc+bw;
            DspFloatType Wl = Fc-bw;
            if(Wu > 1) Wu = 0.995;
            if(Wl < 0.005) Wl = 0.005;

            BiquadSOS bw = ButterworthBandStopSolver(order,Wu,Wl);
            BiquadSOS sos = AnalogBiquadCascade(bw, fc, Sr);

            bc.setCoefficients(sos);
        }
        enum {
            PORT_CUTOFF,
            PORT_BW,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_BW: setBandWidth(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType y = bc.Tick(I);
            return y;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
		void ProcessInplace(size_t n, DspFloatType * out)
        {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(out[i]);
		}
    };
}
