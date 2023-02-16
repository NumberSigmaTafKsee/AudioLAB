#pragma once

#include <Eigen/Eigen>
#include <cmath>
#include <complex>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <future>
#include <vector>
#include <cmath>
#include <chrono>

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;

static double const EAR_Q = 9.26449;				//  Glasberg and Moore Parameters
static double const MIN_BW = 24.7;
static int const ORDER = 1;


namespace Filters::Cochlear
{
    class Timer
    {
    public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count(); }

    private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
    };

    class Biquad
    {
    public:
        Biquad(void);
        Biquad(Eigen::Matrix<double, 5, 1> coeffs);
        ~Biquad(void);
        // Methods
        Eigen::VectorXd process(const Eigen::VectorXd &x, int n);
        Eigen::Matrix<double, 5, 1> coeffs();
    private:
        double a0, a1, a2;
        double b1, b2;
    };

    class CochlearFilterbank
    {
    public:
        // Constructors and destructors
        CochlearFilterbank(double fs, int num_channels, double low_freq);
        ~CochlearFilterbank(void);
        // Methods
        Eigen::MatrixXd process(const Eigen::VectorXd &input, int n);
        Eigen::VectorXd process_channel(const Eigen::VectorXd &input, int n, int ch);
        static Eigen::VectorXd ERBspace(double low_freq, double high_freq, int num_channels);
        // Constants are defined on source file
    private:
        std::vector<std::vector<Biquad> > filters;
        std::vector<std::vector<Biquad> > makeERBFilters(double fs, int num_channels, double low_freq);
    };

    class ModulationFilterBank
    {
    public:
        ModulationFilterBank(double fs, int num_channels, Eigen::VectorXd mf, double q);
        ~ModulationFilterBank(void);
        Eigen::MatrixXd process(const Eigen::VectorXd &input, int n);
        Eigen::VectorXd process_channel(const Eigen::VectorXd &input, int n, int ch);
        static Eigen::VectorXd compute_modulation_cfs(double min_cf, double max_cf, int n);
        // Constants are defined on source file
    private:
        std::vector<Biquad> filters;
        std::vector<Biquad> makeModulationFilters(double fs, int num_channels, Eigen::VectorXd mf, double q);
    };


    Biquad::Biquad(void)
    {
        a0 = 0.0;
        a1 = 0.0;
        a2 = 0.0;
        b1 = 0.0;
        b2 = 0.0;
    }


    Biquad::Biquad(Eigen::Matrix<double, 5, 1> coeffs)
    {
        a0 = coeffs[0];
        a1 = coeffs[1];
        a2 = coeffs[2];
        b1 = coeffs[3];
        b2 = coeffs[4];
    }

    Eigen::Matrix<double, 5, 1> Biquad::coeffs()
    {
        Eigen::Matrix<double, 5, 1> coeffs;
        coeffs << a0, a1, a2, b1, b2;
        return coeffs;
    }

    Biquad::~Biquad(void)
    {
    }






    CochlearFilterbank::CochlearFilterbank(double fs, int num_channels, double low_freq)
    {
        filters = makeERBFilters(fs, num_channels, low_freq);
    }


    CochlearFilterbank::~CochlearFilterbank(void)
    {

    }

    vector< vector<Biquad> >CochlearFilterbank::makeERBFilters(double fs, int num_channels, double low_freq)
    {
        double T = 1.0/fs;
        Eigen::VectorXd ERB(num_channels);
        Eigen::VectorXd B(num_channels);
        Eigen::VectorXd cf = ERBspace(low_freq, fs/2, num_channels);
        double A0 = T;
        double A2 = 0;
        double B0 = 1;
        Eigen::VectorXd B1(num_channels);
        Eigen::VectorXd B2(num_channels);
        Eigen::VectorXd A11(num_channels);
        Eigen::VectorXd A12(num_channels);
        Eigen::VectorXd A13(num_channels);
        Eigen::VectorXd A14(num_channels);
        Eigen::VectorXd gain(num_channels);

        vector< vector<Biquad> > filter_bank = vector< vector<Biquad> >(num_channels);

        complex<double> i (0,1);
        complex<double> aux1, aux2, aux3, aux4, aux5, aux6;

        for (int k=0; k < num_channels; k++)
        {
            ERB[k] = pow(pow((cf[k]/EAR_Q),ORDER) + pow(MIN_BW,ORDER),1/ORDER);
            B[k] = 1.019*2*M_PI*ERB[k];
            B1[k] = -2.0*cos(2*cf[k]*M_PI*T)/exp(B[k]*T);
            B2[k] = exp(-2*B[k]*T);
            A11[k] = -(-B1[k]*T + 2.0*sqrt(3.0+pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
            A12[k] = -(-B1[k]*T - 2.0*sqrt(3.0+pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
            A13[k] = -(-B1[k]*T + 2.0*sqrt(3.0-pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;
            A14[k] = -(-B1[k]*T - 2.0*sqrt(3.0-pow(2.0,1.5))*T*sin(2.0*cf[k]*M_PI*T)/exp(B[k]*T))/2.0;

            aux1 = complex<double>(-2)*exp(complex<double>(4)*i*cf[k]*M_PI*T)*T;
            aux2 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) - sqrt(3.0 - pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
            aux3 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) + sqrt(3.0 - pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
            aux4 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) - sqrt(3.0 + pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
            aux5 = complex<double>(2)*exp(-(B[k]*T) + complex<double>(2)*i*cf[k]*M_PI*T)*T*(cos(2*cf[k]*M_PI*T) + sqrt(3.0 + pow(2.0,1.5))*sin(2*cf[k]*M_PI*T));
            aux6 = pow(-2 / exp(2*B[k]*T) - complex<double>(2)*exp(complex<double>(4)*i*cf[k]*M_PI*T) + complex<double>(2)*(complex<double>(1) + exp(complex<double>(4)*i*cf[k]*M_PI*T))/exp(B[k]*T),4);
            gain[k] = abs(abs((aux1 + aux2)*(aux1 + aux3)*(aux1 + aux4)*(aux1 + aux5)/aux6));
            vector<Biquad> filters = vector<Biquad>(4);
            Eigen::Matrix<double, 5, 1> coeffs1, coeffs2, coeffs3, coeffs4;
            coeffs1 << A0*(1.0/gain[k]), A11[k]*(1.0/gain[k]), A2*(1.0/gain[k]), B1[k], B2[k];
            coeffs2 << A0, A12[k], A2, B1[k], B2[k];
            coeffs3 << A0, A13[k], A2, B1[k], B2[k];
            coeffs4 << A0, A14[k], A2, B1[k], B2[k];
    #ifdef DEBUG
            cout << "Coeffs for channel 1:" << endl;
            cout << "Filter 1" << endl << coeffs1 << endl << endl;
            cout << "Filter 2" << endl << coeffs2 << endl << endl;
            cout << "Filter 3" << endl << coeffs3 << endl << endl;
            cout << "Filter 4" << endl << coeffs4 << endl << endl;
    #endif
            filters[0] = Biquad(coeffs1);
            filters[1] = Biquad(coeffs2);
            filters[2] = Biquad(coeffs3);
            filters[3] = Biquad(coeffs4);
            filter_bank[k] = filters;
        }
        return filter_bank;
    }

    Eigen::VectorXd CochlearFilterbank::ERBspace(double low_freq, double high_freq, int num_channels)
    {
        Eigen::VectorXd cf_array(num_channels);
        double aux = EAR_Q * MIN_BW;
        for (int i=1; i <= num_channels; i++)
        {
            cf_array[i-1] = -(aux) + exp((i)*(-log(high_freq + aux) + log(low_freq + aux))/num_channels) * (high_freq + aux);
        }
        return cf_array;
    }

    Eigen::MatrixXd CochlearFilterbank::process(const Eigen::VectorXd &input, int n)
    {
        Eigen::MatrixXd y(input.rows(), filters.size());
        vector<future<Eigen::VectorXd>> futures;
        for (unsigned int ch=0; ch < filters.size(); ch++)
            futures.push_back(async(launch::async, &CochlearFilterbank::process_channel, this, input, n, ch));
        for (unsigned int ch=0; ch < filters.size(); ch++)
        {
            y.col(ch) = futures[ch].get();
        }
        return y;
    }

    Eigen::VectorXd CochlearFilterbank::process_channel(const Eigen::VectorXd &input, int n, int ch)
    {
            Eigen::VectorXd y1(Eigen::VectorXd::Zero(input.rows()));
            Eigen::VectorXd y2(Eigen::VectorXd::Zero(input.rows()));
            Eigen::VectorXd y3(Eigen::VectorXd::Zero(input.rows()));
            Eigen::VectorXd y4(Eigen::VectorXd::Zero(input.rows()));

            y1 = filters[ch][0].process(input, n);
            y2 = filters[ch][1].process(y1, n);
            y3 = filters[ch][2].process(y2, n);
            y4 = filters[ch][3].process(y3, n);
            return y4;
    }

    ModulationFilterBank::ModulationFilterBank(double fs, int num_channels, Eigen::VectorXd mf, double q)
    {
        filters = makeModulationFilters(fs, num_channels, mf, q);
    }


    ModulationFilterBank::~ModulationFilterBank(void)
    {
    }

    Eigen::VectorXd ModulationFilterBank::compute_modulation_cfs(double min_cf, double max_cf, int n)
    {
        double spacing_factor = pow(max_cf/min_cf,1.0/(n-1));
        Eigen::VectorXd cfs(n);
        cfs[0] = min_cf;
        for (int k = 1; k < n; k++)
            cfs[k] = cfs[k - 1]*spacing_factor;
        return cfs;
    }

    Eigen::MatrixXd ModulationFilterBank::process(const Eigen::VectorXd &input, int n)
    {
        Eigen::MatrixXd y(input.rows(), filters.size());
        vector<future<Eigen::VectorXd>> futures;
        for (unsigned int ch=0; ch < filters.size(); ch++)
            futures.push_back(async(&ModulationFilterBank::process_channel, this, input, n, ch));
        for (unsigned int ch=0; ch < filters.size(); ch++)
        {
            y.col(ch) = futures[ch].get();
        }
        return y;

    }

    Eigen::VectorXd ModulationFilterBank::process_channel(const Eigen::VectorXd &input, int n, int ch)
    {
        return filters[ch].process(input, n);
    }

    std::vector<Biquad> ModulationFilterBank::makeModulationFilters(double fs, int num_channels, Eigen::VectorXd mf, double q)
    {
        vector<Biquad> filters = vector<Biquad>(num_channels);
        for (int k=0; k < num_channels; k++)
        {
        double w0 = 2*M_PI*mf(k)/fs;
            double W0 = tan(w0/2);
            double W02 = pow(W0,2);
            double B0 = W0/q;

            double a0 = 1 + B0 + W02;
            double a1 = 2*W02 - 2;
            double a2 = 1 - B0 + W02;

            Eigen::Matrix<double, 5, 1> f;

            f << B0, 0, -B0, a1, a2;
    #ifdef DEBUG
            cout << "[" << k << "]: " << f << endl;
    #endif
            filters[k] = Biquad(f/a0);
        }
        return filters;
    }
}