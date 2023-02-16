//https://github.com/FelixLaufer/ButterworthLP/blob/master/ButterworthLP.cpp
// MIT License, Copyright (c) 2016 Felix Laufer

#include <vector>
#include <complex>
#include <cmath>

#ifdef DEBUG
    #include <iostream>
    #include <iomanip>
#endif

namespace Filters
{
    /**
    Analog prototype filter is digitalized using bilinear transform and split up into cascaded second-order sections (SOS).
    */
    class ButterworthLP : public FilterProcessor
    {
    public:

        /**
        Generates a Butterworth lowpass filter with a given normalized cutoff frequency and filter order.
        @param normalizedCutoffFrequency   (0, 1)     Normalized cutoff frequency := cuttoffFreq / samplingFreq.
        @param filterOrder                 [1, inf]   Butterworth filter order.
        */
        ButterworthLP(const double normalizedCutoffFrequency, const size_t filterOrder);

        /**
        Generates a Butterworth lowpass filter of a given order with a given cutoff frequency [Hz] in respect of the data sampling frequency [Hz].
        @param samplingFrequency   [1, inf] [Hz]                     Sampling frequency of the data.
        @param cutoffFrequency     [1, samplingFrequency - 1] [Hz]   Cutoff frequency
        @param filterOrder         [1, inf]                          Butterworth filter order.
        */
        ButterworthLP(const double samplingFrequency, const double cutoffFrequency, const size_t filterOrder);

        /**
        Set the filter state to a steady state w.r.t. to the given input value (assuming a constant filter input for infinite time steps in the past;
        see also http://www.emt.tugraz.at/publications/diplomarbeiten/da_hoebenreich/node21.html).
        @param value    Desired steady state output value
        */
        void stepInitialization(const double value);

        /**
        Processes the input value online depending on the current filter state.
        @param input   [-inf, inf]   Input value
        @return        [-inf, inf]   Filter response
        */
        double process(const double input);

        double Tick(double I, double A=1,double X=0, double Y=0) {
            return A*process(I);
        }

        /**
        Processes the input data offline.
        @param *input                                  Ptr to input data.
        @param *output                                 Ptr to output data.
        @param size                    [0, inf]        Length of data.
        @param initialConditionValue   [-inf, inf]     Initializes the filter state to a steady state w.r.t. the given initial value (see stepInitialization).
        @param forwardBackward         {true, false}   Eliminate phase delay by filtering twice: forward and backward.
                                                    Note that forward-backward filtering corresponds to filtering with a 2n-th order filter.
        */
        void filter(const double *input, double *output, const size_t size, const double initialConditionValue = 0.0, const bool forwardBackward = false);

        ~ButterworthLP()
        { }

    private:

        class SOS
        {

        public:

            SOS(const double b0, const double b1, const double b2, const double a1, const double a2, const double gain) :
                _b0(b0), _b1(b1), _b2(b2), _a1(a1), _a2(a2), _gain(gain), _z1(0), _z2(0),
                _preCompStateSpaceOutputVec1(_b1 - _b0*_a1),
                _preCompStateSpaceOutputVec2(_b2 - _b0*_a2)
            {
                #ifdef DEBUG
                    std::cout << std::fixed << std::setprecision(6);
                    size_t w = 10;
                    std::cout <<_b0 << std::setw(w) << _b1 << std::setw(w) << _b2 << std::setw(w) << 1.0 << std::setw(w) << _a1 << std::setw(w) << _a2 << std::endl;
                #endif
            }

            void safeRestoreState(double &z1, double &z2, const bool restore = false);

            void stepInitialization(const double value);

            double process(const double input);

            ~SOS()
            { }

        private:

            const double _b0, _b1, _b2, _a1, _a2, _gain;

            const double _preCompStateSpaceOutputVec1, _preCompStateSpaceOutputVec2;

            double _z1, _z2;

        };

        void addSOS(const SOS sos);

        bool coefficients(const double normalizedCutoffFrequency, const size_t filterOrder);

        size_t _numSOS;

        std::vector<SOS> _sosVec;

    };

    // MIT License, Copyright (c) 2016 Felix Laufer

    
    

    ButterworthLP::ButterworthLP(const double normalizedCutoffFrequency, const size_t filterOrder) :
        FilterProcessor(),
        _numSOS(0)
    {
        if (!coefficients(normalizedCutoffFrequency, filterOrder))
        {
            throw std::domain_error(std::string("Failed to design a filter due to invalid parameters (normalized cutoff frequency and / or filter order) or instability of the resulting digitalized filter."));
        }
    }

    ButterworthLP::ButterworthLP(const double samplingFrequency, const double cutoffFrequency, const size_t filterOrder) :
        ButterworthLP(cutoffFrequency / samplingFrequency, filterOrder)
    { }

    void ButterworthLP::SOS::safeRestoreState(double &z1, double &z2, const bool restore)
    {
        if (restore)
        {
            _z1 = z1;
            _z2 = z2;
        }
        else
        {
            z1 = _z1;
            z2 = _z2;
        }
    }

    void ButterworthLP::SOS::stepInitialization(const double value)
    {
        // Set second-order section state to steady state w.r.t. given step response value
        // http://www.emt.tugraz.at/publications/diplomarbeiten/da_hoebenreich/node21.html
        _z1 = _z2 = value / (1.0 + _a1 + _a2);
    }

    double ButterworthLP::SOS::process(const double input)
    {
        /**
        SOS state space model
        z' = A * z + b * x
        y = c^T * z + d * x
        with A = [-a1  -a2; 1     0]   b = [1; 0]   c = [b1 - b0*a1; b2 - b0*a2]   d = [b0]
        */

        double x = input;
        double y = x;
        double z1_new, z2_new;

        z1_new = -_a1 * _z1 - _a2 * _z2 + x;
        z2_new = _z1;
        y = _preCompStateSpaceOutputVec1 * _z1 + _preCompStateSpaceOutputVec2 * _z2 + _b0 * x;
        _z1 = z1_new;
        _z2 = z2_new;

        // Include SOS gain factor
        y *= _gain;

        return y;
    }

    void ButterworthLP::addSOS(const SOS sos)
    {
        _sosVec.push_back(sos);
        ++_numSOS;
    }

    void ButterworthLP::stepInitialization(const double value)
    {
        double stepResponseValue = value;

        // Propagate step initialization through all second-order sections
        for (size_t i = 0; i < _numSOS; ++i)
        {
            _sosVec[i].stepInitialization(stepResponseValue);
            stepResponseValue = _sosVec[i].process(stepResponseValue);
        }
    }

    double ButterworthLP::process(const double input)
    {
        double x = input;
        double y = x;

        // Cascade all second-order sections s.t. output of SOS i is input for SOS i+1
        for (size_t i = 0; i < _numSOS; ++i)
        {
            y = _sosVec[i].process(x);
            x = y;
        }

        return y;
    }

    void ButterworthLP::filter(const double *input, double *output, const size_t size, const double initialConditionValue, const bool forwardBackward)
    {
        // Save all current SOS states before offline filtering
        std::vector<double> zSaved(_numSOS * 2);
        for (size_t i = 0; i < _numSOS; ++i)
        {
            _sosVec[i].safeRestoreState(zSaved[i], zSaved[i + 1], false);
        }

        // Set initial step response conditions
        stepInitialization(initialConditionValue);

        // Filtering on input data
        for (size_t i = 0; i < size; ++i)
        {
            output[i] = process(input[i]);
        }

        // Additional backward filtering on filtered output data if requested
        if (forwardBackward)
        {
            for (size_t i = size; i > 0;)
            {
                output[--i] = process(output[i]);
            }
        }

        // Restore all SOS states
        for (size_t i = 0; i < _numSOS; ++i)
        {
            _sosVec[i].safeRestoreState(zSaved[i], zSaved[i + 1], true);
        }
    }

    bool ButterworthLP::coefficients(const double normalizedCutoffFrequency, const size_t filterOrder)
    {
        // Assure valid parameters
        if (filterOrder < 1 || normalizedCutoffFrequency <= 0 || normalizedCutoffFrequency > 1)
        {
            return false;
        }

        std::vector<std::complex<double>> poles(filterOrder);

        // Prewarp the analog prototype's cutoff frequency
        double omegaCutoff = 2 * tan(M_PI * normalizedCutoffFrequency);

        double gain = pow(omegaCutoff, filterOrder);
        double initialGain = gain;

        std::complex<double> two(2.0, 0);

        for (size_t i = 0, i_end = (filterOrder + 1) / 2; i < i_end; ++i)
        {
            size_t i2 = 2 * i;

            /** 
            Design the analog prototype Butterworth lowpass filter
            */

            // Generate s-poles of prototype filter
            double phi = (double)(i2 + 1) * M_PI / (2 * filterOrder);
            double real = -sin(phi);
            double imag = cos(phi);

            std::complex<double> pole = std::complex<double>(real, imag);

            /**
            Customize analog prototype filter w.r.t cutoff frequency
            */

            // Scale s-pole with the cutoff frequency
            pole *= omegaCutoff;

            /**
            Digitalize the analog filter
            */

            // Map pole from s-plane to z-plane using bilinear transform
            std::complex<double> s = pole;
            pole = (two + s) / (two - s);

            // Update overall gain in respect of z-pole gain
            gain *= abs((two - s));

            // Ensure z-pole lies in unit circle of z-plane
            if (abs(pole) > 1)
            {
                return false;
            }

            // Add stable z-pole
            poles[i2] = pole;

            // Odd filter order: ignore the second complex conjugate pole
            if (i2 + 1 >= filterOrder)
            {
                break;
            }

            // Do the same as above with the conjugate complex pole
            pole = std::complex<double>(real, -imag);
            pole *= omegaCutoff;
            s = pole;
            pole = (two + s) / (two - s);
            gain *= abs((two - s));
            if (abs(pole) > 1)
            {
                return false;
            }
            poles[i2 + 1] = pole;
        }

        // Distribute the overall gain over all z-poles
        double overallGain = initialGain * (initialGain / gain);
        double distributedPoleGain = pow(overallGain, 1.0 / (double)filterOrder);
        double distributedPolePairGain = distributedPoleGain * distributedPoleGain;

        /**
        Generate second-order sections from conjugate complex z-pole pairs
        */

        for (size_t i = 0, i_end = filterOrder - 1; i < i_end; i += 2)
        {
            addSOS(SOS(1.0, 2.0, 1.0, -(poles[i] + poles[i + 1]).real(), (poles[i] * poles[i + 1]).real(), distributedPolePairGain));
        }

        // Odd filter order: remaining single z-pole requires additional second-order section
        if (filterOrder % 2 == 1)
        {
            addSOS(SOS(1.0, 1.0, 0.0, -poles[filterOrder - 1].real(), 0.0, distributedPoleGain));
        }

        return true;
    }
}