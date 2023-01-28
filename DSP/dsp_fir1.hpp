/*
License: MIT License (http://www.opensource.org/licenses/mit-license.php)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/* (C) 2013-2021 Graeme Hattan & Bernd Porr */

#pragma once

#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <stdexcept>

namespace Casino::Filters
{
    /**
    * Finite impulse response filter. The precision is DspFloatType.
    * It takes as an input a file with coefficients or an DspFloatType
    * array.
    **/
    class Fir1 {
    public:
        /** 
            * Coefficients as a const DspFloatType array. Because the array is const
            * the number of taps is identical to the length of the array.
            * \param _coefficients A const DspFloatType array with the impulse response.
            **/
        template <unsigned nTaps> Fir1(const DspFloatType (&_coefficients)[nTaps]) :
            coefficients(new DspFloatType[nTaps]),
            buffer(new DspFloatType[nTaps]()),
            taps(nTaps) {
            for(unsigned i=0;i<nTaps;i++) {
                coefficients[i] = _coefficients[i];
                buffer[i] = 0;
            }
        }

        /**
        * Coefficients as a C++ vector
        * \param _coefficients is a Vector of doubles.
        **/
        Fir1(std::vector<DspFloatType> _coefficients) {
            initWithVector(_coefficients);
        }

        /**
            * Coefficients as a (non-constant-) DspFloatType array where the length needs to be specified.
        * \param coefficients Coefficients as DspFloatType array.
            * \param number_of_taps Number of taps (needs to match the number of coefficients
            **/
        Fir1(DspFloatType *coefficients, unsigned number_of_taps);

        /** Coefficients as a text file (for example from Python)
        * The number of taps is automatically detected
        * when the taps are kept zero.
            * \param coeffFile Patht to textfile where every line contains one coefficient
            * \param number_of_taps Number of taps (0 = autodetect)
            **/
        Fir1(const char* coeffFile, unsigned number_of_taps = 0);

        /** 
            * Inits all coefficients and the buffer to zero
        * This is useful for adaptive filters where we start with
        * zero valued coefficients.
            **/
        Fir1(unsigned number_of_taps);

        /**
            * Releases the coefficients and buffer.
            **/
        ~Fir1();

        
        /**
            * The actual filter function operation: it receives one sample
            * and returns one sample.
            * \param input The input sample.
            **/
        inline DspFloatType filter(DspFloatType input) {
            const DspFloatType *coeff     = coefficients;
            const DspFloatType *const coeff_end = coefficients + taps;
            
            DspFloatType *buf_val = buffer + offset;
            
            *buf_val = input;
            DspFloatType output_ = 0;
            
            while(buf_val >= buffer)
                output_ += *buf_val-- * *coeff++;
            
            buf_val = buffer + taps-1;
            
            while(coeff < coeff_end)
                output_ += *buf_val-- * *coeff++;
            
            if(++offset >= taps)
                offset = 0;
            
            return output_;
        }


        /**
            * LMS adaptive filter weight update:
        * Every filter coefficient is updated with:
        * w_k(n+1) = w_k(n) + learning_rate * buffer_k(n) * error(n)
            * \param error Is the term error(n), the error which adjusts the FIR conefficients.
            **/
        inline void lms_update(DspFloatType error) {
            DspFloatType *coeff     = coefficients;
            const DspFloatType *coeff_end = coefficients + taps;
        
            DspFloatType *buf_val = buffer + offset;
            
            while(buf_val >= buffer) {
                *coeff++ += *buf_val-- * error * mu;
            }
            
            buf_val = buffer + taps-1;
            
            while(coeff < coeff_end) {
                *coeff++ += *buf_val-- * error * mu;
            }
        }

        /**
            * Setting the learning rate for the adaptive filter.
            * \param _mu The learning rate (i.e. rate of the change by the error signal)
            **/
        void setLearningRate(DspFloatType _mu) {mu = _mu;};

        /**
            * Getting the learning rate for the adaptive filter.
            **/
        DspFloatType getLearningRate() {return mu;};

        /**
            * Resets the buffer (but not the coefficients)
            **/
        void reset();

        /** 
            * Sets all coefficients to zero
            **/
        void zeroCoeff();

        /**
        * Copies the current filter coefficients into a provided array.
        * Useful after an adaptive filter has been trained to query
        * the result of its training.
        * \param coeff_data target where coefficients are copied
        * \param number_of_taps number of doubles to be copied
        * \throws std::out_of_range number_of_taps is less the actual number of taps.
        */
        void getCoeff(DspFloatType* coeff_data, unsigned number_of_taps) const;

        /**
        * Returns the coefficients as a vector
        **/
        std::vector<DspFloatType> getCoeffVector() const {
            return std::vector<DspFloatType>(coefficients,coefficients+taps);
        }

        /**
            * Returns the number of taps.
            **/
        unsigned getTaps() {return taps;};

        /**
            * Returns the power of the of the buffer content:
        * sum_k buffer[k]^2
        * which is needed to implement a normalised LMS algorithm.
            **/
        inline DspFloatType getTapInputPower() {
            DspFloatType *buf_val = buffer;
            
            DspFloatType p = 0;
            
            for(unsigned i = 0; i < taps; i++) {
                p += (*buf_val) * (*buf_val);
                buf_val++;
            }
        
            return p;
        }

    private:
        void initWithVector(std::vector<DspFloatType> _coefficients);
        
        DspFloatType        *coefficients;
        DspFloatType        *buffer;
        unsigned      taps;
        unsigned      offset = 0;
        DspFloatType        mu = 0;
    };


    // give the filter an array of doubles for the coefficients
    Fir1::Fir1(DspFloatType *_coefficients, unsigned number_of_taps) :
        coefficients(new DspFloatType[number_of_taps]),
        buffer(new DspFloatType[number_of_taps]()),
        taps(number_of_taps) {
        for(unsigned int i=0;i<number_of_taps;i++) {
            coefficients[i] = _coefficients[i];
            buffer[i] = 0;
        }
    }

    // init all coefficients and the buffer to zero
    Fir1::Fir1(unsigned number_of_taps) :
        coefficients(new DspFloatType[number_of_taps]),
        buffer(new DspFloatType[number_of_taps]),  
        taps(number_of_taps) {
        zeroCoeff();
        reset();
    }

    void Fir1::initWithVector(std::vector<DspFloatType> _coefficients) {
        coefficients = new DspFloatType[_coefficients.size()];
        buffer = new DspFloatType[_coefficients.size()]();
        taps = ((unsigned int)_coefficients.size());
        for(unsigned long i=0;i<_coefficients.size();i++) {
            coefficients[i] = _coefficients[i];
            buffer[i] = 0;
        }
    }	

    // one coefficient per line
    Fir1::Fir1(const char* coeffFile, unsigned number_of_taps) {

        std::vector<DspFloatType> tmpCoefficients;

        FILE* f=fopen(coeffFile,"rt");
        if (!f) {
            throw std::invalid_argument("Could not open file.");
        }
        for(unsigned int i=0;(i<number_of_taps)||(number_of_taps==0);i++) {
            DspFloatType v = 0;
            int r = fscanf(f,"%lf\n",&v);
            if (r < 1) break;
            tmpCoefficients.push_back(v);
        }
        fclose(f);
        initWithVector(tmpCoefficients);
    }


    Fir1::~Fir1()
    {
        delete[] buffer;
        delete[] coefficients;
    }


    void Fir1::reset()
    {
        memset(buffer, 0, sizeof(DspFloatType)*taps);
        offset = 0;
    }

    void Fir1::zeroCoeff() {
        memset(coefficients, 0, sizeof(DspFloatType)*taps);
        offset = 0;
    }

    void Fir1::getCoeff(DspFloatType* coeff_data, unsigned number_of_taps) const {
        
        if (number_of_taps < taps)
            throw std::out_of_range("Fir1: target of getCoeff: too many weights to copy into target");
    
        memcpy(coeff_data, coefficients, taps * sizeof(DspFloatType));
        if (number_of_taps > taps)
            memset(&coeff_data[taps], 0, (number_of_taps - taps)*sizeof(DspFloatType));
    }
}


