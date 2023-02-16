#include <cmath>
#include <algorithm>

/*
How to use it:
--------------
free_comp.cpp and free_comp.hpp implement a simple C++ class called free_comp.
(People are expected to know what to do with it! If not seek help on a beginner
programming forum!)
The free_comp class implements the following methods:
    void    set_threshold(double value);
            Sets the threshold level; 'value' must be a _positive_ value
            representing the amplitude as a floating point figure which should be
            1.0 at 0dBFS

    void    set_ratio(double value);
            Sets the ratio; 'value' must be in range [0.0; 1.0] with 0.0
            representing a oo:1 ratio, 0.5 a 2:1 ratio; 1.0 a 1:1 ratio and so on

    void    set_attack(double value);
            Sets the attack time; 'value' gives the attack time in _samples_

    void    set_release(double value);
            Sets the release time; 'value' gives the release time in _samples_

    void    set_output(double value);
            Sets the output gain; 'value' represents the gain, where 0dBFS is 1.0
            (see set_threshold())

    void    reset();
            Resets the internal state of the compressor (gain reduction)

    void    process(double *input_left, double *input_right,
                                    double *output_left, double *output_right,
                                    int frames, int skip);
    void    process(double *input_left, double *input_right,
                                    double *output_left, double *output_right,
                                    int frames, int skip);
            Processes a stereo stream of length 'frames' from either two arrays of
            floats or arrays of doubles 'input_left' and 'input_right' then puts
            the processed data in 'output_left' and 'output_right'.
            'input_{left,right}' and 'output_{left,right}' may be the same location
            in which case the algorithm will work in place. '{input,output}_left'
            and '{input,output}_right' can also point to the same data, in which
            case the algorithm works in mono (although if you process a lot of mono
            data it will yield more performance if you modify the source to make the
            algorithm mono in the first place).
            The 'skip' parameter allows for processing of interleaved as well as two
            separate contiguous streams. For two separate streams this value should
            be 1, for interleaved stereo it should be 2 (but it can also have other
            values than that to process specific channels in an interleaved audio
            stream, though if you do that it is highly recommended to study the
            source first to check whether it yields the expected behaviour or not).
*/         

namespace Dynamics
{
        class compressor : StereoFXProcessor
        {

        private:
                double   threshold;
                double   attack, release, envelope_decay;
                double   output;
                double   transfer_A, transfer_B;
                double   env, gain;

        public:
        compressor() : StereoFXProcessor()
        {
                threshold = 1.f;
                attack = release = envelope_decay = 0.f;
                output = 1.f;

                transfer_A = 0.f;
                transfer_B = 1.f;

                env = 0.f;
                gain = 1.f;
        }

        void set_threshold(double value)
        {
                threshold = value;
                transfer_B = output * pow(threshold,-transfer_A);
        }


        void set_ratio(double value)
        {
                transfer_A = value-1.f;
                transfer_B = output * pow(threshold,-transfer_A);
        }


        void set_attack(double value)
        {
                attack = exp(-1.f/value);
        }


        void et_release(double value)
        {
                release = exp(-1.f/value);
                envelope_decay = exp(-4.f/value); /* = exp(-1/(0.25*value)) */
        }


        void set_output(double value)
        {
                output = value;
                transfer_B = output * pow(threshold,-transfer_A);
        }


        void reset()
        {
                env = 0.f; gain = 1.f;
        }


        void process(float *input_left, float *input_right, float *output_left, float *output_right,       int frames)
        {
                double det, transfer_gain;
                for(int i=0; i<frames; i++)
                {
                        det = std::max(std::fabs(input_left[i]),std::fabs(input_right[i]));
                        det += 10e-30f; /* add tiny DC offset (-600dB) to prevent denormals */

                        env = det >= env ? det : det+envelope_decay*(env-det);

                        transfer_gain = env > threshold ? std::pow(env,transfer_A)*transfer_B:output;

                        gain = transfer_gain < gain ?
                                                        transfer_gain+attack *(gain-transfer_gain):
                                                        transfer_gain+release*(gain-transfer_gain);

                        output_left[i] = input_left[i] * gain;
                        output_right[i] = input_right[i] * gain;
                }
        }


        void process(double *input_left, double *input_right,     double *output_left, double *output_right,int frames)
        {
                double det, transfer_gain;
                for(int i=0; i<frames; i++)
                {
                        det = std::max(std::fabs(input_left[i]),std::fabs(input_right[i]));
                        det += 10e-30f; /* add tiny DC offset (-600dB) to prevent denormals */

                        env = det >= env ? det : det+envelope_decay*(env-det);

                        transfer_gain = env > threshold ? std::pow(env,transfer_A)*transfer_B:output;

                        gain = transfer_gain < gain ?
                                                        transfer_gain+attack *(gain-transfer_gain):
                                                        transfer_gain+release*(gain-transfer_gain);

                        output_left[i] = input_left[i] * gain;
                        output_right[i] = input_right[i] * gain;
                }
        }

        void ProcessBlock)size_t n, float ** in, float ** out)
        {
                process(in[0],in[1],out[0],out[1],n);
        }
        void ProcessBlock)size_t n, double ** in, double ** out)
        {
                process(in[0],in[1],out[0],out[1],n);
        }
        };
}