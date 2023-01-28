#include <cmath>
#define max(a,b) (a>b?a:b)

class compressor
{

    private:
            float   threshold;
            float   attack, release, envelope_decay;
            float   output;
            float   transfer_A, transfer_B;
            float   env, gain;

    public:
    compressor()
    {
            threshold = 1.f;
            attack = release = envelope_decay = 0.f;
            output = 1.f;

            transfer_A = 0.f;
            transfer_B = 1.f;

            env = 0.f;
            gain = 1.f;
    }

    void set_threshold(float value)
    {
            threshold = value;
            transfer_B = output * pow(threshold,-transfer_A);
    }


    void set_ratio(float value)
    {
            transfer_A = value-1.f;
            transfer_B = output * pow(threshold,-transfer_A);
    }


    void set_attack(float value)
    {
            attack = exp(-1.f/value);
    }


    void et_release(float value)
    {
            release = exp(-1.f/value);
            envelope_decay = exp(-4.f/value); /* = exp(-1/(0.25*value)) */
    }


    void set_output(float value)
    {
            output = value;
            transfer_B = output * pow(threshold,-transfer_A);
    }


    void reset()
    {
            env = 0.f; gain = 1.f;
    }


    __forceinline void process(float *input_left, float *input_right,float *output_left, float *output_right,       int frames)
    {
            float det, transfer_gain;
            for(int i=0; i<frames; i++)
            {
                    det = max(fabs(input_left[i]),fabs(input_right[i]));
                    det += 10e-30f; /* add tiny DC offset (-600dB) to prevent denormals */

                    env = det >= env ? det : det+envelope_decay*(env-det);

                    transfer_gain = env > threshold ? pow(env,transfer_A)*transfer_B:output;

                    gain = transfer_gain < gain ?
                                                    transfer_gain+attack *(gain-transfer_gain):
                                                    transfer_gain+release*(gain-transfer_gain);

                    output_left[i] = input_left[i] * gain;
                    output_right[i] = input_right[i] * gain;
            }
    }


    __forceinline void process(double *input_left, double *input_right,     double *output_left, double *output_right,int frames)
    {
            double det, transfer_gain;
            for(int i=0; i<frames; i++)
            {
                    det = max(fabs(input_left[i]),fabs(input_right[i]));
                    det += 10e-30f; /* add tiny DC offset (-600dB) to prevent denormals */

                    env = det >= env ? det : det+envelope_decay*(env-det);

                    transfer_gain = env > threshold ? pow(env,transfer_A)*transfer_B:output;

                    gain = transfer_gain < gain ?
                                                    transfer_gain+attack *(gain-transfer_gain):
                                                    transfer_gain+release*(gain-transfer_gain);

                    output_left[i] = input_left[i] * gain;
                    output_right[i] = input_right[i] * gain;
            }
    }

};
