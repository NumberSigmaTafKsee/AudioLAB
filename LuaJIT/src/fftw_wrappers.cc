#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"

using namespace std;

FFTW_R2C_1D_Executor::FFTW_R2C_1D_Executor(int n_real_samples) :
    input_size(n_real_samples),
    input_buffer(fftwf_alloc_real(n_real_samples)),
    output_size(n_real_samples/2 + 1),
    output_buffer((float _Complex*)fftwf_alloc_complex(n_real_samples/2 + 1))
{
    plan = fftwf_plan_dft_r2c_1d(n_real_samples, input_buffer, (float (*)[2])output_buffer, FFTW_ESTIMATE);
}

FFTW_R2C_1D_Executor::~FFTW_R2C_1D_Executor()
{
    fftwf_destroy_plan(plan);
    fftwf_free(input_buffer);
    fftwf_free(output_buffer);
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const float* buffer, int size)
{
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(float)*size);
    memset(&input_buffer[size], 0, sizeof(float)*(input_size - size));
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const vector<float>& vec)
{
    set_input_zeropadded(&vec[0], vec.size());
}

void FFTW_R2C_1D_Executor::execute()
{
    fftwf_execute(plan);
}

vector<float _Complex> FFTW_R2C_1D_Executor::get_output()
{
    return vector<float _Complex>(output_buffer, output_buffer + output_size);
}

FFTW_C2R_1D_Executor::FFTW_C2R_1D_Executor(int n_real_samples) : 
    input_size(n_real_samples/2 + 1),
    input_buffer((float _Complex*)fftwf_alloc_complex(n_real_samples/2 + 1)),
    output_size(n_real_samples),
    output_buffer(fftwf_alloc_real(n_real_samples))
{
    plan = fftwf_plan_dft_c2r_1d(n_real_samples, (float (*)[2])input_buffer, output_buffer, FFTW_ESTIMATE);
}

FFTW_C2R_1D_Executor::~FFTW_C2R_1D_Executor()
{
    fftwf_destroy_plan(plan);
    fftwf_free(input_buffer);
    fftwf_free(output_buffer);
}

void FFTW_C2R_1D_Executor::set_input(const float _Complex* buffer, int size)
{
    assert(size == input_size);
    memcpy(input_buffer, buffer, sizeof(float _Complex)*size);
    memset(&input_buffer[size], 0, sizeof(float _Complex)*(input_size - size));
}

void FFTW_C2R_1D_Executor::set_input(const vector<float _Complex>& vec)
{
    set_input(&vec[0], vec.size());
}

void FFTW_C2R_1D_Executor::execute()
{
    fftwf_execute(plan);
}

vector<float> FFTW_C2R_1D_Executor::get_output()
{
    return vector<float>(output_buffer, output_buffer + output_size);
}
