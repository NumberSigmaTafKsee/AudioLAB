%module fftw_convolve
%{
#include "fftw_wrappers.hh"
#include <vector>

using std::vector;
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector)     std::vector<float>;


%inline %{

   
    // This function computes the discrete convolution of two arrays:
    // result[i] = a[i]*b[0] + a[i-1]*b[1] + ... + a[0]*b[i]
    // a and b can be vectors of different lengths, this function is careful to never
    // exceed the bounds of the vectors.
    vector<float> convolve(const vector<float>& a, const vector<float>& b)
    {
        int n_a = a.size();
        int n_b = b.size();
        vector<float> result(n_a + n_b - 1);

        for (int i = 0; i < n_a + n_b - 1; ++i) {
            float sum = 0.0;
            for (int j = 0; j <= i; ++j) {
                sum += ((j < n_a) && (i-j < n_b)) ? a[j]*b[i-j] : 0.0;
            }
            result[i] = sum;
        }
        return result;
    }

    template <class T>
    vector<T> vector_elementwise_multiply(const vector<T> a, const vector<T> b)
    {
        assert(a.size() == b.size());
        vector<T> result(a.size());
        for (int i = 0; i < result.size(); ++i) {
            result[i] = a[i]*b[i];
        }
        return result;
    }

    // Convolution of real vectors using the Fast Fourier Transform and the convolution theorem.
    // See http://en.wikipedia.org/w/index.php?title=Convolution&oldid=630841165#Fast_convolution_algorithms
    vector<float> fftw_convolve(vector<float>& a, vector<float>& b)
    {
        // Recall that element-wise
        int padded_length = a.size() + b.size() - 1;
        
        // Compute Fourier transform of vector a
        
        FFTW_R2C_1D_Executor fft_a(padded_length);
        fft_a.set_input_zeropadded(a);

        //cout << "a: ";
        //print_float_array(fft_a.input_buffer, fft_a.input_size);

        fft_a.execute();

        //cout << "FFT(a): ";
        //print_complex_array(fft_a.output_buffer, fft_a.output_size);
        //cout << endl;

        // Compute Fourier transform of vector b
        
        FFTW_R2C_1D_Executor fft_b(padded_length);
        fft_b.set_input_zeropadded(b);

        //cout << "b: ";
        //print_float_array(fft_b.input_buffer, fft_b.input_size);

        fft_b.execute();

        //cout << "FFT(b): ";
        //print_complex_array(fft_b.output_buffer, fft_b.output_size);
        //cout << endl;

        // Perform element-wise product of FFT(a) and FFT(b)
        // then compute inverse fourier transform.
        FFTW_C2R_1D_Executor ifft(padded_length);
        assert (ifft.input_size == fft_a.output_size);
        ifft.set_input(vector_elementwise_multiply(fft_a.get_output(), fft_b.get_output()));

        ifft.execute();

        // FFTW returns unnormalized output. To normalize it one must divide each element
        // of the result by the number of elements.
        assert(ifft.output_size == padded_length);
        vector<float> result = ifft.get_output();
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] /= padded_length;
        }

        return result;
    }
%}