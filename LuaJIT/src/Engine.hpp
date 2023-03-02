#pragma once
using namespace std;
#include <iostream>
#include <complex>
#include <vector>

class Engine
{
public:
	Engine(int _size);
	void show_matrix_d();
public:
	vector<complex<double>> fft(vector<complex<double>> input);	
	vector<complex<double>> fixed_fft(vector<complex<double>> input);
	vector<complex<double>> ifft(vector<complex<double>> input);// jeszcze nie ma

private:
	vector<complex<double>> even_odd_decomposer(vector<complex<double>> input);
	vector<complex<double>> sub_fft(vector<complex<double>> input);
public:
	vector<complex<double>> dft(vector<complex<double>> input);
	vector<double> convolution(vector<double> signal, vector<double> filter);
	vector<double> circular_convolution(vector<double> signal, vector<double> filter);


	vector<complex<double>> window_fuction(vector<complex<double>> input);
	
	vector<double> energy_in_bands_log_scale(vector<complex<double>> signal, int bands_per_octave);
	vector<double> band_freq_generator(int octaves, int bands_per_octave);
private:
	const double pi2 = 4 * acos(0.0);
	int size;
	vector<vector<complex<double>>> matrix;
};
