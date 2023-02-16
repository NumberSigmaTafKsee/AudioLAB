using namespace std;
#include "Engine.h"
#include <iostream>
#include <complex>
#include <vector>

	Engine::Engine(int _size)
	{	
		
		size = _size;
		matrix=vector<vector<complex<double>>>(size);
		for (int i = 0; i < size; i++)
		{
			matrix[i] = vector<complex<double>>(size);
			for (int j = 0; j < size; j++)
			{
				matrix[i][j] = polar(1.0, -1* pi2 * i * j/size);
			}
		}
		cout << "engine has been created" << endl;
	}
	vector<complex<double>> Engine:: even_odd_decomposer(vector<complex<double>> input)
	{
		vector<vector<complex<double>>> output = vector<vector<complex<double>>>(2);
		output[0] = vector<complex<double>>(input.size()>>1,0);
		output[1]= vector<complex<double>>(input.size()>>1, 0);
		for (int i = 0; i < (input.size()); i<<1)
		{
			output[0][i] = input[i];
			output[1][i] = input[i + 1];
		}
		if (input.size() > 2)
		{
			output[0] = even_odd_decomposer(output[0]);
			output[1] = even_odd_decomposer(output[1]);
		}
		vector<complex<double>> output2 = vector<complex<double>>(input.size());
		for (int i = 0; i < output[1].size(); i++)
		{
			output2[i] = output[0][i];
			output2[i + output2.size()>>1] = output[1][i];
				
		}
		return output2;
		

	}
	vector<complex<double>> Engine::fixed_fft(vector<complex<double>> input)
	{
		//for length of 8192 = 2^13
		

		vector <complex<double>> output1, output2;

		output2 = even_odd_decomposer(input);		//zamiana 001 100
		int step = 2;
		while (step <= 8192)
		{

			output1 = vector<complex<double>>(step, 0);
			for (int i = 0; i < 8192; i = i + step)	//butterfly step
			{
				for (int j = 0; j < step; j++)
				{
					output1[j] = output2[i + j];
				}
				output1 = sub_fft(output1);
				for (int j = 0; j < step; j++)
				{
					output2[i + j] = output1[j];
				}

			}
			step = step <<1;
		}
		return output2;
	}
	vector<complex<double>> Engine::fft(vector<complex<double>> input)
	{

		int a = input.size();
		int new_size = 1;
		int steps = 0;
		while (a > new_size)
		{
			steps++;
			new_size = new_size * 2;
		}
		
		vector <complex<double>> output1, output2;
		output1 = window_fuction(input);
		//output1 = input;
		output2 = vector<complex<double>>(new_size, 0);	// zero padding
		for (int i = 0; i < input.size(); i++)
		{
			output2[i] = output1[i];
		}
		

		output2 = even_odd_decomposer(output2);		//zamiana 001 100
		int step = 2;	
		while(step <=output2.size())
		{

			output1 = vector<complex<double>>(step, 0);
			for (int i = 0; i < input.size(); i = i + step)	//butterfly step
			{					
				for (int j = 0; j < step; j++)	
				{
					output1[j] = output2[i + j];
				}
				output1 = sub_fft(output1);			
				for (int j = 0; j < output1.size(); j++)
				{
					output2[i + j] = output1[j];
				}

			}
			step = step * 2;
		}
		return output2;
	}
	vector<complex<double>> Engine::ifft(vector<complex<double>> input)
	{
		int a = input.size();
		//tu ma byc uzupelnienie do 2^n
		int b = 0;
		while (a % 2 == 0)
		{
			b++;
			a = a / 2;

		}
		vector <complex<double>> output1, output2;
		output2 = vector<complex<double>>(input.size(), 0);
		output2 = even_odd_decomposer(input);		//zamiana 001 100

		int step = a;
		while (step < input.size())
		{
			output1 = vector<complex<double>>(step, 0);
			for (int i = 0; i < input.size() / step; i = i + step)
			{
				for (int j = 0; j < step; j++)
				{
					output1[j] = output2[i + j];
				}
				output1 = sub_fft(output1);			//butterfly krok
				for (int j = 0; j < output1.size(); j++)
				{
					output2[i + j] = output1[j];
				}
			}
			step = step * 2;
		}
		return output2;
	}
	//do zmiany z fft
	vector<complex<double>> Engine::sub_fft(vector<complex<double>> input)
	{
		int l = input.size() >>1;
		vector<complex<double>> output = vector<complex<double>>(input.size(), 0);
		for (int i = 0; i < l; i++)
		{
			output[i] = input[i] + input[i + l]* polar(1.0, -pi2 * i /(2*l));
			output[i + l] = input[i]+input[i+l] * polar(1.0, (-pi2 * (i+l)) /(2*l));
		}
		return output;
	}
	
	vector<complex<double>> Engine::dft(vector<complex<double>> input)
	{
		vector<complex<double>> complexoutput = vector<complex<double>>(input.size(),0);		
		for (int i = 0; i < size; i++)
		{
			
			for (int j = 0; j < size; j++)
			{					
				complexoutput[i] =  complexoutput[i]+ matrix[i][j] * input[j];// ij ji doesn't matter
			}
			
		}
		
		return complexoutput;

	}	  //coœ jest Ÿle
	void Engine::show_matrix_d()
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				cout << matrix[i][j] << "\t\t";
			}
			cout << endl << endl;
		}
	}
	vector<double> Engine::convolution(vector<double> signal, vector<double> filter)
	{
		int size_filter = filter.size();
		int size_signal = signal.size();
		int shorter = 0;
		if (size_signal > size_filter)
		{
			shorter = size_filter;
			filter.resize(size_signal,0);
		}
		else
		{
			shorter = size_signal;
			signal.resize(size_signal, 0);
		}
		vector<double> output = vector<double>(signal.size(), 0);
		
		for (int i = 0; i < output.size(); i++)
		{
			for (int j = 0; j < i+1; j++)
			{
				output[i] = output[i] + filter[j] * signal[i-j];
			}
			for (int j = i+1; j < output.size(); j++)
			{
				output[i] = output[i] + filter[j] * signal[output.size()+i-j];
			}
			
		}
		return output;
	}
	vector<double> Engine::circular_convolution(vector<double> signal, vector<double> filter)
	{
		int size_filter = filter.size();
		int size_signal = signal.size();
		int shorter = 0;
		if (size_signal > size_filter)
		{
			shorter = size_filter;
		}
		else
		{
			shorter = size_signal;
		}
		int counter = 0;
		vector<double> output = vector<double>(size_signal + size_filter - 1, 0);

		for (int i = 0; i < output.size(); i++)
		{

			if (i < shorter) { counter++; }
			if (i > output.size() - shorter - 1) { counter--; }
			for (int j = 0; j < counter; j++)
			{
				output[i] = output[i] + filter[j] * signal[i - j];
				cout << output[i] << endl;
			}
			cout << endl;
		}
		return output;
	}
	//nieu¿ywane

	vector<double> Engine::energy_in_bands_log_scale(vector<complex<double>> signal_f, int bands_per_octave)
	{
		//powinno byc razy dwa bo tylko pol widma ale i tak bedzie skalowane potem jakos

		vector<double> freq = band_freq_generator(10,bands_per_octave);
	
		vector<double> output = vector<double>(freq.size(),0);
		double step = 44100 / signal_f.size();
		int k = 0;
		for (int i = 0; i < freq.size();i++)
		{
			k=freq[i] / step;
			freq[i] = k; //znalezienie najblizszego k w widmie (indeksu odpowiadajacego f);
			cout << k << "\t";
		}
		cout << endl;
		for (int j = 0; j < freq.size()-1; j++)
		{
			for (int i = freq[j]; i < freq[j+1]; i++)
			{
				output[j] = output[j] + abs(signal_f[i]) * abs(signal_f[i]);
			}

			output[j] =log(output[j]/(freq[j+1]-freq[j]));

		}
		return output;
	}
	vector<complex<double>> Engine::window_fuction(vector<complex<double>> input)
	{
		vector<complex<double>> output = vector<complex<double>>(input.size(), 0);
		for (int i = 0; i < input.size(); i++)
		{
			output[i] = input[i] * (0.54 - 0.46 * cos(pi2 * i / input.size()));
		}
		return output;
	}
	vector<double> Engine::band_freq_generator(int octaves, int bands_per_octave)
	{
		vector<double> output = vector<double>(octaves*bands_per_octave+1, 261.63/16);
		double a = pow(2, 1.0/bands_per_octave);
		for (int i = 0; i < octaves* bands_per_octave; i++)
		{
			output[i + 1] = output[i] * a;
		}
		return output;
	}
