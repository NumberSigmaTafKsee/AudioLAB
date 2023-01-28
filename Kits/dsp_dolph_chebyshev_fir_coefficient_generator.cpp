#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

#define PI 3.1415926535
void Low_Pass(double Fs, double Crossover_point1,int Filter_Length, double *h)
{
	int M = (Filter_Length - 1) / 2;
	double omiga_c = Crossover_point1 * 2.0 * PI / Fs;
	for (int n = 0; n < Filter_Length; n++)
	{
		if (n == M)
		{
			h[n] = omiga_c / PI;
		}
		else
		{
			h[n] = sin(omiga_c*(n-M))/(PI*(n-M));
		}
	}
	return;
}
void High_Pass(double Fs, double Crossover_point1, int Filter_Length, double *h)
{
	int M = (Filter_Length - 1) / 2;
	double omiga_c = Crossover_point1 * 2 * PI / Fs;
	for (int n = 0; n < Filter_Length; n++)
	{
		if (n == M)
		{h[n] = 1.0 - omiga_c / PI;}
		else
		{
			h[n] = -1.0*sin(omiga_c*(n - M)) / (PI*(n - M));
		}
	}
	return;
}
void Band_Pass(double Fs, double Crossover_point1, double Crossover_point2, int Filter_Length, double *h)
{
	int M = (Filter_Length - 1) / 2;
	double omiga_c1 = Crossover_point1 * 2.0 * PI / Fs;
	double omiga_c2 = Crossover_point2 * 2.0 * PI / Fs;
	for (int n = 0; n < Filter_Length; n++)
	{
		if (n == M)
		{
			h[n] = (omiga_c2 - omiga_c1) / PI;
		}
		else
		{
			h[n] = (sin(omiga_c2*(n - M)) / (PI*(n - M))) - (sin(omiga_c1*(n - M)) / (PI*(n - M)));
		}
	}
	return;
}
void Band_Stop(double Fs, double Crossover_point1, double Crossover_point2, int Filter_Length, double *h)
{
	int M = (Filter_Length - 1) / 2;
	double omiga_c1 = Crossover_point1 * 2.0 * PI / Fs;
	double omiga_c2 = Crossover_point2 * 2.0 * PI / Fs;
	for (int n = 0; n < Filter_Length; n++)
	{
		if (n == M)
		{
			h[n] = 1.0 - (omiga_c2- omiga_c1) / PI;
		}
		else
		{
			h[n] = (sin(omiga_c1*(n - M)) / (PI*(n - M))) - (sin(omiga_c2*(n - M)) / (PI*(n - M)));
		}
	}
	return;
}
double Chebyshev_Polynomials(int n, double x)
{
	if (x <= -1.0)
	{
		return pow(1.0, n)*cosh(n*log(sqrt(x*x - 1) - x));
	}
	else if (x >= 1.0)
	{
		return cosh(n*log(x + sqrt(x*x - 1)));
	}
	else
	{
		return cos(n*acos(x));
	}
}
void Get_Filter_Window_Using_Dolph_Chebyshev(double gama, int Filter_Length, double *Filter_Window)
{
	int M = (Filter_Length - 1) / 2;
	double beta = cosh((0.5/M)*log(gama + sqrt(gama*gama - 1)));
	for (int n = 0; n < Filter_Length; n++)
	{
		Filter_Window[n] = 0.0;
		for (int k = 1; k <= M; k++)
		{
			Filter_Window[n] += (Chebyshev_Polynomials(2 * M, beta*cos(k*PI / Filter_Length))*cos(2.0*(n - M)*k*PI / Filter_Length));
		}
		Filter_Window[n] = (gama + Filter_Window[n] * 2.0) / Filter_Length;
	}
	return;
}
void Coefficient_Generator(double Fs, double Crossover_point1, double Crossover_point2, int Filter_Type, int Filter_Length, int *Quantized_Coefficient, double *Original_Coefficient, int &Error_Code)
{
	double *h = NULL,*Filter_Window=NULL;	//Ideal impulse response And Window function
	double gama = 60.0;		//Main lobe amplitude/Side flap amplitude Unit:dB
	int max = (Filter_Length - 1) / 2;
	Error_Code = 65535;

	if ((Fs <= Crossover_point1 * 2) || (Fs <= Crossover_point2 * 2) || (Fs <= 0.0) || (Crossover_point1 <= 0) || (Crossover_point2 <= 0))
	{
		Error_Code = 1;		//Frequency set error.
		return;
	}
	if ((Crossover_point2 <= Crossover_point1) && ((Filter_Type == 2) || (Filter_Type == 3)))
	{
		double temp = Crossover_point2;
		Crossover_point2 = Crossover_point1;
		Crossover_point1 = temp;
	}
	
	if ((0 == Filter_Length % 2) || (Filter_Length < 3))
	{
		Error_Code = 2;		//Filter_Length error.
		return;
	}

	h = (double *)malloc(Filter_Length * sizeof(double));
	if (NULL == h)
	{
		Error_Code = 3;		//No memory available for vector h.
		return;
	}

	Filter_Window = (double *)malloc(Filter_Length * sizeof(double));
	if (NULL == Filter_Window)
	{
		Error_Code = 4;		//No memory available for vector window.
		free(h);
		return;
	}

	switch (Filter_Type)				//Get ideal filter impulse response
	{
	case 0:Low_Pass(Fs, Crossover_point1, Filter_Length, h); break;
	case 1:High_Pass(Fs, Crossover_point1, Filter_Length, h); break;
	case 2:Band_Pass(Fs, Crossover_point1, Crossover_point2, Filter_Length, h); break;
	case 3:Band_Stop(Fs, Crossover_point1, Crossover_point2, Filter_Length, h); break;
	default:
		Error_Code = 5;		//Undefined filter type.
		free(h);
		free(Filter_Window);
		return;
		break;
	}

	Get_Filter_Window_Using_Dolph_Chebyshev(gama, Filter_Length, Filter_Window);	//Get Window function

	Original_Coefficient[max]= h[max] * Filter_Window[max];
	for (int i = 0; i < Filter_Length; i++)			//Get final filter coefficient
	{
		double a, b;
		Original_Coefficient[i] = h[i] * Filter_Window[i];
		a=modf((Original_Coefficient[i]/ Original_Coefficient[max])* 8388607.0, &b);	//24bit,signed,0x7fffff
		Quantized_Coefficient[i] = (int)(b);
	}
	free(h);
	free(Filter_Window);
	Error_Code = 0;		//No error occurs
	return;
}

int main()
{
	ofstream file;
	file.open("out.txt");
	int Error_Code, Quantized_Coefficient[1023] = {0};
	double Original_Coefficient[1023] = { 0.0 };

	double Fs = 48000;
	double Crossover_point1 = 9600, Crossover_point2 = 23999;
	int Filter_Type = 0;
	int Filter_Length = 1023;
	Coefficient_Generator(Fs, Crossover_point1, Crossover_point2, Filter_Type, Filter_Length, Quantized_Coefficient, Original_Coefficient, Error_Code);

	cout << Error_Code;
	system("pause");
	file << '[';
	for (int i = 0; i < 1022; i++)
	{
		file << Quantized_Coefficient[i] << ',';
	}
	file << Quantized_Coefficient[1022] << ']';
	file << '\n' << '\n';
	file << '[';
	for (int i = 0; i < 1022; i++)
	{
		file << Original_Coefficient[i] << ',';
	}
	file << Original_Coefficient[1022] << ']';
	file.close();
	return 0;
}