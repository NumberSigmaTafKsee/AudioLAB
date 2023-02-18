#include <string>
#include <cstdio>

class FIR
{
public:
	FIR(double *coefficients, unsigned number_of_taps);
	~FIR();

	double filter(double input);

	void reset();

	unsigned getTaps() {return taps;};
	
private:
	double        *coefficients;
	double        *buffer;
	unsigned      taps, offset;

};

FIR::FIR(double *coefficients, unsigned number_of_taps) :
   coefficients(coefficients),
   buffer(new double[number_of_taps]()),  
   taps(number_of_taps),
   offset(0)
{}

FIR::~FIR()
{
  delete[] buffer;
  //delete[] coefficients;
}

double FIR::filter(double input)
{
   double *coeff     = coefficients;
   double *coeff_end = coefficients + taps;

   double *buf_val = buffer + offset;

   *buf_val = input;
   double output_ = 0;
	
   while(buf_val >= buffer){
      output_ += *buf_val-- * *coeff++;
   }
	
   buf_val = buffer + taps-1;
	
   while(coeff < coeff_end){
      output_ += *buf_val-- * *coeff++;
   }
	
   if(++offset >= taps){
      offset = 0;
   }
	
   return output_;
}

void FIR::reset()
{
   ::memset(buffer, 0, sizeof(double)*taps);
   offset = 0;
}

/* matlab
%N    = 1;      % Order


Fc1  = 0.35;      % First Cutoff Frequency
Fc2  = 10;      % Second Cutoff Frequency
Fs = 360; % Sampling Rate
flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter
for N = 1:250
% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
dlmwrite(strcat(num2str(N),'.txt'), b);
end

%High Pass Filter%
%{
Fc = 0.9;      % Cutoff Frequency
Fs = 360; % Sampling Rate
flag = 'scale';  % Sampling Flag
Beta = 2;      % Window Parameter
for N = 2:2:250
% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'high', win, flag);
dlmwrite(strcat(num2str(N),'.txt'), b);
end
%}

%Low Pass Filter%
%{
Fc = 2;      % Cutoff Frequency
Fs = 360; % Sampling Rate
flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter
for N = 2:2:250
% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
dlmwrite(strcat(num2str(N),'.txt'), b);
end
%}
*/