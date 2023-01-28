function [n_t_sf, t_sf, fs_sf] = compute_novelty_sf(x_t, t, fs, win_size, hop_size)

%   Name: Karan Pareek
%   Student ID: kp2218
%   ------------------
%
%   Compute spectral flux novelty function.
%
%   Parameters
%   ----------
%   x_t : 1 x T array
%       time domain signal
%   t : 1 x T array
%       time points in seconds
%   fs : int
%       sample rate of x t (samples per second)
%   win_size : int
%       window size (in samples)
%   hop_size : int
%       hop size (in samples)
%
%   Returns
%   -------
%   n_t_sf : 1 x L array
%       novelty function
%   t_sf : 1 x L array
%       time points in seconds
%   fs_sf : float
%       sample rate of novelty function

%% Input Checking

if win_size < 0
    error('Window Length cannot be less than 0');
elseif hop_size >= win_size
    error('Invalid Overlap Length');
end

%% Buffer Function

% Since the computation of the spectral flux occurs in the frequency
% domain, the function will first compute the STFT of the input audio
% signal based on the window and hop length and then proceed with the
% spectral flux calculation.

% We define a variable 'x_buf' that calculates the buffer function for our
% signal 'x' with a length of 'winLength' and 'overlapLength'.

x_buf = buffer(x_t,win_size,hop_size);
n_buf = size(x_buf,2);

%% Windowing

W = window(@hamming, win_size);

% We create a new matrix W that has horizontal length 'n_buf' such that each
% element of 'x_buf' is now multiplied with the window matrix that give us 
% the final windowed signal for conducting the STFT.

W = repmat(W,1,n_buf);
x_win = x_buf .* W; 

%% Computing the FFT of the signal

nfft = win_size;
Fx = fft(x_win,nfft);

% To get the desired form of the STFT result, we first take the absolute
% value of the first half of the function, double it in length and divide
% it by the NFFT parameter to get the final value

Fx = abs(Fx(1:nfft/2 + 1,:));

%% Spectral Flux

% Now that the STFT of the input file has been calculated, the function
% will proceed with the calculation of the spectral flux.

% Calculating the difference between successive columns
specDiff = diff(Fx,1,2);

% Summing all elements from 1 to N/2
specDiff = specDiff((1:win_size/2+1),:);
n_t_sf = sum(specDiff);

% Performing Half wave rectificatio such that the H(x+|x|) = max(0,x)
for n = 1:length(n_t_sf)
    if n_t_sf(n) < 0
        n_t_sf(n) = 0;
    end
end

% Deriving the value of the Spectral Flux (Half Wave Rectification)
n_t_sf = 2/win_size * n_t_sf;

% Defining the time array and the sample rate of the novelty function
t_sf = t(1):(hop_size/fs):(length(n_t_sf)-1)/(fs/hop_size);
fs_sf = fs/hop_size;

end