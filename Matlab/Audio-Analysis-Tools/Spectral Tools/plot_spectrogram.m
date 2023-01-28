function [] = plot_spectrogram(x_t, win_size, hop_size, win_type, fs, nfft)

%   Name: Karan Pareek
%   Student ID: kp2218
%   ------------------
%
%   Calculate and plot the spectrogram of a time-domain audio signal.
%
%   Parameters
%   ----------
%   x_t : 1 x T array
%       time domain signal 
%   win_size : int
%       window size (in samples) 
%   hop_size : int
%       hop size (in samples) 
%   win_type : str
%       window type (one of \texttt{`rect', `hamm', `black'}) 
%   fs : int
%     sample rate (samples per second)
%   nfft : int
%     fft length (in samples)
% 
%   Returns
%   ------- 
%   None

%% Input Checking

if win_size < 0
    error('Window Length cannot be less than 0');
elseif ~ischar(win_type)
    error('Invalid Window Type');
elseif hop_size >= win_size
    error('Invalid Overlap Length');
end

%% Buffer Function

% We define a variable 'x_buf' that calculates the buffer function for our
% signal 'x' with a length of 'winLength' and 'overlapLength'.

x_buf = buffer(x_t,win_size,hop_size);
n_buf = size(x_buf,2);

%% Windowing

if strcmp(win_type, 'rect')
    W = window(@rectwin,win_size);
elseif strcmp(win_type, 'hamm')
    W = window(@hamming,win_size);
elseif strcmp(win_type, 'black')
    W = window(@blackman,win_size);
end

% We create a new matrix W that has horizontal length 'n_buf' such that each
% element of 'x_buf' is now multiplied with the window matrix that give us 
% the final windowed signal for conducting the STFT.

W = repmat(W,1,n_buf); % makes multiple repitions of W with length 'n_buf'
x_win = x_buf .* W; % signal input for the FFT function

%% NFFT and Zero Padding

% We will now define certain conditions to determine the correct value of
% the NFFT or fftlength. For any value of fftlength that is less than the
% windowLength, we ignore it as set it as equal values by default. If
% fftlength exceeds the windowLength, we will Zero Pad the signal by an
% amount equal to the difference between the two.

if nfft <= win_size
    nfft = win_size;
elseif nfft > win_size
    padAmount = nfft - win_size;
    padMat = zeros(padAmount,n_buf); % constructing the zero pad matrix
    x_Pad = [x_buf; padMat]; % concatenating the two matrices
    nfft = size(x_Pad,1); % compute new NFFT
end

%% Computing the FFT of the signal

Fx = fft(x_win,nfft);

% To get the desired form of the STFT result, we first take the absolute
% value of the first half of the function, double it in length and divide
% it by the NFFT parameter to get the final value

Fx = abs(Fx(1:nfft/2 + 1,:))*(2/nfft);
y = 20 * log10(Fx); % Converting the output to the dB scale 

%% Input Parameters for Spectrogram

% In order to plot the spectrogram of the input signal using 'imagesc', we
% need to first define the x and y axis. These are expressed as the Time
% and Frequency Resolution respectively.

freqVect = 0:(fs)/nfft:(fs)/2;
timeVect = 0:1/fs:(length(y)-1)/(fs);

%% Plotting the Spectrogram

imagesc(timeVect,freqVect,y);
c = colorbar;
title('Spectrogram - Using the FFT Function');
set(gca,'YDir','normal') % Reverses the axes
c.Label.String = 'Power/Frequency (dB/Hz)'; % Naming the colorbar
xlabel('Time (secs)');
ylabel('Frequency (kHz)');
axis tight;
    
end