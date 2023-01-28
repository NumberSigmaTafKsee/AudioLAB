function [mfccs, fs_mfcc] = compute_mfccs(filepath, win_size, hop_size, ... 
                                min_freq, max_freq, num_mel_filts, n_dct)
%   Name: Karan Pareek
%   Student ID: kp2218
%
%   Compute MFCCs from audio file.
%   
%   Parameters
%    ----------
%    filepath : string
%        path to .wav file 
%    win size : int
%     spectrogram window size (samples) 
%    hop size : int
%        spectrogram hop size (samples) 
%    min freq : float
%        minimum frequency in Mel filterbank (Hz) 
%    max freq : float
%        maximum frequency in Mel filterbank (Hz) 
%    num mel filts: int
%        number of Mel filters 
%    n dct: int
%        number of DCT coefficients
%        
%    Returns
%     -------
%    mfccs : n dct-1 x NT array
%        MFCC matrix (NT is number spectrogram frames)
%    fs mfcc : int
%        sample rate of MFCC matrix (samples/sec)                   

%% Initialization

% Reading an audio file                           
[x_t,fs] = audioread(filepath);

% Defining the parameters for the spectrogram
noverlap = win_size - hop_size;
nfft = win_size;

%% Buffer Function

% We define a variable 'x_buf' that calculates the buffer function for our
% signal 'x' with a length of 'winLength' and 'overlapLength'.

x_buf = buffer(x_t,win_size,noverlap);
n_buf = size(x_buf,2);

%% Windowing

W = window(@hamming,win_size);

% We create a new matrix W that has horizontal length 'n_buf' such that each
% element of 'x_buf' is now multiplied with the window matrix that give us 
% the final windowed signal for conducting the STFT.

W = repmat(W,1,n_buf); % makes multiple repitions of W with length 'n_buf'
x_win = x_buf .* W; % signal input for the FFT function

%% Computing the FFT of the signal

Fx = fft(x_win,nfft);

% To get the desired form of the STFT result, we first take the absolute
% value of the first half of the function and divide it by the NFFT 
% parameter to get the final value

Fx = abs(Fx(1:nfft/2 + 1,:))/nfft;
fs_mfcc = fs/hop_size; % Sample rate of MFCC

%% Mel Filterbank

% Given the min and max frequency values, we can derive triangular shaped
% mel filters to calculate the MFCC of the input signal

mel_min = hz2mel(min_freq);
mel_max = hz2mel(max_freq);

% Given the number of filters and the fact that the filters will be placed
% quidistant to each other in the Mel scale, we can derive the range of Mel
% frequency that we want to work with.

mel_filters = linspace(mel_min,mel_max,(num_mel_filts+2));
    
% Converting back to the frequency scale
freq_scale = mel2hz(mel_filters);

%% Filter Construction

% Creating a [M x (N/2+1)] matrix
spectrum_len = size(Fx,1);

% Blank matrix for the filterbank
filterbank = zeros(num_mel_filts,spectrum_len);  

% converting frequency to nearest bin
filts = floor((nfft+1)*freq_scale/fs);

% Here, the Mel Filterbank is constructed using the index values derived
for m = 2:(num_mel_filts+1)
    
    f_left = filts(m-1);
    f_center = filts(m);
    f_right = filts(m+1);
    
    for k = f_left:f_center
        filterbank(m-1,k) =+ (k-filts(m-1))/(filts(m)-filts(m-1));
    end
    
    for k = f_center:f_right
        filterbank(m-1,k) =+ (filts(m+1)-k)/(filts(m+1)-filts(m));
    end
   
end

% normalize each Mel filter to sum to 1
for l = 1:size(filterbank,1)
    filterbank(l,:) = filterbank(l,:)/sum(filterbank(l,:));
end

% The Mel spectrum can now be derived by multiplying the Filter matrix with
% the DFT of the signal
mel_spectrum = filterbank * Fx;

%% Mel Power Spectrum

% Doing a dB conversion
mel_power_spectrum = 20*log10(mel_spectrum);

%% DCT

% Calculating the DCT of n_dct elements of each column
cols = size(mel_power_spectrum,2);
dct_matrix = zeros(n_dct,size(mel_power_spectrum,2));

% Calculating the DCT of the mel_spectrum matrix
for n = 1:cols
    dct_matrix(:,n) =+ dct(mel_power_spectrum(1:n_dct,n));
end

% Removing the first DC component
mfccs = dct_matrix(2:end,:); % Removing the first DC coeffiecient

% Normalizing the MFCC matrix using global maximization
delta = (mfccs - min(mfccs(:)));
mfccs = delta ./ max(delta(:));

end
