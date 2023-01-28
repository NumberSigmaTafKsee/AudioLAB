function [n_t_le, t_le, fs_le] = compute_novelty_ef(x_t, t, fs, win_size, hop_size)

%   Name: Karan Pareek
%   Student ID: kp2218
%   ------------------
%
%   Compute envelope follower novelty function.
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
%   n_t_le : 1 x L array
%       novelty function
%   t_le : 1 x L array
%       time points in seconds
%   fs_le : float
%       sample rate of novelty function

%% Buffer Function

% We define a variable 'x_buf' that calculates the buffer function for our
% signal 'x' with a length of 'winLength' and 'overlapLength'.

x_buf = buffer(x_t,win_size,hop_size);
n_buf = size(x_buf,2);

%% Windowing

% We create a new matrix W that has horizontal length 'n_buf' such that each
% element of 'x_buf' is now multiplied with the window matrix that give us 
% the final windowed signal for conducting the STFT.

W = window(@hamming, win_size); 
W = repmat(W,1,n_buf);

%% Log Energy Derivative

% Now that the basic premise has been established, the function will now
% square the elements of each section of the envelope follower, multiply by
% the window function and sum each section.

x_buf = abs(x_buf); % |x(n + mh)|
x_win = x_buf .* W; % |x(n + mh)| * w(n)

col = size(x_win,2); % No. of columns in the buffer
blank = zeros(1,col); % Creating a blank matrix

% Summing all elements of each column of x_win
for c = 1:col
    blank(:,c) = blank(:,c) + sum(x_win(:,c));
end

% Deriving the local energy
n_t_le = 1/win_size * blank;
n_t_le = n_t_le/max(abs(n_t_le)); % Maximizing the novelty function

% Defining the time array and the sample rate of the novelty function
t_le = t(1):(hop_size/fs):(length(n_t_le)-1)/(fs/hop_size);
fs_le = fs/hop_size;

end