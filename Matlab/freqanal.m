% simplpnfa.m - matlab program for frequency analysis
%               of the simplest lowpass filter:
%
%                    y(n) = x(n)+x(n-1)}
%
%               the way people do it in practice.

B = [1,1]; % filter feedforward coefficients
A = 1;     % filter feedback coefficients

N=128;     % FFT size = number of COMPLEX sinusoids
fs = 1;    % sampling rate in Hz (arbitrary)

Bzp = [B, zeros(1,N-length(B))]; % zero-pad for the FFT

H = fft(Bzp);   % length(Bzp) should be a power of 2

if length(A)>1  % we're not using this here
  Azp = [A,zeros(1,N-length(A))]; % but show it anyway.
  % [Should guard against fft(Azp)==0 for some index]
  H = H ./ fft(A,N); % denominator from feedback coeffs
end

% Discard the frequency-response samples at
% negative frequencies, but keep the samples at
% dc and fs/2:

nnfi = (1:N/2+1);     % nonnegative-frequency indices
Hnnf = H(nnfi);       % lose negative-frequency samples
nnfb = nnfi-1;        % corresponding bin numbers
f = nnfb*fs/N;        % frequency axis in Hz
gains = abs(Hnnf);    % amplitude response
phases = angle(Hnnf); % phase response

plotfile = 'simplpnfa.eps';
swanalmainplot;    % final plots and comparison to theory
