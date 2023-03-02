function h = golay_response(fnamebase)
% h = golay_response(fnamebase)
% 
% fnamebase:      The base filename for the measurement data
% h:              Estimated impulse response
%
% golay_response uses the golay code excitations stored in
% 'golayA.wav' and 'golayB.wav' to estimate the impulse response of
% the measured system, which is returned in h and written to disk.
% 
% RealSimPLE Project
% Edgar Berdahl, 6/10
%
% e.g. golay_response('hpf');
%      refers to measurement data stored in 'hpfRespA.wav' and 'hpfRespB.wav.'
%      For future reference, the linear response term is scaled and
%      written to the file hpfImpResp.wav.



% Load signals
a = wavread('golayA.wav');
b = wavread('golayB.wav');
L = length(a);
[respA, fs] = wavread(sprintf('%sRespA.wav',fnamebase));
[respB, fs] = wavread(sprintf('%sRespB.wav',fnamebase));


% Compute the impulse response h
h = fftfilt(a(L:-1:1),respA) + fftfilt(b(L:-1:1),respB);
h = h / (2*L);


% Remove the delay due to the Golay codes
h = h(length(a):length(h));


% Write a scaled version to disk
wavwrite(h/max(abs(h))*0.9999,fs,sprintf('%sImpResp.wav',fnamebase));



% Plot the impulse response
figure(1)
plot(([1:length(h)])/fs,h)
xlabel('Time [sec]')
ylabel('Impulse Response')
grid on



% Plot the magnitude response.
figure(2)
Fh = fft(h);
N = length(h);
semilogx(linspace(0,fs/2,N/2+1),20*log10(abs(Fh(1:N/2+1))))
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on



% Plot the minimum phase portion of the response
figure(3)
Fhminphase = mps(Fh);    % min phase version
semilogx(linspace(0,fs/2,N/2+1),unwrap(angle(Fhminphase(1:N/2+1))))
xlabel('Frequency [Hz]')
ylabel('Angle [radians]')
grid on