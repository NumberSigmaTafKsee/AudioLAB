function h = sinesweeps_response(fnamebase, range, imptime)
% h = sinesweeps_response(fnamebase, range, imptime)
% 
% fnamebase:      The base filename for the measurement data
% range:          Limit on the dynamic range of the inverse filter
% imptime:        The length of the impulse response in seconds
%
% sinesweeps_response uses the swept sine excitation stored in
% 'sinesweeps.wav' to calculate an inverse filter using the FFT.
% The inverse filter is applied cyclically by dividing in the 
% frequency domain.  For this reason, the dynamic range of the 
% inverse filter is limited by range to help avoid divide by zero errors.  
% The impulse response time is finally used to separate out the linear
% term of the response from the nonlinear terms.
% 
% RealSimPLE Project
% Edgar Berdahl, 6/10
%
% e.g. sinesweeps_response('nonlinear2', 100, 0.4);
%      refers to measurement data stored in nonlinear2Resp.wav
%      For future reference, the linear response term is scaled and
%      written to the file nonlinear2ImpResp.wav.




% Load signals
sdbl = wavread('sinesweeps.wav');
L = length(sdbl);
scycl = sdbl(L/2+1:L);

[rdbl, fs] = wavread(sprintf('%sResp.wav',fnamebase));
rcycl = rdbl(L/2+1:L);


% Here we limit the dynamic range of the excitation (and equivalently
% the inverse filter) and display the result.
figure(1)
Fs = fft(scycl);
N = length(scycl);
hold off
plot(linspace(0,fs/2,N/2+1),20*log10(abs(Fs(1:N/2+1))))
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on

ind = find(abs(Fs) < max(abs(Fs))/range);
Fs(ind) = ones(size(ind)) * max(abs(Fs))/range;
h = ifft(fft(rcycl)./Fs);
h = real(h);
hold on
plot(linspace(0,fs/2,N/2+1),20*log10(abs(Fs(1:N/2+1))),'r')
legend('Actual Excitation','Dynamic Range-Limited Excitation')
wavwrite(h/max(abs(h))*0.9999,fs,sprintf('%sImpResp.wav',fnamebase));




% Plot the full response.
figure(2)
plot([1:length(h)]/fs,h)
xlabel('Time [sec]')
ylabel('Full-length Response')
grid on




% Plot the log of the absolute value of the response to
% help show the presence of any nonlinear terms.
figure(3)
dBh = log(abs(h));
plot([1:length(h)]/fs,max(dBh, min(dBh)*ones(size(dBh))))
xlabel('Time [sec]')
ylabel('LOG(|Full-length Response|)')
grid on




% Plot the linear term of the response only.
figure(4)
h = h(1:floor(imptime*fs));
plot([1:length(h)]/fs,h)
xlabel('Time [sec]')
ylabel('Linear Term Of Response')
grid on




% Plot the magnitude spectrum of the linear response term.
figure(5)
Fh = fft(h);
N = length(h);
semilogx(linspace(0,fs/2,N/2+1),20*log10(abs(Fh(1:N/2+1))))
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on




% Plot the minimum phase portion of the linear response term.
figure(6)
Fhminphase = mps(Fh);    % min phase version
semilogx(linspace(0,fs/2,N/2+1),unwrap(angle(Fhminphase(1:N/2+1))))
xlabel('Frequency [Hz]')
ylabel('Angle [radians]')
grid on