% Windowed Frequency Sampling FIR Graphic Equalizer Design
% Designing a Programmable Graphic Equalizer that runs on matlab, using a
% digital filter to perfrom the spectral enhancments
% Fsample=44.1kHz
% Center Frequencies occur at 62.5Hz, 125 Hz, 250 Hz, 500 Hz, 1kHz, 2kHz,
% 4kHz, 8kHz, 16kHz

%FIRST PART
%Uncomment these to get  values for "awful" test case
EQdBsettings=[12,-6,-12,-12,-6,12,-6,12,-6];
Dk_delays_msec=[250 400 520 660 750 1220];
alphak_gains=[0.7 0.6 0.5 0.33 0.2 0.8];
equalize_and_reverb_wavfile('Zarathustra.wav',EQdBsettings,Dk_delays_msec,alphak_gains,'Zarawful.wav');
%stop uncommenting for 'awful' settings
%If magnitude/phase response of FIR plot needed, uncomment the part in the
%function FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum) (Line 217
%through 234)

%SECOND PART
%uncomment these to get 'awesome' settings
%EQdBsettings = [9, 6, 3, 0, 0, 0, 0, 3, 6];
%Dk_delays_msec = [2 100 150 160 170 180];
%alphak_gains=[0.5 0.35 0.33 0.3 0.2 0.1];
%equalize_and_reverb_wavfile('Zarathustra.wav', EQdBsettings, Dk_delays_msec, alphak_gains, 'Zarawesome.wav');
%stop uncommented for 'awesome' settings
%If magnitude/phase response of FIR plot needed, uncomment the part in the
%function FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum) (Line 217
%through 234)

clear all

function hn = equalizer(equalizer_dB)
%function hn = equalizer(equalizer_dB)
%Inputs:
%equalizer_dB = equalizer dB settings needed to process
%Outputs
%hn = output impulse response with equalizer settings

%Select the frequency sample spacing, given to be 44.1kHz
Fs=44100;
M=707; %44100/62.5=705.6 and round up to nearest odd so use 707

equalizer_db_values=equalizer_dB;
%Given frequency settings for equalizer
equalizer_freqs = [62.5 125 250 500 1000 2000 4000 8000 16000];
%Proper spacing for frequency sampling goes from 0 to 1 with no samples at
%F=1
equalizer_spacing = 0:1/M:(M-1)/M;
%Convert into Digital Frequency
i=[];
F_equalizer = equalizer_freqs/Fs;

%Determine which exact frequency sample value indexes are closest to the
%given band center frequencies
%i will give you the indices that correspond to the closest frequencies
%given the FIR fitler with M=707
%i= 2 3 5 9 17 33 65 129 258
for freq=F_equalizer
    [x,index]=min(abs((freq-equalizer_spacing)));
    i=[i index];
end

HF=zeros(1,(M+1)/2);

%This will create a sampling array with length M all of 1's
sampling_array=ones(1,M);
a=0;

%this will set equalizer values to the corresponding indices in the
%sampling array
for z=i
    a=a+1;
    sampling_array(z)=equalizer_db_values(a);
end

%Set samples before first bin to same value as center of first bin
first_index=i(1)-1;
set_first_bin=equalizer_db_values(1);
sampling_array(1:first_index)=ones(1,first_index)*set_first_bin;

%Set samples after last bin to same values as center of last bin
%Only take into account F=0 to 0.5
last_index=i(9)+1;
%i(9)=258 + 1 = 259 so start setting values at this point
set_last_bin=equalizer_db_values(9);
at_center=(M+1)/2;
sampling_array(last_index:(at_center))=ones(1,at_center-i(9))*set_last_bin;
%This will set all values at index 259 to the last equalizer setting given

% Linear Interpolation of dB values
% Will set in between values using a linear slope
% y=mx+b
%i= 2 3 5 9 17 33 65 129 258
for find_slope = 2:9
    w=0;
    in=find_slope-1;
    for num=i(in)+1:i(find_slope)-1
        w=w+1;
        var = w*((equalizer_db_values(find_slope)-equalizer_db_values(in)) / (i(find_slope)-i(in)));
        sampling_array(num)=equalizer_db_values(find_slope-1)+ var;
    end
end

%Take into consideration Nyquist for F=0.5 to 1, should be even symmetric
for nyquist=(1:at_center-1)
    sampling_array(at_center+nyquist)=sampling_array(at_center+1-nyquist);
end

%Convert dB into linear magnitudes
%20logx = dB values
%10^dB/20=linear
linear_freq_array=10.^(sampling_array./20);

%FIR Filter by Freq Sampling
hn=FIR_Filter_By_Freq_Sample(linear_freq_array,1);

%Filter with Turkey
hn = hn .* tukeywin(M).';
H_filter = fft(hn, 2^12);

%Plot Frequency plot with Turkey window added
F=0:1/length(H_filter):1-1/length(H_filter);
figure(6)
grid on
plot(F, 20*log10(abs(H_filter)))
ylabel('Magnitude Response in dB')
xlabel('Digital Frequency F')
title('Equalizer Magnitude Response')
%Make axis go from F=0 to 1
axis([0 1 -12 12])
end

%This function will add delays needed
function [echo_filter_hn]=echo_filter(Dk_delays_msec,alphak_gains,Fsample)
%function [echo_filter_hn]=echo_filter(Dk_delays_msec,alphak_gains,Fsample)
%Inputs:
%Dk_delays_msec: delay times in an array
%alphak_gains: gain factors corresponding to each delay time
%Fsample = sampling frequency rate
%Outputs:
%echo_filter_hn = echo unit sample response hn with the input parameters
delays=round(Dk_delays_msec/1000*Fsample);
echo_filter_hn=zeros(1,max(delays));
echo_filter_hn(1)=1;
echo_filter_hn(delays)=alphak_gains;
end

function Processed_wav=equalize_and_reverb_wavfile(inwavfilename, EQdBsettings, Dk_delays_msec, alphak_gains, outwavfilename)
[xn,Fs]=audioread(inwavfilename);
%function Processed_wav=equalize_and_reverb_wavfile(inwavfilename, EQdBsettings, Dk_delays_msec, alphak_gains, outwavfilename)
%[xn,Fs]=audioread(inwavfilename);
%Inputs:
%inwavfilename=name of the input waveform that needs to be processed with
%equalizer
%EQdBsettings = settings of equalizer
%Dk_delays_msec = delays that are assigned for audio processing
%alphak_gains = array of atteniation for each of the delays
%outwavfilename = new output audio that has been
%processed/equalized/reverbed will go here

%If there are two audio channels,Y will be a 2-dimensional array
%One channel's data will be read into Y(:,1), and the 2nd channel data will
%be in Y(:,2))
%Y is of type double and are normalized values between -1 and 1
hn_equalizer=equalizer(EQdBsettings);
hn_delay=echo_filter(Dk_delays_msec, alphak_gains,Fs);
hn_integrated=fftconv(hn_equalizer,hn_delay);

%Apply the echo and Eq to both channels of a stereo (2-channel) song
%This will check the width of the audio file
%If it's 1 then its mono. For stereo, you convolve hn with the left
%Recombine them into a Nx2 Matrix
%Mono - audio convolved with hn
length_xn = size(xn);
if length_xn(2) == 1 %Represents one mono channel
%Convolving FFT xn and hn
yn = fftconv(xn.', hn_integrated).';
yn = yn / max(abs(yn)); %prevent clipping
elseif length_xn(2) == 2 %Represents stereo (2-channel)
one_x = xn(:,1).';
two_x = xn(:,2).';
one_channel = fftconv(one_x, hn_integrated).';
two_channel = fftconv(two_x, hn_integrated).';
yn = [one_channel / max(abs(one_channel)), two_channel/ max(abs(two_channel))];
end
Processed_wav = yn;
%yn is outputted as the processed wave
audiowrite(outwavfilename, yn, Fs);
%Write output wav back to outwavfilename to be viewed later
end

function [hn,H_F,F]=FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum)
%function [hn,HF,F]=FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum)
%Inputs:
% hn - impulse response of filter (same length as HF_mag_samples)
% HF - complex frequency response of filter (estimated H(F) values found by FFT or freqz)
% F – digital frequency values corresponding to the estimated H(F)values
%Outputs:
% HF_mag_samples – H[k] Magnitude response samples for desired filter
% figurenum - Figure # to plot frequency responses

%Create a frequency vector with length M
length_M = length(HF_mag_samples);
n = 0:length_M-1;

%Create different phase values
phase_diff = exp(-j* pi * n * (length_M-1)/length_M);
H_k = HF_mag_samples .* phase_diff;
hn = real(ifft(H_k));
%Calcualte fft of hn function
H_F = fft(hn, 2^12);
F = 0:1/length(H_F):1-1/length(H_F);

%Uncomment to verify magnitude and phase response of FIR filter
%Plot magnitdue response to verify results
%figure(figurenum)
%subplot(211)
%stem(n/length_M, HF_mag_samples)
%hold on
%plot(F, abs(H_F))
%grid on
%title('Frequency Response of FIR Filter')
%xlabel('Digital Frequency F (cycles/sample)')
%ylabel('Magnitude Response')
%hold off
%Plot phase response to verify results
%subplot(212)
%stem(n/length_M, angle(phase_diff)/pi)
%hold on
%plot(F, angle(H_F)/pi)
%grid on
%xlabel('Digital Frequency F (cycles/sample)')
%ylabel('Phase/pi Response')
end

function yn = fftconv(xn, hn)
%function yn=fftconv(xn,hn)
%This function will perform fast linear time-domain convolution of two
%arbitrary, finite-length sequences, using frequency domain spectral
%multiplication and Fast Fourier Transforms (FFT)
%Inputs:
%xn = xn with length M in time domain
%hn=unit sample response with length N in time domain
%Outputs:
%yn= Output sequence in the time domain
%Zero Padding xn and hn
	length_x = length(xn);
	length_h = length(hn);
	total = length_x + length_h;
	m = 2^nextpow2(total);
	padding_x = [xn zeros(1, m-length_x)];
	padding_h = [hn zeros(1, m-length_h)];

%taking the FFT of the zero padded xn and hn
	x_fft = fft(padding_x);
	h_fft = fft(padding_h);
	yfft = x_fft .* h_fft;
    %take inverse FFT to get y back in time domain
	yn = ifft(yfft);
end
