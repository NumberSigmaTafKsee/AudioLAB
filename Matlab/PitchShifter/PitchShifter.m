clc
clf
clear all

% Read an audio file
soundFile = 'F:\Fall 2020\ECE 484\rooster.mp3';
[input1, Fs] = audioread(soundFile);

if ~contains(soundFile, 'tone')
    if contains(soundFile, '.mp3')
        outputSoundFile = strrep(soundFile, '.mp3', '-shifted.wav');
    elseif contains(soundFile, '.wav')
        outputSoundFile = strrep(soundFile, '.wav', '-shifted.wav');
    end
end

numChannels = size(input1,2);

% Set up parameters
hopSizeA = 32;
frameSize = 2048;
figure(1)
plot(abs(fft(input1, Fs)))
title("Frequency Magnitude Plot for 440 Hz Input");
xlabel("Frequency (Hz - 1)");
ylabel("Magnitude")

Ratio = 1.06;
[num, dem] = rat(Ratio);
hopSizeS = int64(Ratio * hopSizeA);

numHops = ceil(size(input1,1)/hopSizeA);
magsA = zeros(frameSize);
phasesA = zeros(frameSize);

magsS = zeros(frameSize);
phasesS = zeros(frameSize);
framesS = zeros(frameSize);

% Pad the input array
padArray = zeros(frameSize, numChannels);
input = [input1;padArray];
output = zeros(3*size(input, 1), numChannels);

% Hann window
w = 0.5*(1-cos(2*pi*(0:frameSize-1)'/(frameSize)));
wS = w;

% Stacking windows with hop ratio < 0.25 will become large and start
% clipping.
% Pre-calculate this amount so we can scale the window accordingly.
testWindowA = zeros(2 * frameSize, 1);
numHopPerFrameA = frameSize/hopSizeA;
for hop = 1:numHopPerFrameA
    testStartA = (hop - 1) * hopSizeA + 1;
    testEndA = testStartA + frameSize - 1;
    testWindowA(testStartA:testEndA) = testWindowA(testStartA:testEndA) + w;
end
maxA = max(testWindowA);
w = w/maxA;

testWindowS = zeros(2 * frameSize, 1);
numHopPerFrameS = frameSize/hopSizeS;
for hop = 1:numHopPerFrameS
    testStartS = (hop - 1) * hopSizeS + 1;
    testEndS = testStartS + frameSize - 1;
    testWindowS(testStartS:testEndS) = testWindowS(testStartS:testEndS) + wS;
end
maxS = max(testWindowS);
wS = wS/maxS;


for channel = 1:numChannels
    for hopA = 1:numHops
        
        % Perform STFT Analysis step
        frameStartA = (hopA - 1) * hopSizeA + 1;
        frameEndA = frameStartA + frameSize - 1;
        window = w .* input(frameStartA:frameEndA, channel);
        window = circshift(window, frameSize/2);
        framefft = fft(window);
        magsA = abs(framefft).';
        phasesA = atan2(imag(framefft), real(framefft)).';
        
        %Do Phase stuff
        omegakh = 2 * pi * (0:frameSize-1) * hopSizeA / frameSize;
        if hopA == 1
          magsS = magsA;
          phasesS = phasesA;
        else
          magsS = magsA;
          phaseDiff = phasesA - prevPhasesA;
          test1 = phaseDiff - omegakh;
          test1 = test1 - round(test1/(2*pi))*2*pi;
          test1 = (test1 + omegakh) * Ratio;
          phasesS = test1 + prevPhasesS;
        end
        prevPhasesA = phasesA;
        prevPhasesS = phasesS;
      
        framesS = magsS .* cos(phasesS) + 1i * magsS .* sin(phasesS);

        % Synthesis step
        frameStartS = (hopA - 1) * hopSizeS + 1;
        frameEndS = frameStartS + frameSize - 1;
        frameSifft = real(ifft(framesS));
        windowS = circshift(frameSifft, frameSize/2);
                
        output(frameStartS:frameEndS, channel) = output(frameStartS:frameEndS, channel) + windowS.';
    end
end
output1 = resample(output, dem, num);
if ~contains(soundFile, 'tone')
    audiowrite(outputSoundFile, output1, Fs)
end
figure(2)
plot(abs(fft(output1, Fs)))
title("Frequency Magnitude Plot for Hop Size = 32, Frame Size = 2048, Semitone Up");
xlabel("Frequency (Hz - 1)");
ylabel("Magnitude")
sound(output1, Fs)

%memes = input1;
%memes(:,1) = output1(1:size(input1,1), 2);
%sound(memes, Fs)
