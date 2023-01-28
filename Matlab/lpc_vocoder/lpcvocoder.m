
clc
clear
close all

[audio,fs] =audioread('file4.wav');
fprintf('Sampling frequency is %d Hz\n',fs);

p=12; % lpc filter order
fr_size=0.03;    %speech signal is divided into frames of 30ms

frame_length=round(fs*fr_size); %sample length of a frame


%% ENCODINDG %%
lpc_coeffs=[];
gain=[];
for frame_no=1:frame_length:(length(audio)-frame_length)
    frame=audio(frame_no:frame_no+frame_length-1).*hamming(331);
    autocorrelation=xcorr(frame);
    
    % pitch calculation
    [pk,ind]=findpeaks(autocorrelation,'MinPeakDistance',65); % this value can be played for different pitch periods
    pitch_period(frame_no:(frame_no+frame_length-1))=1/mean(fs./diff(ind)); % pitch period detection
    
    % voiced/unvoiced detection according to zero crossing rate
    zcr= sum(abs(diff(frame>0)))/length(frame);
    %disp(zcr);
    if(zcr<0.3)% threshold is found by experimenting
        voiced_unvoiced(frame_no:(frame_no+frame_length-1))=1;    % voiced/unvoiced count
    else
        voiced_unvoiced(frame_no:(frame_no+frame_length-1))=0;
    end
    
    
    % levinson-durbin algorithm to extract filter coefficients and gain
    a_matrix= zeros(p+1,p+1);% initializing coefficient matrix
    R=autocorrelation(((length(autocorrelation)+1)/2):end); % taking only nonnegative side of the autocorrelation vector
    
    %LEVINSON - DURBIN ALGORITHM starts
    JJ(1)=R(1); % order 1
    for l=2:p+1
        %step1
        temp_sum=0;
        for i=2:(l-1)
            temp_sum=temp_sum+a_matrix(i,(l-1)).*R(l-i+1);
        end
        k_l=(R(l)+temp_sum)./JJ(l-1);
        
        %step2
        a_matrix(l,l)= -k_l;
        a_matrix(1,l)=1;
        for i=2:(l-1)
            a_matrix(i,l)=a_matrix(i,(l-1))-k_l.*a_matrix((l-i+1),(l-1));
        end
        
        %step3
        JJ(l)=JJ(l-1).*(1-k_l^2);
    end
    
    a=a_matrix((1:l),l)'; % lpc filter coefficients are found
    
    estimated_frame=filter([0 -a(2:end)],1,frame);
    error=frame-estimated_frame;  %% calculating the error between estimated frame and original frame
    
    lpc_coeffs=[lpc_coeffs a];% concatenating lpc coefficients of each frame
    %gain calculation
    if(voiced_unvoiced(frame_no)==0) % unvoiced case, equation (8) in the report
        temp_gain= sum(error(1:length(error)).^2);
        gain= [gain sqrt(temp_gain)];
        
    else                             % voiced case, equation (9) in the report, modified the denominator part as it does not perform well
        N_over_T=floor(frame_length/pitch_period(frame_no));
        %denominator=frame_length;
        temp_gain=(pitch_period(frame_no).*sum(error(1:N_over_T*pitch_period(frame_no)).^2)/2);
        gain= [gain sqrt(temp_gain)];
    end
end

%% DECODING %%
count=1;
for frame_no=1:frame_length:(length(audio)-frame_length)
    if(voiced_unvoiced(frame_no)==1) % if voiced send impulse train
        excitation=zeros(frame_length,1);
        excitation(1:fs*pitch_period(frame_no):length(excitation))=1;
        
    else % if unvoiced send random noise
        excitation=randn(1,frame_length);
    end
    
    synthesized_frame=filter(1,lpc_coeffs((count-1)*(p+1)+1:(count)*(p+1)),excitation.*gain(count));
    count=count+1;
    synthesized_audio(frame_no:frame_no+frame_length-1)=synthesized_frame;
end

%% TESTING AND PLOTTING %%

soundsc(audio,fs); %original speech
pause(3);
soundsc(synthesized_audio,fs); %synthesized speech


figure,plot(audio);
title('Input Speech Signal');
figure,plot(synthesized_audio);
title('Decoded Speech Signal');