Fs=1000;
 x = chirp([0:1/Fs:2],0,2,500);  # freq. sweep from 0-500 over 2 sec.
 step=ceil(20*Fs/1000);    # one spectral slice every 20 ms
 window=ceil(100*Fs/1000); # 100 ms data window

 ## test of automatic plot
 [S, f, t] = specgram(x);
 specgram(x, 2^nextpow2(window), Fs, window, window-step);
