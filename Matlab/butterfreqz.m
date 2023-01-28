fc = 300;
fs = 1000;

[b,a] = butter(6,fc/(fs/2));

freqz(b,a,[],fs)

subplot(2,1,1)
ylim([-100 20])
