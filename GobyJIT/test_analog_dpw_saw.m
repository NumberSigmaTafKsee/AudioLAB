Analog;
audio_dc_block_filter;

dc = audio_dc_block_filter.DCBlock(20/44100)

x = Analog.DPWSaw(44100)
x.setFrequency(440)

v = zeros(1,256);
% it has to warm up a bit for the DC to settle
for j=1:10
for i=1:256
    v(i) = x.Tick();
end
end
plot(v);
pause;
