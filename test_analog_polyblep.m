Analog;
audio_dc_block_filter;

dc = audio_dc_block_filter.DCBlock(20/44100)

x = Analog.PolyBLEP(44100,Analog.PolyBLEP_SAWTOOTH)
x.setFrequency(440);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();    
end
end
plot(v);
pause;

x = Analog.PolyBLEP(44100,Analog.PolyBLEP_SQUARE)
x.setFrequency(440);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();    
end
end
plot(v);
pause;

x = Analog.PolyBLEP(44100,Analog.PolyBLEP_TRIANGLE)
x.setFrequency(440);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();    
end
end
plot(v);
pause;

