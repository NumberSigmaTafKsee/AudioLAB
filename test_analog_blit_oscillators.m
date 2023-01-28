Analog;
audio_dc_block_filter;

dc = audio_dc_block_filter.DCBlock(10/44100)

x = Analog.BlitSaw()
x.setFrequency(4000);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();    
end
end
plot(v);
pause;

x = Analog.BlitSquare()
x.setFrequency(4000);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();    
end
end
plot(v);
pause;

x = Analog.BlitTriangle()
x.setFrequency(440);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();     
end
end
plot(v);
pause;

