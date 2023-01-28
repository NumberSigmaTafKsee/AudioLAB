va_blit_oscillators;
audio_dc_block_filter;

dc = audio_dc_block_filter.DCBlock(20/44100)

x = va_blit_oscillators.BlitSaw();    
x.setFrequency(1000);
v = zeros(1,256);
% it has to warm up a bit for the DC to settle
for j=1:10
for i=1:256
    v(i) = x.Tick();

end
end
plot(v);
pause;


x = va_blit_oscillators.BlitSquare();
x.setFrequency(1000);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();

end
end
plot(v);
pause;


% this triangle is from integrating the square
% It's not stable 
x = va_blit_oscillators.BlitTriangle();
x.setFrequency(440);
v = zeros(1,256);
for j=1:10
for i=1:256
    v(i) = x.Tick();        
end
end
plot(v);
pause;

%{
sv = Analog.AnalogSVF()
sv.setCutoff(5000)
sv.setQ(0.5)
v = zeros(1,256);
for i=1:256
    v(i) = sv.Tick(x.Tick()+1.5);
end
plot(v);
pause;
}%