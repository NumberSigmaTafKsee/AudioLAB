va_blit_oscillators;
o = zeros(1,256);
x = va_blit_oscillators.BlitSaw();
% it needs to warm up because there is a dc offset in the start
for j=1:10
for i=1:256
    o(i) = x.Tick();
end
end
plot(o);
pause;