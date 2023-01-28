function [ux] = deri_fft(u, Lx)
% The input function vector, u, should be valued at a uniform grid
% Since fft needs the function vector to be periodic, u should be periodic
% Lx is the length of the grid
N = length(u); % The length of u should be even. For the odd case, it has a litte difference
kx = 2*pi/Lx*[0:N/2-1,0,-N/2+1:-1];
ux = ifft(1i*kx.*fft(u));
ux = real(ux);
end

