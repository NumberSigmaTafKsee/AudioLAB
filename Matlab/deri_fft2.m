function [ux, uy] = deri_fft2(u, Lx, Ly)
% The input function matrix, u, should be valued at a uniform grid
% u(i,j) = u(x_j, y_i), it is compatible with [X,Y] = meshgrid(x,y) in Matlab
% Since fft2 needs the function matrix to be periodic, u should be periodic
% Lx is the length of the grid on the x axis; Ly is the length of the grid on the y axis
[M,N] = size(u); % The length of u should be even. For the odd case, it has a litte difference
kx = 2*pi/Lx*ones(N,1)*[0:M/2-1,0,-M/2+1:-1];
ky = 2*pi/Ly*[0:N/2-1,0,-N/2+1:-1]'*ones(1,M);
ux = ifft2(1i*kx.*fft2(u)); uy = ifft2(1i*ky.*fft2(u));
ux = real(ux); uy = real(uy);
end
