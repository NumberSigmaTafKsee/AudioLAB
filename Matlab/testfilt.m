% testfilt.m - test of seven filtering realizations
% 
%    ytr = transpfilt - frequency-shifted transposed            
%    ynl = nlattfilt  - frequency-shifted normalized lattice
%    yst = statefilt  - frequency-shifted optimum state-space
%    ydf = df2filt    - frequency-shifted direct-form-II
%    yca = cascfilt   - cascade of 4th order sections
%    yde = decoupfilt - cascade of decoupled realizations 
%    yfi = filter     - full-length coefficients 

% -------------------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058, USA
%
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% -------------------------------------------------------------------------


N=5; type=0; 
w0 = 0.1*pi; Dw = 0.2*pi;
G0 = 0; G = 12; GB = 9; Gs = 0.1; 

[Be,Ae,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type,Gs);

[gamma,d] = dir2latt(Bh,Ah);           % normalized lattice coefficients
[A,B,C,D] = dir2state(Bh,Ah);          % state-space coefficients
[gamma,dec] = dir2decoup(Bh,Ah);       % decoupled realizations

b = cas2dir(Be); a = cas2dir(Ae);      % convolve all numerator and denominator 4th order sections

seed = 2005; randn('state',seed); x = randn(1000,1);

ytr = transpfilt(Bh,Ah,w0,x);
ynl = nlattfilt(gamma,d,w0,x);
yst = statefilt(A,B,C,D,w0,x);
ydf = df2filt(Bh,Ah,w0,x);
yca = cascfilt(Be,Ae,x);
yde = decoupfilt(gamma,dec,w0,x);
yfi = filter(b,a,x);

error1 = norm(ytr-ynl) + norm(ytr-yst) + norm(ytr-ydf) + norm(ytr-yde);
error2 = norm(ytr-yca);
error3 = norm(ytr-yfi);


fprintf('\noutput    function      description \n');
fprintf('---------------------------------------------------------------------\n');
fprintf(' ytr      transpfilt    frequency-shifted transposed  \n');      
fprintf(' ynl      nlattfilt     frequency-shifted normalized lattice \n');
fprintf(' yst      statefilt     frequency-shifted optimum state-space \n');
fprintf(' ydf      df2filt       frequency-shifted direct-form-II \n');
fprintf(' yde      decoupfilt    frequency-shifted decoupled realizations \n');
fprintf(' yca      cascfilt      cascade of 4th order sections \n');
fprintf(' yfi      filter        use full-length coefficients (least accurate)\n\n');
fprintf('                                                                   error\n');
fprintf('                                                                  --------\n');
fprintf('norm(ytr-ynl) + norm(ytr-yst) + norm(ytr-ydf) + norm(ytr-yde)  =  %5.1d\n', error1);
fprintf('                                                norm(ytr-yca)  =  %5.1d\n', error2);
fprintf('                                                norm(ytr-yfi)  =  %5.1d\n', error3);
