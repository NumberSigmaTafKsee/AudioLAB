% hpeqex0.m - example output from hpeq.m
%  
% Butterworth equalizer - boost
% Butterworth equalizer - cut
% Elliptic equalizer 
% Chebyshev-2 lowpass shelving filter 
% Elliptic highpass shelving filter 
% Chebyshev-1 bandpass filter
% Elliptic bandpass filter
% Elliptic bandstop filter
% Chebyshev-2 bandstop filter 
% Elliptic lowpass filter 
% Chebyshev-2 highpass filter 

% --------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% --------------------------------------------------------------

clear all;

% -----------------------------------------------------------------------------

disp(' ');
disp('Butterworth boost example:');
disp('--------------------------');
disp('N = 5; fs=40; f0 = 5; Df = 4; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;'); 
disp('G0 = 0; G = 12; GB = 9; type = 0;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],20*log10([GB,GB]),''b.'');');

N = 5; fs=40; f0 = 5; Df = 4; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs; 
G0 = 0; G = 12; GB = 9; type = 0;

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type)

f = linspace(0,20,1001); w = 2*pi*f/fs;

H = 20*log10(abs(fresp(B,A,w)));

% H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;

plot(f,H,'r-', [f1,f2],([GB,GB]),'b.');
ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
title('Butterworth Boost {\it N} = 5');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;

% -----------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Butterworth cut example:');
disp('------------------------');
disp('N = 5; fs=40; f0 = 5; Df = 4; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;'); 
disp('G0 = 0; G = -12; GB = -9; type = 0;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],[GB,GB],''b.'');');

N = 5; fs=40; f0 = 5; Df = 4; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs; 
G0 = 0; G = -12; GB = -9; type = 0;

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type)

f = linspace(0,20,1001); w = 2*pi*f/fs;

H = 20*log10(abs(fresp(B,A,w)));

% H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;

plot(f,H,'r-', [f1,f2],[GB,GB],'b.');
ylim([-14 14]); ytick(-12:3:12);
xlim([0,20]); xtick(0:2:20);
title('Butterworth Cut {\it N} = 5');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;

% -----------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Elliptic example:');
disp('-----------------');
disp('N = 4; fs=40; f0 = 5; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;')
disp('G0 = 0; G = 12; GB = 11.99; type = 3; Gs = 0.01; tol=eps;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],[GB,GB],''b.'');');

N = 4; fs=40; f0 = 5; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs; 
G0 = 0; G = 12; GB = 11.99; type = 3; Gs = 0.01; tol=eps;

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = 20*log10(abs(fresp(B,A,w)));

% H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;

plot(f,H,'r-', [f1,f2],[GB,GB],'b.');
ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
title('Elliptic {\it N} = 4');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;

% -----------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Chebyshev-2 lowpass shelving filter:');
disp('------------------------------------');
disp('N = 4; fs=40; f0 = 0; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;')
disp('G0 = 0; G = 12; GB = 0.01; type = 2;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],[GB,GB],''b.'');');

N = 4; fs=40; f0 = 0; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs; 
G0 = 0; G = 12; GB = 0.01; type = 2;

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type)

H = 20*log10(abs(fresp(B,A,w)));

% H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;

plot(f,H,'r-', [f1,f2],[GB,GB],'b.');
ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
title('Chebyshev-2 Lowpass Shelf {\it N} = 4');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Elliptic highpass shelving filter:');
disp('----------------------------------');
disp('N = 5; fs=40; f0 = fs/2; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;')
disp('G0 = 0; G = 12; GB = 11.99; type = 3; Gs = 0.01; tol=eps;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],[GB,GB],''b.'');');

N = 5; fs=40; f0 = fs/2; Df = 2; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs; 
G0 = 0; G = 12; GB = 11.99; type = 3; Gs = 0.01; tol=eps;

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = 20*log10(abs(fresp(B,A,w)));

% H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;

plot(f,H,'r-', [f1,f2],[GB,GB],'b.');
ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
title('Elliptic Highpass Shelf {\it N} = 5');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Ordinary Chebyshev-1 bandpass filter:');
disp('-------------------------------------');
disp('fs=40; f1=4; f2=6; f1s=3.5; f2s=6.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;');
disp('G0 = -Inf; G = 0; GB = -0.10; type = 3; Gs = -30; tol=eps;');
disp('[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [f1,f2],10.^([GB,GB]/20),''b.'', [f1s,f2s],10.^([Gs,Gs]/20),''b.'');');

fs=40; f1=4; f2=6; f1s=3.5; f2s=6.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;
G0 = -Inf; G = 0; GB = -0.1; type = 1; Gs = -30; tol=eps;

[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [f1,f2],10.^([GB,GB]/20),'b.', [f1s,f2s],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Chebyshev-1 Bandpass');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Ordinary elliptic bandpass filter:');
disp('----------------------------------');
disp('fs=40; f1=4; f2=6; f1s=3.5; f2s=6.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;');
disp('G0 = -Inf; G = 0; GB = -0.10; type = 3; Gs = -30; tol=eps;');
disp('[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [f1,f2],10.^([GB,GB]/20),''b.'', [f1s,f2s],10.^([Gs,Gs]/20),''b.'');');

fs=40; f1=4; f2=6; f1s=3.5; f2s=6.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;
G0 = -Inf; G = 0; GB = -0.10; type = 3; Gs = -30; tol=eps;

[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [f1,f2],10.^([GB,GB]/20),'b.', [f1s,f2s],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Elliptic Bandpass');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Ordinary elliptic bandstop filter:');
disp('----------------------------------');
disp('fs=40; f1=4; f2=6; f1s=3.5; f2s=6.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;');
disp('G0 = 0; G = -Inf; GB = -0.10; type = 3; Gs = -30; tol=eps;');
disp('[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [f1,f2],10.^([GB,GB]/20),''b.'', [f1s,f2s],10.^([Gs,Gs]/20),''b.'');');

fs=40; f1=12; f2=14; f1s=11.5; f2s=14.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;
G0 = 0; G = -Inf; GB = -30; type = 3; Gs = -0.1; tol=eps;

[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [f1,f2],10.^([GB,GB]/20),'b.', [f1s,f2s],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Elliptic Bandstop');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Ordinary Chebyshev-2 bandstop filter:');
disp('-------------------------------------');
disp('fs=40; f1=12; f2=14; f1s=11.5; f2s=14.5; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;');
disp('G0 = 0; G = -Inf; GB = -30; type = 3; Gs = -0.1; tol=eps;');
disp('[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [f1,f2],10.^([GB,GB]/20),''b.'', [f1s,f2s],10.^([Gs,Gs]/20),''b.'');');

fs=40; f1=11.5; f2=14.5; f1s=12; f2s=14; w1=2*pi*f1/fs; w2=2*pi*f2/fs; w1s=2*pi*f1s/fs; w2s=2*pi*f2s/fs;
G0 = 0; G = -Inf; GB = -0.1; type = 2; Gs = -30; 

[w0,Dw] = bandedge(w1,w2,1); Dws = w2s-w1s;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [f1,f2],10.^([GB,GB]/20),'b.', [f1s,f2s],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Chebyshev-2 Bandstop');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Elliptic lowpass filter with specs:');
disp('------------------------------------');
disp('fs=40; f2=4; f2s=4.5; w0=0; Dw=2*pi*f2/fs; Dws=2*pi*f2s/fs;');
disp('G0 = -Inf; G = 0; GB = -0.10; type = 1; Gs = -30; tol=eps;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [0,f2],10.^([GB,GB]/20),''b.'', [0,f2s],10.^([Gs,Gs]/20),''b.'');');

fs=40; f2=4; f2s=4.5; w0=0; Dw=2*pi*f2/fs; Dws=2*pi*f2s/fs;
G0 = -Inf; G = 0; GB = -0.10; type = 3; Gs = -30; tol=eps;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [0,f2],10.^([GB,GB]/20),'b.', [0,f2s],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Elliptic Lowpass');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% ------------------------------------------------------------------------------------

disp(' '); disp('type <RET> to continue...'); pause;

disp(' ');
disp('Chebyshev-2 highpass filter with specs:');
disp('---------------------------------------');
disp('fs=40; f1=14; f1s=16; w0=pi; Dw=2*pi*(20-f1)/fs; Dws=2*pi*(20-f1s)/fs;');
disp('G0 = -Inf; G = 0; GB = -0.30; type = 3; Gs = -0.1; tol=eps;');
disp('Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs; H = abs(fresp(B,A,w));');
disp('plot(f,H,''r-'', [f1,20],10.^([GB,GB]/20),''b.'', [f1s,20],10.^([Gs,Gs]/20),''b.'');');

fs=40; f1=14; f1s=16; w0=pi; Dw=2*pi*(20-f1)/fs; Dws=2*pi*(20-f1s)/fs;
G0 = -Inf; G = 0; GB = -30; type = 2; Gs = -0.1; tol=eps;

Nexact = hpeqord(G0, G, GB, Gs, Dw, Dws, type), N=ceil(Nexact);

[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

H = abs(fresp(B,A,w));

% H = abs(fresp(Bh,Ah,w,w0));  % alternative computation of frequency response

plot(f,H,'r-', [f1,20],10.^([GB,GB]/20),'b.', [f1s,20],10.^([Gs,Gs]/20),'b.');
ylim([0,1.05]); ytick(0:0.2:1);
xlim([0,20]); xtick(0:2:20);
title('Chebyshev-2 Highpass');
xlabel('{\it f}  (kHz)'); ylabel('absolute units');
grid;


% -------------------------------------------------

%disp(' '); 
disp('type <RET> to finish...'); pause;
close all;




