% bandedge_a.m - bandedge frequencies for analog equalizer designs
%
% Usage: [f1,f2] = bandedge_a(f0,Df,dir);
%    
%        [f1,f2] = bandedge_a(f0,Df);          equivalent to dir=0
%        [f1,f2] = bandedge_a(f0,Df,0);        calculate f1,f2 from f0,Df
%        [f0,Df] = bandedge_a(f1,f2,1);        calculate f0,Df from f1,f2
%
% f0  = center frequency in Hz
% Df  = bandwidth in Hz
% dir = 0,1, direction of calculation
%
% f1,f2 = left and right bandedge frequencies in Hz 
%
% notes: dir=0 case computes f1,f2 from f0,Df as follows:
%            f1 =-Df/2 + sqrt(f0^2 + Df^2/4)
%            f2 = Df/2 + sqrt(f0^2 + Df^2/4)
%        dir=1 case computes f0,Df from f1,f2 as follows:
%            Df = f2-f1
%            f0 = sqrt(f1*f2)
%
%        LP case: f0=0   ==> f1=0 and f2=Df=cutoff freqency measured from DC
%        HP case: f0=Inf ==> f1=Df=cutoff frequency measured from DC and f2=Inf

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

function [f1,f2] = bandedge_a(f0,Df,dir);

if nargin==0, help bandedge_a; return; end
if nargin==2, dir=0; end

if dir==0,
   f1 = -Df/2 + sqrt(f0^2 + Df^2/4);
   f2 =  Df/2 + sqrt(f0^2 + Df^2/4);
else
   f1 = sqrt(f0*Df);
   f2 = Df - f0;
end

if f0==Inf & dir==0,
   f1 = Df;
   f2 = Inf;
end

if Df==Inf & dir==1,
   f1 = Inf;
   f2 = f0;
end   





