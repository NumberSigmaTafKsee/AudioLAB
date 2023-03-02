function [z, p, k, d] = ibilinear(zd, pd, kd, fs, fp, fp1)
%IBILINEAR Inverse bilinear transformation with optional frequency
%   prewarping.
%   [Z,P,K] = IBILINEAR(Zd,Pd,Kd,Fs) converts the z-domain transfer
%   function specified by Zd, Pd, and Kd to a s-transform discrete
%   equivalent obtained from the inverse bilinear transformation:
%
%      H(s) = H(z) |
%                  | z = (s+2*Fs)/(s-2*Fs)
%
%   where column vectors Zd and Pd specify the zeros and poles, scalar
%   Kd specifies the gain, and Fs is the sample frequency in Hz.
%   [NUM,DEN] = IBILINEAR(NUMd,DENd,Fs), where NUM and DEN are
%   row vectors containing numerator and denominator transfer
%   function coefficients, NUM(z)/DEN(z), in descending powers of
%   z, transforms to s-transform coefficients NUM(s)/DEN(s).
%   [A,B,C,D] = IBILINEAR(Ad,Bd,Cd,Dd,Fs) is a state-space version.
%   Each of the above three forms of IBILINEAR accepts an optional
%   additional input argument that specifies prewarping. For example,
%   [Z,P,K] = IBILINEAR(Zd,Pd,Kd,Fs,Fp) applies prewarping before
%   the inverse bilinear transformation so that the frequency responses
%   before and after mapping match exactly at frequency point Fp
%   (match point Fp is specified in Hz).
%
%   See also BILINEAR, IMPINVAR.
%   Author(s): J.N. Little, 4-28-87
%   	   J.N. Little, 5-5-87, revised
%   Copyright (c) 1988-98 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1998/06/03 14:41:57 $
%   Gene Franklin, Stanford Univ., motivated the state-space
%   approach to the bilinear transformation.
%
%	 Adapted/plagiarised from BILINEAR.m by Paul Eccles, 2 April 2005
[mn,nn] = size(zd);
[md,nd] = size(pd);
if (nd == 1 & nn < 2) & nargout ~= 4	% In zero-pole-gain form
	if mn > md
		error('Numerator cannot be higher order than denominator.')
	end
	if nargin == 5		% Prewarp
		fp = 2*pi*fp;
		fs = fp/tan(fp/fs/2);
	else
		fs = 2*fs;
   end

   % prune extra zeroes. There must be a better way of doing this.
	for i = size(zd):-1:1,
   	if (zd(i) < -0.99) & (zd(i) > -1.01)
         zd = zd(1:(i-1));		% Prune final element.
      else
         break;					% Leave the 'for' loop.
   	end
	end
	p = fs*(pd-1)./(pd+1);
	z = fs*(zd-1)./(zd+1);
	k = kd*prod((fs)-p)./prod((fs)-z);
elseif (md == 1 & mn == 1) | nargout == 4 %
	if nargout == 4		% State-space case
		ad = zd; bd = pd; cd = kd; dd = fs; fs = fp;
		error(abcdchk(ad,bd,cd,dd));
		if nargin == 6			% Prewarp
			fp = fp1;		% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		end
	else			% Transfer function case
		if nn > nd
			error('Numerator cannot be higher order than denominator.')
		end
		num = zd; den = pd;		% Decode arguments
		if nargin == 4			% Prewarp
			fp = fs; fs = k;	% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		else
			fs = kd;			% Decode arguments
		end
		% Put num(s)/den(s) in state-space canonical form.
		[ad,bd,cd,dd] = tf2ss(num,den);
	end
	% Now do state-space version of bilinear transformation:
	t = 1/fs;
   r = sqrt(t);
   I = eye(size(ad));
   a = 2*fs*(ad-I)/(ad+I);
   a1 = I-a/(2*fs);
   b = a1*bd*r/t;
   c = cd*a1/r;
   d = dd - c/a1*b*t/2;

   if nargout == 4				% state-space form result
		z = a; p = b; k = c;
	else								% transfer function form result
		[z,p] = ss2tf(a,b,c,d);
	end
else
	error('First two arguments must have the same orientation.')
end

