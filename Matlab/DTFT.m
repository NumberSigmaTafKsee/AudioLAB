function [X, w] = DTFT(x, n, range, num)
% x and n are representation of signal in time domain
% range = [lbw, rbw] of w (omega)
% num is number of equispaced points from lbw to rbw 
lbw = range(1);
rbw = range(2);
stepw = (rbw - lbw)/(num - 1);
w = lbw:stepw:rbw;

%// Old tech
% This will work for case [0, *pi], [-*pi, +*pi] but not this one [-*pi, +**pi]
%if (lbw == 0) 
%    k = 0:1:(num - 1);
%elseif (lbw < 0)
%    k = -((num - 1)/2):1:((num - 1)/2);
%end
%
%M = k(length(k)); 
%w = (rbw / M) * k;
%
%e = exp(1); % Define e. I like this.
%X = x * (e .^ (-1i*pi/M)) .^ (n'*k); 
%//> End Old techh

X = x * exp(-1i * n' * w);
end