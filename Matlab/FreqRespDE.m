function [H, w] = FreqRespDE(x_coeff, y_coeff, range, num)
% Frequency Respone from Difference Equation
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
%//> End Old tech

m = 0:(length(x_coeff) - 1);
l = 0:(length(y_coeff) - 1);

numerator = x_coeff * exp(-1i * m' * w);
denominator =  y_coeff * exp(-1i * l' * w);

H = numerator ./ denominator;
end