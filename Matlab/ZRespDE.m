function [H, w, r] = ZRespDE(x_coeff, y_coeff, rangeW, numW, rangeR, numR)
% Z Respone from Difference Equation
% rangeW = [lbw, rbw] of w (omega)
% numW is number of equispaced points from lbw to rbw
% rangeR = [lbr, rbr] of r (radius)
% numR is number of equispaced points from lbr to rbr
lbw = rangeW(1);
rbw = rangeW(2);
stepw = (rbw - lbw)/(numW - 1);
w = lbw:stepw:rbw;

%//> Old Tech
% This will work for case [0, *pi], [-*pi, +*pi] but not this one [-*pi, +**pi]
%if (lbw == 0) 
%    k = 0:1:(numW - 1);
%elseif (lbw < 0)
%    k = -((numW - 1)/2):1:((numW - 1)/2);
%end
%
%M = k(length(k)); 
%w = (rbw / M) * k;
%//> End Old tech

m = 0:(length(x_coeff) - 1);
l = 0:(length(y_coeff) - 1);

lbr = rangeR(1);
rbr = rangeR(2);
stepr = (rbr - lbr)/(numR - 1);
r = lbr:stepr:rbr;

%// Old tech
%H = zeros(length(r), length(k));
%for r_i = 1:1:length(r)
%	numerator = (x_coeff .* ((r(r_i)) .^ (-m))) * exp(-1i * m' * w);
%	denominator = (y_coeff .* ((r(r_i)) .^ (-l))) * exp(-1i * l' * w);
%	H(r_i,:) = numerator ./ denominator;
%end
%//> End Old tech

H = zeros(length(r), length(w));
for r_i = 1:1:length(r)
	numerator = (x_coeff .* ((r(r_i)) .^ (-1 .* m))) * exp(-1i * m' * w);
	denominator = (y_coeff .* ((r(r_i)) .^ (-1 .* l))) * exp(-1i * l' * w);
	H(r_i,:) = numerator ./ denominator;
end
end