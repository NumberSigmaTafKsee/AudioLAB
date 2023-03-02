function [X, w, r] = DTZT(x, n, rangeW, numW, rangeR, numR)
% x and n are representation of signal in time domain
% rangeW = [lbw, rbw] of w (omega)
% numW is number of equispaced points from lbw to rbw 
% rangeR = [lbr, rbr] of r (radius)
% numR is number of equispaced points from lbr to rbr
% DTZT is more general than DTFT

lbw = rangeW(1);
rbw = rangeW(2);
stepw = (rbw - lbw)/(numW - 1);
w = lbw:stepw:rbw;

%// Old tech
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

lbr = rangeR(1);
rbr = rangeR(2);
stepr = (rbr - lbr)/(numR - 1);
r = lbr:stepr:rbr;

%// Old tech
%X = zeros(length(r), length(k));
%for r_i = 1:1:length(r)
%  X(r_i,:) = (x .* (r(r_i)) .^ (-n) ) * ((1 * e) .^ (-1i*pi/M)) .^ (n'*k);
%end
%//> End old tech

X = zeros(length(r), length(w));
for r_i = 1:1:length(r)
  X(r_i,:) = (x .* ((r(r_i)) .^ (-1 .* n))) * exp(-1i * n' * w);
end
end