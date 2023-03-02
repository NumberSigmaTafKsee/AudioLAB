function [x, n] = CexpSig(a, b, lb, rb)
% lb is left bound
% rb is right bound
% x(n) = exp((a+bj)*n);
n = [lb : rb];
x = exp((a + b .* 1i) .* n);
end