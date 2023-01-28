function [x, n] = SinSig(A, w, p, lb , rb)
% A is amplitude
% w is omega
% p is phase
% lb is left bound
% rb is right bound
n = [lb : rb];
x = A .* sin(w .* n + p);
end