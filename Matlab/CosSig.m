function [x, n] = CosSig(A, w, p, lb , rb)
% A is amplitude
% w is omega
% p is phase
% lb is left bound
% rb is right bound
n = [lb : rb];
x = A .* cos(w .* n + p);
end