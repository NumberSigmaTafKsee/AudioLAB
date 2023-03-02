function [x, n] = RectSig(n0, n1, lb, rb)
% n0 is where RectSignal start giving 1
% n1 is where RectSignal end giving 1
% lb is left bound
% rb is right bound
n = [lb : rb];
x = [double((n - n0) >= 0 & (n1 - n) >= 0)];
end