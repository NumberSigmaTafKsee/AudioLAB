function [x, n] = RandSig(lb, rb)
% lb is left bound
% rb is right bound
n = [lb : rb];
x = randn(1, rb - lb + 1);
end