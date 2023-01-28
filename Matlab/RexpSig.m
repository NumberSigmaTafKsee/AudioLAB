function [x, n] = RexpSig(a, lb, rb)
% a is the base
% lb is left bound
% rb is right bound
% x(n) = a^n from lb to rb
n = [lb : rb];
x = a .^ n;
end