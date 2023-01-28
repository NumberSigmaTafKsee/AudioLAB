function [x, n] = UnitStepSig(n0, lb, rb)
% n0 is where UnitStepSignal start giving 1
% lb is left bound
% rb is right bound
n = [lb : rb];
x = [double((n - n0) >= 0)];
end