function [xa, t] = ReconstructThis(x, n, Ts, rangeT, deltaT, mode)
% x and n are 2 last outputs of SsampleThis()
% Ts is sampling interval. Fs = 1/Ts.
% rangeT = [lbt, rbt]. t runs from lbt from rbt.
% deltaT is dt.
% if mode == 1 then use default formula
% if mode == 2 then use cubic spline interpolation
nTs = n * Ts;
Fs = 1/Ts;
t = rangeT(1):deltaT:rangeT(2);
if mode == 1
    xa = x * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
end
if mode == 2
    xa = spline(nTs, x, s);
end
end