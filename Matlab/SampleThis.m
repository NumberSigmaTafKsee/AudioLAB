function [X, w, x, n] = SampleThis(analog, Ts, rangeN, numN, rangeW, numW)
% analog is an anonymous function which is an analog signal aka
% continous-time sginal.
% Ts is sampling interval. Fs = 1/Ts.
% rangeN = [lbn, rbn]. n runs from lbn to rbn.
% numN is number of n.
% rangeW = [lbw, rbw]. w runs from lbn to rbn.
% numW is number of w.

lbn = rangeN(1);
rbn = rangeN(2);
stepn = (rbn - lbn)/(numN - 1);
n = lbn:stepn:rbn;

x = analog(n, Ts);

[X, w] = DTFT(x, n, rangeW, numW);
end