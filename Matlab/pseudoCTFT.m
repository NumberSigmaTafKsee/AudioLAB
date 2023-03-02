function [Xa, w] = pseudoCTFT(analog, rangeT, deltaT, rangeW, numW)
% analog is an anonymous function which is an analog signal aka
% continous-time sginal.
% rangeT = [lbt, rbt]. t runs from lbt from rbt.
% deltaT is dt.
% rangeW = [lbw, rbw]. w runs from lbn to rbn.
% numW is number of w.
t = rangeT(1):deltaT:rangeT(2);
xa = analog(t);

lbw = rangeW(1);
rbw = rangeW(2);
stepw = (rbw - lbw)/(numW - 1);
w = lbw:stepw:rbw;

Xa = xa * exp(-1i*t'*w) * deltaT;
end