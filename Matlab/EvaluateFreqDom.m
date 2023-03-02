function [magX, angX, realX, imX] = EvaluateFreqDom(X)
% The function's name is long but shotert than these lines.
magX = abs(X);
angX = angle(X);
realX = real(X);
imX = imag(X);
end