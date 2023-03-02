function [] = PlotEvaluateFreqDom(X, w, range)
% The function's name is long but shotert than these lines.
rbw = range(2);

[magX, angX, realX, imX] = EvaluateFreqDom(X);
subplot(2,2,1);
plot(w/rbw, magX);
xlabel('Frequency in pi units'); title('Magnitude Part'); ylabel('Magnitude');

subplot(2,2,3);
plot(w/rbw, angX/rbw);
xlabel('Frequency in pi units'); title('Angle Part'); ylabel('Radians');

subplot(2,2,2);
plot(w/rbw, realX);
xlabel('Frequency in pi units'); title('Real Part'); ylabel('Real');

subplot(2,2,4);
plot(w/rbw, imX);
xlabel('Frequency in pi units'); title('Imaginary Part'); ylabel('Imaginary');
end