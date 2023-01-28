[x, n] = RexpSig(0.5, 0, 500);
range = [0, pi];
num = 501;
[X, w] = DTFT(x, n, range, num);
[magX, angX, realX, imX] = EvaluateFreqDom(X);
PlotEvaluateFreqDom(X, w, range);
