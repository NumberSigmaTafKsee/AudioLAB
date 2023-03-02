x_coeff = [0.0181, 0.0543, 0.0543, 0.0181];
y_coeff = [1.0000, -1.7600, 1.1829, -0.2781];
range = [0, pi];
num = 501;
[H, w] = FreqRespDE(x_coeff, y_coeff, range, num);
PlotEvaluateFreqDom(H, w, range);
