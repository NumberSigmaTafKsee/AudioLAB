function [x, n] = Crosscorrelate(x1, n1, x2, n2)
[x2, n2] = Fold(x2, n2);
[x, n] = Convole(x1, n1, x2, n2);
end