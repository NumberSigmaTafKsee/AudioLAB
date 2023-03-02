function [x, n] = Fold(x1, n1)
x = fliplr(x1);
n = -fliplr(n1);
end