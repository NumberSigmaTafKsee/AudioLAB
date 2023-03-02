function [x_odd, x_even] = OddEvenSyn(x)
x_odd = 0.5 .* (x - fliplr(x));
x_even = 0.5 .* (x + fliplr(x));
end