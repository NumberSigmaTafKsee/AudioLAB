function [x, n] = Convole(x1, n1, x2, n2)
% This function requires two sets of data with fixed length.
% In other word, you must define the range of n for each signal before
% using this function.
lb = n1(1) + n2(1);
rb = n1(length(n1)) + n2(length(n2));

n = [lb : rb];

x = conv(x1, x2);
end