function [x, n] = InScale(x1, n1, a)
% WARNING: third parameter aka a only recieves positive integer

% If you want to use negative integer, you should use the function fliplr() for x1 and n1
% then use abs(that negative integer) then InScale()

% Or you can use Fold()

% How about real number? Well if this line is still existed that mean I
% have not implemented it yet

% Indicate the right bound of n
i = 0;
while true
    i = i + 1;
    if ((i + 1) * a > n1(length(n1))) 
        break;
    end
end
rb = i;

% Indicate the left bound of n
i = 0;
while true
    i = i - 1;
    if ((i - 1) * a < n1(1))
        break;
    end
end
lb = i;

n = [lb : rb];
x = zeros(1, length(n));

% Get the value from x1 to input into x
for i = 1:length(x)
    x(i) = GetValueAt(x1, n1, n(i) * a);
end
end