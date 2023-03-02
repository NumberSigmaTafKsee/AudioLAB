function [x, n] = Add(x1, n1, x2, n2)
% x1 and n1 are illumination of first signal
% x2 and n2 are illumination of second signal

n = [min(min(n1),min(n2)) : max(max(n1),max(n2))];

y1 = zeros(1,length(n)); 
y2 = y1;

y1((n>=min(n1))&(n<=max(n1))==1) = x1;
y2((n>=min(n2))&(n<=max(n2))==1) = x2;

x = y1 + y2;
end