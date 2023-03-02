% sigav.m - signal averaging
%
% y = sigav(D, N, x)
%
% D = length of each period
% N = number of periods
% x = row vector of length at least ND (doesn't check it)
% y = length-D row vector containing the averaged period
% It averages the first N blocks in x

function y = sigav(D, N, x)

y = 0;

for i=0:N-1,
       y = y + x((i*D+1) : (i+1)*D);       % accumulate i-th period
end

y = y / N;
