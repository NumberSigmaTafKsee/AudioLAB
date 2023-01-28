% ddiff.m - ideal lowpass differentiator FIR filter
%
% h = ddiff(wc, N) = row vector
%
% N = 2M+1 = filter length (odd)
% wc in rads/sample

function h = ddiff(wc, N)

M = (N-1)/2;

for k = -M:M,
  if k == 0,
    h(k+M+1) = 0;
  else
    h(k+M+1) = wc * cos(wc * k) / (pi * k) - sin(wc * k) / (pi * k^2);
  end
end
