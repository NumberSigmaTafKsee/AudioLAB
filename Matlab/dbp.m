% dbp.m - ideal bandpass FIR filter
%
% h = dbp(wa, wb, N) = row vector
%
% N = 2M+1 = filter length (odd)
% wa, wb = cutoff frequencies in [rads/sample]

function h = dbp(wa, wb, N)

M = (N-1)/2;

for k = -M:M,
   if k == 0,
      h(k+M+1) = (wb - wa) / pi;
   else
      h(k+M+1) = sin(wb * k) / (pi * k) - sin(wa * k) / (pi * k);
   end
end
