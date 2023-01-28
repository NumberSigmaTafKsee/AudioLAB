function freqplot(fdata, ydata, symbol, ttl, xlab, ylab)
% FREQPLOT - Plot a function of time, Matlab version.

  if nargin<6, ylab=''; end
  if nargin<5, xlab='Frequency (Hz)'; end
  if nargin<4, ttl=''; end
  if nargin<3, symbol=''; end
  if nargin<2, fdata=0:length(ydata)-1; end

  plot(fdata,ydata,symbol); grid;
  if ttl, title(ttl); end
  if ylab, ylabel(ylab); end
  xlabel(xlab);


