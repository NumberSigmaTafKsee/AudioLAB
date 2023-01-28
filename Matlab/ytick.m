% ytick.m - set ytick marks
%
% Usage: ytick(tick)
%
% tick = vector of ticks
%
% examples: ytick([0 2 4 6 8 10]);
%           ytick(0:2:10);
%

% S. J. Orfanidis - 1999

function ytick(tick)

if nargin==0, help ytick; return; end

set(gca, 'ytick', tick);

