% xtick.m - set xtick marks
%
% Usage: xtick(tick)
%
% tick = vector of ticks
%
% examples: xtick([0 2 4 6 8 10]);
%           xtick(0:2:10);
%

% S. J. Orfanidis - 1999

function xtick(tick)

if nargin==0, help xtick; return; end

set(gca, 'xtick', tick);

