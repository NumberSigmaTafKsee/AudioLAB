function saveplot(filename)
% SAVEPLOT - Save current plot to disk in a PostScript file.
%            This version is compatible only with Matlab.

cmd = ['print -deps ',filename]; % for color, use '-depsc'
disp(cmd); eval(cmd);
