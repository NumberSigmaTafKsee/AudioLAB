% swanalmainplot.m
% Compare measured and theoretical frequency response.
% This script is invoked by swanalmainplot.m and family,
% and requires context set up by the caller.

figure(N+1); % figure number is arbitary

subplot(2,1,1);
ttl = 'Amplitude Response';
freqplot(f,gains,'*k',ttl,'Frequency (Hz)','Gain');
tar = 2*cos(pi*f/fs); % theoretical amplitude response
hold on; freqplot(f,tar,'-k'); hold off;
text(-0.08,mean(ylim),'(a)');

subplot(2,1,2);
ttl = 'Phase Response';
tpr = -pi*f/fs; % theoretical phase response
pscl = 1/(2*pi);% convert radian phase shift to cycles
freqplot(f,tpr*pscl,'-k',ttl,'Frequency (cycles)',...
	'Phase shift (cycles)');
hold on; freqplot(f,phases*pscl,'*k'); hold off;
text(-0.08,mean(ylim),'(b)');
saveplot(plotfile); % set by caller
