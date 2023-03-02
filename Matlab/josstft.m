x = 1:1:1024
nframes = int32(1024/128)
N=1024
M=128
R=128
w=hanning(M)
Xtwz = zeros(N,nframes); % pre-allocate STFT output array
M = length(w);           % M = window length, N = FFT length
zp = zeros(N-M,1);       % zero padding (to be inserted)
xoff = 0;                % current offset in input signal x
Mo2 = int32((M-1)/2);    % Assume M odd for simplicity here
for m=1:nframes
  xt = x(xoff+1:xoff+M); % extract frame of input data
  xtw = w .* xt;         % apply window to current frame
  xtwz = [xtw(Mo2+1:M,1); zp; xtw(1:Mo2,1)]; % windowed, zero padded
  xtwz(:,m) = fft(xtwz); % STFT for frame m
  xoff = xoff + R;       % advance in-pointer by hop-size R
end
plot(xtwz)
