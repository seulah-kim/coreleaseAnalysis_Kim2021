function out = preprocess(in)

% Low pass filter at 2kHz
Fs=10000;
y=lowpass(in',2000,Fs);
% Savitsky-golay filter
yT=sgolayfilt(y,5,27); % polynomial order of 5 and framelength of 27
% Median filter using 0.5ms window
out=movmedian(yT',6,2);