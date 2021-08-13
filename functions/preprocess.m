function out = preprocess(in)
% This function accepts an input matrix (observations x time), pre-processes, and then outputs a matrix that is filtered in 3-steps as described under Methods in Kim et al. 2021

% Low pass filter at 2kHz
Fs=10000;
y=lowpass(in',2000,Fs);
% Savitsky-golay filter
yT=sgolayfilt(y,5,27); % polynomial order of 5 and framelength of 27
% Median filter using 0.5ms window
out=movmedian(yT',6,2);
