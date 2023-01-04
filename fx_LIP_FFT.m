function [ Frequency, Yfft, Frequency_half, Yfft_half ] = fx_LIP_FFT(MyData,Fs,Lpad)

MyData = MyData - mean(MyData);
     L = length(MyData);

% Window options. c is correction factor.
% MyWindow =  tukeywin(L,0.25)'; c = 2.28571;
MyWindow = ones(size(MyData)); c = 1;

MyData = MyData.*MyWindow;
[ Frequency, Yfft, Frequency_half, Yfft_half ] = myFFTc( MyData, Fs, Lpad);
Yfft      = Yfft*c;
Yfft_half = Yfft_half*c;

end

