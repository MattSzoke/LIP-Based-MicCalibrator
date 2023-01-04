function [ Frequency, Yfft, Frequency_half, Yfft_half ] = myFFTc( y, Fsampling, Lpad)
% FFT transformaTHOR: One and two sided
% Written by Matt Szoke, m.szoke@vt.edu
% Use as:
% [ Frequency, Yfft ] = myFFTc( y, Fsampling, Lpad)
% y: signal to FFT
% Lpad: padding length, in number of timesteps

y  = y - mean(y);
dT = 1/Fsampling; % time step size
Ly = length(y);

if Lpad > Ly
	f  = linspace(0.0, 1.0/(2.0*dT), Lpad/2); % frequency vector
	Y = fft(y,Lpad)/Ly;
	L = Lpad;
else
	L = Ly;
	f  = linspace(0.0, 1.0/(2.0*dT), L/2); % frequency vector
	Y = fft(y,L)/L;
end

    
Yfft_half        = 2*abs(Y(1:(L/2)));
Frequency_half   = f ;

Yfft        = Y;
Frequency   = [-fliplr(f), f(2:end)] ;

end
