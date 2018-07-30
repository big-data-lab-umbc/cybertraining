function y=periodogram2(x)
% Author: Wang Xianju
% April 1th, 2002
% This function computes the 2-D periodogram based on Fourier Method.
   
x=fft2(x);
x=fftshift(x);

x=x.*conj(x);

x = x;

y=x;