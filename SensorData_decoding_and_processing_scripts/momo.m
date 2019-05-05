function f=momo(a,n,ns)
%  To generate Morlet Wavelet in Fourier Domain
nu=a*ns*(0:n/2)'/n;
f=exp(-1/sqrt(2)*((nu-1)/.220636).^2);
