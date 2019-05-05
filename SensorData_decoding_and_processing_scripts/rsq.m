function [ Out ] =rsq(observation,prediction)
%R2=rsq(x,y)   ...where y = observation, x = modelled observation
% Calculates the R-squared correlation coefficient by:
%
%Rsq =   1   -  sum( (observations(i) -  model(i))^2)
%               --------------------------------------------
%               (Length(observations)-1) * (var(observations)
%
%
%        i.e. Rsq = 1-SSR/SST   (in statistical terminology)
%
% Brian Scanlon, School of Physics, NUI Galway, Ireland
% 06/05/2015

if nargin<2
    error('Not enough input arguments, xdata and y data required');
elseif nargin>2
    error('Too many input arguments, xdata and y data only required');
end

SSR=sum((observation - prediction).^2);
SST=var(observation)*(length(observation)-1);

Out=  1 - SSR/SST;
end

