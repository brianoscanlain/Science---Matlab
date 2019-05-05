%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automated Whitecap Extraction image processing algorithm.
% 
% For algorithm descroption see:
% Callaghan and White, (2009), Automated Processing of Sea Surface Images
% for the Determination of Whitecap Coverage, Vol. 26, pp.383-394
%
% Please contact Adrian Callaghan before using this code.
% callaghan.adrian@gmail.com
%
% Disclaimer:
% This code has not been rigorously tested and may contain bugs.
% All queries should be directed to callaghan.adrian@gmail.com
%
% This code version has been specifically written to handle 5 Mega Pixel
% images and may not run correctly with images of lower resolution.
%
% Adrian Callaghan 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculates the derivative of a vector.

function deriv = calcDeriv(yData)

[nx ny] = size(yData);
if ny > 1
    yData = tData';
end

deltaX = 0.01;
deriv = diff(yData)./deltaX;
deriv = [0;deriv];      %add 0 to compensate for diff function

