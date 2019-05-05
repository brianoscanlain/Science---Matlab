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

function [smoothed] = mySmooth_Deriv(origPerInc)

[cx cy]= size(origPerInc);
if cx < cy
    origPerInc = fliplr(origPerInc);
else
    origPerInc = flipud(origPerInc);
end
%Get the size of the perInc Vector
len = length(origPerInc);
smoothed = zeros(len,1);
smoothed(1) = origPerInc(1);
smoothed(2) = origPerInc(2);
% smoothed(3) = mean(origPerInc(2:4));    %Added with weighting 3:1
smoothed(len-1) = origPerInc(len-1);
smoothed(len) = origPerInc(len);

for dd=3:len-2
    smoothed(dd) = sum(origPerInc(dd-2:dd+1))/4;
end
% smoothed(end-1) = mean(smoothed(end-2:end));
[cx cy]= size(smoothed);
if cx < cy
    smoothed = fliplr(smoothed);
else
    smoothed = flipud(smoothed);
end