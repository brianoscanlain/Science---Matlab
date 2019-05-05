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
% 
% Function to make PIP.  

function [threshold,PIP] =  makePip_BS(image,minPix)

abandon = 0;
maxPix = (floor(max(max(image))*100))/100;

runThresh = maxPix;     %Have a running threshold to track progress until endThresh
threshCount = 0;        %Initialise the counter
whitePix = [];
while runThresh > minPix
    runThresh = runThresh - .005;
    %Begin counting the iterations
    %Is it necessary to have this many iterations if the threshold cutoff
    %is 0.9-0.95?
    threshCount = threshCount+1;            
    %Populate the threshold vector
    threshold(threshCount,1) = runThresh;   
    %Get all pixels with intensities greater than threshold value
    imBW =  image>threshold(threshCount);       
    %Calculates the number of white pixels at each new threshold
    whitePix(threshCount,1) = length(find(imBW == 1));    
end

%Calculate PIP
PIP(1:length(threshold)-1,1) = (diff(whitePix) *100)./whitePix(2:end,1);
PIP = [1;PIP];

%Some housekeeping left over from original AWE
if PIP(end) > 2*PIP(end-1,1)
    PIP = PIP(1:end-1);
    threshold = threshold(1:end-1,1);
end
