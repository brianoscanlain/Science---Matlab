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
%calculate percentage whitecap coverage

function W = calcW_BS(editImage,thresh,showImages)
%Script to calculate W from the image given its threshold
dispImage = editImage;
dispImage(dispImage < thresh) = 0;
[nr nc] = size(dispImage);
numWhite = length(find(dispImage ~= 0)); W = numWhite*100 ./(nr*nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eliminating possible images that process but have no whitecaps and
%therefore change the overall values when averaged together
%if W < 0.002
 %   dispImage = 0;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if showImages
    figure;
    imshow(dispImage); title(['Thresh : ' num2str(thresh) '. W = ' num2str(W) ' %']);
end