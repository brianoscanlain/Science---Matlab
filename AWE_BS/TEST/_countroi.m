function [ROI_count,roiDX] = countroi( editImage,threshold )
%UNTITLED2 Summary of this function goes here
%   Uses region props to count the number of regions per binary image using
%   the values in the threshold array. The roiDX is the 1st derivative of
%   the

hh=length(threshold);
ROI_count=zeros(hh,1);

for i=1:hh
dispImage=editImage>threshold(i,1);
ROI_count(i,1)=length(regionprops(dispImage,'PixelIdxList'));
end

% Now time for smoothing and filtering out any peaks.





%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%            Smooth now
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

roiDX=[0;diff(ROI_count)./ROI_count(2:end,1)];
roiDX=flipud(roiDX);  %smoothed::

% Two 4 point average running means ::
for ii=3:hh-2
   roiDX(ii,1)=sum(roiDX(ii-2:ii+1))/4; 
end
roiDX(end-1,1)=mean(roiDX(end-2:end));

for ii=3:hh-2
   roiDX(ii,1)=sum(roiDX(ii-2:ii+1))/4; 
end
roiDX(end-1,1)=mean(roiDX(end-2:end));
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if roiDX(end-1,1)<roiDX(end,1) && roiDX(end-1,1)<roiDX(end-2)
    roiDX(end-1,1)=mean([roiDX(end,1);roiDX(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if roiDX(end-1,1)>roiDX(end,1) && roiDX(end-1,1)>roiDX(end-2)
    roiDX(end-1,1)=mean([roiDX(end,1);roiDX(end-2,1)]);
end

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$                Filter to remove any singularities        $$$$$$$$$$
%$$$$$$                                                          $$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for ii=1:hh-3
    if (roiDX(ii,1)<0 && roiDX(ii+1,1)>0 && roiDX(ii+2,1)>0 && ...
            roiDX(ii+3,1)<0) || (roiDX(ii,1)>0 && roiDX(ii+1,1)<0 &&...
            roiDX(ii+3,1)<0 && roiDX(ii+3,1)>0)
       %a two value peak/singularity 
       roiDX(ii+1:ii+2,1)=[mean(roiDX([ii,ii+3],1));mean(roiDX([ii...
           ,ii+3],1))];
    end
end
for iii=2:hh-1
   
    if (roiDX(iii,1)<0 && roiDX(iii-1,1)>0 && roiDX(iii+1,1)>0) || ...
            (roiDX(iii,1)>0 && roiDX(iii-1,1)<0 && roiDX(iii+1,1)<0)
        % 1 singularity
        roiDX(iii+1,1)=mean(roiDX([iii,iii+2],1));
    end
end