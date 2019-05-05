function [LocalMax,LocalMin,peakChange] = AWE_BS_deriv_analysis(...
    sdx,s2dx,thres,deriv,filtT,AC,threshLM)
%UNTITLED2 Summary of this function goes here
%   Finds the peak local maximum and minimum of sdx and the remaining maxs
%   and mins in the array. Max and mins before the local max are not
%   recorded
%Globals
global I_max   cutOff;
% sdx=PIP_sdx;
% s2dx=PIP_s2dx;
% thres=threshold;
% deriv=1;

% sdx=PIP_s2dx;
% s2dx=PIP_s2dx;
% thres=threshold;
% deriv=1;

MaxCount=0;
MinCount=0;
peakChange=0;
LocalMax=[];
LocalMin=[];
%cutOff=sum(isnan(sdx))/2; %now defined in global variables of master file.

% Find all maximums and minimums %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AC==1
for i=4:length(thres)-1
    if (sdx(i,1)==sdx(i-1,1) || sdx(i,1) > sdx(i-1,1)) && sdx(i,1) > ...
            sdx(i+1,1) && mean(abs(s2dx(i-1:i+1,:)))>filtT %&& mean((sdx(i-1:i+1,:)))>filtT
            %sdx(i+1,1) && mean(abs(s2dx(i-1:i+1,:)))>0.1   % '.>0' for 0.01 resolution, '>0.01' for 0.005 resolution
        %the mean(abs.. above can be used to filter out any peaks close to
        %zero.
        MaxCount=MaxCount+1;
        LocalMax(MaxCount,[1 2])=[thres(i,1),sdx(i)];
    end
    if (sdx(i,1) == sdx(i-1,1) || sdx(i,1) < sdx(i-1,1)) && sdx(i) < sdx(i+1)
        MinCount=MinCount+1;
        LocalMin(MinCount,[1 2])=[thres(i),sdx(i)];
    end
end
elseif AC==0
    for i=4:length(thres)-1
    if (sdx(i,1)==sdx(i-1,1) || sdx(i,1) > sdx(i-1,1)) && sdx(i,1) > ...
            sdx(i+1,1) && mean(abs(s2dx(i-1:i+1,:)))>filtT && mean((sdx(i-1:i+1,:)))>threshLM
            %sdx(i+1,1) && mean(abs(s2dx(i-1:i+1,:)))>0.1   % '.>0' for 0.01 resolution, '>0.01' for 0.005 resolution
        %the mean(abs.. above can be used to filter out any peaks close to
        %zero.
        MaxCount=MaxCount+1;
        LocalMax(MaxCount,[1 2])=[thres(i,1),sdx(i)];
    end
    if (sdx(i,1) == sdx(i-1,1) || sdx(i,1) < sdx(i-1,1)) && sdx(i) < sdx(i+1)
        MinCount=MinCount+1;
        LocalMin(MinCount,[1 2])=[thres(i),sdx(i)];
    end
end 
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cutout data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(LocalMin) && ~isempty(LocalMax)
    LocalMin(LocalMin(:,1)>230,:)=[];
    LocalMax(LocalMax(:,1)>230,:)=[];

if deriv==1
[maxPeakVal,maxPeakPos]=max(LocalMax(:,2));
LocalMax(LocalMax(:,1)<LocalMax(maxPeakPos,1),:)=[];   %remove peaks before max peak

if length(LocalMax(:,1))>2
    if abs(LocalMax(3,1) - LocalMax(1,1))>0.2
        peakChange=1; %theres exists a large separation between peaks here
    end
end
end
end
% 
% [minPeakVal,minPeakPos]=min(Peaks(:,1));
% Troughs(Troughs(:,1) < (Peaks(minPeakPos,1)),:)=[];

LocalMax=flipud(LocalMax);
LocalMin=flipud(LocalMin);
end

