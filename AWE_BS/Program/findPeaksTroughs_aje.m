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
%Script to find maxima and minima of a series of points

function [allPeaks,allTroughs,peakChange] = findPeaksTroughs_aje(timeSeries,threshold,secDeriv,deriv,res)
%Threshold is descending from max to min
lenTS = length(timeSeries);
peakCount = 0;
troughCount = 0;
peakChange = 0;
allPeaks = [];
allTroughs = [];
%This is to do with the resolution that is chosen for the PIP generation
%Only introduce the threshold for the identification of first derivative
%peaks for the finer resolution as there are commonly more peaks at this
%resolution.  Basically is setting a noise floor for the "strength" of the
%potential peak.  It must be greater than the value of secondPeakThresh.
if res == 0.005
      secondPeakThresh = 0.01;
elseif res == 0.01
    secondPeakThresh = 0;
end
                           
 for uu = 4:lenTS-1
    %Find time-series' peaks
    if (timeSeries(uu) == timeSeries(uu-1) || timeSeries(uu) > timeSeries(uu-1))...
            && timeSeries(uu) > timeSeries(uu+1) && mean(abs(secDeriv(uu-1:uu+1))) > secondPeakThresh
        peakCount = peakCount+1;
        allPeaks(peakCount,[1 2]) = [threshold(uu),timeSeries(uu)];
    end
    %Find time-series' troughs
    if (timeSeries(uu) == timeSeries(uu-1) || timeSeries(uu) < timeSeries(uu-1))...
            && timeSeries(uu) < timeSeries(uu+1)
        troughCount = troughCount+1;
        allTroughs(troughCount,[1 2]) = [threshold(uu),timeSeries(uu)];
    end
end

if ~isempty(allPeaks) && ~isempty(allTroughs)
    %Cut out unnecessary values above thresholds of 0.95 and below 0.5
    allPeaks(allPeaks(:,1) > .97,:) = [];
    %allPeaks(allPeaks(:,1) < 0.5,:) = [];
    allTroughs(allTroughs(:,1) > .97,:) = [];
    %allTroughs(allTroughs(:,1) < .5,:) = [];
    
    
    if deriv == 1
         [maxPeakVal,maxPeakPos] = max(allPeaks(:,2)) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Uncomment this to get rid of pictures that have likely no
        %%whitecaps (maxPeakVal < 0.5 suggests no whitecaps)
        %if maxPeakVal < 0.5 && maxPeakPos(allPeaks(:,2)) < 1
         %  allPeaks = [];
        %end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %[maxPeakVal,maxPeakPos] = max(allPeaks(:,2));
        allPeaks(allPeaks(:,1) < allPeaks(maxPeakPos,1),:) = [];
        %If there are 2 or more peaks after the maximum peak, then only take
        %the last two peaks
        
        if length(allPeaks(:,1)) > 2
            if abs(allPeaks(3,1) - allPeaks(1,1)) > 0.2
                peakChange = 1; %To do with cases where there is a very large seperation between peaks
            end
            %        allPeaks = allPeaks(1:2,:);
            %        [maxPeakVal,maxPeakPos] = max(allPeaks(:,2));
        end
        
        %     [minPeakVal,minPeakPos] = min(allPeaks(:,1));
        %     allTroughs(allTroughs(:,1) < (allPeaks(minPeakPos,1)),:) = [];
    end
    allPeaks = flipud(allPeaks);
    allTroughs = flipud(allTroughs);
end

% figure,
% plot(threshold,(timeSeries),'-k');
% hold on
% plot(allPeaks(:,1),allPeaks(:,2),'r.','markersize',20);
% plot(allTroughs(:,1),allTroughs(:,2),'b.','markersize',20);