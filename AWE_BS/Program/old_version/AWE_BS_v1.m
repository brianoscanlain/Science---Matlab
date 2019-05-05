% Function that finds the most suitable threshold to distinguidh whitecaps,
% using derivative analysis. 
% 

function [PIP,smoothed_2,threshold,...
    GreyIm,processCode,opt,BGSlope,finalThresh,maxPeak,flag]=...
    AWE_BS(imName,imDir,normPix)

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Globals
global  columnCrop rowCrop;

showImages=0;
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&           Allocate Memory                &&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                        Load Image and find W_c and PIP
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
raw=im2double(imread(fullfile(imDir,imName))); 
GreyIm = rgb2gray(raw(columnCrop(1):end-columnCrop(2),rowCrop(1):end-rowCrop(2),:));
BlueIm = raw(columnCrop(1):end-columnCrop(2),rowCrop(1):end-rowCrop(2),3);


       if max(max(raw))>0.3

%find the slope::
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
meanGrey=mean(GreyIm,2);  
[pGrey,~]=polyfit([1:length(meanGrey)]',meanGrey,1);


%meanGreen=mean(GreenIm,2);
%[pGreen,~]=polyfit([1:length(meanGreen)]',meanGreen,1);

meanBlue=mean(BlueIm,2);
[pBlue,~]=polyfit([1:length(meanBlue)]',meanBlue,1);

[~,  minSlopeIndex]=min(abs([pGrey(1) pBlue(1)]));
if minSlopeIndex==1
    %greyscale is most suited
    opt=1;
    Pn=pGrey;        
elseif minSlopeIndex==2

    %blue channel is the least noisest, use this::
    GreyIm=BlueIm;
    opt=3;
    Pn=pBlue;
end


clear BlueIm GreenIm blueSlope greenSlope greySlope mean* p*
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%The following adjusts the columns pixels accoring to the slope of the mean
%pixel value for each row.
if normPix==1
   if Pn < 0
       BGSlope=Pn;
        backG = [1:size(GreyIm,1)]'.*Pn;
        GreyIm_Lin = GreyIm - repmat(backG,1,size(GreyIm,2));
         GreyIm_Lin(GreyIm_Lin>1) = 1;
        if max(max(GreyIm_Lin)) > 1
            GreyIm_Lin =  GreyIm_Lin - (max(max(GreyIm_Lin))-1);
            GreyIm_Lin(GreyIm_Lin < 0) = 0;
        end
        
        GreyIm = GreyIm_Lin;
    else
        %slope is positive, do nothing.
        BGSlope=0;
    end
end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%       Make PIP and threshold arrays::
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   [threshold,PIP] = makePip_BS(GreyIm,0.3);
   smCoef = 0.1;
   res=0.005;
   possWCap=1;
   flag=NaN;
   
[~,smoothed_2,firstDeriv,secondDeriv,~,~,threshold,PIP] = derivAnalysis_BS(threshold,PIP,smCoef);

            
            %findPeaksTroughs returns a 2-column array.
            %Col1 is threshold
            %Col2 is value at that threshold
            [firstPeaks,firstTroughs,peakChange1] = findPeaksTroughs_aje(firstDeriv,threshold,secondDeriv,1,res);
            %Now if the firstPeaks is empty, then assume image has no whitecaps
            %and skip all the decision tree.
            
            if ~isempty(firstPeaks)
                [secondPeaks,secondTroughs,peakChange2] = findPeaksTroughs_aje(secondDeriv,threshold,secondDeriv,2,res);
                %The lowest threshold in firstPeaks corresponds to the largest
                %positive peak of the first derivative of the PIP. This value is
                %used as the baseline below which all other potential thresholds
                %are deemed useless.
                %Now find the neg2pos zero crossings of the first derivative in
                %decreasing threshold sense above the minimum threshold.
                
                if length(firstPeaks(1,1)) > 1
                    indexNP = find(threshold > firstPeaks(1,1) & threshold < firstPeaks(end,1));
                else
                    indexNP = find(threshold > firstPeaks(1,1));
                end
                firstNegPosCross = [];
                if ~isempty(indexNP)
                    if indexNP(1) == 1
                        indexNP2 = find(firstDeriv(1:indexNP(end)) < 0);
                    else
                        indexNP2 = find(firstDeriv(indexNP(1):indexNP(end)) < 0)+indexNP(1)-1;
                    end
                    if ~isempty(indexNP2)
                        firstNegPosCross = [threshold(indexNP2(end)),firstDeriv(indexNP2(end))];
                    end
                end
               
                %Now find the first positive to negative crossing after the
                %negative to positive crossing
                firstPosNegCross = indexNP2(find(diff(indexNP2)>1)+1);
                if ~isempty(firstPosNegCross)
                    firstPosNegCross = [threshold(firstPosNegCross(end)),firstDeriv(firstPosNegCross(end))];
                elseif isempty(firstPosNegCross);
                    %allPeaks = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [maxFDVal maxFDPos] = max(firstDeriv);
                if maxFDPos < 5
                   %Then probably no whitecaps
                   processCode = 0;
                   finalThresh(1,1:2) = [1,0];
                   possWCap = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                secondPeaks(secondPeaks(:,1) < min(firstPeaks(:,1)),:) = [];
                if ~isempty(secondPeaks) && ~isempty(secondTroughs)
                    firstTroughs(firstTroughs(:,1) < min(secondPeaks(:,1)),:) = [];
                    secondTroughs(secondTroughs(:,1) < min(secondPeaks(:,1)),:) = [];
                    
                else
                    %If there are no second Peaks then more than likely there are
                    %no whitecaps so empty firstTroughs.
                    firstTroughs = [];
                    secondTroughs = [];
                end
                

                
                %Assume that a second maximum peak will be found.  This peak
                %typically has a much larger threshold than the final chosen
                %threshold.It somewhat represents the whitecap contribution
                %without streaks.
                flag = 0;
                %Now choose the threshold
                if ~isempty(firstPeaks) && possWCap == 1 && ~isempty(firstTroughs)
                    
                    finalThresh(1,1) = min(firstTroughs(:,1));
                    finalThresh(1,2) = calcW_BS(GreyIm,finalThresh(1,1),showImages);
                    derivVal = firstDeriv(threshold == finalThresh(1,1));
                    processCode = 1;
                    if derivVal < 0 && ~isempty(firstNegPosCross) && min(firstNegPosCross(:,1)) < finalThresh(1,1)
                        start = bsearch(threshold,.9);
                        finish = bsearch(threshold,min(firstNegPosCross(:,1)));
                        posLen = length(find(firstDeriv(start:finish,1)>= 0));
                        negLen = length(find(firstDeriv(start:finish,1)< 0));
                        if negLen >= posLen
                            %Then choose the neg2pos zero crossing of the first
                            %derivative as the threshold
                            finalThresh(1,1) = min(firstNegPosCross(:,1));
                            finalThresh(1,2) = calcW_BS(GreyIm,finalThresh(1,1),showImages);
                            processCode = 2;
                        end
                    end
                else
                    %No whitecaps
                    processCode = 0;
                    finalThresh(1,1:2) = [1,0];
                end 
                
                if isempty(firstTroughs) && ~isempty(firstPeaks) && possWCap == 1 ...
                        && ~isempty(secondPeaks) && processCode == 0
                    %Added for CE image. 07/12/2010 AC.
                    %Scenario when there are no troughs in the first
                    %derivative, but there is a well-defined peak in the
                    %first derivative.  The reason that there are no
                    %troughs in the first derivative is probably that the
                    %image is bright so Region 1 in the PIP does not
                    %contain many points.
                    finalThresh(1,1) = max(secondPeaks(:,1));
                    finalThresh(1,2) = calcW_BS(GreyIm,finalThresh(1,1),showImages);
                    processCode = 3;
                end
                
                %Now provide the max peak as an alternative threshold that could
                %exclude streaks
                if length(firstPeaks(:,1)) > 1 && processCode ~= 0 %%%%%%%%%%%%%%%%%%%%%%%%%%
                    maxPeak(1,1) = firstPeaks(end,1);
                    maxPeak(1,2) = calcW_BS(GreyIm,maxPeak(1,1),showImages);
                else
                    %No max peak
                    maxPeak(1:2) = [1,0];
                    if processCode ~= 0
                        %Flag is set to 1 if no additional maximum peak is found,
                        %but a threshold is found.
                        %Otherwise if it is found it's value is 0.
                        flag = 1;
                    end
                end
               
            else
                %Here if findPeaksTroughs returned an empty firstPeaks array
                %This means that there were no whitecaps

                finalThresh = [1,0];
                maxPeak(1:2) = [1,0];
                processCode = 10;
            end
        else
            %we are here if the image wasn't read properly or image was too
            %dark indicating a nighttime image
            finalThresh = [NaN,NaN];
            maxPeak = [NaN,NaN];
            flag = NaN;
            opt = NaN;
            processCode = NaN;
            disp('Image Not Read Properly Or else too dark');
        end

