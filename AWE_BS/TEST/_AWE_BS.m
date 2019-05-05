function [pCode,finalThresh,smoothen,PIP_sdx,PIP_s2dx,PIP_s3dx,...
    opt,maxPeak,editImage,threshold,BGSlope] = AWE_BS(imName,imDir,normPix,...
    divAdd,filtT,AC,threshLM)
%Finds the most suitable threshold for a given image, using derivative
%analysis. 
%   imName=image file name.
%   imDir=directory path of image file
%   normPix=the case where the rows of the image are corrected/normalized
%   divAdd= a case which adds a negligible value onto the W_c array when
%       calculating the PIP. This reduces the chance of dividing by one and 
%       creating spikes in the tail end of the resulting PIP
%   filtT= the value which is used as a threshold when distinguishing the
%   local Max and Min vaules of the PIP.
%
%  This function calls on AWE_BS_deriv.m and AWE_BS_deriv_analysis.m

%initial values are set::
pCode = 0;
finalThresh(1,1:2)=[255 0];
        

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Find derivatives, local max and mins too.
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

[PIP_sdx,PIP_s2dx,PIP_s3dx,PIP_s4dx,smoothen,threshold,editImage,peak,opt,...
    BGSlope]=AWE_BS_derivT(imName,imDir,normPix,divAdd);

%Find the Local mins and maxs::

     [LocalMax,LocalMin,peakChange] = AWE_BS_deriv_analysis(...
         PIP_sdx,PIP_s2dx,threshold,1,filtT,AC,threshLM);
     
     
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     
%First filter for images, if its intensity is not over 100, then its too
%dark.
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if max(max(editImage)) <= 100
        rgb=editImage>255;
        editImage = [];
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Purge the LocalMax and Mins below the peak value. The peak value will be
% unchanged from normalisation as a result from the initial PIP calculation
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if ~isempty(LocalMax)
LocalMax=LocalMax(LocalMax(:,1)>peak(1,1),:);
end
if ~isempty(LocalMin)
LocalMin=LocalMin(LocalMin(:,1)>peak(1,1),:);
end




        if ~isempty(editImage)
           if ~isempty(LocalMax)
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Here the results show that whitecap pixels exist in the PIP. A threshold
%to separate background pixels and whitecap pixels is now carried out::
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Get max and mins of 2nd derivatives::
[LocalMax2,LocalMin2,peakChange2]=AWE_BS_deriv_analysis(PIP_s2dx,...
    PIP_s2dx,threshold,2,filtT,AC,threshLM);
           
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Purge the LocalMax and Mins below the peak value. The peak value will be
% unchanged from normalisation as a result from the initial PIP calculation
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if ~isempty(LocalMax2)
LocalMax2=LocalMax2(LocalMax2(:,1)>peak(1,1),:);
end
if ~isempty(LocalMin2)
LocalMin2=LocalMin2(LocalMin2(:,1)>peak(1,1),:);
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% First attempt at finding the threshold. 
%
%Firslty, the 'peak' found is used to equivalently find the threshold
%values of the negative to positive crossing of the 1st derivative (N2P),
%then the P2N of the 2nd derivative, then the N2P of the 3rd derivative
%(N2P2) and finally the finalP2N of the 4th derivative.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

 [N2P,P2N,N2P2, finalP2N] = n2pcrossover(LocalMax,threshold,PIP_sdx,PIP_s2dx,PIP_s3dx,PIP_s4dx,peak );
if ~isempty(N2P)
   pCode=1;
   finalThresh(1,1)=(N2P);
   if ~isempty(P2N)
       pCode=2;
       finalThresh(1,1)=threshold(P2N);
       if ~isempty(N2P2)
           pCode=3;
           finalThresh(1,1)=threshold(N2P2);
           if ~isempty(finalP2N)
               pCode=4;
              finalThresh(1,1)=threshold(finalP2N);
           end
       end
   end
else    %no N2P, P2N, N2P2, finalP2N found::
    pCode=-1;
    finalThresh(1,1)=255;
end
finalThresh(1,2) = sum(sum(editImage>finalThresh(1,1)))*100./(size(editImage,1)*size(editImage,2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%7
% second  threshold using the second peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
                %Now provide the max peak as an alternative threshold that could
                %exclude streaks
                if length(LocalMax(:,1)) > 1 && pCode ~= 0  
                    dispImage = editImage;
                    dispImage(dispImage < (LocalMax(2,1))) = 0;
                    [nr nc] = size(dispImage);
                    numWhite = length(find(dispImage ~= 0));
                    maxPeak(1,1)=LocalMax(2,1);
                    maxPeak(1,2) = numWhite*100 ./(nr*nc);
                    clear dispImage
                    
                else
                    %No max peak
                    maxPeak(1,1:2) = [255,0];

                end
                
                
            else
                %Here if findPeaksTroughs returned an empty LocalMax array
                %This means that there were no whitecaps
                
                finalThresh(1,1:2) = [255,0];
                maxPeak(1,1:2) = [255,0];
                pCode = 10;
            end
            else
            %we are here if the image wasn't read properly or image was too
            %dark indicating a nighttime image
            finalThresh(1,1:2) = [NaN,NaN];
            maxPeak(1,1:2) = [NaN,NaN];
            flag = NaN;
            opt = NaN;
            pCode = NaN;

        end 


        


end

