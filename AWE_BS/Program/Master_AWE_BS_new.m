 %Set two things:: workDir (on line 4) and outDir (on line 9)

clear all; close all; clc;
%workDir='D:\Work\Whitecap\ASIP_dep3';   %set location of images::
% workDir='C:\Users\Administrator\Desktop\sea2sky\soap';
workDir='C:\Users\Brian\Desktop\TomBell';
cd(workDir); 
folders=dir([workDir '\S*']);
folders(([folders(:).isdir]==0))=[];    %only save folders (directories)

%set up the output directory;
% outDir='C:\Users\Administrator\Desktop\sea2sky\SOAP_Results\Whitecap_results';
outDir='C:\Users\Brian\Desktop\TomBell\SOAPresults';


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%              Set the processing to developer mode
%
%  The developer mode allows you to view the programs threshold values, and
%  adjust them by using the following keys::
%
%  up    arrow key  == threshold increases by 1
%  down  arrow key  == threshold decreases by 1
%  left  arrow key  == threshold decreases by 15
%  right arrow key  == threshold increases by 15
%
%  the keypad 'zero' was hit, set t=1 and break
%  the 'n' key hit, toggle normPix
%    
%
%  return key continues to the final figure
%
%  The final figure shows the PIP, its 1st, 2nd, 3rd & 4th derivatives. It 
%  also shows the local max and min points on the 1st and 2nd derivatives, 
%  the old and new threshold values
%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
normPix=1;     % normalise the vertical columns of the image using the slope 
%normPix=0;

tuning = 1;    %tuning on               %manual tuning of the threshold
derivShow=0;   %show the plot of the PIP, its derivatives and the thresholds
WaWb_analysis=1; %separation of the 
writing=1;     %save files etc
lightSave=1;  %doesn't save tiff images, and the rgb 
ROI_analysis=0;
step=1;       %The interval between sucessive images to be processed
minThreshold=0.2;    %Minimum threshold. Can be set to set a minimum limit for T

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Global Variable Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
global rowCrop;   %sampling used for the butterworth cuttoff. This will 
%rowCrop=[100 50];                   %sharpen the step changes in PIP
rowCrop=[5 5];

global columnCrop;       %take 5 pixels off the original image
%columnCrop=[25 25];
columnCrop=[5 5];       
       

for i=1:length(folders)

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
files=dir([workDir '\' folders(i).name '\*.jpg']);         %find images
cd(outDir);
if ~exist(folders(i).name, 'dir')                          %set output path
   mkdir(folders(i).name);
end
cd([outDir '/' folders(i).name]);    
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if exist('WC.txt')==2                          %if WC.txt exists, resume.
    WC=load('WC.txt');
    start=WC(end,1)+step;                      %reads WC, and resumes proc.
    clear WC;
else
    start=1;                                   %Starts at 1st image    
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   for ii=start:step:length(files)
   %ii=xxx(iii);
   %for ii= 26901:50:34350%:2:length(files)
   
 clc;       
disp([num2str(ii) ' of ' num2str(length(files))]);
disp([num2str(ceil(length(files(ii:step:end)))) ' images left to process']);
imName=files(ii).name;
imDir=[workDir '\' folders(i).name];

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                Find PIP, threshold values etc using AWE::
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for d=1:10
[PIP,smoothed_2,threshold,editImage,processCode,opt,BGSlope,finalThresh,...
    maxPeak,flag]=AWE_BS(imName,imDir,normPix,minThreshold);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                Full on Developer Noob mode Engaged!
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if tuning ==1 && ~isnan(processCode)
            statuss=1;
            deleteBtn=0;
            t=finalThresh(1,1);
 
            while statuss==1 
                  %Plot the original and BW image
                  dispImage=editImage>t;
                  subplot(1,2,1),h1=imshow(editImage);
                  title([imName(1:end-7) ' # = ' num2str(ii) ' , Code = ' num2str(processCode),' opt= ' num2str(opt) ' normPix = ' num2str(normPix)]);
                  subplot(1,2,2);
                  h2=imshow(dispImage);
                  title(['threshold = ' num2str(t)...
                      ' ,   W(%) = ' num2str(sum(sum(dispImage)).*100./...
                      (size(editImage,1)*size(editImage,2))) ...
                      '  w/ streaks t = ' num2str(maxPeak(1,1))]);

set(gcf,'OuterPosition', [20 250 1910 980])

                  
                  %Wait for a button to be pressed
                  while waitforbuttonpress~=1; end
                  adj=double(get(gcf,'CurrentCharacter'));
                  if adj == 31 && t>0    %         down arror key
                     t=t-.005;
                  elseif adj == 30 && t<1    %     up arrow key
                     t=t+.005;
                  elseif adj == 29     %  right arror key
                      if t<.9
                      t=t+.1;
                      else
                          t=1;
                      end
                  elseif adj == 28 && t>.1     %   left arrow key
                     t=t-0.1;   
                     
                  elseif adj == 13 || adj ==32  %hit return/space key to finish the loop
                     statuss=0; close all; break;
                     
                  elseif adj==48   %the keypad 'zero' was hit, set t=1 and break
                     dispImage=editImage>1;
                     statuss=0; close all; break;
                     
                  elseif adj==127   %the 'Delete' was hit, break and don't write and skip to next image
                     processCode=0/0; deleteBtn=1; 
                     statuss=0; close all; break;
                     
                  elseif adj==110  %the 'n' key hit, toggle normPix
                      
                      if normPix==1
                          normPix=0;
                      else
                          normPix=1;
                      end
                      close all; break;
                      
                      
                  end        
            end

 end
 
 if statuss==0   %break d loop
     break;
 end

end    %d loop  (used in case normPix is changed during tuning mode:
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$ developer Mode ended %%%%%%%%%%%%%%%%%%%%%%%%



%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           %Stage A and stage B analysis
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Here I will use my old program as a function distinguish the whitecap
%pixels as either stage-A or stage-B::
if ~isnan(processCode) && WaWb_analysis==1
    close all
[Wa,Wb,W,PixIdx_A,PixIdx_B,rgb,fill,deleteBtn]=Wa_Wb_program(editImage,dispImage);


else 
    W=0; Wa=0; Wb=0; PixIdx_A=[]; PixIdx_B=[]; 
end

%write the W values to a 'WC.txt' file 
if writing==1 && deleteBtn==0
fid = fopen('WC.txt', 'a');
fprintf(fid, '%d %s        %6.6f          %12.6f          %18.6f          \r\n', ii, files(ii).name(1:end-4),Wa,Wb,W);
fclose(fid);
%save the regionprops PixIdx data to file
  if lightSave==1
      save(files(ii).name(1:end-4),'PixIdx_A', 'PixIdx_B','t','rowCrop','columnCrop','fill');
  else  
      %save rgb image
      imwrite(rgb, [files(ii).name(1:end-4) '.tiff']);
      save(files(ii).name(1:end-4),'PixIdx_A', 'PixIdx_B','rgb','t','rowCrop','columnCrop','fill');
  end
  end

clear Wa Wb W PixIdx* rgb
    end
    
    

  
    
    
end
    
