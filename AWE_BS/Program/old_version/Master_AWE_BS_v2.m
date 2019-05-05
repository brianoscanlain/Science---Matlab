%Set two things:: workDir (on line 4) and outDir (on line 9)

clear all; close all; clc;
%workDir='D:\Work\Whitecap\ASIP_dep3';   %set location of images::
workDir='D:\Work\Whitecap\Extensive_analysis_2014';
cd(workDir);
folders=dir([workDir '\K*']);
folders(([folders(:).isdir]==0))=[];    %only save folders (directories)

%set up the output directory;
outDir='D:\Work\Whitecap\Extensive_analysis_2014\results';


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
ROI_analysis=0;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Global Variable Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
global rowCrop;   %sampling used for the butterworth cuttoff. This will 
%rowCrop=[100 50];                   %sharpen the step changes in PIP
rowCrop=[5 5];

global columnCrop;       %take 5 pixels off the original image
%columnCrop=[25 25];
columnCrop=[5 5];       
       

for i=1%:length(folders)
     files=dir([workDir '\' folders(i).name '\*.jpg']);
     cd(outDir);

    if ~exist(folders(i).name, 'dir')
       mkdir(folders(i).name);
    end
    cd([outDir '/' folders(i).name]);
  
    
%&&&&&&&&&&&&&&&qloyment 3, ii will start with the following:
%    1,6,11,16,21,26,31,36,41,46
%        ^^                                (^ is the current marker)
%    for 7am to 9am::   
%    1,6,11,16,21,26,31,36,41,46
%                             ^^    (^^ currently processed)
%    now do the next batch::

%
%    This will achieve the desired 0.2Hz processing 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%26946
% xxx=[34011:34350];
% logic=ones(1,length(xxx));
% logic(1:50:end)=0;
% logic(6:50:end)=0;
% logic(11:50:end)=0;
% logic(16:50:end)=0;
% logic(21:50:end)=0;
% logic(26:50:end)=0;
% logic(31:50:end)=0;
% logic(36:50:end)=0;
% logic(41:50:end)=0;
% logic(46:50:end)=0;
% xxx=xxx(logic==1);
% disp([num2str(sum(logic)) ' images left to process']);

   for ii=1:2:length(files)
   %ii=xxx(iii);
   %for ii= 26901:50:34350%:2:length(files)
   
 clc;       
disp([num2str(ii) ' of ' num2str(length(files))]);
disp([num2str(length(files(1:2:end))-ii/2) ' images left to process']);

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                Find PIP, threshold values etc using AWE::
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
imName=files(ii).name;
imDir=[workDir '\' folders(i).name];
[PIP,smoothed_2,threshold,editImage,processCode,opt,BGSlope,finalThresh,...
    maxPeak,flag]=AWE_BS(imName,imDir,normPix);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                Full on Developer Noob mode Engaged!
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if tuning ==1 && ~isnan(processCode)
            statuss=1;
            t=finalThresh(1,1);
 
            while statuss==1 
                  %Plot the original and BW image
                  dispImage=editImage>t;
                  subplot(1,2,1),h1=imshow(editImage);
                  title(['# = ' num2str(ii) ' , Code = ' num2str(processCode),' opt= ' num2str(opt) ' normPix = ' num2str(normPix)]);
                  subplot(1,2,2);
                  h2=imshow(dispImage);
                  title(['threshold = ' num2str(t)...
                      ' ,   W(%) = ' num2str(sum(sum(dispImage)).*100./...
                      (size(editImage,1)*size(editImage,2))) ...
                      '  w/ streaks t = ' num2str(maxPeak(1,1))]);

set(gcf,'OuterPosition', [20 20 1920 1080])

                  
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
                  elseif adj == 13   %hit return key to finish the loop
                     statuss==0; close all; break;
                  end        
            end

end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$ developer Mode ended %%%%%%%%%%%%%%%%%%%%%%%%



%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           %Stage A and stage B analysis
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Here I will use my old program as a function distinguish the whitecap
%pixels as either stage-A or stage-B::
if ~isnan(processCode) && WaWb_analysis==1
    close all
[Wa,Wb,W,PixIdx_A,PixIdx_B,rgb,fill]=Wa_Wb_program(editImage,dispImage);


else 
    W=0; Wa=0; Wb=0; PixIdx_A=[]; PixIdx_B=[]; 
end

%write the W values to a 'WC.txt' file 
if writing==1 
fid = fopen('WC.txt', 'a');
fprintf(fid, '%d %s        %6.6f          %12.6f          %18.6f          \r\n', ii, files(ii).name(1:end-4),Wa,Wb,W);
fclose(fid);
%save the regionprops PixIdx data to file
imwrite(rgb, [files(ii).name(1:end-4) '.tiff']);
%save rgb image
save(files(ii).name(1:end-4),'PixIdx_A', 'PixIdx_B','rgb','t','rowCrop','columnCrop','fill');
end

clear Wa Wb W PixIdx* rgb
    end
    
    

    
    
    
end
    
