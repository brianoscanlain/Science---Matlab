%Set two things:: workDir (on line 4) and outDir (on line 9)

clear all; close all; clc;
% workDir='D:\Work\Whitecap\ASIP_dep3';   %set location of images::
workDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images';
cd(workDir);
folders=dir([workDir '\_*']);
folders(([folders(:).isdir]==0))=[];    %only save folders (directories)

%set up the output directory;
outDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\Results';



normPix=1;     % normalise the vertical columns of the image using the slope 
ROI_analysis=0;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Global Variable Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
global ButterSample;   %sampling used for the butterworth cuttoff. This will 
ButterSample=0.1;                    %sharpen the step changes in PIP

global rowCrop;       %take 5 pixels off the original image
rowCrop=[100 50];

global columnCrop
columnCrop=[25 25];

global I_max;          %The maximum intensity value of each image, note; that
I_max=253;                    %matlab indices start at 1, so instead of the normal
                    %[0:255] range, it will be succeeded with [1:256] 

global cutOff;       %Crops the final derivative data and replaces outliers 
cutOff=9;                   %with Nan values
                   


for i=1%1:length(folders)
     files=dir([workDir '\' folders(i).name '\*.jpg']);
     cd(outDir);
    if ~exist(folders(i).name, 'dir')
       mkdir(folders(i).name);
    end
    cd([outDir '/' folders(i).name]);
  
    %Allocate memory::
    finalThreshold=ones(length(files),20)./0;
    Option=ones(length(files),20)./0;
    tunedT=ones(length(files),1)*255;
    finalW=zeros(length(files),20);
    tunedW=zeros(length(files),1);
    
for ii=1:length(files)    

    disp([num2str(ii) ' of ' num2str(length(files))])
imName=files(ii).name;
imDir=[workDir '\' folders(i).name '\'];

divAdd=1; filtT=0.1; AC=1; threshLM=1;

[processCode,finalThresh,smoothen1,~,~,~,...
    opt,maxPeak,editImage,threshold,BGSlope] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,1)=finalThresh(1,1);
Option(ii,1)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,1)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));

divAdd=1; filtT=0.8;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,2)=finalThresh(1,1);
Option(ii,2)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,2)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));

divAdd=1; filtT=0.06;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,3)=finalThresh(1,1);
Option(ii,3)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,3)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.04;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,4)=finalThresh(1,1);
Option(ii,4)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,4)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.02;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,5)=finalThresh(1,1);
Option(ii,5)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,5)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
%turn off divAdd
divAdd=0; filtT=0.1;

[processCode,finalThresh,smoothen0,~,~,~,...
    ~,~,~,threshold0,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,6)=finalThresh(1,1);
Option(ii,6)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,6)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
 filtT=0.08;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,7)=finalThresh(1,1);
Option(ii,7)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,7)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
 filtT=0.06;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,8)=finalThresh(1,1);
Option(ii,8)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,8)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.04;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,9)=finalThresh(1,1);
Option(ii,9)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,9)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.02;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,10)=finalThresh(1,1);
Option(ii,10)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,10)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
%Use the extra addition to filter local Maxes further::
AC=0;
divAdd=1; filtT=0.1;

[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,11)=finalThresh(1,1);
Option(ii,11)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,11)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
 filtT=0.08;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,12)=finalThresh(1,1);
Option(ii,12)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,12)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
 filtT=0.06;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,13)=finalThresh(1,1);
Option(ii,13)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,13)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.04;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,14)=finalThresh(1,1);
Option(ii,14)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,14)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
filtT=0.02;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,15)=finalThresh(1,1);
Option(ii,15)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,15)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
%Use the extra addition to filter local Maxes further again::
AC=0;
divAdd=0; filtT=0.6; threshLM=-0.5;

[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,16)=finalThresh(1,1);
Option(ii,16)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,16)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
threshLM=-0.2;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,17)=finalThresh(1,1);
Option(ii,17)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,17)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
threshLM=0;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,18)=finalThresh(1,1);
Option(ii,18)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,18)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
threshLM=0.2;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,19)=finalThresh(1,1);
Option(ii,19)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,19)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));
threshLM=0.5;
[processCode,finalThresh,~,~,~,~,...
    ~,~,~,~,~] = AWE_BS(imName,imDir,normPix,divAdd,filtT,AC,threshLM);
finalThreshold(ii,20)=finalThresh(1,1);
Option(ii,20)=processCode;
dispImage=editImage>finalThresh(1,1);
finalW(ii,20)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                Full on Developer Noob mode Engaged!
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
                     t=t-1;
                  elseif adj == 30 && t<255    %     up arrow key
                     t=t+1;
                  elseif adj == 29     %  right arror key
                      if t<241
                      t=t+15;
                      else
                          t=255;
                      end
                  elseif adj == 28 && t>14     %   left arrow key
                     t=t-15;   
                  elseif adj == 13   %hit return key to finish the loop
                     statuss==0; close all; break;
                  end        
            end
            
%get some info about the background pixels::
BGPix_std(ii,1)= std(double(editImage(~dispImage)));
BGPix_mean(ii,1)= mean(double(editImage(~dispImage)));    

 plot(threshold,smoothen1,'k');
set(gcf,'OuterPosition', [20 20 1920 1080])
hold on;
plot((threshold0),smoothen0,'r')

plot(finalThreshold(ii,1),smoothen1(bsearch(threshold,finalThreshold(ii,1)),1),'*r')
plot(finalThreshold(ii,2),smoothen1(bsearch(threshold,finalThreshold(ii,2)),1),'*g')
plot(finalThreshold(ii,3),smoothen1(bsearch(threshold,finalThreshold(ii,3)),1),'*b')
plot(finalThreshold(ii,4),smoothen1(bsearch(threshold,finalThreshold(ii,4)),1),'*c')
plot(finalThreshold(ii,5),smoothen1(bsearch(threshold,finalThreshold(ii,5)),1),'*m')
plot(finalThreshold(ii,6),smoothen1(bsearch(threshold,finalThreshold(ii,6)),1),'sr')
plot(finalThreshold(ii,7),smoothen1(bsearch(threshold,finalThreshold(ii,7)),1),'sg')
plot(finalThreshold(ii,8),smoothen1(bsearch(threshold,finalThreshold(ii,8)),1),'sb')
plot(finalThreshold(ii,9),smoothen1(bsearch(threshold,finalThreshold(ii,9)),1),'sc')
plot(finalThreshold(ii,10),smoothen1(bsearch(threshold,finalThreshold(ii,10)),1),'sm')
plot(finalThreshold(ii,11),smoothen1(bsearch(threshold,finalThreshold(ii,11)),1),'vr')
plot(finalThreshold(ii,12),smoothen1(bsearch(threshold,finalThreshold(ii,12)),1),'vg')
plot(finalThreshold(ii,13),smoothen1(bsearch(threshold,finalThreshold(ii,13)),1),'vb')
plot(finalThreshold(ii,14),smoothen1(bsearch(threshold,finalThreshold(ii,14)),1),'vc')
plot(finalThreshold(ii,15),smoothen1(bsearch(threshold,finalThreshold(ii,15)),1),'vm')
plot(finalThreshold(ii,16),smoothen1(bsearch(threshold,finalThreshold(ii,16)),1),'^r')
plot(finalThreshold(ii,17),smoothen1(bsearch(threshold,finalThreshold(ii,17)),1),'^g')
plot(finalThreshold(ii,18),smoothen1(bsearch(threshold,finalThreshold(ii,18)),1),'^b')
plot(finalThreshold(ii,19),smoothen1(bsearch(threshold,finalThreshold(ii,19)),1),'^c')
plot(finalThreshold(ii,20),smoothen1(bsearch(threshold,finalThreshold(ii,20)),1),'^m')

%update threshold values            
finalThresh(1,:)=[t sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2))]; %update the choices
tunedT(ii,1)=t; 
dispImage=editImage>t;
tunedW(ii,1)=sum(sum(dispImage)).*100./(size(editImage,...
                1)*size(editImage,2));

%plot the new threshold
plot(finalThresh(1,1),smoothen1(bsearch(threshold,finalThresh(1,1)),1),'ok')
%legend and titles
xlabel('thresholds','fontsize',15); legend('PIP1','PIP0','.05 1','.02 1'...
    ,'.01 1','.005 1','0 1','0.05 0','0.02 0', '0.01 0','0.005 0','0 0',...
    '0.05 BS','0.02 BS','0.01 BS','0.005 BS','0 BS','0.5 0 BS -0.5',...
    '0.5 0 BS -0.2','0.5 0 BS 0','0.5 0 BS 0.2','0.5 0 BS 0.2',...
    '0.5 0 BS 0.5','tuned t');

title(num2str(finalThreshold(ii,:)));
set(gca,'Xlim',[0 255]);
while waitforbuttonpress~=1; end       

    close all;
    
end

save results finalThreshold finalW tunedT Option BG* tunedW;

end
    
