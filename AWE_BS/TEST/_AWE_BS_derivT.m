function [PIP_sdx,PIP_s2dx,PIP_s3dx,PIP_s4dx,smoothen2,threshold,...
    editImage,peak,opt,BGSlope]=AWE_BS_derivT(imName,imDir,normPix,divAdd)

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Globals
global I_max columnCrop ButterSample rowCrop;
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&           Allocate Memory                &&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
peak(1,1:2)=[0 0];
% res=0;
W_c=zeros(I_max+3,1);
PIP=zeros(I_max+3,1);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                        Load Image and find W_c and PIP
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
raw=im2double(imread(fullfile(imDir,imName))); 
        editImage = rgb2gray(image(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2),:));
    elseif opt == 3
        editImage = image(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2),3);
    end
    


%find the slope::
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
meanGrey=mean(GreyIm,2);  
[pGrey,sGrey]=polyfit([1:length(meanGrey)]',meanGrey,1);


meanGreen=mean(GreenIm,2);
[pGreen,sGreen]=polyfit([1:length(meanGreen)]',meanGreen,1);

meanBlue=mean(BlueIm,2);
[pBlue,sBlue]=polyfit([1:length(meanBlue)]',meanBlue,1);

[~,  minSlopeIndex]=min(abs([pGrey(1) pGreen(1) pBlue(1)]));
if minSlopeIndex==1
    %greyscale is most suited
    opt=1;
    Pn=pGrey;
        editImage = rgb2gray(origImage);     %Convert to grayscale
        editImage = im2double(editImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2)));
        
elseif minSlopeIndex==2
    %green channel of the image is the least noisiest, use this::
      GreyIm=GreenIm;
      opt=2;
      Pn=pGreen;
elseif minSlopeIndex==3
    %blue channel is the least noisest, use this::
    GreyIm=BlueIm;
    opt=3;
    Pn=pBlue;
end

    if opt == 1

        %colEditImage = im2double(origImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2),:));
    elseif opt == 2%Option 2 ..... Green
        %Convert to double precision and crop
        editImage = origImage(:,:,2);       %Convert to green
        editImage = im2double(editImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2)));
        %colEditImage =im2double(origImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2),:));
    elseif opt == 3
        %Option 3 ..... Blue
        %Convert to double precision and crop
        editImage = origImage(:,:,3);       %Convert to blue
        editImage = im2double(editImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2)));
        %colEditImage = im2double(origImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2),:));
    end

BGSlope=Pn(1); %record the slope for export
clear BlueIm GreenIm blueSlope greenSlope greySlope mean* p*
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%The following adjusts the columns pixels accoring to the slope of the mean
%pixel value for each row.
if normPix==1
    scale=polyval(Pn,[1:size(GreyIm,1)]);
    X=zeros(size(GreyIm,1),size(GreyIm,2));
    centre=scale(round(median(1:size(GreyIm,1))));
    
    for u=1:size(GreyIm,1)
        X(u,:)=scale(u)/centre; %create a matrix of coefficients for normalisation
    end
    
    oldIm=double(GreyIm);  %have to convert to double 
    newIm=oldIm./X;      %adjust the 
    GreyIm=uint8(newIm);
    clear newIm oldIm n X scale centre
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


threshold(1:I_max+3,1)=fliplr([0:I_max+2]);   %all potential thresholds / intensities
%Find the W_c values (as similarly used in Sugihara 2007). These values
%calculate the number of whitecap pixels for all possible threshold values
for i=1:I_max+3
    W_c(i,1)=sum(sum(GreyIm>=256-i));     %pixel index reversed:::
end

W_c(:,1)=W_c(:,1).*100./W_c(256,1);   %convert it to percentage of total pixels

%Calculate the PIP (reference Callaghan 2009)
PIP(2:I_max+3,1)=(diff(W_c)*100);%./W_c(2:end,1);
PIP(1,1)=NaN;

if divAdd==1
%add a minimal fraction onto Wc to stop the PIP reaching 100% values due to
%divisions by zero a few lines below::
W_c=min(W_c(W_c>0))+W_c;
end

%find the highest PIP value:
[peak(1,2) peak(1,1)]=max(PIP);
peak(1,1)=threshold(peak(1,1));
PIP(2:I_max+3,1)=PIP(2:I_max+3,1)./W_c(2:end,1);
clear W_c;

%get rid of any NaN's in the PIP
hh=length(PIP);
for iii=1:hh
    if  isnan(PIP(hh+1-iii,1))
        threshold(hh+1-iii)=[];
        PIP(hh+1-iii)=[];
    end
end
hh=length(PIP);

PIP_dx=zeros(hh,1);PIP_sdx=zeros(hh,1);
PIP_2dx=zeros(hh,1);PIP_s2dx=zeros(hh,1);                   
PIP_3dx=zeros(hh,1);PIP_s3dx=zeros(hh,1);
PIP_4dx=zeros(hh,1);PIP_s4dx=zeros(hh,1);
smooth=zeros(hh+20,1);smoothen=zeros(hh,1);


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                  Filter and smooth the PIP::
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Filter number 1:: 
[but1,but2]=butter(4,ButterSample);
PIP=[flipud(PIP(1:20));PIP];
smooth=filtfilt(but1,but2,PIP);
smooth(1:20)=[];
if  max(threshold)>254
    smooth(1:3)=[];
    threshold(1:3)=[];   %top values removed
    hh=length(smooth);
end
clear PIP;

%Smooth number 1:: 
%using a running mean, with a window of length=4.
smooth=flipud(smooth);
smoothen=zeros(hh,1);
smoothen([1,2,hh-1,hh],1)=smooth([1,2,hh-1,hh],1);
smoothen(3)=mean(smooth(2:4));                 %added with weighting 3:1 :P
for ii=3:hh-2
   smoothen(ii,1)=sum(smooth(ii-2:ii+1))/4; 
end
smoothen(end-1,1)=mean(smoothen(end-2:end));
smooth=flipud(smooth); smoothen=flipud(smoothen);
%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if smoothen(end-1,1)<smoothen(end,1) && smoothen(end-1,1)<smoothen(end-2)
    smoothen(end-1,1)=mean([smoothen(end,1);smoothen(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if smoothen(end-1,1)>smoothen(end,1) && smoothen(end-1,1)>smoothen(end-2)
    smoothen(end-1,1)=mean([smoothen(end,1);smoothen(end-2,1)]);
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% %Filter number 2::
[but1,but2]=butter(6,0.1);
smooth2=filtfilt(but1,but2,smoothen);
%Smooth number 2::
smooth2=flipud(smooth2);
smoothen2=zeros(length(smooth2),1);
smoothen2([1,2,hh-1,hh],1)=smooth2([1,2,hh-1,hh],1);
for ii=3:hh-2
   smoothen2(ii,1)=sum(smooth2(ii-2:ii+1))/4; 
end
smoothen2(end-1,1)=mean(smoothen2(end-2:end));
smooth2=flipud(smooth2); smoothen2=flipud(smoothen2);
%Rare occurence when we get a dip in the smoothen2 at the end from this
%smoothening process, we correct it:: :P
if smoothen2(end-1,1)<smoothen2(end,1) && smoothen2(end-1,1)<smoothen2(end-2)
    smoothen2(end-1,1)=mean([smoothen2(end,1);smoothen2(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if smoothen2(end-1,1)>smoothen2(end,1) && smoothen2(end-1,1)>smoothen2(end-2)
    smoothen2(end-1,1)=mean([smoothen2(end,1);smoothen2(end-2,1)]);
end

% Test all these smoothing and filtering, and see if they are actually necessary. 
%Was he just being over-cautious here or does it actually make difference
%
% plot(PIP,'k')
% hold
% plot(smooth,'g')
% plot(smoothen,'r')
% plot(smooth2,'b')
% plot(smoothen2,'c')
% legend('PIP','filt','filt,smooth','f,s,f','f,s,f,s')



%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$      1st DERIVATION OF PIP                          %%%%%%%%%%%%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Now find the 1st derivative of PIP (smoothed and filtered)::
PIP_dx=[0;diff(smoothen2)];
PIP_sdx=flipud(PIP_dx);  %smoothed::   The data is flipped before smoothed
for ii=3:hh-2
   PIP_sdx(ii,1)=sum(PIP_dx(ii-2:ii+1))/4; 
end
PIP_sdx(end-1,1)=mean(PIP_sdx(end-2:end));
%PIP_sdx=flipud(PIP_sdx); %data refilipped
%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if PIP_sdx(end-1,1)<PIP_sdx(end,1) && PIP_sdx(end-1,1)<PIP_sdx(end-2)
    PIP_sdx(end-1,1)=mean([PIP_sdx(end,1);PIP_sdx(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if PIP_sdx(end-1,1)>PIP_sdx(end,1) && PIP_sdx(end-1,1)>PIP_sdx(end-2)
    PIP_sdx(end-1,1)=mean([PIP_sdx(end,1);PIP_sdx(end-2,1)]);
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear PIP_dx;


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$                Filter to remove any singularities        $$$$$$$$$$
%$$$$$$                where d2()/dx will be Nan                 $$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for ii=1:hh-3
    if (PIP_sdx(ii,1)<0 && PIP_sdx(ii+1,1)>0 && PIP_sdx(ii+2,1)>0 && ...
            PIP_sdx(ii+3,1)<0) || (PIP_sdx(ii,1)>0 && PIP_sdx(ii+1,1)<0 &&...
            PIP_sdx(ii+3,1)<0 && PIP_sdx(ii+3,1)>0)
       %a two value peak/singularity 
       PIP_sdx(ii+1:ii+2,1)=[mean(PIP_sdx([ii,ii+3],1));mean(PIP_sdx([ii...
           ,ii+3],1))];
    end
end
for iii=2:hh-1
   
    if (PIP_sdx(iii,1)<0 && PIP_sdx(iii-1,1)>0 && PIP_sdx(iii+1,1)>0) || ...
            (PIP_sdx(iii,1)>0 && PIP_sdx(iii-1,1)<0 && PIP_sdx(iii+1,1)<0)
        % 1 singularity
        PIP_sdx(iii+1,1)=mean(PIP_sdx([iii,iii+2],1));
    end
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&           2nd Derivative of PIP    &&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Now find the 2st derivative of PIP (smoothed)::
PIP_2dx=[0;diff(PIP_sdx)];
PIP_s2dx=flipud(PIP_2dx);  %smoothed::
for ii=3:hh-2
   PIP_s2dx(ii,1)=sum(PIP_2dx(ii-2:ii+1))/4; 
end
PIP_s2dx(end-1,1)=mean(PIP_s2dx(end-2:end));

%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if PIP_s2dx(end-1,1)<PIP_s2dx(end,1) && PIP_s2dx(end-1,1)<PIP_s2dx(end-2)
    PIP_s2dx(end-1,1)=mean([PIP_s2dx(end,1);PIP_s2dx(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if PIP_s2dx(end-1,1)>PIP_s2dx(end,1) && PIP_s2dx(end-1,1)>PIP_s2dx(end-2)
    PIP_s2dx(end-1,1)=mean([PIP_s2dx(end,1);PIP_s2dx(end-2,1)]);
end

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear PIP_2dx;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$                Filter to remove any singularities        $$$$$$$$$$
%$$$$$$                where d3()/dx will be Nan                 $$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for ii=1:hh-3
    if (PIP_s2dx(ii,1)<0 && PIP_s2dx(ii+1,1)>0 && PIP_s2dx(ii+2,1)>0 && ...
            PIP_s2dx(ii+3,1)<0) || (PIP_s2dx(ii,1)>0 && PIP_s2dx(ii+1,1)<0 &&...
            PIP_s2dx(ii+3,1)<0 && PIP_s2dx(ii+3,1)>0)
       %a two value peak/singularity 
       PIP_s2dx(ii+1:ii+2,1)=[mean(PIP_s2dx([ii,ii+3],1));mean(PIP_s2dx([ii...
           ,ii+3],1))];
    end
end
for iii=2:hh-1
   
    if (PIP_s2dx(iii,1)<0 && PIP_s2dx(iii-1,1)>0 && PIP_s2dx(iii+1,1)>0) || ...
            (PIP_s2dx(iii,1)>0 && PIP_s2dx(iii-1,1)<0 && PIP_s2dx(iii+1,1)<0)
        % 1 singularity
        PIP_s2dx(iii+1,1)=mean(PIP_s2dx([iii,iii+2],1));
    end
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&           3rd Derivative of PIP    &&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Now find the 2st derivative of PIP (smoothed)::
PIP_3dx=[0;diff(PIP_s2dx)];
PIP_s3dx=flipud(PIP_3dx);  %smoothed::
for ii=3:hh-2
   PIP_s3dx(ii,1)=sum(PIP_3dx(ii-2:ii+1))/4; 
end
PIP_s3dx(end-1,1)=mean(PIP_s3dx(end-2:end));

%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if PIP_s3dx(end-1,1)<PIP_s3dx(end,1) && PIP_s3dx(end-1,1)<PIP_s3dx(end-2)
    PIP_s3dx(end-1,1)=mean([PIP_s3dx(end,1);PIP_s3dx(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if PIP_s3dx(end-1,1)>PIP_s3dx(end,1) && PIP_s3dx(end-1,1)>PIP_s3dx(end-2)
    PIP_s3dx(end-1,1)=mean([PIP_s3dx(end,1);PIP_s3dx(end-2,1)]);
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear PIP_3dx;


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&           4th Derivative of PIP    &&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Now find the 2st derivative of PIP (smoothed)::
PIP_4dx=[0;diff(PIP_s3dx)];
PIP_s4dx=flipud(PIP_4dx);  %smoothed::
for ii=3:hh-2
   PIP_s4dx(ii,1)=sum(PIP_4dx(ii-2:ii+1))/4; 
end
PIP_s4dx(end-1,1)=mean(PIP_s4dx(end-2:end));

%Rare occurence when we get a dip in the smoothen at the end from this
%smoothening process, we correct it:: :P
if PIP_s4dx(end-1,1)<PIP_s4dx(end,1) && PIP_s4dx(end-1,1)<PIP_s4dx(end-2)
    PIP_s4dx(end-1,1)=mean([PIP_s4dx(end,1);PIP_s4dx(end-2,1)]);
end
%and vice versa if the opposite occurrs:
if PIP_s4dx(end-1,1)>PIP_s4dx(end,1) && PIP_s4dx(end-1,1)>PIP_s4dx(end-2)
    PIP_s4dx(end-1,1)=mean([PIP_s4dx(end,1);PIP_s4dx(end-2,1)]);
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear PIP_4dx;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&           Normalise and crop data  &&&&&&&&&&&&&&&&&&&&&&&&
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Normalise and crop data
smooth2=smooth2./max(abs(smooth));
smoothen2=smoothen2./max(abs(smoothen));
PIP_sdx=PIP_sdx./max(abs(PIP_sdx));
PIP_s2dx=PIP_s2dx./max(abs(PIP_s2dx));
PIP_s3dx=PIP_s3dx./max(abs(PIP_s3dx));
PIP_s4dx=PIP_s4dx./max(abs(PIP_s4dx));

smooth2=smooth2(13:end-13,1);
smoothen2=smoothen2(13:end-13,1);
PIP_sdx=PIP_sdx(13:end-13,1);
PIP_s2dx=PIP_s2dx(13:end-13,1);
PIP_s3dx=PIP_s3dx(13:end-13,1);
PIP_s4dx=PIP_s4dx(13:end-13,1);
threshold=threshold(13:end-13,1);
editImage=GreyIm;

Pek=bsearch(threshold(:,1),peak(1,1));

% if PIP_sdx(Pek(1,1)-3,1)<0 && PIP_sdx(Pek(1,1)+3,1)>0 && PIP_s2dx(Pek(1,1),1)<0
%     %then the peak(1,1) represents the peak of the PIP.
%     res=1;
%     peak=Pek;
% end
% ind=[1:cutOff,I_max-cutOff+1:I_max]; finalVal= [];%0/0;0/0;0/0;0/0;0/0;0/0;];
% smooth(ind,1)=[];
% smoothen(ind,1)=[];
% PIP_sdx(ind,1)=[];
% PIP_s2dx(ind,1)=[];
% PIP_s3dx(ind,1)=[];
% PIP_s4dx(ind,1)=[];
% threshold(ind,1)=[];
% editImage=GreyIm;
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear finalVal ind;
clear finalVal ind;


