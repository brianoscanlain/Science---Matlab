%Script to obtain the greyscale values of each image
 clear all;
 
cd('C:\Users\Brian\Dropbox\knorr_k\whitecap\BWARD_request_2013\Raw_data\Raw_images');
raw=pwd;
cd('C:\Users\Brian\Dropbox\knorr_k\whitecap\BWARD_request_2013\Raw_data');
proc=pwd;



folders=dir('P*');
folders([folders.isdir]==0)=[];

for i=1%:length(folders);
    
    cd(folders(i).name)
    
    %load threshold data
    %threshName=[folders(i).name(1:9) '-' folders(i).name(11:14) '_QC2.mat'];
   % load(threshName);
    %cd('processed_colour');
    files=dir('*.tiff');
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
   %Allocate memory 
   stageA_Pix_intensities=zeros(length(files),256);
   stageB_Pix_intensities=zeros(length(files),256);
   threshold=zeros(length(files),1);
   name = repmat(char(0),length(files),13);
   ROI=struct('stageA',{},'stageB',{});
   rowCrop=[ones(length(files),1) zeros(length(files),1)];
   columnCrop=[ones(length(files),1) zeros(length(files),1)];
   %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
    
        for ii=1:length(files)
            cd([proc filesep folders(i).name]);
%             cd ..; cd('processed_colour');
            x_proc=im2double(imread(files(ii).name));
            if size(x_proc,3) <3
               x_proc(:,:,2)=x_proc(:,:,1);
               x_proc(:,:,3)=x_proc(:,:,1); 
            end
            
            x2=x_proc(:,:,2);   %green colour (for identifying stage B)
            x3=x_proc(:,:,3);   %blue colour (for identifying stage A)
            

            stageA=[];
            stageB=[];
            clc;
            disp([num2str((ii/length(files))*100) '% ...' num2str(i) '/' num2str(length(folders))]);
            

c3=bwconncomp(x3);
c2=bwconncomp(x2);
stageA=regionprops(c3,'PixelIdxList');
stageB=regionprops(c2,'PixelIdxList');


ROI(ii).stageA=stageA;
ROI(ii).stageB=stageB;

            
            
             
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            %load raw image
 
            cd(raw);
            xx=im2double(imread([files(ii).name(1:13) '.jpg']));
            
            %allow for the two cropping options used in this study::
            if size(x_proc,1)==951
               columnCrop(ii,:)=[5 5];
            elseif size(x_proc,1)==911
               columnCrop(ii,:)=[25 25]; 
            end
            
            if size(x2,2)==1271
               rowCrop(ii,:)=[5 5];
            elseif size(x2,2)==1131
                rowCrop(ii,:)=[100 50];
            end

            %crop the raw image to the same as the others::

            xx = rgb2gray(xx(columnCrop(ii,1):end-columnCrop(ii,2),...
            rowCrop(ii,1):end-rowCrop(ii,2),:));
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            %find the most suitable threshold used::
            WPix=[];
            if ~isempty(stageA)
                for j=1:length(stageA)
                    WPix=[WPix;stageA(j).PixelIdxList];
                end
            end
            if ~isempty(stageB)
                for k=1:length(stageB)
                    WPix=[WPix;stageB(k).PixelIdxList];
                end
            end
            
            
            if ~isempty(WPix)
                intensities=xx(WPix);
            else
                intensities=[1];
            end
               
            threshold(ii,1:2)=[min(intensities) length(WPix)*100/(size(xx,1)*size(xx,2))];
            
            
            
%             %find AWE threshold for each image:
%             for tt=1:length(QC2)
%                 if strcmpi(QC2(tt).name,raw.name)>0
%                    
%                     threshold(ii)=QC2(tt).finalValues(1,1);  
%                 end
%             end
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            
            %save the name of each image
            name(ii,1:13)=files(ii).name(1:13);
    dte(ii,1)=datenum(2011,str2num(name(ii,1:2)),str2num(name(ii,3:4)),...
    str2num(name(ii,5:6)),str2num(name(ii,7:8)),str2num(name(ii,9:10)));
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
        end %end of file loop
    
ROI=ROI';        
        %Save the data
        cd('C:\Users\Brian\Desktop\Whitecap_Project\SOAP_Results\thresholds');
        save(folders(i).name,'threshold','name','dte','ROI','rowCrop',...
            'columnCrop'); 
        
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% [dte,srrt]=sort(dte);
% threshold=threshold(srrt);
% name=name(srrt,:);
% ROI=ROI(srrt);
% rowCrop=rowCrop(srrt);
% columnCrop=columnCrop(srrt);
% %    
% ind=dte>=datenum(2011,07,03,06,59,13)&dte<datenum(2011,07,03,09,00,38);
% thresh=threshold(ind);
% datenumber=dte(ind);
% whitecapRegions=ROI(ind);
% RowCrop=rowCrop(ind);
% ColumnCrop=columnCrop(ind);
% name=name(ind,:);
% clear srrt ind;
% cd('C:\Users\Brian\Dropbox\knorr_k\whitecap\BWARD_request_2013\Raw_data\Image_data')
% newDir=pwd;
% cd ..; cd ..;
% load('W_data_deployment_3.mat');
% %sort the W data in ascending order, then crop it to the interval::
% [dte_W,srrt]=sort(dte_W);
% Wa=Wa(srrt);
% Wb=Wb(srrt);
% W=W(srrt);
% ind=dte_W>=datenum(2011,07,03,06,59,13)&dte_W<datenum(2011,07,03,09,00,38);
% dte_W=dte_W(ind);
% W=W(ind);
% Wa=Wa(ind);
% Wb=Wb(ind);
% 
%    ROI(length(Wa)+1:end)=[];
%    Name = repmat(char(0),length(Wa),13);
%    t=zeros(length(Wa),2)./0;
%    rowCrop=[ones(length(Wa),1) zeros(length(Wa),1)];
%    columnCrop=[ones(length(Wa),1) zeros(length(Wa),1)];
% for iiii=1:length(Wa)
%       ROI(iiii).stageA=[];
%    ROI(iiii).stageB=[]; 
% end
%    
% for k=1:length(datenumber)
%    bs=bsearch(dte_W,datenumber(k,1));
%     Name(bs,:)=name(k,:);
%     t(bs)=thresh(k,:);
%     rowCrop(bs,:)=RowCrop(k,:);
%     columnCrop(bs,:)=ColumnCrop(k,:);
%     ROI(bs).stageA=whitecapRegions(k).stageA;
%     ROI(bs).stageB=whitecapRegions(k).stageB;
%     
% end
% 
% t=t(:,1);
% 
% cd(newDir);
% save Bward_data W Wa Wb dte_W t Name ROI rowCrop columnCrop;
% 
% %Move images from folder to other folder::
% origDir='D:\Work\Whitecap\ASIP_dep3\raw_images';
% movDir='C:\Users\Brian\Dropbox\knorr_k\whitecap\BWARD_request_2013\Raw_data\Raw_images'
% cd('C:\Users\Brian\Dropbox\knorr_k\whitecap\BWARD_request_2013\Raw_data\Processed_images');
% file=dir('*.tiff');
% 
% for  u=1:length(file)
%    cd(origDir);
%    xx=imread([file(u).name(1:end-5) '.jpg']);
%    cd(movDir); 
%    imwrite(xx, [file(u).name(1:end-5) '.jpg']);
%    
% end 

end %end of folder loop