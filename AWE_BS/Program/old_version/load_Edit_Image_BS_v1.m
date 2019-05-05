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
% Load and crop image.
% Convert to grayscale or take blue or green colour channel only
%

function [editImage,imageNumber,origImage,colEditImage,BGSlope,rowIndex,...
    colIndex] = load_Edit_Image_aje(imagePath,opt)


s=(imagePath(end-16:end-4));

imageNumber =datenum([2011,str2num(s(1:2)),str2num(s(3:4)),str2num(s(5:6))...
    ,str2num(s(7:8)),str2num(s(9:10))]); 

% imageDate = str2num(imagePath(end-21:end-14));
%Load image
%error check the reading of the file
badImage = 0;

% %Check file size
% info = [];
% info = imfinfo(imagePath);
% if info.FileSize < 138000
%     %Then the image is probably too dark to analyse
%     badImage = 1;
% end

%Now make sure the image can be read
try
    origImage = imread(imagePath);
catch
    %Here if the image cannot be read.
    badImage = 1;
end

if badImage == 0
    %Image size
    [imR imC] = size(origImage(:,:,1));
    
    %Parts of image to be analysed
    rowIndex = [100 imR-50];
    colIndex = [25 imC-25];
    
    %
    if opt == 1
        %Option 1.....grayscale
        %Convert to double precision and crop
        editImage = rgb2gray(origImage);     %Convert to grayscale
        editImage = im2double(editImage(rowIndex(1):rowIndex(2),colIndex(1):colIndex(2)));
        
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
    colEditImage = [];
    
%Added correction::
[numRows numCols] = size(editImage);    
%Normalise the image using the BGSlope::
meanBG=mean(editImage,2);  
[BGSlope,~]=polyfit([1:length(meanBG)]',meanBG,1);
BGSlope=BGSlope(1);
     
m = BGSlope;   
    %Remove Slope
    if m < 0
        backG = [1:numRows]'.*m;
        editImage_Lin = editImage - repmat(backG,1,numCols);
        editImage_Lin(editImage_Lin>1) = 1;
        if max(max(editImage_Lin)) > 1
            editImage_Lin =  editImage_Lin - (max(max(editImage_Lin))-1);
            editImage_Lin(editImage_Lin < 0) = 0;
        end
        
        fixImage = editImage_Lin;
    else
        fixImage = editImage;
    end
    
editImage=fixImage;

else
    %error reading in the image
    origImage = [];
    colEditImage = [];
    editImage = [];
    BGSlope=[];
end

