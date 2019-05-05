function [Wa,Wb,W,PixIdx_A,PixIdx_B,rgb,fill,deleteBtn] = Wa_Wb_program(origImage,BW)
%[Wa,Wb,W,stageA_idx,stageB_idx,colourImage]=Wa_Wb_program(editImage,dispImage,finalThres(1,1));
%   Detailed explanation goes here
%
%
%Taken from the old program::
%    ANAYLSIS INSTRUCTIONS
%--------------------------------------------------------------------------
%1. Run .matfile, if 'all files are now processed' is displayed, please
% redefine the current directory such that it contains the relevant images.
%
%
%2. Two images are displayed; here you have the chance to proceed with
% processing that image or to skip it. Compare the raw and processed images
% for defects. These include reflections from the sun, birds, etc (anything
% that has white pixel values and is neither Stage A nor Stage B). 
%        
%         Click of the mouse == continue  
%           Keyboard Key ==  skip image
%
%
%3. if continue; The two images are now displayed with greater
% magnification. Here it is possible to identify the white areas as either
% Stage A or B by comparing the two images. There are four different
% comands that are executable from here.
%
%       <Q key> == crops the full image
%       <del key> == restarts the image           (press key twice in total)
%       <return key> == exits the cropping loop
%       <any other key> == user defined cropping is small areas
%
%
% Cropping an Image
%------------------
%Click and drag (after hitting <any other key> to define the cropping box.
%Once an image is cropped, the new cropped image is displayed. The User has
%two options;
%
%    Left mouse button click == defines white pixels in box as Stage A
%    Right mouse button click == defines white pixels in box as Stage B
%
%    Keep repeating this until there are no white pixels showing. 
%
% Personally, I crop the stage A regions using <any other key> (using the 
% left click mouse button to define them as Stage A regions) until I am
% left with Stage B. Then I use the <Q key> to crop and count the remaining
% pixels. I then right click to define them as Stage B.
%
%
%Once done analysing the image, The program sums the total white pixel
%counts from the cropping and compares them to a known value (calcWcap). 
%A threshold can be set, I have it set at 0.1% difference. So if error is
%lower than the threshold, the results are printed to a WC.txt file and the
%rgb image is saved to a 'processed' folder. 
%   If the error exceeds the threshold, a error is displayed, showing the
%error % value. After considering the magnitude of the error, the user
%then has two options; 
%               <del key>  == repeat image
%               <return key> == accept error, print results
%
%Source of errors;
%-----------------
%While using the cropping tool, if a perimeter of the rectangle goes
%through a white pixel patch, uncertanty can arise. I will have to
%nvestigate this further.
%
%
%Pausing and resuming;
%--------------------
%Each image is defined a 'count' number. If the variable 'count' doesn't 
%exist, it is set to a value 1. If you process to count= 100, for example,
%and you stop AWE_mod6.m, all the user has to do is either keep 
%matlabs variable 'count' value saved to memory (keep matlab running) or
%specify count = 101 at the command window before running AWE_mod6.m again
%to resume.
finished=0;
a=0;                   %initial stage A pixel count
b=0;                      %initial stage B pixel count
wpixcount= sum(sum(BW));
PixIdx_A=[];
PixIdx_B=[];
fill=[];
deleteBtn=0;
BW_orig=BW;

if wpixcount==0       %empty image, return zero values...
    Wa=0; Wb=0; W=0;
    rgb=BW;
    finished=2;
else                  %whitecap present, prepare rgb image
    BWsize=size(BW);
    numPixels=BWsize(1)*BWsize(2);
    numWhitePix=sum(sum(BW));
    calcWcap = numWhitePix.*100./numPixels;
    rgb=im2uint8(BW);
    rgb(:,:,2)=rgb(:,:,1);
    rgb(:,:,3)=rgb(:,:,1);
end

while finished==0
            close all;
%%%%%%%%%%%%%%%%%%%%%%%%Position windows%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            figure,imshow(origImage);
            scnsize = get(0,'ScreenSize');
            position1 = get(1,'Position');
            outerpos1 = get(1,'OuterPosition');
            borders1 = outerpos1 - position1;
            edge = -borders1(1)/2;
            pos1 = [edge,scnsize(4)*(1/3),scnsize(3)/2-edge,scnsize(4)*2/3];
            pos2 = [scnsize(3)/2 + edge, pos1(2), pos1(3), pos1(4)];

            set(1,'OuterPosition',pos1) ; 
           
            figure,imshow(rgb);
           
            set(2,'OuterPosition',pos2) ;
            clc;
            disp(wpixcount);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       while waitforbuttonpress~=1; end
            adj=double(get(gcf,'CurrentCharacter'));
             if   adj==13           %if <ret key> is hit, the cropping loop 
                                    %ends, keeping the 'a' and 'b' values.
                  
                  break
                  
             elseif adj ==127           %if <del key> is hit, the loop is 
                  a=zeros(1,1);   %forced to end, sets 'a' and b'' to    
                  b=zeros(1,1);   %zero values. This stops fprintf
                  deleteBtn=1;
                  break                 %writing the values to file.  
             elseif adj==113            %press Q key to cropp remaining
                 close all;
                 imshow(BW);
            
                 
                calcWcap1 = length(find(BW == 1));
                %8888888888888888888888888888888888888888888888888888888888
                %PixIdx=regionprops(BW,'PixelIdxList');
           
%                 if isempty(PixIdx)
%                     PixIdx=[];
%                 end
                %888888888888888888888888888888888888888888888888888888888
                
             while waitforbuttonpress~=0; end
                mouseside=get(gcf,'SelectionType');
                abm=sum(mouseside);
                STR1=[110   111   114   109    97   108];   %normal in ASCII
                au=sum(STR1);
                STR2=[97   108   116];                      %alt in ASCII
                bu=sum(STR2);
                
                if abm==au              %left mouse click (stage A)
               % PixIdx_A=[PixIdx_A;PixIdx];
                for m=1:BWsize(1)           %scan all pixels
                    for n=1:BWsize(2)       %and if one is white, colour it
                        if rgb(m,n)==255;    %blue
                           rgb(m,n,1)=0;  
                           rgb(m,n,2)=0;
                           rgb(m,n,3)=255;
                        end
                    end
                end
                
                
                a=a+calcWcap1;
                wpixcount=wpixcount-calcWcap1;
                    
                elseif abm==bu          %right mouse click (stage B)
               % PixIdx_B=[PixIdx_B;PixIdx];
                for m=1:BWsize(1)           %scan all pixels
                    for n=1:BWsize(2)       %and if one is white, colour it
                        if rgb(m,n)==255;    %green
                           rgb(m,n,1)=0;  
                           rgb(m,n,3)=0;
                           rgb(m,n,2)=255;
                        end
                    end
                end                    
                    b=b+calcWcap1;
                    wpixcount=wpixcount-calcWcap1;
                    pause(.2);%stop the issue of right clicking n desktop
                %%%% Case for unwanted white pixels caused by seagulls %%%%
                elseif abm == 648   % If the middle button is hit, white is blacked out
                    for m=1:BWsize(1)           %scan all pixels
                    for n=1:BWsize(2)       %and if one is white, colour it
                        if rgb(m,n)==255;    %green
                           rgb(m,n,1:3)=0;
                        end
                    end
                    end
                    %%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
             break                      
             else                       %image recropping loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                 h=imrect(gca, []);
                 pos = getPosition(h);
                 posX=round(pos-0.5);
                
             x1=posX(1,1);
             x2=posX(1,2);
             x3=posX(1,3);
             x4=posX(1,4);
             if x1<1
                 x3=x3+x1;
                 x1=1; end
             if x2<1
                 x4=x4+x2;
                 x2=1; end
             if x3>BWsize(2)-x1             %boundary conditions for 
                 x3=BWsize(2)-x1; end        %coloring in rgb image
             if x4>BWsize(1)-x2
                 x4=BWsize(1)-x2; end
             posY=[x1 x2 x3 x4];   
                 BW1=imcrop(BW, posY);
                 rgb1=imcrop(rgb, posY); 
            
            %alternative method is to use imrect, get value of the coords,
            %imcrop the image using the 
            clc;
            close all;
            imshow(rgb1);
            
                 
                calcWcap1 = length(find(BW1 == 1));
                
                %8888888888888888888888888888888888888888888888888888888888
%                 PixIdx=regionprops(BW1,'PixelIdxList');
                
               
                %8888888888888888888888888888888888888888888888888888888888
                
             while waitforbuttonpress~=0; end
                mouseside=get(gcf,'SelectionType');
                abm=sum(mouseside);
                STR1=[110   111   114   109    97   108];   %normal in ASCII
                au=sum(STR1);
                STR2=[97   108   116];                      %alt in ASCII
                bu=sum(STR2);
                
                if abm==au          %left mouse click (Stage A)
                %PixIdx_A=[PixIdx_A;PixIdx];
                    for m=x2:x2+x4           %scan all pixels
                    for n=x1:x1+x3       %and if one is white, colour it
                        if rgb(m,n)==255;    %blue
                           rgb(m,n,1)=0;  
                           rgb(m,n,2)=0;
                           rgb(m,n,3)=255;
                        
                        end
                    end
                end           
              
                    a=a+calcWcap1;
                    wpixcount=wpixcount-calcWcap1;
                    
                elseif abm==bu      %right mouse click (Stage B)
                %PixIdx_B=[PixIdx_B;PixIdx];    
                for m=x2:x2+x4           %scan all pixels
                    for n=x1:x1+x3       %and if one is white, colour it
                        if rgb(m,n)==255;    %green
                           rgb(m,n,1)=0;  
                           rgb(m,n,2)=255;
                           rgb(m,n,3)=0;
                        end
                    end
                end  
                    b=b+calcWcap1;
                    wpixcount=wpixcount-calcWcap1;
                    pause(.2);
                    %%%% Case for unwanted white pixels caused by seagulls %%%%
                elseif abm == 648   % If the middle button is hit, white is blacked out
                    for m=x2:x2+x4           %scan all pixels
                    for n=x1:x1+x3       %and if one is white, colour it
                        if rgb(m,n)==255;    %green
                           rgb(m,n,1:3)=0;
                        end
                    end
                    end
                    %%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
                    calcWcap=calcWcap-(calcWcap1*100/numPixels);
                    wpixcount=wpixcount-calcWcap1;
                    
                end
             
             disp(wpixcount);
             
             BW((x2):(x2+x4),(x1):(x1+x3))=0;
             
        
             end
        
end %while loop, only finished if all pixels are done
        
if finished~=2
Wa=a*100/numPixels;
Wb=b*100/numPixels;
W=Wa+Wb;
end






if W>0
 %correct method for obtaining Pixel Index for whitecap regions::
x2=rgb(:,:,2);
x3=rgb(:,:,3);
c3=bwconncomp(x3);
c2=bwconncomp(x2);
PixIdx_A=regionprops(c3,'PixelIdxList');
PixIdx_B=regionprops(c2,'PixelIdxList');
%get the Wa and Wb pixel indices::
da=[]; db=[];
if ~isempty(PixIdx_A)
    for l=1:length(PixIdx_A)
        da=[da;PixIdx_A(l).PixelIdxList];
    end
    da=da(:);
end

if ~isempty(PixIdx_B)
    for ll=1:length(PixIdx_B)
        db=[db;PixIdx_B(ll).PixelIdxList];
    end
    db=db(:);
end

%Find the fill value, if appropriate::
x2(BW_orig)=255;
x2(da)=0;
x2(db)=0;
c1=bwconncomp(x2);
fill=regionprops(c1,'PixelIdxList');
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

imshow(rgb);
pause(0.3);
end
close all;


end  

