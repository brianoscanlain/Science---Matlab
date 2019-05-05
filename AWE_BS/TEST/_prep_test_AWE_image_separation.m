%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Training routine for the AWE algorithm. The following images were found to
%have issues::
%
%1. the absence of whitecap images, where AWE finds whitecaps. opt=10 fails
%
%2. the AWE finds a very close threshold opt=4,3,2 work
%
%3. Threshold is very high, W underestimated. opt=4,3,2 fail
%
%4. Threshold is low, W is overestimated opt=4,3,2 fail

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Objectives, run the AWE algorithm again, this time record as much
%information about each image, (derivative, image intensity patterns, image
%size, Regions of interest counts, mean and std background pixel intensity,
% mean and std whitecap pixel intensity, blue, green, gray slopes)
%
%Evaluate this data and identify patterns of the above mentioned values for
%the 4 cases of images. 
%
%Then, try to apply corrections to the AWE code to help reduce these errors
%from occurring again.

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%The following data represents arrays of indices which identify the image
%files to load from the dataset labelled ASIP_deployment_2 from the Knorr
%cruise 2011.

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Threshold looks good
good=[101;151;301;601;701;1001;1501;1651;1751;1801;6;106;256;606;906;...
1556;1606;1656;1756;];

%should be empty image, but not
null_fail=[];

%empty image, but has whitecaps
null_fail2=[];

%T set high
T_high=[1;51;201;251;401;501;801;901;951;1101;1151;1201;1401;1551;56;...
156;306;356;406;806;956;1406;1456;1706;1806];

%T set low
T_low=[451;651;456;1006];

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%set output directories::
lowDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\low_t';
goodDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\good_t';
highDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\high_t';
null_fail2Dir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\presence_null';
null_failDir='C:\Users\Brian\Desktop\Whitecap_Project\AWE_test_images\null_presence';

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%set input diectory and scan all jpg files::
imageDir='D:\Work\Whitecap\ASIP_dep4';
cd(imageDir);
files=dir('*.jpg');
totalIm=length(null_fail)+length(null_fail2)+length(T_high)+length(T_low);

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%put good images into directory::
for i=1:length(good)
    copyfile(files(good(i)).name,goodDir)
    clc;
    disp([num2str((i*100)/totalIm) '% complete']);
end

%put null_presence into Directory::
for j=1:length(null_fail)
     copyfile(files(null_fail(j)).name,null_failDir)
    clc;
    disp([num2str(((i+j)*100)/totalIm) '% complete']);
end

%put presence_null into directory::
for k=1:length(null_fail2)
     copyfile(files(null_fail2(k)).name,null_fail2Dir)
    clc;
    disp([num2str(((i+k+j)*100)/totalIm) '% complete']);
end

%put images with a resulting high T in directory::
for m=1:length(T_high)
   copyfile(files(T_high(m)).name,highDir); 
       clc;
    disp([num2str(((i+k+j+m)*100)/totalIm) '% complete']);
end

%put images with a resulting low T in directory::
for n=1:length(T_low)
    copyfile(files(T_low(n)).name,lowDir);
        clc;
    disp([num2str(((i+k+j+m+n)*100)/totalIm) '% complete']);
end

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
