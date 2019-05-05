%Script to run the AWE algorithm, for the selected test images.
%
% 5 image image classes are defined;
%
% Good= Current AWE Threshold is good
% NP=   Null image resulting in white regions after current AWE processing
% PN=   Whitecap image being registered as Null after current AWE proc.
% High= Current AWE threshold is set high.
% Low=  Current AWE threshold is set low.


%This script loads the images, finds the AWE threshold, and allows the user
%to manually tune the most suitable one. Then, as much information about
%the processing is saved to a master file. The objective is to try and
%locate the culprit within the algorithm which has trouble distinguishing
%between the types of errors listed above (obviously 'Good' is not an
%error, and is here to give insight into when the AWE works very well).
%
%
%
% Straight away, I believe the NP error to ocurr due to the very
% sensitively set peak threshold within the *deriv.mat file. maybe
% adjusting this to another constant, or perhaps a dependant variable may
% reduce or eliminate this error. Each adjustment can be tested using the
% 'Good' and 'NP' values.
%
% When seagulls exist, and appear brighter than the observed whitecaps,
% they seem to play havoc with the AWE processing. I think this effect
% cannot be eliminated. If there exists a big difference between whitecap
% PI and seagull PI, then the AWE T will not be set accurately.
%
% Another problem which seems to induce the NP error, is early in the
% morning. This is due to strong PI gradients across the length of the
% image, making front facing wave facelets in the forefront of images 
% appear as whitecaps. Playing around with the PI normalizer could
% potentially redce this effect
%
%
% Another idea is to look at the 'opt' values of each image and see when
% opt 2,3 & 4 spread light on the errors stated above.
%
%
% Sofar, it seems fit to record the following:
% AWE-T, tuned-T, slope_of_blue/green/gray_profile, 
%
%