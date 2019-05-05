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
%Perform the smoothing and derivative analysis

function [smoothed,smoothed_2,smFirstDeriv,smSecondDeriv,smThirdDeriv,smFourthDeriv,threshold,PIP] = derivAnalysis_BS(threshold,PIP,smCoef)

% fN = fftAWE(PIP);
% fN2 = fftAWE(detrend(PIP));

%Begin butterworth filtering and smoothing
[b,a] = butter(4,smCoef);      %0.5 times the nyquist frequency
PIP = [flipud(PIP(1:20));PIP];
smoothed = filtfilt(b,a,PIP);
smoothed(1:20) = [];
PIP(1:20) = [];
%Now cut it here
startPos = bsearch(threshold,0.99);
smoothed(1:startPos) = [];
threshold(1:startPos) = [];
PIP(1:startPos) = [];

%Implement running average smoothing
% smoothed = mySmooth_2(smoothed);

if smoothed(end-1) < smoothed(end) && smoothed(end-1) < smoothed(end-2)
    smoothed(end-1) = mean([smoothed(end);smoothed(end-2)]);
end
if smoothed(end-1) > smoothed(end) && smoothed(end-1) > smoothed(end-2)
    smoothed(end-1) = mean([smoothed(end);smoothed(end-2)]);
end
%Implement further running average smoothing
%[b,a] = butter(6,0.1);
%smoothed_2 = filtfilt(b,a,smoothed);
% smoothed_2 = mySmooth_2(smoothed);
% 
% if smoothed_2(end-1) < smoothed_2(end) && smoothed_2(end-1) < smoothed_2(end-2)
%     smoothed_2(end-1) = mean([smoothed_2(end);smoothed_2(end-2)]);
% end
% if smoothed_2(end-1) > smoothed_2(end) && smoothed_2(end-1) > smoothed_2(end-2)
%     smoothed_2(end-1) = mean([smoothed_2(end);smoothed_2(end-2)]);
% end
smoothed_2 = smoothed;
firstDeriv = calcDeriv(smoothed_2);

%Now smooth the first derivative and get the second derivative
%of that
smFirstDeriv = firstDeriv;
smFirstDeriv = mySmooth_Deriv(firstDeriv);


% hold on;plot(smFirstDeriv,threshold,'-bo');
%Here I have to try to eliminate any single negative first derivatives as
%they will cause problems later on

%Remove any pos2neg changes of 2 elements or less
for to = 1:length(smFirstDeriv)-3
    if smFirstDeriv(to) <0 && smFirstDeriv(to+1) >0 && smFirstDeriv(to+2) >0 && smFirstDeriv(to+3) <0
        meanVal = mean([smFirstDeriv(to) smFirstDeriv(to+3)]);
        smFirstDeriv(to+1) = meanVal;
        smFirstDeriv(to+2) = meanVal;
    elseif smFirstDeriv(to) >0 && smFirstDeriv(to+1) <0 && smFirstDeriv(to+2) <0 && smFirstDeriv(to+3) >0
        meanVal = mean([smFirstDeriv(to) smFirstDeriv(to+3)]);
        smFirstDeriv(to+1) = meanVal;
        smFirstDeriv(to+2) = meanVal;
    end
    
end

for ty = 2:length(smFirstDeriv)-1
    if smFirstDeriv(ty) < 0 && smFirstDeriv(ty-1) > 0 && smFirstDeriv(ty+1) > 0
        smFirstDeriv(ty) = mean([smFirstDeriv(ty-1);smFirstDeriv(ty+1)]);
    end
    if smFirstDeriv(ty) > 0 && smFirstDeriv(ty-1) < 0 && smFirstDeriv(ty+1) < 0
        smFirstDeriv(ty) = mean([smFirstDeriv(ty-1);smFirstDeriv(ty+1)]);
    end
end
secondDeriv = calcDeriv(smFirstDeriv);
%Now smooth the second Derivative in order to get the third derivative
smSecondDeriv = secondDeriv;
% smSecondDeriv = mySmooth_Deriv(secondDeriv);
%Remove any pos2neg changes of 2 elements or less
for to = 1:length(smSecondDeriv)-3
    if smSecondDeriv(to) >0 && smSecondDeriv(to+1) <0 && smSecondDeriv(to+2) <0 && smSecondDeriv(to+3) >0
        meanVal = mean([smSecondDeriv(to) smSecondDeriv(to+3)]);
        smSecondDeriv(to+1) = meanVal;
        smSecondDeriv(to+2) = meanVal;
    end
end
%Here I have to try to eliminate any single negative elements
for ty = 2:length(smSecondDeriv)-1
    if smSecondDeriv(ty) < 0 && smSecondDeriv(ty-1) > 0 && smSecondDeriv(ty+1) > 0
        smSecondDeriv(ty) = mean([smSecondDeriv(ty-1);smSecondDeriv(ty+1)]);
    end
end

thirdDeriv = calcDeriv(smSecondDeriv);

smThirdDeriv = mySmooth_Deriv(thirdDeriv);

fourthDeriv = calcDeriv(smThirdDeriv);

smFourthDeriv = mySmooth_Deriv(fourthDeriv);

%Now normalise the derivatives
%%%This is the section that might need changing: normalise in a different
%%%way?????
smFirstDeriv = smFirstDeriv(3:end-3)./max(smFirstDeriv(3:end-3));
smSecondDeriv = smSecondDeriv(3:end-3)./max(smSecondDeriv(3:end-3));
smThirdDeriv = smThirdDeriv(3:end-3)./max(smThirdDeriv(3:end-3));
threshold = threshold(3:end-3);
PIP = PIP(3:end-3);
smoothed = smoothed(3:end-3);
smoothed_2 = smoothed_2(3:end-3);

% %Chop the beginning and end of PIP,threshold,smoothed,smoothed_2,
% %smFirstDeriv,smSecondDeriv,smThirdDeriv,smFourthDeriv
% smoothed = smoothed(10:end);%./max(abs(smoothed));
% smoothed_2 = smoothed_2(10:end);%./max(abs(smoothed_2));
% smFirstDeriv = smFirstDeriv(10:end)./max(abs(smFirstDeriv(10:end)));
% smSecondDeriv = smSecondDeriv(10:end)./max(abs(smSecondDeriv(10:end)));
% smThirdDeriv = smThirdDeriv(10:end)./max(abs(smThirdDeriv(10:end)));
% smFourthDeriv = smFourthDeriv(10:end)./max(abs(smFourthDeriv(10:end)));
% threshold = threshold(10:end);
% PIP = PIP(10:end);
    