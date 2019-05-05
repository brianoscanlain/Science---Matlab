function [ output_args ] = hist2d_Filt(x,y,N)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

debug=1;
if debug==1
    x=bow.symeo.distance;
    y=bow.symeo.velocity;
    N=100;
end

%Settings:
ZScore=4; %Span the matrix out +/-Zscore*STD(x|y)

%=========================
%Create population matrix:
%=========================
numR = length(x);
stdX = std(x); meanX = mean(x);
stdY = std(y); meanY = mean(y);
%Define Boundaries of Matrix:
leftEdge = max(meanX-stdX*ZScore,min(x));
riteEdge = min(meanX+stdX*ZScore,max(x));
botEdge = max(meanY-stdY*stdTimes,min(y));
topEdge = min(meanY+stdY*stdTimes,max(y));
%Create matrix:
Hist2d=zeros(N,N);
%Define Axes of Marix:
HaxisX=linspace(leftEdge,riteEdge,N);
HaxisY=linspace(botEdge,riteEdge,N);
%=========================

%Main Loop
%=========================
for i=1:2 %
    
   for j=1:numR
       
    
    
   end
end



%=========================
end

