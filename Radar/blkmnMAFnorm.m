function [Xnorm,MvAvg] = blkmnMAFnorm(X,n,PaddingType)
%blkmnMAFnorm is a Moving Average Filter (MAF), using a Blackman Window.
%   Inputs:
%             X  =  signal (1 dimensional)   
%             n  =  window size
%   PaddingType  =  'none'     (nopadding at all, we will loose data on ends of X!)
%                   'zeros'    (pads the start and end with zeros)
%                   'real'     (pads the start and end with reflected X data
%                   
%   Outputs:
%         Xnorm  =  Normalised X data (X - MvAvg)
%         MvAvg  =  Moving average found for the data.
%
%
%  Brian Scanlon, NUIG, March 2018
%--------------------------------------------------------------------------

%-------------------------------------
%Debugging preamble (ignore when calling real function():
%-------------------------------------
% nfo5=hdf5info('F:\Dropbox (AirSea Laboratory)\fuelwave\radars\Wellenzahl\23Feb_Salthill_10degrees_8metres\20180223_134713_rec.hdf');   %20 degree
% data2=hdf5read(nfo5.GroupHierarchy.Datasets(2));
% data3=hdf5read(nfo5.GroupHierarchy.Datasets(3));
% data4=hdf5read(nfo5.GroupHierarchy.Datasets(4));
% data5=hdf5read(nfo5.GroupHierarchy.Datasets(5));
% data6=hdf5read(nfo5.GroupHierarchy.Datasets(6));
% data7=hdf5read(nfo5.GroupHierarchy.Datasets(7));
% x1=data3(:,770);
% n=16;
% PaddingType='real';
% %Calculate the FFT:
% yI=fft(x1);
% yI=abs(yI/length(x1));
% yI=yI(1:length(x1)/2);
% yI(2:end-1)=2*yI(2:end-1);
% X=yI;
% c=3e8;
% Bsweep=abs(data5(70)-data4(70))*1e9; %Bandwidth of the sweep (in Frequency)
% Tsweep=data7(70)/1000;   %convert Tsweep from ms to seconds
% f=(c/2)*Tsweep/Bsweep.*(0:length(yI)-1)*1000;
%---------------------


%-------------------------------------
%Prepare the X data:
%-------------------------------------
%Measure the length of the X array:
nLen=length(X);
%Ensure the X array has a columnwise expansion:
if size(X,1)>size(X,2)
    X=X';
    swap=1; %take note if we swapped, so we can undo the swapping for Argouts
else
    swap=0;
end

%pad the beginning and end of the X array with zeros:
if strcmp(PaddingType,'zero')
Xpadded=[zeros(1,n/2), X, zeros(1,n/2)];
elseif strcmp(PaddingType,'real')
Xpadded=[fliplr(X(1:n/2)), X, fliplr(X(end-n/2+1:end))]; 
elseif strcmp(PaddingType,'none') || isempty(PaddingType)
Xpadded=X;    
end
nLenPad=length(Xpadded);
MovingAvg=zeros(1,nLen);
%-------------------------------------


%-------------------------------------
%Create the window function:
%-------------------------------------
alpha=0.16;
a0=(1-alpha)/2;
a1=1/2;
a2=alpha/2;
BlkWin = a0 - a1*cos(2*pi*(0:n-1)/(n-1))  + a2*cos(4*pi*(0:n-1)/(n-1));
BaseWin=zeros(1,nLenPad);
%-------------------------------------



%--------------------------------------------------------
%Loop for each point to calculate the new Moving average
%--------------------------------------------------------
for i=1:nLenPad-n
    i;
    iWin=BaseWin;
    iWin(i:i+n-1)=BlkWin;
    MovingAvg(i)=mean(iWin.*Xpadded); %overwrite the MvAvg for each iteration
end
%-------------------------------------    

% plot(MovingAvg)
% figure
% plot(MovingAvg*mean(X)/mean(MovingAvg))
% hold on
% plot(X);
% figure
%Vargout:
Xnorm=X-(MovingAvg*mean(X)/mean(MovingAvg));
MvAvg=(MovingAvg*mean(X)/mean(MovingAvg));
    
end






