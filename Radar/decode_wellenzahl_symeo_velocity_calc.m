%Script for decoding the SPI data:
clear all; clc; close all;
%The HDF file has eight datasets enclosed:
%  1. FMCW nested structure
%            - Logical array [1 x 2048 array] ...==1 if FMCW option is used
%  2. Data_imag [2048 x 2048]
%  3. Data_real [2048 x 2048]
%  4. FrequencyStartGHz [1 x 2048]
%  5. FrequencyStopGHz [1 x 2048]
%  6. localtime [1 x 2048]
%  7. sweeptime [1 x 2048]
%  8. time [1 x 2048]
%
% Brian Scanlon, Galway, 27/02/2018
%-----------------------------------------
c=299792458; %speed of light
alpha=0.16; %Blackman filter coefficients
a0=(1-alpha)/2; %Blackman filter coefficients
a1=1/2;%Blackman filter coefficients
a2=alpha/2;%Blackman filter coefficients
ZoomFactor=10;
minRangePeak=1.5;
PLL_mode=1; %This must be set in the radar's config.ini!

%Sweeptime matrix [ sweeptimeidx,PLLrampmode ]:
swptm=[3.84e-1; 6.4e-1; 1.536; 2.688; 5.376; 13.312]*1e-3;
swptm(:,2)=swptm(:,1)*2;
%ADC Baseband sample rate in Hz [Sweeptimeidx, PLLrampmode ]:  %This is the
%frequency sampling rate, not the frequency step of the chirp!
BB_ADC_sample_Rate=[5.333e6;3.2e6;1.333e6;7.6190e5;3.8095e5;1.5385e5];
BB_ADC_sample_Rate(:,2)=BB_ADC_sample_Rate(:,1)./2;


%Load the data:
nfo5=hdf5info('sti3_pll1_Corridor_const_spd_gain2.hdf');   %Triangle modulation & velocity test
FMCW=hdf5read(nfo5.GroupHierarchy.Datasets(1));
I=hdf5read(nfo5.GroupHierarchy.Datasets(2));
Q=hdf5read(nfo5.GroupHierarchy.Datasets(3));
FqStrt=hdf5read(nfo5.GroupHierarchy.Datasets(4)).*1e9;
FqStp=hdf5read(nfo5.GroupHierarchy.Datasets(5)).*1e9;
localtime=hdf5read(nfo5.GroupHierarchy.Datasets(6));
Tsweep_idx=hdf5read(nfo5.GroupHierarchy.Datasets(7));
time=hdf5read(nfo5.GroupHierarchy.Datasets(8));
%Timestamp processed:
dte=datenum(datenum(1970,1,1,0,0,0)+(time(:)/86400));


%Calculate various properties:
Bsweep=median(FqStp)-median(FqStrt);
Tsweep=swptm(mode(Tsweep_idx)+1,PLL_mode+1); %plus one addition as Matlab starts at one, not zero!
BBadc_FS=BB_ADC_sample_Rate(mode(Tsweep_idx)+1,PLL_mode+1); %plus one addition as Matlab starts at one, not zero!
Chirp_Hz=1/median(diff(time));  %Sample rate of chirps
CentreFq=median((FqStp + FqStrt)/2);
FqStep=Bsweep/(length(time));
TsweepStep=Tsweep/size(I,2);
ChirpSlope=Tsweep/Bsweep;



%======================================
%Reconstruction and Blackman filtering:
%======================================
%Reconstruct the I & Q data to form the full signal:
ADC = sqrt(I.^2 + Q.^2) .* cos( atan2(I,Q) ); %row = distance, column = velocity
%ADC=ADC(:,1650:1670); % Take a handful of chirps
%Detrend the data:
ADC = detrend(ADC,0);
if PLL_mode==1  %Split up the two slopes and apply separate Blackman filters!
    ADC_uc=ADC(1:floor(size(ADC,1)/2),:); %up chirp
    ADC_dc=ADC(floor(size(ADC,1)/2)+1:end,:); % down chirp
    nuc=size(ADC_uc,1);
    ndc=size(ADC_dc,1);
    BlkWinU = (a0 - a1*cos(2*pi*(0:nuc-1)/(nuc-1))  + a2*cos(4*pi*(0:nuc-1)/(nuc-1)))';
    BlkWinD = (a0 - a1*cos(2*pi*(0:ndc-1)/(ndc-1))  + a2*cos(4*pi*(0:ndc-1)/(ndc-1)))';
    ADC_uc=ADC_uc.*(repmat(BlkWinU,1,size(ADC_uc,2)));
    ADC_dc=ADC_dc.*(repmat(BlkWinD,1,size(ADC_dc,2)));
else
    %Calculate the blackman window filtered ADC data:
    n=size(ADC,1);
    BlkWin = (a0 - a1*cos(2*pi*(0:n-1)/(n-1))  + a2*cos(4*pi*(0:n-1)/(n-1)))';
    %ADC=ADC.*(WinD*WinV);
    ADC = ADC.*(repmat(BlkWin,1,size(ADC,2)));
end
%======================================



%======================================
%Compute FFT
%======================================
%calculate the 2-D FFT:   %NB:: Do not use FFT2, it doesn't calculate the
%correct result!, why exactly, I don't know, but I guess that it isn't giving chirp-independent fft's, but instead must be taking into account information from neighbouring chirps!
if PLL_mode==1
    %Up chirp padding:
    padsU=zeros(size(ADC_uc,1)*ZoomFactor,size(ADC_uc,2));%FFT padding:
    padsU(1:size(ADC_uc,1),1:size(ADC_uc,2))=ADC_uc;%FFT padding:
    ADC_uc=padsU;%FFT padding:
    %Down chirp padding:
    padsD=zeros(size(ADC_dc,1)*ZoomFactor,size(ADC_dc,2));%FFT padding:
    padsD(1:size(ADC_dc,1),1:size(ADC_dc,2))=ADC_dc;%FFT padding:
    ADC_dc=padsD;%FFT padding:
    for  i=1:size(ADC_uc,2)
        levelsU(:,i)=fft(ADC_uc(:,i));
        levelsD(:,i)=fft(ADC_dc(:,i));
    end
    levelsU=(levelsU(1:size(levelsU,1)/2,:));
    phasorU=angle(levelsU);
    levelsU_cmplx=(levelsU);
    levelsU=abs(levelsU);
    levelsU(2:end-1,:)=levelsU(2:end-1,:).*2;
    levelsD=(levelsD(1:size(levelsD,1)/2,:));
    phasorD=angle(levelsD);
    levelsD_cmplx=(levelsD);
    levelsD=abs(levelsD);
    levelsD(2:end-1,:)=levelsD(2:end-1,:).*2; 
else
    %padd the data:
    pads=zeros(size(ADC,1)*ZoomFactor,size(ADC,2));%FFT padding:
    pads(1:size(ADC,1),1:size(ADC,2))=ADC;%FFT padding:
    ADC=pads;%FFT padding:
    for  i=1:size(ADC,2)
        levels(:,i)=fft(ADC(:,i));
    end
    levels=abs(levels(1:size(levels,1)/2,:));
    levels(2:end-1,:)=levels(2:end-1,:).*2;
end
%======================================







%======================================
% Calculate the Distance (and velocity)
%======================================
if PLL_mode==1
    %Calculate the axes:
    fDistU=(0:size(levelsU,1)-1) ./ (size(levelsU,1)*ZoomFactor) *BBadc_FS;
    fDistD=(0:size(levelsD,1)-1) ./ (size(levelsD,1)*ZoomFactor) *BBadc_FS;
    RangeUC = fDistU * (ChirpSlope*c/2);
    RangeDC = fDistD * (ChirpSlope*c/2);
    %Filter instrument noise range (below the set minRange)
    minRangeIndxU=min(find(RangeUC>minRangePeak));
    levelsU(1:minRangeIndxU,:)=0;
    minRangeIndxD=min(find(RangeDC>minRangePeak));
    levelsD(1:minRangeIndxD,:)=0;
    %Allocate memory:
    RangeCalcU=zeros(size(levelsU,2),3)./0;
    RangeCalcD=zeros(size(levelsD,2),3)./0;
    fBeatU=zeros(size(levelsU,2),3)./0;
    fBeatD=zeros(size(levelsD,2),3)./0;
    PeakIndU=zeros(size(levelsU,2),3)./0;
    PeakIndD=zeros(size(levelsD,2),3)./0;
    PeakLvlU=zeros(size(levelsU,2),3)./0;
    PeakLvlD=zeros(size(levelsD,2),3)./0;
    phaseVU=zeros(size(levelsU,2),3)./0;
    phaseVD=zeros(size(levelsD,2),3)./0;
    %Perform computations:
    for x=1:size(levelsU,2) %cycle through each chirp
        [pkLvlsU,pksU]=findpeaks(levelsU(:,x),'MinPeakDistance',10); %Have a minimum spacing between peaks
        [pkLvlsD,pksD]=findpeaks(levelsD(:,x),'MinPeakDistance',10); %Have a minimum spacing between peaks
        [~,srtTransU]=sort(pkLvlsU,'descend'); %find the translation index map to rank peaks by level (descending)
        [~,srtTransD]=sort(pkLvlsD,'descend'); %find the translation index map to rank peaks by level (descending)
        pksU=pksU(srtTransU); %sort peak indices (descending)
        pksD=pksD(srtTransD); %sort peak indices (descending)
        pkLvlsU=pkLvlsU(srtTransU);
        pkLvlsD=pkLvlsD(srtTransD);
        PeakIndU(x,1:3)=pksU(1:3);
        PeakIndD(x,1:3)=pksD(1:3);
        fBeatU(x,1:3)=fDistU(pksU(1:3));
        fBeatD(x,1:3)=fDistD(pksD(1:3));
        Velocity(x,1)=(fBeatD(x,1)-fBeatU(x,1)) /4/CentreFq*c;
        Range(x,1)=(fBeatD(x,1)+fBeatU(x,1))/4*ChirpSlope*c;
        %PeakLvlU(x,1:3)=pkLvlsU(1:3);
        %PeakLvlD(x,1:3)=pkLvlsD(1:3);
        %phaseVU(x,1:3)=phasorU(pksU(1:3));
        %phaseVD(x,1:3)=phasorD(pksD(1:3));
        %RangeCalcD(x,1:3)=RangeDC(pksD(1:3));
        %RangeCalcU(x,1:3)=RangeUC(pksU(1:3));
        
    end
else
%Calculate the axes:
fDist=(0:size(levels,1)-1) ./ (size(levels,1)*ZoomFactor) *BBadc_FS;
%fVel=(0:size(levels,2)-1) ./ (size(levels,2)*ZoomFactor) *Chirp_Hz;
%fVel = fVel - fVel(size(ADC,2)*ZoomFactor/2); %shift the velocity to +/- zero
%Convert these frequency scales to Range and velocity scales:
Range = fDist * (ChirpSlope*c/2);
%Velocity = fVel * (c / (2 * CentreFq - (Bsweep/2)));
%Filter instrument noise below minRange:
minRangeIndx=min(find(Range>minRangePeak));
levels(1:minRangeIndx,:)=0;
%Find the fBeat for each chirp:
RangeD=zeros(size(levels,2),3)./0;
fBeat=zeros(size(levels,2),3)./0;
PeakInd=zeros(size(levels,2),3)./0;
PeakLvl=zeros(size(levels,2),3)./0;
phaseV=zeros(size(levels,2),3)./0;
for x=1:size(levels,2) %cycle through each chirp
    [pkLvls,pks]=findpeaks(levels(:,x),'MinPeakDistance',10); %Have a minimum spacing between peaks
    [~,srtTrans]=sort(pkLvls,'descend'); %find the translation index map to rank peaks by level (descending)
    pks=pks(srtTrans); %sort peak indices (descending)
    pkLvls=pkLvls(srtTrans);
    PeakInd(x,1:3)=pks(1:3);
    fBeat(x,1:3)=fDist(pks(1:3));
    PeakLvl(x,1:3)=pkLvls(1:3);
%    phaseV(x,1:3)=phasor(pks(1:3));
    RangeD(x,1:3)=Range(pks(1:3));
end
%==
end
%======================================


% fft2=abs(fft(fBeat(:,1)));  
% fft2=fft2[1:int(fft2.shape[0]/2)]
% fft2_=np.fft.fftshift(fft2)
%plot(Velocity,fBeat(:,1))
plot(Range(:,1))
figure
plot(Velocity(:,1))
median(diff(dte))*86400