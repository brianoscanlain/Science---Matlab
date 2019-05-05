%Script for decoding the SPI data:

%The HDF file has eight datasets enclosed:
%  1. FMCW nested structure
%            - Logical array [1 x 2048 array] ...==1 if FMCW option is used
%  2. Data_imag [2048 x 2048]
%  3. Data_real [2048 x 2048]
%  4. FrequencyStartGHz [1 x 2048]
%  5. FrequencyStopGHz [1 x 2048]
%  6. localtime [1 x 2048]
%  7. sweeptime [1 x 2048]   %sweeptime index, not the actual sweeptime!!
%  8. time [1 x 2048]
%
% Brian Scanlon, Galway, 27/02/2018
%-----------------------------------------
c=299792458; %speed of light
alpha=0.16; %Blackman filter coefficients
a0=(1-alpha)/2; %Blackman filter coefficients
a1=1/2;%Blackman filter coefficients
a2=alpha/2;%Blackman filter coefficients


nfo5=hdf5info('sti0_pll1_gain4.hdf');   %Triangle modulation & velocity test

%Identify the file and create a hdf handle:
% nfo5=hdf5info('28Mar18_1.5_23m_2.hdf');
% nfo5=hdf5info('28Mar18_0_23m.hdf');
%nfo5=hdf5info('20180607_155706_rec.hdf');   %Triangle modulation & velocity test
%nfo5=hdf5info('triangle_roof_23_03_18.hdf');   %Triangle modulation & velocity test
% nfo5=hdf5info('10Degrees.hdf');   %10 degree
% nfo5=hdf5info('20Degrees.hdf');   %20 degree
% nfo5=hdf5info('25Degrees.hdf');   %25 degree
% nfo5=hdf5info('30Degrees.hdf');   %30 degree
% nfo5=hdf5info('35Degrees.hdf');   %35 degree
%Extract the data (I already know there are 8 datasets,
%if uncertain, use >>{nfo5.GroupHierarchy.Datasets.Name})
data1=hdf5read(nfo5.GroupHierarchy.Datasets(1));
data2=hdf5read(nfo5.GroupHierarchy.Datasets(2));
data3=hdf5read(nfo5.GroupHierarchy.Datasets(3));
data4=hdf5read(nfo5.GroupHierarchy.Datasets(4));
data5=hdf5read(nfo5.GroupHierarchy.Datasets(5));
data6=hdf5read(nfo5.GroupHierarchy.Datasets(6));
data7=hdf5read(nfo5.GroupHierarchy.Datasets(7));
data8=hdf5read(nfo5.GroupHierarchy.Datasets(8));

%Timestamp processed:
dte=datenum(datenum(1970,1,1,0,0,0)+(data8(:)/86400));



for i=1:length(data8)
    %Reconstruct original waveform
     A=sqrt( data2(:,i).^2 + data3(:,i).^2); 
     P=atan2(data2(:,i) , data3(:,i) );
     x=A.*cos(P); %reconstruct real part only (omitt "+iSin(P)")
    % x=A.*(cos(P)+ sqrt(-1)*sin(P)); %reconstruct real and imag part 

    %Normalise the x data using the Blkmnn filter:    
    n=size(x,1);
    BlkWin = (a0 - a1*cos(2*pi*(0:n-1)/(n-1))  + a2*cos(4*pi*(0:n-1)/(n-1)))';
    x=BlkWin.*detrend(x,0);
   
    %find the FFT of x, and plot it:
    y=fft(x);
    w=angle(y);
    y=abs(y/length(x)); %not needed if you use rfft()
    y=y(1:length(x)/2+1); % " "
    w=w(1:length(x)/2+1); % " "
    y(2:end-1)=2*y(2:end-1); % " "
    w(2:end-1)=2*w(2:end-1); % " "
    %convert the normalised fft x-axis to Distance: (in loop incase settings change during operation)

    Bsweep=abs(data5(i)-data4(i))*1e9; %Bandwidth of the sweep (in Frequency)
    Tsweep=double(data7(i))./1000;   %convert Tsweep from ms to seconds
    normFq=(0:length(y)-1); %corresponding fft frequency [normalised] (0: Nyq_Fq)
    TsweepStep=Tsweep/(length(y)-1);
    actualFq=normFq  / ((length(y)-1)*TsweepStep);  % (0:Nyq_fq*DeltaFreq), DeltaFreq= fftlength/Tsweep[seconds];
    chirpSlope=Tsweep/Bsweep;  %Used to convert values in frequency space to time space.
    %Range = c*tau/2:
    Range=(c/2)*chirpSlope.*actualFq; %Range [m] (here, chirpSlope*actualFq = time delay tau)
    
    %Filter peaks f<1.5 (instrument errors)
    y(Range<1.5)=0;
    
    %Find the top three peaks:
    [DistancePeaks(i,1:3),LevelPeaks(i,1:3),wPeaks(i,1:3)]=WellenzahlLocMax(y,Range,w);
end



%select the maximum-level peaks:
LevelPeaks(DistancePeaks<1.5)=0;
[~,LvlMax]=max(LevelPeaks,[],2);
LvlMax=[LvlMax==1 LvlMax==2 LvlMax==3]; %convert it to a logical




%Figure Distance timeseries:
%===========================
figure
plot((dte-dte(1))*86400,DistancePeaks(LvlMax),'.')
ylabel('distance [m]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[0 6])
%Figure Distance vs Level:
%========================
figure;
plot(DistancePeaks,LevelPeaks,'.')
xlabel('distance [m]','Fontname','Times','Fontsize',14);
ylabel('Amplitude of return signal','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[0 1e3])
hleg3=legend('1st max','2nd max','3rd max');
%========================




%Filter out spikes:
DistancePeaks(abs(diff(DistancePeaks))>0.2)=0./0;
www=1:length(DistancePeaks);
DistancePeaks=interp1(www(~isnan(DistancePeaks(:,1))),DistancePeaks(~isnan(DistancePeaks(:,1)),1),www,'linear','extrap');
fprintf('\n\n Hsig =  %03.2f [m]\n\n',4*sqrt(std(DistancePeaks)))


%Convert the distance back into beat frequency
fBeat=DistancePeaks(LvlMax)*Bsweep/Tsweep/c*2/1000;

%Clean up the fBeat signal:
fBeat(diff(fBeat)> 3*mad(diff(fBeat)) | diff(fBeat)<-3*mad(diff(fBeat)))=0./0;
fBeat=interp1(find(~isnan(fBeat)),fBeat(~isnan(fBeat)),1:length(fBeat),'linear','extrap')';
fBeat=DaFixer(1:2048,fBeat,'greater',median(fBeat)+2*mad(fBeat));

K=c/122e9; %Wavelength of the Radar signal:
Velocity= [0;diff(fBeat*1000)*K/2];


%Figure Velocity timeseries:
%===========================
figure
plot((dte-dte(1))*86400,Velocity,'.')
ylabel('Velocity [m/s]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[-15 15])


figure
plot((dte-dte(1))*86400,Velocity,'.')
hold on;
plot((dte-dte(1))*86400,runmean(Velocity,10),'.')
ylabel('Velocity [m/s]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[-15 15])
hleg2=legend('raw Velocity','runmean(velocity,10)')



%Next, let's integrate and find the displacement from this velocity:
Vel=runmean(Velocity,10);
DeltaT=median((diff(dte)))*86400; %time step between samples [seconds]
figure;
plot((dte-dte(1))*86400,cumsum(Vel).*DeltaT,'-','linewidth',1.4)
hold on;
plot((dte-dte(1))*86400,detrend(DistancePeaks(LvlMax),0),'.')
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[-0.5 2.5])
hleg3=legend('Velocity-integrated','measured')
ylabel('Displacement [m]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);












fBeat(diff(fBeat)> 3*mad(diff(fBeat)) | diff(fBeat)<-3*mad(diff(fBeat)))=0./0;
fBeat=interp1(find(~isnan(fBeat)),fBeat(~isnan(fBeat)),1:length(fBeat),'linear','extrap');
%
%fBeat(fBeat>median(fBeat)+0.5*mad(fBeat) | fBeat<median(fBeat)-0.5.*mad(fBeat))=0./0;
fBeat=DaFixer(1:2048,fBeat,'greater',median(fBeat)+2*mad(fBeat));

%Low pass filter:
fc = 0.05;
fs = 1; %1/median(diff(dte*86400));

[b,a] = butter(6,fc/(fs/2));
%freqz(b,a)
fBeatF_0_5 = filter(b,a,fBeat);


figure;
plot((dte-dte(1))*86400,fBeat,'.')
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
ylabel('Beat frequency [Hz]','Fontname','Times','Fontsize',14);


yB=fft(fBeat);
veloc=angle((yB));
yB=abs(yB/length(fBeat));

yB=yB(1:length(fBeat)/2+1);
yB(2:end-1)=2*yB(2:end-1);
veloc=angle(yB);
yB=abs(yB/length(fBeat));





