%Script for decoding the SPI data:

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

%Identify the file and create a hdf handle:

nfo5=hdf5info('triangle_roof_23_03_18.hdf');   %Triangle modulation & velocity test

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
    timestep=(diff(data6-data6(1)));
    x1=data2(:,i);
    x2=data3(:,i);
    
    %Normalise the X1 and X2 data using the Blkmnn filter:
    alpha=0.16;
    a0=(1-alpha)/2;
    a1=1/2;
    a2=alpha/2;
    n=size(x2,1);
    BlkWin = (a0 - a1*cos(2*pi*(0:n-1)/(n-1))  + a2*cos(4*pi*(0:n-1)/(n-1)))';
    x1f=BlkWin.*x1;
    x2f=BlkWin.*x2;
    
    %find the FFT of x, and plot it:
    yI=fft(x1f);
    yI=abs(yI/length(x1));
    yI=yI(1:length(x1)/2+1);
    yI(2:end-1)=2*yI(2:end-1);
    yQ=fft(x2f);
    yQ=abs(yQ/length(x2));
    yQ=yQ(1:length(x2)/2+1);
    yQ(2:end-1)=2*yQ(2:end-1);
    
    c=3e8;
    Bsweep=abs(data5(i)-data4(i))*1e9; %Bandwidth of the sweep (in Frequency)
    Tsweep=data7(i)/1000;   %convert Tsweep from ms to seconds
    f=(c/2)*Tsweep/Bsweep.*(0:length(yI)-1)*1000; %Range [m] (/1000 to convert from Nyquist to real Fq)
    %Find the top three peaks:
    [Sdist(i,1:3),Slvl(i,1:3)]=WellenzahlLocMax(nanmean([yQ yI],2),f);    
end


%Apply a 1.5 metre threshold:
Slvl(Sdist<1.5)=0;
%Identify the highest level peaks above 1.5 metre threshold:
[~,LvlMax]=max(Slvl,[],2);
LvlMax=[LvlMax==1 LvlMax==2 LvlMax==3]; %convert it to a logical
Distance=Sdist(LvlMax); 


%Figure of Distance timeseries:
%===========================
figure
plot((dte-dte(1))*86400,detrend(Distance),'.')
ylabel('distance [m]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[-5 5])



%Convert the distance back into beat frequency
fBeat=Sdist(LvlMax)*Bsweep/Tsweep/c*2/1000;
K=c/122e9; %Wavelength of the Radar signal:
Velocity= [0;diff(fBeat*1000)*K/2];
VelFilt=runmean(Velocity,10);
%Plot Velocity:
figure
plot((dte-dte(1))*86400,detrend(VelFilt),'.')


%let's try Integrating the velocity and see what we get:
Fs=2048/((dte(end)-dte(1))*86400);
SigmaV=VelFilt*(1/Fs);
Dis4mV=SigmaV.*0;
for u=76:size(VelFilt,1)
Dis4mV(u)=sum(SigmaV(76:u));
end
figure
plot((dte-dte(1))*86400,detrend(Distance),'.')
hold on;
plot((dte-dte(1))*86400,detrend(Dis4mV),'.')
ylabel('distance [m]','Fontname','Times','Fontsize',14);
xlabel('time elapsed [sec]','Fontname','Times','Fontsize',14);
set(gca,'yscale','linear','fontsize',14,'FontName',...
    'Times','linewidth',1.2,'ylim',[-5 5])











fBeat(diff(fBeat)> 3*mad(diff(fBeat)) | diff(fBeat)<-3*mad(diff(fBeat)))=0./0;
%
%fBeat(fBeat>median(fBeat)+0.5*mad(fBeat) | fBeat<median(fBeat)-0.5.*mad(fBeat))=0./0;
fBeat=DaFixer(1:2048,fBeat,'greater',median(fBeat)+mad(fBeat));

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
yB=abs(yB/length(fBeat));
yB=yI(1:length(fBeat)/2+1);
yB(2:end-1)=2*yB(2:end-1);





