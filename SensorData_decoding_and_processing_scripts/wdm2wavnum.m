function [wdm2D] = wdm2wavnum(Wavenumber,Direction,Frequency,Amplitude,ElevationVariance,Options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Debug
% Wavenumber=data.wdm.Wavenumber;
% Direction=data.wdm.Direction;
% Frequency=data.wdm.Frequency;
% Amplitude=data.wdm.Amplitude;
% Options=data.wdm.options;
% ElevationVariance=mean(var(data.wdm.elevation));
%Preamble

nv=Options.nv; %nv = Number of voices in Wavelet transformation.
WDM_NormFactor=1.03565;
Kmax=7; %set the maximum wavenumber, usually= fix(2*pi/diameter of array)
%TODO!! write a Kmax calculator, using sensorLocMat &
%data.tag.vn100.ship.dx & .dy data.
%
%Prepare the Wavenumber data.
Kfac=120/Kmax;
K=round(Wavenumber*Kfac); %put wavenumbers in terms of
ik=find(K > Kmax*Kfac); %find wavenumbers greater than K-maximum
K(ik)=(Kmax*Kfac+1)*ones(size(ik)); %Set K(K>Kmax)=Kmax
%Prepare Wavenumber arrays:
Kn=(1:Kfac*Kmax+1)./Kfac;
KK=zeros(360,121); %new wavenumber array with individual 1 degree direction bins
df=Frequency*log(2)/nv; %Frequency step sizes

%Make sure the direction range is (0 360] i.e.  0<DIR<=360
Direction=mod360(Direction);
%Direction=mod(Direction,360); Direction(Direction==0)=360;
ii=find(Direction == 0); Direction(ii)=Direction(ii)+360;




%Run a loop to calculate the directional content of the wavenumbers
for j = 1:size(Direction,2)%Kmax*Kfac+1
    ik=[];
    ik=find(K==j); %find the wavenumber per j-th bin
    if ~isempty(ik)
        for di = 1:360
            id=[];
            id=find(Direction(ik)==di); %find the direction indices for di-th bin
            KK(di,j)=KK(di,j) + sum(Amplitude(ik(id)).^2)/nv*WDM_NormFactor/size(Amplitude,1);
        end
    end
end

KK=(KK./(ones(360,1)*Kn))*Kfac*180/pi;
corr=ElevationVariance/(sum(sum(KK).*Kn)/Kfac*pi/180); 
KK=KK.*(corr);
Hs=4*sqrt(sum(sum(KK).*Kn)/Kfac*pi/180); 


%Calculate the Power spectral density:
E=zeros(360,size(Frequency,2))./0;
kMEAN=zeros(360,size(Frequency,2))./0;
kMAX=zeros(360,size(Frequency,2))./0;
kAApMEAN=zeros(360,size(Frequency,2))./0;
for f=1:length(Frequency)
    for di=1:360
        idi=find(Direction(:,f)==di);
        E(di,f)=sum(Amplitude(idi,f).^2)/nv*WDM_NormFactor/size(Amplitude,1);
        kMEAN(di,f)=nanmean(Wavenumber(idi,f));
        kAApMEAN(di,f)=sum((Amplitude(idi,f).*conj(Amplitude(idi,f))).*Wavenumber(idi,f))...
            /sum(Amplitude(idi,f).*conj(Amplitude(idi,f)));
        if ~isempty(idi)
            [Emax,Dirmax]=nanmax(Amplitude(idi,f).*conj(Amplitude(idi,f)));
            kMAX(di,f) = Wavenumber(idi(Dirmax),f);
        else
            kMAX(di,f) = 0;
        end
    end 
end
%Convert to spectral density as follows:
E=(E./(ones(360,1)*(df.*Frequency)))*180/pi;
corr=ElevationVariance/(sum(sum(E).*df.*Frequency)*pi/180);
Hs=4*sqrt(sum(sum(E).*df.*Frequency)*pi/180);
% Set Nan values of kMEAN & kAApMEAN to zero:
kMEAN(isnan(kMEAN))=0;
kAApMEAN(isnan(kAApMEAN))=0;
%;
figure;
plot(Frequency,sum(E(:,:)).*Frequency*pi/180);
%Plot the frequency Directional data:
figure;
contour(1:360,Frequency,E'.*(Frequency'.^4 *ones(1,360)));

%Plot the wave numbers:
figure;
loglog(Kn,sum(KK).*Kn*pi/180,'.-'); grid on
xlabel('wavenumber [m^{-1}]')
ylabel('spectral density [m^3]')
%Energy spectrum contour plot
figure
contour(1:360,Kn,KK'.*(Kn'.^1*ones(1,360)));grid
title('Energy spectrum')
ylabel('wavenumber [m^{-1}]')
xlabel('direction [degrees]')
legend('energy weighted average K','K max','K average')
xlabel('Frequency [degrees]')
%Plot KAApMEAN, kMAX, kMEAN
figure;
plot(Frequency,mean(kAApMEAN(:,:)));
hold on
plot(Frequency,mean(kMAX(:,:)));
plot(Frequency,mean(kMEAN(:,:)));
title('Wavenumber Information per Frequency Bin')
legend('energy weighted average K','K max','K average')
xlabel('Frequency [degrees]')
%Output:
wdm2D.K_max=kMAX;
wdm2D.K_mean=kMEAN;
%wdm2D.K=Kn;
wdm2D.F_spectrum=E;
wdm2D.F=Frequency;
wdm2D.Direction=(1:360);
end


function x=mod360(x)% Changed to make mod360(360) = 360 not 0.
%  function x = mod360(x)
%  Given array of angles x [deg], returns x in range 1-360. i.e. mod(x,360).
x = mod(x,360);
x(x==0) = 360;
end




