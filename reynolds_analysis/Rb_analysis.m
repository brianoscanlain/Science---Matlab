%Analysis of the reynolds numbers mentioned by Zhao and toba 2001, Sugihara
%2007, Godjin and Murphy 2010.
%
%
%
%


%Load data
%--------
cd('C:\Users\Brian\Desktop\Whitecap_Project\KnorrData');
load merged_metdata_knorr201.mat

%remove the Nan's
data_met(28365:28370,30)=linspace(data_met(28365,30),data_met(28370,30),6);
data_met(28365:28370,5 )=linspace(data_met(28365,2 ),data_met(28370,2 ),6);

% Run a runmean on the data 60 minute window for each 1 minute value
metDTE=data_met(:,1);              %relative datenum
PRES=data_met(:,5);                %pressure  
SSS=data_met(:,30);                %Sea surface Salinity
clear data_met header

%Remove the spikes in the Salinity data
SSS=DaFixer(metDTE,SSS,'less',32,[734693 734696]);



cd('C:\Users\Brian\Desktop\Whitecap_Project\Wave_data\');
load('knorrwavedata.mat', 'brian_hs', 'brian_dte_filt','brian_Tp','ec_dte','ec_hs','brian_wa','ec_wa');

cd('C:\Users\Brian\Desktop\Whitecap_Project\Images\');
load('data.mat', 'dte','tws','Wa','airt','sst','W','Wb');

brian_Tp=brian_Tp';
Hs=zeros(207,1);
interval=(1/144)*ones(207,1);
Hf=zeros(207,1);
Va=zeros(207,1);
Vw=zeros(207,1);
sss=zeros(207,1);
pres=zeros(207,1);
Rb=zeros(207,1);
Rh1=zeros(207,1);
Rh2=zeros(207,1);
wavAng=zeros(207,1);

[cdn,u]=cdntc(tws,10,airt);
ustar=u./(u.^(1/2)); %This is wrong, ==cdn*(u)^2

for i=1:length(dte)
    
    if dte(i)<max(brian_dte_filt)+(1/24)   % make sure the dte is within range + 1 hour
        
       bin=0;
        
        
     while  sum(bin)<1  
        bin=brian_dte_filt>dte(i)-interval(i,1) & brian_dte_filt<=dte(i)+interval(i,1);
        interval(i,1)=interval(i,1) + (1/144);
     
     end
     
     Hs(i,1)=nanmean(brian_hs(bin));
     Hf(i,1)=(2*pi)/nanmean(brian_Tp(bin));  
     wavAng(i,1)=nanmean(brian_wa(bin));
     
     elseif dte(i)<max(ec_dte)+(1/24)
        bin=0;
        
        while  sum(bin)<1  
        bin=ec_dte>dte(i)-interval(i,1) & ec_dte<=dte(i)+interval(i,1);
        interval(i,1)=interval(i,1) + (1/144);
     
        end
      Hs(i,1)=nanmean(ec_hs(bin));
      wavAng(i,1)=nanmean(ec_wa(bin));
      Hf(i,1)=nan;
      
    else
        Hs(i,1)=nan;
        Hf(i,1)=nan;
        wavAng(i,1)=nan;
        
    end
    
    
    
    %calculate the sss and pres values:
    metind=metDTE<=dte(i)&metDTE>dte(i)-1/144;
    sss(i,1)=nanmean(SSS(metind));
    pres(i,1)=nanmean(PRES(metind));
   
    %Calculate the kinematic viscosities 
    Va(i,1)=viscair(airt(i,1));
    Vw(i,1)=sw_visc(sss(i,1),sst(i,1),pres(i,1));
    
   

   %Rb(i,1)= ((ustar(i,1))^2)/(Va(i,1)*Hf(i,1));
   Rh1(i,1)=ustar(i,1)*Hs(i,1)/Va(i,1);
   Rh2(i,1)=ustar(i,1)*Hs(i,1)/Vw(i,1);
   Rb(i,1)=((ustar(i,1))^2)/(Va(i,1)*Hf(i,1));
end


clear V* SSS PRES Hf* met* p s* i cdn* u* tws int* bin brian*


%do a fit for the three wave breaking parameters::
%
%  1.  Rb
%  2.  Rh (using kinematic air viscosity)
%  3.  Rh (using kinematic water viscosity)
%
%      W = a ( parameter ) ^(b)
%
%  A nested loop is run to find the highest values of Rsquared



%Remove the nan data:
ind=~isnan(Rh1);
ind(1)=0;


a1=linspace( 1e-7,1e-2,1000);

b1=linspace( 0.5,2,1000);
    
denomA=(var(Wa(ind))*(length(Wa(ind))-1));
denomB=(var(Wb(ind))*(length(Wb(ind))-1));
denomC=(var(W(ind))*(length(W(ind))-1));


Rsq_A(1:1000,1:1000)=0;
Rsq_B(1:1000,1:1000)=0;
Rsq(1:1000,1:1000)=0;



for iii=1:1000
    
    for ii=1:1000
       
    Rsq_A(iii,ii)= 1- (sum(    (Wa(ind) - a1(iii).*((Rh1(ind))./(1000)).^(b1(ii))).^2 ) ) /denomA;
    Rsq_B(iii,ii)= 1- (sum(    (Wb(ind) - a1(iii).*((Rh1(ind))./(1000)).^(b1(ii))).^2 ) ) /denomB;
    Rsq(iii,ii)= 1- (sum(    (W(ind) - a1(iii).*((Rh1(ind))./(1000)).^(b1(ii))).^2 ) ) /denomC;
   
   
    end
    
    
end
close all;
plot(b1,max(Rsq_A),'b',b1,max(Rsq_B),'g',b1,max(Rsq),'k');
set(gca,'YLim',[-1 1]);
set(gca,'XLim',[b1(1) b1(end)]);
title('b coefficients');
[Ca,Ia]=max(max(Rsq_A));
[Cb,Ib]=max(max(Rsq_B));
[Co,Io]=max(max(Rsq));

hleg=legend(['Wa (b=' num2str(b1(Ia)) ' )  R^2= ' num2str(Ca)],['Wb (b='...
    num2str(b1(Ib)) ' )  R^2= ' num2str(Cb)],['W (b=' num2str(b1(Io))...
    ' )  R^2= ' num2str(Co)]);
xlabel('b values'); ylabel('R^2 coefficients');

figure;
plot(a1./1000,max(Rsq_A'),'b',a1./1000,max(Rsq_B'),'g',a1./1000,max(Rsq'),'k');
set(gca,'YLim',[-1 1]);
set(gca,'XLim',[a1(1)/1000 a1(end)/1000]);
title('a coefficients');
[Ca,Ia]=max(max(Rsq_A'));
[Cb,Ib]=max(max(Rsq_B'));
[Co,Io]=max(max(Rsq'));

hleg=legend(['Wa (a=' num2str(a1(Ia)/1000) ' )  R^2= ' num2str(Ca)],['Wb (a=' ...
    num2str(a1(Ib)/1000) ' )  R^2= ' num2str(Cb)],['W (a=' num2str(a1(Io)/1000)...
    ' )  R^2= ' num2str(Co)]);
xlabel('a values'); ylabel('R^2 coefficients');
