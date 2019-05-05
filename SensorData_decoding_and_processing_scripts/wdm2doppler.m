function [ output_args ] = wdm2doppler( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Debug:
K_spec=data.WDM2D.K_mean;
Frequency=data.WDM2D.F;
Uship=-mean(data.gps.gprmc.speedOverGround_knots);   %Negative ship speed to maintain same convention as Fenton and Mckee 1990
Direction=data.WDM2D.Direction;
kMEAN=data.WDM2D.K_mean;
kMAX=data.WDM2D.K_max;
E_spec=data.WDM2D.F_spectrum;

kMEAN16b=data.WDM2D16b.K_mean;
kMAX16b=data.WDM2D16b.K_max;
E_spec16b=data.WDM2D16b.F_spectrum;
%---------
%Preamble:
%---------
g = 9.81; % m/s^2 gravity
knots2ms=0.514444;
%---------

%We calculate the three Frequencies of interest:
%
%   |---------------*--------------|  (i-th Frequency bin)
%  fLow(i)    Frequency(i)      fUp(i)
%
%Calculate the Frequency bin ranges, and  upper and lower bounds:
[fUp, fLow, dfInfo] = freqPMdf(Frequency);
%Calculate the Frequency spacing for each bin:
dFreq=fUp-fLow;
%
% plot(Frequency,Frequency,'rs')
% hold on
% plot(fUp,fUp,'go')
% plot(fLow,fLow,'b*')
%
%Convert the WDM wavenumbers to Frequency space (assuming Infinite depth)
Kf_spec=sqrt(g.*K_spec)/(2*pi);



%----------------------
%Prepare K-plus and K-minus conversion matrices:
%----------------------
UshipDir=(knots2ms*Uship*cos(Direction*pi/180)); %Calculate the directional-dependant speed:
%Calculate the opposite direction to "Direction" vector
DirOpp = mod(Direction-180,360);
DirOpp(DirOpp==0)=360;
%----------------------
BlueShift=(((4*pi^2)*4*Frequency'.^2)/g) ./ (1 - sqrt( 1- (2*pi)*4/g*Frequency'.*UshipDir  )   ).^2;
RedShift=(((4*pi^2)*4*Frequency'.^2)/g) ./ (1 + sqrt( 1- (2*pi)*4/g*Frequency'.*UshipDir  )   ).^2;
K_0=((2*pi*ones(length(Frequency),length(Direction)).*Frequency').^2)/g;
%Calculate Fcrit: (point at which the (complex U>0) K- & k+ curves meet).
Fcrit=abs(sqrt(g./(2*abs(UshipDir).^2*g)-g./(2*abs(UshipDir)).^2.*abs(UshipDir))/2/pi);
%----------------------


% Send the Wavenumber and Energy components from directional bins [89,90.5],
% [90.5,91], [269 270.5] & [270.5 271] to their nearest bins (i.e. [88]
% [92] [268] [272]. This is because these angles are extremes and it is
% extremely difficult to accurately deduce if the waves are following,
% approaching, and also if they should be blue shifted or red shifted.
%
% According to my understanding, Blue shifting waves, which were initially
% defined as approaching the ship, is a consequence of the waves "following"
% the ship, thus, Blue shift should change the direction of the propogation
% of the waves, as well as relocate the Energy and wavenumber content at a
% particular observed Frequency to the Doppler-corrected Frequency
% (estimate of the Frequency of the waves in a stationary frame of
% reference).
%
% Conversely, Red shifting waves is the expected resolution for waves
% approaching the ship (so long as they require a Doppler correction in the
% first place). Thus the initial direction, which was assigned by the
% wavelet transform algorithm, is correct and doe not need to be changed. As
% per the Blue shift, the red shift should relocate the Energy and wavenumber
% content at a particular observed Frequency to the Doppler-corrected
% Frequency (estimate of the Frequency of the waves in a stationary frame of
% reference).
E_spec([88 92 268 272],:)=[(E_spec(88,:) + E_spec(89,:) + 0.5*E_spec(90,:));...
    (E_spec(92,:) + E_spec(91,:) + 0.5*E_spec(90,:));...
    (E_spec(268,:) + E_spec(269,:) + 0.5*E_spec(270,:));...
    (E_spec(272,:) + E_spec(271,:) + 0.5*E_spec(270,:))];
K_spec([88 92 268 272],:)=[(K_spec(88,:) + K_spec(89,:) + 0.5*K_spec(90,:));...
    (K_spec(92,:) + K_spec(91,:) + 0.5*K_spec(90,:));...
    (K_spec(268,:) + K_spec(269,:) + 0.5*K_spec(270,:));...
    (K_spec(272,:) + K_spec(271,:) + 0.5*K_spec(270,:))];
%I noticed that I can distinguish between BS+ RS+ BS- RS- in Frequency
%space using an additional tuning parameter. It shifts the Frequencies to a
%higher Frequencies. This has no theoretical basis, it will need to be
%investigated.
DecisionTuningParam=sqrt(g/4/pi/pi);
Kf_spec=sqrt(g/4/pi/pi).*K_spec.^(0.5)./DecisionTuningParam;
%Prepare new Doppler-corrected matrices:
F_dopSft=ones(size(Direction,2),1).*Frequency;
Dir_dopSft=ones(1,size(Frequency,2)).*Direction';
E_dopShft=zeros(size(E_spec,1),size(E_spec,2));
K_dopShft=zeros(size(K_spec,1),size(K_spec,2));
dirOp=mod(Direction+180,360); dirOp(dirOp==0)=360; %Opposite direction vector


%----------------------
%Decision nested loops:
%----------------------
for dir=1:180
    % We now have to decide which Doppler shift to apply. It can be
    % BlueShift(theta), RedShift(theta), BlueShift(theta+180), or RedShift(theta+180).
    %    For each directional bin, we will assess the wavenumber (e.g. using
    %    kMAX or kMEAN across the frequency range.
    
    if dir==89 || dir==90 || dir==91
        %Do nothing, just skip these bins due to the waves being orthogonal
        %to the direction of the ship. We transfer the K and energy content
        %from the bins to the nearest 'acceptable' neighbors so no
        %information should be lost.
    else
        %Determine the most suitable Doppler correction solution for each
        %of the Frequency bins.
        figure
        plot(Frequency,Kf_spec(dir,:))
        hold on
        plot(Frequency,Kf_spec(dirOp(dir),:))
        plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,dir)))
        plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,dir)))
        plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,dir+180)))
        plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,dir+180)))
        set(gca,'YLim',[0 0.3])
        
        for f=1:length(Frequency)
            if Kf_spec(dir,f)>0 || Kf_spec(dirOp(dir),f)>0 %K(di,f) > zero, Doppler correct the content:
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                %                  Freq < Fcrit
                %==========================================================
                if Frequency(f)<Fcrit(dir) % For Fq below Fcrit:
                    BSneg=sqrt(g/4/pi/pi).*sqrt(RedShift(f,dirOp(dir))); %compute the BS- threshold
                    RSneg=sqrt(g/4/pi/pi).*sqrt(BlueShift(f,dirOp(dir))); %In freq space
                    %Analyse the content for Direction(dir) first:
                    %=============================================
                    if Kf_spec(dir,f)>0
                        if Kf_spec(dir,f)>BSneg %K > BS-
                            %Not RS+
                            %check the 3 remaining scenarios:
                            if Kf_spec(dir,f)>mean([RSneg BSneg]) % K>Kcrit
                                %Not RS+ or BS-
                                %check the remaining two
                                if Kf_spec(dir,f)>RSneg          %K>RS-
                                    %not RS+, BS- or RS-... must be BS+!
                                    %--------------
                                    %BS+ correction
                                    %--------------
                                    F_dopSft(dir,f)=sqrt(BlueShift(f,dir)).*sqrt(g/4/pi/pi);
                                    Dir_dopSft(dir,f)=dirOp(dir);
                                    %--------------
                                    
                                else                           %Kcrit < K < RS-
                                    %--------------
                                    %RS- correction
                                    %--------------
                                    F_dopSft(dir,f)=sqrt(RedShift(f,dirOp(dir))).*sqrt(g/4/pi/pi);
                                    Dir_dopSft(dir,f)=dirOp(dir);
                                    %--------------
                                end
                            else% BS- < K < Kcrit
                                %--------------
                                %BS- correction
                                %--------------
                                F_dopSft(dir,f)=sqrt(BlueShift(f,dirOp(dir))).*sqrt(g/4/pi/pi);
                                %--------------
                                
                            end
                        else %  0< K < BS-,
                            %assign RS+ correction
                            %--------------
                            %RS+ correction
                            %--------------
                            F_dopSft(dir,f)=sqrt(RedShift(f,dir)).*sqrt(g/4/pi/pi);
                            %--------------
                        end
                    end
                    %=============================================
                    %Analyse the content for dirOp(dir) now: (waves
                    %considered to be opposite direction:
                    %=============================================
                    if Kf_spec(dirOp(dir),f)>0
                        if Kf_spec(dirOp(dir),f)>BSneg %K > BS-
                            %Not RS+
                            %check the 3 remaining scenarios:
                            if Kf_spec(dirOp(dir),f)>mean([RSneg BSneg]) % K>Kcrit
                                %Not RS+ or BS-
                                %check the remaining two
                                if Kf_spec(dirOp(dir),f)>RSneg          %K>RS-
                                    %not RS+, BS- or RS-... must be BS+!
                                    %--------------
                                    %BS+ correction
                                    %--------------
                                    F_dopSft(dirOp(dir),f)=sqrt(BlueShift(f,dir)).*sqrt(g/4/pi/pi);
                                    Dir_dopSft(dirOp(dir),f)=dir;
                                    %--------------
                                    
                                else                           %Kcrit < K < RS-
                                    %--------------
                                    %RS- correction
                                    %--------------
                                    %......
                                    F_dopSft(dirOp(dir),f)=sqrt(RedShift(f,dirOp(dir))).*sqrt(g/4/pi/pi);
                                    Dir_dopSft(dirOp(dir),f)=dir;
                                    %--------------
                                end
                            else% BS- < K < Kcrit
                                %--------------
                                %BS- correction
                                %--------------
                                F_dopSft(dirOp(dir),f)=sqrt(BlueShift(f,dirOp(dir))).*sqrt(g/4/pi/pi);
                                %--------------
                                
                            end
                        else %  0< K < BS-,
                            %assign RS+ correction
                            %--------------
                            %RS+ correction
                            %--------------
                            F_dopSft(dirOp(dir),f)=sqrt(RedShift(f,dir)).*sqrt(g/4/pi/pi);
                            %--------------
                        end
                    end
                    
                    %=============================================
                    
                else                                   %For Fq above Fcrit:
                    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    %                  Freq < Fcrit
                    %==========================================================
                    %There can only be two real solutions: RS+ or BS+
                    RSpos=sqrt(g/4/pi/pi).*sqrt(BlueShift(f,dir)); %calculate the Red shift value at f
                    BSpos=sqrt(g/4/pi/pi).*sqrt(RedShift(f,dir)); %calculate the Blue shift value at f
                    %=============================================
                    %Analyse the content for Direction(dir) first:
                    %=============================================
                    if Kf_spec(dir,f)>0
                        if Kf_spec(dir,f)>=mean([RSpos BSpos]) % k>= middle point between RS and BS
                            %--------------
                            %BS+ correction
                            %--------------
                            F_dopSft(dir,f)=sqrt(BlueShift(f,dir)).*sqrt(g/4/pi/pi);
                            Dir_dopSft(dir,f)=dirOp(dir);
                            %--------------
                        else
                            %--------------
                            %RS+ correction
                            %--------------
                            F_dopSft(dir,f)=sqrt(RedShift(f,dir)).*sqrt(g/4/pi/pi);
                            %--------------
                        end
                    end
                    %=============================================
                    %Analyse the content for dirOp(dir) now: (waves
                    %considered to be opposite direction:
                    %=============================================
                    if Kf_spec(dirOp(dir),f)>0 %Do the same for the 
                        if Kf_spec(dirOp(dir),f)>=mean([RSpos BSpos]) % k>= middle point between RS and BS
                            %--------------
                            %BS+ correction
                            %--------------
                            F_dopSft(dirOp(dir),f)=sqrt(BlueShift(f,dir)).*sqrt(g/4/pi/pi);
                            Dir_dopSft(dirOp(dir),f)=dir;
                            %--------------
                        else
                            %--------------
                            %RS+ correction
                            %--------------
                            F_dopSft(dirOp(dir),f)=sqrt(RedShift(f,dir)).*sqrt(g/4/pi/pi);
                            %--------------
                        end
                        
                    end
                end
                
            else %K(di,f) ==0, no K data recorded for this bin, skip to next
            end
            
        end
    end
end











%Plot the RedShift and BlueShift solutions:
figure
plot(Frequency,(g/4/pi/pi).*sqrt(RedShift(:,[271:360 1:89])))
title(['RedShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
figure
plot(Frequency,(g/4/pi/pi).*sqrt(BlueShift(:,[271:360 1:89])))
title(['BlueShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');




%Plot the solutions for theta=...:
theta=250;
spread=10;
figure
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[theta])),'linewidth',1.2) %This represents waves coming toward the ship
hold on;
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship, Cp<<Uship (ship overtaking waves)
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship Cp
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[theta])),'linewidth',1.2)
eval(['legend(''Red shift ' num2str(theta) '^\circ'',''Red shift ' num2str(mod(theta+180,360)) '^\circ'',''Blue shift ' num2str(mod(theta+180,360)) ' ^\circ'',''Blue Shift' num2str(theta) '^\circ'');']);
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
eval(['title(''Solution for \theta=' num2str(theta) '^\circ'',''fontsize'',14,''fontname'',''times'');']);
DirSpread=mod([theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
DirSpreadOpp=mod(180+[theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
%plot(Frequency,mean(kMAX(DirSpread,:)));
%plot(Frequency,mean(kAApMEAN(DirSpread,:)));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN(DirSpread,:))));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN(DirSpreadOpp,:))));



%Plot the solutions for theta=...:
theta=250;
spread=10;
figure
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[theta])),'linewidth',1.2) %This represents waves coming toward the ship
hold on;
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship, Cp<<Uship (ship overtaking waves)
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship Cp
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[theta])),'linewidth',1.2)
eval(['legend(''Red shift ' num2str(theta) '^\circ'',''Red shift ' num2str(mod(theta+180,360)) '^\circ'',''Blue shift ' num2str(mod(theta+180,360)) ' ^\circ'',''Blue Shift' num2str(theta) '^\circ'');']);
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
eval(['title(''Solution for \theta=' num2str(theta) '^\circ'',''fontsize'',14,''fontname'',''times'');']);
DirSpread=mod([theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
DirSpreadOpp=mod(180+[theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
%plot(Frequency,mean(kMAX(DirSpread,:)));
%plot(Frequency,mean(kAApMEAN(DirSpread,:)));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN16b(DirSpread,:))));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN16b(DirSpreadOpp,:))));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMAX16b(DirSpread,:))));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMAX16b(DirSpreadOpp,:))));

%Plot the solutions for theta=...: %this is wrong but interesting. the
%kmean is divided by a value.
theta=260;
spread=10;
figure
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[theta])),'linewidth',1.2) %This represents waves coming toward the ship
hold on;
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(RedShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship, Cp<<Uship (ship overtaking waves)
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[mod(theta+180,360)])),'linewidth',1.2) %This represents waves following the ship Cp
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(BlueShift(:,[theta])),'linewidth',1.2)
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
%eval(['title(''Solution for \theta=' num2str(theta) '^\circ'',''fontsize'',14,''fontname'',''times'');']);
DirSpread=mod([theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
DirSpreadOpp=mod(180+[theta-spread:theta+spread],360); DirSpread(DirSpread==0)=360;
%plot(Frequency,mean(kMAX(DirSpread,:)));
%plot(Frequency,mean(kMAX(DirSpreadOpp,:)));
%plot(Frequency,mean(kAApMEAN(DirSpread,:)));
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN(DirSpread,:)))./sqrt(g/4/pi/pi)   ,'k--' );
plot(Frequency,sqrt(g/4/pi/pi).*sqrt(mean(kMEAN(DirSpreadOpp,:)))./sqrt(g/4/pi/pi), 'k-*'   );
eval(['title(''Solution for \theta=' num2str(theta) '^\circ'',''fontsize'',14,''fontname'',''times'');']);
eval(['legend(''Red shift ' num2str(theta) '^\circ'',''Red shift '...
    num2str(mod(theta+180,360)) '^\circ'',''Blue shift ' ...
    num2str(mod(theta+180,360)) ' ^\circ'',''Blue Shift' ...
    num2str(theta) '^\circ'',''F(K_{mean}(\theta))'', ''F(K_{mean}(\theta+180))'');']);





figure;
contour([1:360],Frequency,(kMEAN(:,:)'));
set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
title('K_{MEAN} - Wavenumbers per Direction per Frequency','fontsize',14,'fontname','times');
ylabel('Frequency [Hz]','fontsize',14,'fontname','times');
xlabel('Direction [degrees]','fontsize',14,'fontname','times');
hcb=colorbar;
hL = ylabel(hcb,'Wavenumber [m^{-1}]','fontsize',14,'fontname','times');
set(hL,'Rotation',90);


end  % Main Function end


function [K_corrected,E_corrected] = RedShiftFcn [RedShift]
%Subfunction to perform the Doppler Red shift. Red shift considers waves
%approaching the ship.
F_corrected=RedShift(f,dir)
Redshift(
K_corrected=K_spec
K_corrected=E_spec(Freq_obs
end



function [] = BlueShifter []
%Subfunction to perform the Doppler Blue shift
end

































%         %Old coded solution:
%         % K- solution
%         RedShift2 = 1./(2*UshipDir.^2).*( (g-4*pi*UshipDir.*Frequency') -sqrt(g*(g-4*2*pi*UshipDir.*Frequency')) );
%         % K+ solution
%         BlueShift2 =  1./(2*UshipDir.^2).*( (g-4*pi*UshipDir.*Frequency') +sqrt(g*(g-4*2*pi*UshipDir.*Frequency')) );
%
%
%         figure
%         plot(Frequency,(g/4/pi/pi).*sqrt(BlueShift2(:,[271:360 1:89])))
%         title(['BlueShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
%         set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
%         xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
%         ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
%         figure
%         plot(Frequency,(g/4/pi/pi).*sqrt(RedShift2(:,[271:360 1:89])))
%         title(['RedShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
%         set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
%         xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
%         ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
%
%
%         figure
%         plot(Frequency,(g/4/pi/pi).*sqrt(BlueShift2(:,[271:360 1:89 91:269])))
%         title(['BlueShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
%         set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
%         xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
%         ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
%         figure
%         plot(Frequency,(g/4/pi/pi).*sqrt(RedShift2(:,[271:360 1:89 91:269])))
%         title(['RedShift \pm 89^\circ, U=' num2str(Uship) '  m/s'])
%         set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
%         xlabel('F_{obs} [Hz]','fontsize',14,'fontname','times');
%         ylabel('F_{waves} [Hz]','fontsize',14,'fontname','times');
%
%
%         % K- solution
%         RedShift2 = 1./(2*UshipDir.^2).*( (g-4*pi*UshipDir.*Frequency') -sqrt(g*(g-4*2*pi*UshipDir.*Frequency')) );
%         % K+ solution
%         BlueShift2 =  1./(2*UshipDir.^2).*( (g-4*pi*UshipDir.*Frequency') +sqrt(g*(g-4*2*pi*UshipDir.*Frequency')) );
%         %K solutions for the lower and upper bin limits:
%         RedShiftPdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fLow') -sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         RedShiftMdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fLow') -sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         BlueShiftPdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fUp') + sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         BlueShiftMdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fLow') + sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         % K0 on station solution
%         K_0=((2*pi*ones(length(Frequency),length(Direction)).*Frequency').^2)/g;
%         % Calculate shifted Frequencies and shifted df (bin size values)
%         F_RedShift = sqrt(g*RedShift)/2/pi; % dopler shifted frequecy using (-) solution %  dopler shifted frequecy using (-) solution
%         F_BlueShift = sqrt(g*BlueShift )/2/pi; % dopler shifted frequecy using (-) solution %
%         %Correct the orthogonal wave singularities (spikes in RedShift and BlueShift due
%         %to the denominator U -> -> 0 m/s. This occurs for waves at orthogonal
%         %angles to the ship's heading, when the observed speed between the ship and
%         %waves become zero. We correct this by applying the K_0 relationship to
%         %RedShift and BlueShift when U<0.1 m/s:
%         RedShift(:,(abs(UshipDir)<0.1)')=K_0(:,(abs(UshipDir)<0.1)');
%         BlueShift(:,(abs(UshipDir)<0.1)')=K_0(:,(abs(UshipDir)<0.1)');
%
%         df_RedShift=(sqrt(g*RedShiftPdf2D)/2/pi - sqrt(g*RedShiftMdf)/2/pi);	% doppler shifted frequency bin width corresponding to df.
%         df_BlueShift=(sqrt(g*BlueShiftPdf)/2/pi - sqrt(g*BlueShiftMdf)/2/pi);	% doppler shifted frequency bin width corresponding to df.
%         %Calculate The critical Frequency:
%         SMx.fcrit2D = (sqrt(g./(2*abs(SMx.u2D)).^2*g)-g./(2*abs(SMx.u2D)).^2.*abs(SMx.u2D))/2/pi;
%         %
%         %Correct the RedShift spikes cause by U -> -> 0 m/s
%
%         RedShiftMdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fLow') -sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         BlueShiftPdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fUp') + sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%         BlueShiftMdf = 1./(2*UshipDir.^2).*((g-4*pi*UshipDir.*fLow') + sqrt(g*(g-8*pi*UshipDir.*Frequency')) );
%
%
%         %----------------------
%
%


