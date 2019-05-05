clc; clear all; close all;
addpath('scripts')
%The script processes tar.gz files placed into the /Data/raw/ directory.
% bow171101_16.tar.gz (example tar file)

%==========================================================================
%                        Untar Routine
%==========================================================================
%Preamble
cleanup=1; %set to one if you want to clean up
string_pattern= '*.tar.gz';
tarDir= ['.' filesep 'Data' filesep 'raw' filesep];
% make the target directory if non existent
if  ~exist('./Data/mat/','dir') 
    mkdir('./Data/mat/'); 
end	
%Scan for files in raw directory:
fwave=dir([tarDir string_pattern]);
%disp(char(fwave.name));
disp(['I found ' num2str(length(fwave)) ' ' string_pattern ' files '] );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:size(fwave,1)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   disp(['Untaring ' fwave(ii).name]);
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   untar( [tarDir fwave(ii).name] , tarDir(1:end-1)); % don't use the '/' at the end
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if cleanup==1
	disp(['Removing ' fwave(ii).name]);
   	delete([tarDir fwave(ii).name]); % remove the tar file
   end
end
%clean workspace:
clear fwave;
%==========================================================================




%==========================================================================
%                   Decoding routine
%==========================================================================
%Preamble:
portSymeoTag='.symeo1';
bowSymeoTag='.symeo2';
stbdSymeoTag='.symeo3';
portVNtag='.VN1004';
bowVNtag='.VN1005';
stbdVNtag='.VN1006';
matDir=['.' filesep 'Data' filesep 'mat' filesep];
searchTag='.symeo1';
gpsTag='.GPS8';
%scan for file names:
fname=dir([tarDir '*' searchTag]);

for i=1:length(fname)
    %remove the 'symeo1' from the name:
    fname(i).name=fname(i).name(1:end-length(searchTag));
    %--------------------------
    %Symeo Decoding subroutine:
    %--------------------------
    disp(['Decoding '  fname(i).name '* Symeo files now...']);
    eval(['!python "' pwd filesep 'symeo2mat.py" ' pwd tarDir(2:end)  ' *symeo*']);
    %-------------------------- 
    %VN100 Decoding subroutine:
    %--------------------------
    disp(['Decoding '  fname(i).name '* VN100 files now...']);
    eval(['!python "' pwd filesep 'vn2mat.py" ' pwd tarDir(2:end)  ' *vn*']);
    %move mat files to mat folder:
    movefile([pwd tarDir(2:end) fname(i).name '*.mat'],[pwd matDir(2:end)]);
    %--------------------------
    %GPS decoding subroutine:
    %--------------------------
    % GPS ASCII filtering:
    %---------------------
    fid=fopen([tarDir fname(i).name gpsTag],'r');
    str=textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    delete([tarDir fname(i).name gpsTag]);
    %remove non-ascii values:
    str=char(str{1});
    %write filtered text to file:
    fig=fopen([tarDir fname(i).name gpsTag],'w+'); 
    for ln=1:size(str,1)
        sentence=str(ln,str(ln,:) >= char(31) & str(ln,:) <= char(127) | str(ln,:)==char(09));
        if sum(~isspace(sentence))~=0
            fprintf(fig,'%s\r\n',sentence);
        end
    end   
    fclose(fig);
    %--------------------
    %Read the GPS data:
    gps=gpsToMat([pwd tarDir(2:end) fname(i).name gpsTag]);
    %Load in the symeo and VN100 data:
    %port:
    port=load([pwd matDir(2:end) fname(i).name portSymeoTag '.mat']);
    port.VN100=load([pwd matDir(2:end) fname(i).name portVNtag '.mat']);
    port.VN100=port.VN100.vn100;
    %bow:
    bow=load([pwd matDir(2:end) fname(i).name bowSymeoTag '.mat']);
    bow.VN100=load([pwd matDir(2:end) fname(i).name bowVNtag '.mat']);
    bow.VN100=bow.VN100.vn100;
    %starboard (stbd)
    stbd=load([pwd matDir(2:end) fname(i).name stbdSymeoTag '.mat']);
    stbd.VN100=load([pwd matDir(2:end) fname(i).name stbdVNtag '.mat']);
    stbd.VN100=stbd.VN100.vn100;
    %Formulate good data name and save relevant data:  (yyyymmddTHHMM.mat)
    save([pwd matDir(2:end) '20' fname(i).name(1:6) 'T' fname(i).name(8:11)...
        '.mat'],'gps','port','bow','stbd');
    if cleanup==1
        %cleanup folder:
        delete([pwd matDir(2:end) fname(i).name portSymeoTag '.mat']);
        delete([pwd matDir(2:end) fname(i).name portVNtag '.mat']);
        delete([pwd matDir(2:end) fname(i).name portVNtag '.pckl*']);
        delete([pwd matDir(2:end) fname(i).name bowSymeoTag '.mat']);
        delete([pwd matDir(2:end) fname(i).name bowVNtag '.mat']);
        delete([pwd matDir(2:end) fname(i).name bowVNtag '.pckl*']);
        delete([pwd matDir(2:end) fname(i).name stbdSymeoTag '.mat']);
        delete([pwd matDir(2:end) fname(i).name stbdVNtag '.mat']);
        delete([pwd matDir(2:end) fname(i).name stbdVNtag '.pckl*']);
    end
    %clean workspace:
    clear gps* bow* stbd* port* i* se* str* ln fi* fname ans;
end
%==========================================================================






%==========================================================================
%                   Data Filtering routine
%==========================================================================
matDir=['.' filesep 'Data' filesep 'mat' filesep];
passFiltNum=4;
runMeanWindowSize=35;  %setting it to 27, so that we have a window size of 1sec worth of data
zScoreM=1; %outlier threshold for runmean filter
zScoreDiff=1; %outlier threshold for spatial differential filter

plotting=1;
%Generating common timestamp series
fname=dir([pwd matDir(2:end) '*.mat']);
for i=1:length(fname)
   data=load([pwd matDir(2:end) fname(i).name]); 
   %Report on the quality of the raw data  (make a function for doing this)
   %append the report (string output) to the data struct.
   %data.rawDataQuality.portSymeo.
   
   %remove the zero & nan readings and replace them with interpreted results:
   for tag={'port','bow','stbd'}
       %interpolation prepassing: Initial interp of zero readings and Nan's values:
       %^^ use 'extrap' option incase there are Nan's at start or end of data.
       %SYMEO:
      eval(['data.' char(tag) '.symeo.distance_filt=interp1(find(data.'...
       char(tag) '.symeo.distance>0.1 & isnan(data.' char(tag) '.symeo.distance)==0 & '...
       '[0 (diff(data.' char(tag) '.symeo.distance))]<std(diff(data.' char(tag) ...
       '.symeo.distance))) ,data.' char(tag) '.symeo.distance(data.'...
       char(tag) '.symeo.distance>0.1 & [0 (diff(data.' char(tag) ...
       '.symeo.distance))]<std(diff(data.' char(tag) '.symeo.distance))'...
       ' & ~isnan(data.' char(tag) '.symeo.distance)),(1:length(data.'...
       char(tag) '.symeo.distance)),''linear'',''extrap'');']); 
       %Redo the symeo and VN100 moxatimes:  (I noticed that they are a bit
       %sporadic, and this throws off any interpolating that is done based
       %on these time series).
       eval(['data.' char(tag) '.VN100.moxatime=linspace(data.' char(tag)...
           '.VN100.moxatime(1),data.' char(tag) '.VN100.moxatime(end),length(data.'...
           char(tag) '.VN100.moxatime));']);
       eval(['data.' char(tag) '.symeo.moxatime=linspace(data.' char(tag)...
           '.symeo.moxatime(1),data.' char(tag) '.symeo.moxatime(end),length(data.'...
           char(tag) '.symeo.moxatime));']);
       %
       %VN100 (Filter Nan's):
       for axes=1:3
            %roll pitch yaw:
            eval(['data.' char(tag) '.VN100.ypr_filt(' num2str(axes) ',:)=interp1(find(isnan(data.'...
                char(tag) '.VN100.ypr(' num2str(axes) ',:))==0),data.' char(tag) ...
                '.VN100.ypr(' num2str(axes) ',~isnan(data.' char(tag) '.VN100.ypr(' num2str(axes) ',:))),(1:length(data.'...
                char(tag) '.VN100.ypr(' num2str(axes) ',:))),''linear'',''extrap'');']);  
            %angular velocities:
            eval(['data.' char(tag) '.VN100.ang_filt(' num2str(axes) ',:)=interp1(find(isnan(data.'...
                char(tag) '.VN100.ang(' num2str(axes) ',:))==0),data.' char(tag) ...
                '.VN100.ang(' num2str(axes) ',~isnan(data.' char(tag) '.VN100.ang(' num2str(axes) ',:))),(1:length(data.'...
                char(tag) '.VN100.ang(' num2str(axes) ',:))),''linear'',''extrap'');']);  
            %acceleration:
            eval(['data.' char(tag) '.VN100.acc_filt(' num2str(axes) ',:)=interp1(find(isnan(data.'...
                char(tag) '.VN100.acc(' num2str(axes) ',:))==0),data.' char(tag) ...
                '.VN100.acc(' num2str(axes) ',~isnan(data.' char(tag) '.VN100.acc(' num2str(axes) ',:))),(1:length(data.'...
                char(tag) '.VN100.acc(' num2str(axes) ',:))),''linear'',''extrap'');']);  
       end

       
      if plotting==1
          eval(['figure; plot(data.' char(tag) '.symeo.distance); hold on;']);
          eval(['plot(data.' char(tag) '.symeo.distance_filt);']);
      end
       
       for ii=1:passFiltNum
           %------------------
           %  SYMEO
           %------------------
            %filter out the first derivitive spikes:
            %  (diff(distance) >zScore * std(diff(distance)) == Nan)
            eval(['data.' char(tag) '.symeo.distance_filt([0 (diff(data.' ...
                char(tag) '.symeo.distance))]>=zScoreDiff*std(diff(data.' char(tag) ...
                '.symeo.distance_filt)))=0./0;']);
            %interp Nan's
            eval(['data.' char(tag) '.symeo.distance_filt=interp1(find('...
                'isnan(data.' char(tag) '.symeo.distance_filt)==0),data.' ...
                char(tag) '.symeo.distance_filt(~isnan(data.' char(tag)...
                '.symeo.distance_filt)),(1:length(data.' char(tag) ...
                '.symeo.distance_filt)),''linear'',''extrap'');']); 
            %Compute runmean
            eval(['dist_m=runmean(data.' char(tag) '.symeo.distance_filt ,27);']);
            %set Nan values if  |raw-runmean| > std (runmean)  
            eval(['data.' char(tag) '.symeo.distance_filt((abs(data.' char(tag) ...
            '.symeo.distance_filt-dist_m))>zScoreM*std(runmean( data.' char(tag)... 
            '.symeo.distance_filt,' num2str(runMeanWindowSize) ')))=0./0;']);        
            %interpolate Nan values:
            eval(['data.' char(tag) '.symeo.distance_filt=interp1(find('...
                'isnan(data.' char(tag) '.symeo.distance_filt)==0),data.' ...
                char(tag) '.symeo.distance_filt(~isnan(data.' char(tag)...
                '.symeo.distance_filt)),(1:length(data.' char(tag) ...
                '.symeo.distance_filt)),''linear'',''extrap'');']); 
            
            
            
            if plotting==1
                 eval(['plot(data.' char(tag) '.symeo.distance_filt);']);
            end
       end
%        figure
%        plot(data.port.symeo.moxatime,data.port.symeo.distance_filt);
%        hold on
%        plot(data.port.symeo.moxatime,data.port.symeo.distance_filt,'b.','markersize',5);
%        ccdd=((abs([0 (diff(data.port.symeo.distance))])>=1*std(diff(data.port.symeo.distance_filt))));
%        plot(data.port.symeo.moxatime(ccdd),data.port.symeo.distance_filt(ccdd),'g.','markersize',5)

      
       if plotting==1
           h=legend('raw','remove zeros & Nans','1st filt pass',...
               '2nd filt pass','3rd filt pass', '4th filt pass');
           set(h,'fontname','times','fontsize',12);
           set(gca,'fontname','times','fontsize',12,'linewidth',1.20...
               ,'xticklabel',{},'YLim',[0 30]);
           eval(['title(''Filtering the symeo signal (' char(tag) ...
               ' 20071112T0000.mat)'',''fontname'',''Times'',''fontsize'''...
               ',15);']);
           ylabel('distance (m)','fontname','Times','fontsize',15)
       end
   
   end
   save([pwd matDir(2:end) fname(i).name],'-struct','data'); 
   
end
%==========================================================================
%==========================================================================






%==========================================================================
%       Ship motion correction:
%==========================================================================
matDir=['.' filesep 'Data' filesep 'mat' filesep];
knots2ms=0.514444;
%The translation matrix represents the translation required to point the
%x-axis towards the bow! (not the absolute position of the sensors x axis, instead it would be the negative of this value)
translation_matrix=[0 0 0; 0 0 0; 90 0 180] .* (pi/180); %columns(port bow stbd)  rows(roll pitch yaw)
%for CV:
%translation_matrix=[0 0 0; 0 0 0; -37.2 0 127.2] .* (pi/180); %columns(port bow stbd)  rows(roll pitch yaw)
%z-axis to propogate upwards.
L = [0 0 0; 0 0 0; 0 0 0];	% [port(x,y,z); bow(x,y,z); stbd(x,y,z)]
                            % for now assume symeo is located @ imu, (L=[0,0,0;...])
%Bandpass filter design:
fdegree=4;
Flow=0.03;%Hz
Fhigh=7;%Hz
%
fname=dir([pwd matDir(2:end) '*.mat']);
for i=1:length(fname)
   data=load([pwd matDir(2:end) fname(i).name]); 
   %translate the VN coords to RHS ship coordinate system
   ind=1;
   for tag={'port' 'bow' 'stbd'}      
       
        %Angular rate data (rotation direction is reversed for y & z axes):
        %Note, the rotation of this data results in a rotation of the
        %ship.pitch value (where Y-axis now points to ship's port).
        eval(['data.' char(tag) '.VN100.acc_trans=trans([data.' char(tag) '.VN100.acc_filt(1,:)'...
            '; -data.' char(tag) '.VN100.acc_filt(2,:)'...
            '; -data.' char(tag) '.VN100.acc_filt(3,:)'...
            '], translation_matrix(:,' num2str(ind) ') , 0);']);
       %Acceleration data:
        eval(['data.' char(tag) '.VN100.ang_trans=trans([data.' char(tag) '.VN100.ang_filt(1,:)'...
            '; -data.' char(tag) '.VN100.ang_filt(2,:)'...
            '; -data.' char(tag) '.VN100.ang_filt(3,:)'...
            '], translation_matrix(:,' num2str(ind) ') , 0);']);

        ind=ind+1;
        %---------
        %
        %calculate the heave (at sensor location):
        eval(['[data.' char(tag) '.VN100.ship.pitch, data.' char(tag) ...
           '.VN100.ship.roll, data.' char(tag) ...
           '.VN100.ship.surge, data.' char(tag) ...
           '.VN100.ship.sway, data.' char(tag) ...
           '.VN100.ship.heave] = heave( data.' char(tag) ...
           '.VN100.acc_trans, data.' char(tag) ...
           '.VN100.ang_trans ,length(data.' char(tag) '.VN100.moxatime)'...
           '/((data.bow.VN100.moxatime(end)-'...
           'data.bow.VN100.moxatime(1))*86400)); ']);
       %---------
       %
       %Interpolate the heave, roll and pitch:
       eval(['data.' char(tag) '.VN100.ship.heave=interp1(data.' char(tag)...
           '.VN100.moxatime, data.' char(tag) '.VN100.ship.heave, data.'...
           char(tag) '.symeo.moxatime,''pchip'' , ''extrap'')'';']);
       eval(['data.' char(tag) '.VN100.ship.roll=interp1(data.' char(tag)...
           '.VN100.moxatime, data.' char(tag) '.VN100.ship.roll, data.'...
           char(tag) '.symeo.moxatime,''pchip'' , ''extrap'')'';']);
       eval(['data.' char(tag) '.VN100.ship.pitch=interp1(data.' char(tag)...
           '.VN100.moxatime, data.' char(tag) '.VN100.ship.pitch, data.'...
           char(tag) '.symeo.moxatime,''pchip'' , ''extrap'')'';']);
       %---------
       %
       %Calculate the Dx and Dy distance offsets at water level (ship F.O.R.)
       eval(['[data.' char(tag) '.VN100.ship.dx, data.' char(tag) '.VN100.'...
           'ship.dy]=DxDy(data.' char(tag) '.symeo.distance_filt,data.'...
           char(tag) '.VN100.ship.roll.*(pi/180), data.' char(tag) ...
           '.VN100.ship.pitch.*(pi/180));']);
       %---------
       %
       %translate symeo to negative distance, add the heave value
       eval(['data.' char(tag) '.symeo.distance_motCor=((-data.' char(tag)...
           '.symeo.distance_filt''.*(cos(data.' char(tag) '.VN100.ship.pitch'...
           '.*(pi/180)).*cos(data.' char(tag) '.VN100.ship.roll.*(pi/180))))'...
           ' + data.' char(tag) '.VN100.ship.heave)'';']);
       %---------
       %apply bandpass filter on the distance, dx & dy data:
       eval(['FS=length(data.' char(tag) '.symeo.moxatime)/((data.' char(tag) '.symeo.'...
           'moxatime(end)-data.' char(tag) '.symeo.moxatime(1))*86400);']); %calculate the sampling Fq
       [wavB,wavA]=butter(fdegree,[Flow Fhigh]/(FS/2));
       %filt distance:
       eval(['data.' char(tag) '.symeo.distance_motCorF=filtfilt(wavB,wavA,data.'...
           char(tag) '.symeo.distance_motCor) + mean(data.' char(tag) ...
           '.symeo.distance_motCor);']);
       %filt dx:
       eval(['data.' char(tag) '.VN100.ship.dxF=filtfilt(wavB,wavA,data.'...
           char(tag) '.VN100.ship.dx) + mean(data.' char(tag) ...
           '.VN100.ship.dx);']);
       %filt dy:
       eval(['data.' char(tag) '.VN100.ship.dyF=filtfilt(wavB,wavA,data.'...
           char(tag) '.VN100.ship.dy) + mean(data.' char(tag) ...
           '.VN100.ship.dy);']);
       
       %ship course/heading from differentiating GPS locations (We need to measure
       %this for future, as we need ship heading in the event of ship drift
       %due to current, and if the ship is stationary).
      eval(['data.' char(tag) '.VN100.ship.heading=ang360b(180/pi*atan2('...
           'interp1(data.gps.gprmc.moxaDTE,sin(data.gps.gprmc.CourseMadeGood_true/180*pi)'...
           ',data.' char(tag) '.symeo.moxatime,''pchip'',''extrap''),interp1(data'...
           '.gps.gprmc.moxaDTE,cos(data.gps.gprmc.CourseMadeGood_true/180*pi),'...
           'data.' char(tag) '.symeo.moxatime,''pchip'',''extrap'')))'';']);
   end
   
   save([pwd matDir(2:end) fname(i).name],'-struct','data');
end


%==========================================================================
%Interpolate to achieve Common symeo timeseries:
%==========================================================================
shrt={'port' 'bow' 'stbd'};
[~,i]=min([length(data.port.symeo.distance_motCorF); ...
    length(data.bow.symeo.distance_motCorF);
    length(data.stbd.symeo.distance_motCorF)]);
shrt=(shrt(i));  %select the shortest timestamp

matDir=['.' filesep 'Data' filesep 'mat' filesep];
passFiltNum=4;
runMeanWindowSize=27;  %setting it to 27, so that we have a window size of 1sec worth of data
plotting=1;
%Generating common timestamp series
fname=dir([pwd matDir(2:end) '*.mat']);
for i=1:length(fname)
   data=load([pwd matDir(2:end) fname(i).name]); 
   %Report on the quality of the raw data  (make a function for doing this)
   %append the report (string output) to the data struct.
   %data.rawDataQuality.portSymeo.
   
   %remove the zero & nan readings and replace them with interpreted results:
   for tag={'port','bow','stbd'}
       if strcmp(char(tag),char(shrt))~=1
           %interpolation prepassing: Initial interp of zero readings and Nan's values:
           %^^ use 'extrap' option incase there are Nan's at start or end of data.
           %SYMEO:
           eval(['data.' char(tag) '.symeo.distance_motCorF=interp1(data.' char(tag)...
               '.symeo.moxatime, data.' char(tag) '.symeo.distance_motCorF, data.'...
               char(shrt) '.symeo.moxatime,''pchip'' , ''extrap'');']);
           eval(['data.' char(tag) '.symeo.distance_filt=interp1(data.' char(tag)...
               '.symeo.moxatime, data.' char(tag) '.symeo.distance_filt, data.'...
               char(shrt) '.symeo.moxatime,''pchip'' , ''extrap'');']);
           eval(['data.' char(tag) '.VN100.ship.dxF=interp1(data.' char(tag)...
               '.symeo.moxatime, data.' char(tag) '.VN100.ship.dxF, data.'...
               char(shrt) '.symeo.moxatime,''pchip'' , ''extrap'')'';']);
           eval(['data.' char(tag) '.VN100.ship.dyF=interp1(data.' char(tag)...
               '.symeo.moxatime, data.' char(tag) '.VN100.ship.dyF, data.'...
               char(shrt) '.symeo.moxatime,''pchip'' , ''extrap'')'';']);
           eval(['data.' char(tag) '.symeo.moxatime= data.'...
               char(shrt) '.symeo.moxatime;']);
       end
   end
       save([pwd matDir(2:end) fname(i).name],'-struct','data');
end
%==========================================================================





%==========================================================================
% 1-D spectral analysis:
%==========================================================================
matDir=['.' filesep 'Data' filesep 'mat' filesep];
%FFT settings:
nfft=2048;


fname=dir([pwd matDir(2:end) '*.mat']);
for i=1:length(fname)
    data=load([pwd matDir(2:end) fname(i).name]);  %load data file
    for tag={'port' 'bow' 'stbd'}
        eval(['FS=length(data.' char(tag) '.symeo.moxatime)/((data.' char(tag) '.symeo.'...
            'moxatime(end)-data.' char(tag) '.symeo.moxatime(1))*86400);']); %calculate the sampling Fq
        %----------------------------------------
        %Compute the FFT and spectrum parameters:
        %----------------------------------------
        %compute 1-D spectrum
        eval(['[data.wav1D.' char(tag) '.Spec, data.wav1D.' char(tag)...
            '.fvec]=pwelch(detrend(data.' char(tag) '.symeo.distance_motCorF,0)'...
            ',hanning(nfft),nfft/2,nfft,FS);']);
        %-----------
        %moments
        %-----------
        %compute zero moment:
        eval(['m0 = (trapz(data.wav1D.' char(tag)...
            '.fvec,data.wav1D.' char(tag) '.Spec));']);
        %compute 1st moment:
        eval(['m1 = (trapz(data.wav1D.' char(tag)...
            '.fvec,(data.wav1D.' char(tag) '.fvec).*data.wav1D.' char(tag) '.Spec));']);
        %compute 2nd moment:
        eval(['m2 = (trapz(data.wav1D.' char(tag)...
            '.fvec,(data.wav1D.' char(tag) '.fvec.^2).*data.wav1D.' char(tag) '.Spec));']);
        %comupte 4th moment:
        eval(['m4 = (trapz(data.wav1D.' char(tag)...
            '.fvec,(data.wav1D.' char(tag) '.fvec.^4).*data.wav1D.' char(tag) '.Spec));']);
        %comupte -1st (negative) moment:
        eval(['mN1 = (trapz(data.wav1D.' char(tag)...
            '.fvec,(data.wav1D.' char(tag) '.fvec.^(-1)).*data.wav1D.' char(tag) '.Spec));']);
        %comupte -1st (negative) moment:
        eval(['mN2 = (trapz(data.wav1D.' char(tag)...
            '.fvec,(data.wav1D.' char(tag) '.fvec.^(-2)).*data.wav1D.' char(tag) '.Spec));']);
        %-----------
        %Parameters:
        %-----------
        %integral sig wave height: (significant wave height from 0th moment)
        eval(['data.wav1D.' char(tag) '.Hm0 = 4*sqrt(trapz(data.wav1D.' char(tag)...
            '.fvec,data.wav1D.' char(tag) '.Spec));']);
        %Sig wave height: (significant wave height from 4*sqrt(var(elevation)))
        eval(['data.wav1D.' char(tag) '.Hsig = 4*sqrt(nanstd(data.' char(tag)...
            '.symeo.distance_motCorF));']);
        %Period Tm1 =m0/m1
        eval(['data.wav1D.' char(tag) '.Tm1 = m0/m1;']);
        %Period Tm2 =sqrt(m0/m2)
        eval(['data.wav1D.' char(tag) '.Tm2 = sqrt(m0/m2);']);
        %Peak Period Tp
        eval(['data.wav1D.' char(tag) '.Tp = 1/data.wav1D.' char(tag) '.fvec('...
            'data.wav1D.' char(tag) '.Spec==max(data.wav1D.' char(tag) '.Spec));']);
        %Peak Period Tpc
        eval(['data.wav1D.' char(tag) '.Tpc = mN2 * m1 / (m0^2);']);
        %-----------
        % calulate Tz0 and periods based on 1st and 2nd moment
        % for simplicity we do this for all members of wav_spec,
        % might be interesting to see how the motion correction influences these parameters
    end
    save([pwd matDir(2:end) fname(i).name],'-struct','data');
end
%==========================================================================











%==========================================================================
% 2-D spectral analysis:
%==========================================================================
% Set parameters required
% X and Y sensor locations in ship RHS coordinates [columns = x y] [rows = port; bow; stbd]
SensorLocationMatrix=[-0.77, 3.22; 0 ,0; -0.77, -3.22]; %FSG
%SensorLocationMatrix=[-1.87 1.095; 0,0 ; -1.87,-1.095]; %Celtic Voyager!!
CartLocMat=[-SensorLocationMatrix(:,2), ... %translate to cartesian coordinates
    SensorLocationMatrix(:,1)]; % 1st, 2nd & 3rd columns: Port, Bow & Stbd
% ns=4;       %
% n=4096;     % sample window
% lf=.0625;   % Lower Fq limit
% hf=1;       % Upper Fq limit
% nv=4;       % Number of voices (divisions per octave)
% %change to log-based freq limits:
% lp=log(lf)/log(2);lp=floor(lp);
% hp=log(hf)/log(2);hp=ceil(hp);



%scan for files:
matDir=['.' filesep 'Data' filesep 'mat' filesep];
fname=dir([pwd matDir(2:end) '*.mat']);
%
for i=1:length(fname)
    %
    data=load([pwd matDir(2:end) fname(i).name]);  %load data file
    %----------------------------------------------------------------------
    %Calculate the polar coordinates of the sensor locations:
    [data.wdm.theta,data.wdm.radius,data.wdm.origin]=polarCoords([...
        -data.port.VN100.ship.dyF+CartLocMat(1,1), ...
        data.port.VN100.ship.dxF+CartLocMat(1,2)],...
        [-data.bow.VN100.ship.dyF+CartLocMat(2,1),...
        data.bow.VN100.ship.dxF+CartLocMat(2,2)],...   %1st column is Port (x,y)
        [-data.stbd.VN100.ship.dyF+CartLocMat(3,1), ...%2nd  column is Bow (x,y)
        data.stbd.VN100.ship.dxF+CartLocMat(3,2)]);    %3rd column is Stbd (x,y)
    %Calculate the pair data:
    [data.wdm.pair] = RelativePolarCoords(data.wdm.theta, data.wdm.radius); 
    % Indexing for calculating pair-specific  values;
    % 1st, 2nd & 3rd columns: [ port - bow  |  port - stbd  |  bow  - stbd ]
    %-----------------------------------------------------------------------
    data.wdm.elevation=[data.port.symeo.distance_motCorF' data.bow.symeo.distance_motCorF'...
        data.stbd.symeo.distance_motCorF'];
    data.wdm.options.nLen=4096;
    data.wdm.options.lp=-5; %in octaves , equivalent to 2^(lp) Hz (lp=log(0.0625)/log(2);lp=floor(lp))
    data.wdm.options.hp=-1; %in octaves, equivalent to 2^(hp) Hz   (hp=log(1)/log(2);hp=ceil(hp))
    data.wdm.options.ns=27.8; %Symeo sampling freq
    data.wdm.options.nv=4; %number of voices per octave
    data.wdm.pair=data.wdm.pair;%pair; %run RelativePolarCoords.m from Proc_data.m
    %---------------
    % Wavelet analysis
    %=================
    [data.wdm.Amplitude,data.wdm.Wavenumber,data.wdm.Direction,...
        data.wdm.DirectionUnCorrected,data.wdm.Frequency]=wdmBS(...
        data.wdm.elevation,data.wdm.pair,data.wdm.options);
    %================
    % Normalize spatial and temporal coWavelet to Wavenumbers and Frequency
    wdm2wavnum(data.wdm.Wavenumber,data.wdm.Direction,data.wdm.Frequency,...
        data.wdm.Amplitude,mean(var(data.wdm.elevation)),data.wdm.options);
    %Save data:
    save([pwd matDir(2:end) fname(i).name],'-struct','data');
end
    
%==========================================================================




