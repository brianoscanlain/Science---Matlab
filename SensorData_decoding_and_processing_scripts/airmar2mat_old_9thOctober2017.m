function [airmar]  =  airmar2mat(filepath)
% TODO this decoding script is not FINISHED!!!!!! 
% TODO think about nicer variable naming ...
% TODO decoding can be quite slow ...
%
%airmar2mat reads an NMEA sentences of AIRMAR logging data into a MATLAB 
%structure array.Airmar2mat requires the filepath inclusive of filename as
%an input variable. The function output is a structure array:
%
%  airmar.time.moxa
%             .GPS
%
%        .ship.lat   (GPGGA)
%             .lon   (GPGGA)
%             .TrueShipCourse  (GPVTG)  (heading considering drift of ship etc.)
%             .MagShipCourse   (GPVTG)  (compas reading of heading)
%             .UShipKMH  (GPVTG)
%             .UShipKNOT (GPVTG)                   
%             .altitude  (GPGGA)                
%             .d_YAW_dt  (TIROT)        (temporal rate of change of heading)
%
%        .MET.pres     (WIMDA, bar)     (air pressure)
%            .airT     (WIMDA)          (air temperature)
%            .RWS      (WIMDA)          (Relative wind speed in Knots)
%            .RWSknot  (WIMWV)          (Relative wind speed, assumes ship is stationary)
%            .TWD                       (True wind direction)
%            .TWS                       (True wind speed)
%
%
%
 
 
%temporary allocations:
%======================
%sl:: filepath='F:\Dropbox (AirSea Laboratory)\ultrawave\raw\voyager_2017_01\raw_svpak\170107_160000.00.GPS7';
script_name = 'airmar2mat';
%sl:: temp_dir = './temp/'; mymkdir(temp_dir);
%sl:: tempAIRMmat = [temp_dir 'tempAIRM.mat'];
%~~~
 
%Open file:
%==========
fid=fopen([filepath],'rt');
%sl:: AIRMtext=textscan(fid,'%s',Inf,'Delimiter','\n'); % scan each line as string
%sl:: for some reason 'Inf' does not work anymor "error: textscan: REPEAT = inf is too large"
%sl:: use -1 instead which reads till the END
AIRMtext=textscan(fid,'%s',-1,'Delimiter','\n'); % scan each line as string
fclose(fid);
AIRMtext=AIRMtext{1};   %go one level into the cell
%~~~
 
%Parse timestamps and recorded data fields:
%==========================================
[moxa_date, AIRM_text] = strtok(AIRMtext);   %Extract the Moxa timestamps
[syntaxID,~] = strtok(AIRM_text,',');    %Extract the NMEA syntax stamps (i.e. $GPGGA)
moxa_date=char(moxa_date); %convertt to string
syntaxID=char(syntaxID);  syntaxID=syntaxID(:,2:7);%convert to string
 
moxa_date = moxa_date(:,1:16);  %trim any excess 
moxa_date(:,3:end+2)=moxa_date;   %convert the moxa string to conventional
moxa_date(:,1)='2';               %yyyymmddHHMMSS.FFF  format
moxa_date(:,2)='0';               %
%~~~
airmar_raw.moxa_datenum = datenum(moxa_date,'yyyymmddHHMMSS.FFF' );
 
%Find a list of $GPGGA entries:
airmar_raw.gpgga_index=strmatch('$GPGGA',syntaxID(:,1:6));
airmar_raw.wimda_index=strmatch('$WIMDA',syntaxID(:,1:6));
airmar_raw.wimwv_index=strmatch('$WIMWV',syntaxID(:,1:6));
airmar_raw.yxxdr_index=strmatch('$YXXDR',syntaxID(:,1:6));
airmar_raw.gpvtg_index=strmatch('$GPVTG',syntaxID(:,1:6));
airmar_raw.gpzda_index=strmatch('$GPZDA',syntaxID(:,1:6));
airmar_raw.tirot_index=strmatch('$TIROT',syntaxID(:,1:6));
 
 
 
%run the scripts
% airmar.gpgga=gpggaread(AIRM_text(gpgga_index));
% alternatively decode GPGGA singal with the strok based GPXread.m
[airmar_raw.gpgga, airmar_raw.gpgga_index] = GPXread(AIRMtext, '$GPGGA');
%
airmar_raw.wimda=wimdaread(AIRM_text(airmar_raw.wimda_index));
airmar_raw.wimwv=wimwvread(AIRM_text(airmar_raw.wimwv_index));
airmar_raw.tirot=tirotread(AIRM_text(airmar_raw.tirot_index));
%airmar.gpzda=gpzdaread(AIRM_text(gpzda_index));
[airmar_raw.gpzda, airmar_raw.tgpzda_index] = GPXread(AIRMtext, '$GPZDA');
airmar_raw.yxxdr=yxxdrread(AIRM_text(airmar_raw.yxxdr_index));
%
[airmar_raw.gpvtg, airmar_raw.gpvtg_index] = GPXread(AIRMtext, '$GPVTG');
 
 
%%%%%%%%%%%%%%%%%%
% -> use moxa time stamp to interpolate all data on GPS time stamp (e.g. from GPZDA)
% note some times we have long gaps of nan data, need to linearly interpolate the GPS time first.
% use "nearest" for interpolation as this will work also for course and wind direction
%
% interpolate/copy data to airmar structure for output
airmar = airmar_raw.gpzda; % start with GPZDA dte(moxa time), gdte(gps time), and lz(time zone)
% 
airmar.gdte_WAS_NAN = isnan(airmar.gdte); % note nans in GPS data!
airmar.gdte = inpaint_nans(airmar.gdte); % fix gaps in the gps time stamp % TODO MAY NEED TO BE MORE SOFISTICATED!!!!
 
% GPGGA: dte, tod, fixQ, lat, lon, alt, Nsat, age_dgps, hor_pos_dilution
for varbl = {'lat', 'lon', 'alt', 'fixQ', 'Nsat'}
airmar.(char(varbl)) = interp1(airmar_raw.moxa_datenum(airmar_raw.gpgga_index), airmar_raw.gpgga.(char(varbl)), airmar.dte,'nearest'); % ship speed over ground from GPS
end
 
% GPVTG ship speed and course over ground (TrueShipCourse)
airmar.sog_knots = interp1(airmar_raw.moxa_datenum(airmar_raw.gpvtg_index), airmar_raw.gpvtg.sog, airmar.dte,'nearest'); % ship speed over ground from GPS
airmar.cog = interp1(airmar_raw.moxa_datenum(airmar_raw.gpvtg_index), airmar_raw.gpvtg.ttmg, airmar.dte,'nearest'); %TrueShipCourse from GPS
 
 
% WIMWV wind speed and direction
airmar.RWS_knots = interp1(airmar_raw.moxa_datenum(airmar_raw.wimwv_index), airmar_raw.wimwv.RWS_knots, airmar.dte,'nearest');
airmar.RWD = interp1(airmar_raw.moxa_datenum(airmar_raw.wimwv_index), airmar_raw.wimwv.RWD, airmar.dte,'nearest');
 
airmar.AirTemperatureCelsius = interp1(airmar_raw.moxa_datenum(airmar_raw.wimda_index), airmar_raw.wimda.AirTemperature, airmar.dte,'nearest');
airmar.AirPressureBAR = interp1(airmar_raw.moxa_datenum(airmar_raw.wimda_index), airmar_raw.wimda.pressueBar, airmar.dte,'nearest');
 
 
% fix NaNs at first or last entry (from interpolation)
for varbl = {'RWS_knots', 'RWD', 'AirTemperatureCelsius', 'AirPressureBAR'}
if isnan(airmar.(char(varbl))(1))
airmar.(char(varbl))(1)=airmar.(char(varbl))(2);
end
if isnan(airmar.(char(varbl))(end))
airmar.(char(varbl))(end)=airmar.(char(varbl))(end-1);
end
 
end
 
%next step is to scan the data for Nan entries and perhaps run some
%despiking routines.
 
 
%2.Check to see if the data are all the same lengths. In this example, the
%GPZDA data is one sentence shorter than the rest. 
 
%3. Test compatibility with matlab
