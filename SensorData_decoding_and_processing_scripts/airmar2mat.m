function [airmar]  =  airmar2mat(filepath)
% Script updated:
% > Now runs on internal subfunctions and doesnt call GPXread().m 
% > Takes ~4 seconds to run 1-hour-long airmar file.
% > Names and structures may have changed slightly (see below description)
% > Interpolates GPS Nan entries
% > Extrapolates missing data (missing sentences at end of file) 
%
%airmar2mat reads an NMEA sentences of AIRMAR logging data into a MATLAB 
%structure array.Airmar2mat requires the filepath inclusive of filename as
%an input variable. The function output is a structure array:
%
%  airmar.time.moxa        
%             .gdte        (GPZDA)      (GPS UTC time)
%             .gdte_local  (GPZDA)      (GPS local time)
%
%        .ship.lat   (GPGGA)
%             .lon   (GPGGA)
%             .TrueShipCourse  (GPVTG)  (heading considering drift of ship etc.)
%             .MagShipCourse   (GPVTG)  (compas reading of heading)
%             .sog_kmh  (GPVTG)
%             .sog_knots (GPVTG)                   
%             .altitude  (GPGGA)                
%             .d_YAW_dt  (TIROT)        (temporal rate of change of heading)
%
%        .MET.pres     (WIMDA, bar)     (air pressure)
%            .airT     (WIMDA)          (air temperature)
%            .RWS      (WIMDA)          (Relative wind speed in Knots)
%            .RWSknot  (WIMWV)          (Relative wind speed, assumes ship is stationary)
%            .TWD                       (True wind direction)
%
% =========================================================================
tic
%======================
%temporary allocations:
%======================
script_name = 'airmar2mat';
%~~~
%==========
%Open file:
%==========
fid=fopen([filepath],'rt');
%AIRMtext=textscan(fid,'%s',-1,'Delimiter','\n'); % scan each line as string 
%The option -1 reads no data in octave: errMsg: "textscan: N = 0, no data read"
AIRMtext=textscan(fid,'%s','Delimiter','\n'); % scan each line as string
fclose(fid);
AIRMtext=AIRMtext{1};   %go one level into the cell
%~~~
%==========================================
%Parse timestamps and recorded data fields:
%==========================================
[moxa_date, AIRMtext] = strtok(AIRMtext);   %Extract the Moxa timestamps
[syntaxID,~] = strtok(AIRMtext,',');    %Extract the NMEA syntax stamps (i.e. $GPGGA)
moxa_date=char(moxa_date); %convertt to string
syntaxID=char(syntaxID);  syntaxID=syntaxID(:,2:7);%convert to string
moxa_date = moxa_date(:,1:16);  %trim any excess 
moxa_date(:,3:end+2)=moxa_date;   %convert the moxa string to conventional
moxa_date(:,1)='2';               %yyyymmddHHMMSS.FFF  format
moxa_date(:,2)='0';               %
moxa_dte = datenum(moxa_date,'yyyymmddHHMMSS.FFF' );
airmar.time.moxa=moxa_dte;
%~~~
disp('- finished importing file -');toc; 
%=================================================
% Identify the sentences and all the subfunctions:
%=================================================
%Find a list of $GPGGA entries:
airmar_raw.gpgga_index=strmatch('$GPGGA',syntaxID(:,1:6));
airmar_raw.wimda_index=strmatch('$WIMDA',syntaxID(:,1:6));
airmar_raw.wimwv_index=strmatch('$WIMWV',syntaxID(:,1:6));
airmar_raw.yxxdr_index=strmatch('$YXXDR',syntaxID(:,1:6));
airmar_raw.gpvtg_index=strmatch('$GPVTG',syntaxID(:,1:6));
airmar_raw.gpzda_index=strmatch('$GPZDA',syntaxID(:,1:6));
airmar_raw.tirot_index=strmatch('$TIROT',syntaxID(:,1:6));
%run the scripts
disp('- running sub functions to parse sentence info:: -');toc; 
airmar_raw.gpgga = gpggaread(AIRMtext(airmar_raw.gpgga_index));
airmar_raw.wimda = wimdaread(AIRMtext(airmar_raw.wimda_index));
airmar_raw.wimwv = wimwvread(AIRMtext(airmar_raw.wimwv_index));
airmar_raw.tirot = tirotread(AIRMtext(airmar_raw.tirot_index));
airmar_raw.gpzda = gpzdaread(AIRMtext(airmar_raw.gpzda_index));
airmar_raw.yxxdr = yxxdrread(AIRMtext(airmar_raw.yxxdr_index));
airmar_raw.gpvtg = gpvtgread(AIRMtext(airmar_raw.gpvtg_index));
%%%%%%%%%%%%%%%%%%
% -> use moxa time stamp to interpolate all data on GPS time stamp (e.g. from GPZDA)
% note some times we have long gaps of nan data, need to linearly interpolate the GPS time first.
% use "nearest" for interpolation as this will work also for course and wind direction

%==========================================================================
%             Deal with varying length data:
%==========================================================================
disp('- Interpolating data and tidying up -');toc; 
%Theoretical length
airmar.time.moxa=airmar.time.moxa(1:8:end); %airmar prints 8-sentence chunks
Len=length(airmar.time.moxa);
%Calculate lengths:
Lgpg = length(airmar_raw.gpgga_index);
Lwmd = length(airmar_raw.wimda_index);
Lwmv = length(airmar_raw.wimwv_index);
Ltir = length(airmar_raw.tirot_index);
Lgpz = length(airmar_raw.gpzda_index);   
Lyxx = length(airmar_raw.yxxdr_index)/2; %should be twice the length
Lgpv = length(airmar_raw.gpvtg_index);
% test lengths with Theoretical, if wrong, add number of Nan's and Extrapolate:
if Lgpz<Len || sum(isnan(airmar_raw.gpzda.gdte))>0
    %ignore Nan indexed values:
    ind=(1:Lgpz);
    airmar_raw.gpzda.gdte=interp1(ind(~isnan(airmar_raw.gpzda.gdte)),...
        airmar_raw.gpzda.gdte(~isnan(airmar_raw.gpzda.gdte)),(1:Len));
    airmar_raw.gpzda.gdte_local=interp1(ind(~isnan(airmar_raw.gpzda.gdte)),...
        airmar_raw.gpzda.gdte_local(~isnan(ind(~isnan(airmar_raw.gpzda.gdte_local)))),(1:Len));
end
if Lwmd<Len
    airmar_raw.wimda.AirTemperature=interp1((1:Lwmd),airmar_raw.wimda.AirTemperature ,(1:Len),'linear','extrap');
    airmar_raw.wimda.pressueBar=interp1((1:Lwmd),airmar_raw.wimda.pressueBar,(1:Len),'linear','extrap');
    airmar_raw.wimda.pressueInches=interp1((1:Lwmd),airmar_raw.wimda.pressueInches,(1:Len),'linear','extrap');
    airmar_raw.wimda.RWD=interp1((1:Lwmd),airmar_raw.wimda.RWD,(1:Len),'linear','extrap');
    airmar_raw.wimda.TWD=interp1((1:Lwmd),airmar_raw.wimda.TWD,(1:Len),'linear','extrap');
    airmar_raw.wimda.TWS=interp1((1:Lwmd),airmar_raw.wimda.TWS,(1:Len),'linear','extrap');
    airmar_raw.wimda.TWS_knots=interp1((1:Lwmd),airmar_raw.wimda.TWS_knots,(1:Len),'linear','extrap');
end
if Lwmv<Len
    airmar_raw.wimwv.RWD=interp1((1:Lwmv),airmar_raw.wimwv.RWD,(1:Len),'linear','extrap');
    airmar_raw.wimwv.RWS_knots=interp1((1:Lwmv),airmar_raw.wimwv.RWS_knots,(1:Len),'linear','extrap');
end
if Ltir<Len
    airmar_raw.tirot.ShipRateOfTurn=interp1((1:Ltir),airmar_raw.tirot.ShipRateOfTurn,(1:Len),'linear','extrap');
end
if Lgpg<Len
    airmar_raw.gpgga.altitude=interp1((1:Lgpg),airmar_raw.gpgga.altitude,(1:Len),'linear','extrap');
    airmar_raw.gpgga.latitude=interp1((1:Lgpg),airmar_raw.gpgga.latitude,(1:Len),'linear','extrap');
    airmar_raw.gpgga.longitude=interp1((1:Lgpg),airmar_raw.gpgga.longitude,(1:Len),'linear','extrap');
    airmar_raw.gpgga.GPSquality=interp1((1:Lgpg),airmar_raw.gpgga.GPSquality,(1:Len),'linear','extrap');
end
if Lyxx<Len
    airmar_raw.yxxdr.A.WindChillTemperature=interp1((1:Lyxx),...
        airmar_raw.yxxdr.A.WindChillTemperature,(1:Len),'linear','extrap');
end
if Lgpv<Len
    airmar_raw.gpvtg.RelativeShipHeading=interp1((1:Lgpv),airmar_raw.gpvtg.RelativeShipHeading,...
        (1:Len),'linear','extrap');
    airmar_raw.gpvtg.ShipSpeedOverGround_KMh=interp1((1:Lgpv),airmar_raw.gpvtg.ShipSpeedOverGround_KMh,...
        (1:Len),'linear','extrap');
    airmar_raw.gpvtg.ShipSpeedOverGround_KNOTS=interp1((1:Lgpv),airmar_raw.gpvtg.ShipSpeedOverGround_KNOTS,...
        (1:Len),'linear','extrap');
    airmar_raw.gpvtg.TrueShipHeading=interp1((1:Lgpv),airmar_raw.gpvtg.TrueShipHeading,...
        (1:Len),'linear','extrap');
end



%====================================================
%Place the relevant information in the output struct:
%====================================================
airmar.time.gdte = airmar_raw.gpzda.gdte;
airmar.time.gdte_local = airmar_raw.gpzda.gdte_local;
airmar.ship.lat = airmar_raw.gpgga.latitude;
airmar.ship.lon = airmar_raw.gpgga.longitude;
airmar.ship.altitude = airmar_raw.gpgga.altitude;
airmar.ship.TrueShipCourse = airmar_raw.gpvtg.TrueShipHeading;
airmar.ship.SOG_kmh = airmar_raw.gpvtg.ShipSpeedOverGround_KMh;
airmar.ship.sog_knots = airmar_raw.gpvtg.ShipSpeedOverGround_KNOTS;
airmar.ship.d_YAW_dt = airmar_raw.gpgga;
airmar.MET.pres = airmar_raw.wimda.pressueBar;
airmar.MET.airT = airmar_raw.wimda.AirTemperature;
airmar.MET.RWD = airmar_raw.wimda.RWD;
airmar.MET.TWD = airmar_raw.wimda.TWD;
airmar.MET.RWS = airmar_raw.wimda.TWS;
airmar.MET.RWS_knots = airmar_raw.wimda.TWS_knots;
end








%==========================================================================
%                 SUB FUNCTIONS
%==========================================================================



function [TIROT] = tirotread( sentences)
%tirotread reads in NMEA $TIROT-syntaxed sentences and returns a structured
%array of the results.
%
%=======
% $TIROT,<1>,<2>*hh<CR><LF>
%=======
%<1> Signed rate of turn, degrees per minute, to the nearest 0.1 degree. Negative 
%    values indicate the bow is turning to port. 
%<2> Status: A = Data Valid; V = Data Invalid. 
%==============================================
%
%
if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end
for i=1:size(sentences,1)
    %XX(i,:)=textscan(sentences(i,:),'%[$TIROT] %f %c \r\n','delimiter',',');
    XX(i,:)=textscan(sentences(i,:),'%s %f %s','delimiter',',');
end
TIROT.ShipRateOfTurn=([XX{:,2}]');
end









function [GPVTG] = gpvtgread( sentences)
%gpvtgread reads in NMEA $GPVG-syntaxed sentences and return a data
%structure
%
%=======
% $GPVTG,<1>,<2>,<3>,<4>,<5>,<6>,<7>,<8>,<9>*hh<CR><LF>
%=======
%<1> Course Over Ground, degrees True, to the nearest 0.1 degree 
% <2> T = True 
% <3> Course Over Ground, degrees Magnetic, to the nearest 0.1 degree 
% <4> M = Magnetic 
% <5> Speed Over Ground, knots, to the nearest 0.1 knot 
% <6> N = Knots 
% <7> Speed Over Ground, km/hr, to the nearest 0.1 km/hr 
% <8> K = km/hr 
% <9> Mode indicator: 
%      nmea_options =['$GPGGA'
%                     '$GPGLL'
%                     '$GPGSA'
%                     '$GPGSV'
%                     '$GPRMC'
%                     '$GPVTG'
%                     '$GPZDA'
%                     '$SDDBS'];
% (The only values transmitted by the WX Series WeatherStation Sensor for the Mode indicator are A, D, and N. )
%
%

if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end

for i=1:size(sentences,1)
    XX(i,:)=textscan(sentences(i,:),'%s %f %s %f %s %f %s %f %s %s','delimiter',',');
end
GPVTG.TrueShipHeading=([XX{:,2}]');
GPVTG.RelativeShipHeading=([XX{:,4}]');
GPVTG.ShipSpeedOverGround_KNOTS=([XX{:,6}]');
GPVTG.ShipSpeedOverGround_KMh=([XX{:,8}]');
end







function [GPZDA] = gpzdaread( sentences)
%gpvtgread reads in NMEA $GPVG-syntaxed sentences and return a data
%structure
%
%=======
% $GPZDA,<1>,<2>,<3>,<4>,<5>,<6>*hh<CR><LF>
%=======
% <1> UTC time of day, in the form hhmmss 
% <2> UTC day, 01 to 31 
% <3> UTC month, 01 to 12 
% <4> UTC year (four digits, e.g. 2006) 
% <5> Local time zone hours, 00 to +/-13 hrs 
% <6> Local time zone minutes, 00 to +59 
if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end

for i=1:size(sentences,1)
    %XX(i,:)=textscan(sentences(i,:),'%[$GPZDA] %s %f %f %f %f %f \r\n','delimiter',',');
    XX(i,:)=textscan(sentences(i,:),'%s %s %f %f %f %f %f','delimiter',',');
end
hhmmss=char([XX{:,2}]');
%UTC time:
GPZDA.gdte=datenum(([XX{:,5}]'),([XX{:,4}]'),([XX{:,3}]'),...
    str2num(hhmmss(:,1:2)),str2num(hhmmss(:,3:4)),...
    str2num(hhmmss(:,5:6)));
%Local time:
GPZDA.gdte_local=GPZDA.gdte+(([XX{:,6}]')./24)+(([XX{:,7}]')./1440);
end







function [ YXXDR] = yxxdrread(sentences )
%yxxdrread reads NMEA $YXXDR-syntaxed sentences and returns data structure
%array.
%=======
% Three versions of $YXXDR Syntax, type A, B, & C:
% Type A:
%========
% $YXXDR ,<1>, <2>, <3>, <4>,<5>, <6>, <7>, <8>,<9>, <10>,<11>,<12>,<13>,<14>,<15>,<16>*hh<CR><LF>
% <1> C=temperature
% <2> Calculated "relative" wind chill temperature, degrees Celsius, to the nearest 0.1 degree
% <3> C = degrees C
% <4> WCHR (ID indicating relative wind chill)
% <5> C = temperature
% <6> Calculated "theoretical" wind chill temperature, degrees Celsius, to the nearest 0.1 degree
% <7> C = degrees C
% <8> WCHT (ID indicating theoretical wind chill)
% <9> C = temperatureC = temperature
% <10> Calculated heat index, degrees Celsius, to the nearest 0.1 degree
% <11> C = degrees C
% <12> HINX (ID indicating heat index) 
% <13>  P = pressure 
% <14> Actual measured barometric pressure, or "station pressure", bars, to the nearest 0.001 bar 
% <15>  B = bars %  Set up a list of valid NMEA strings
%
%
%Type B:
%========
%$YXXDR,<1>, <2>, <3>, <4>,<5>, <6>, <7>, <8>*hh<CR><LF> 
%======
% <1>  A = angular displacement 
% <2> Pitch: oscillation of vessel about its latitudinal axis.  Bow moving up is positive.  Value reported to the nearest 0.1 degree. 
% <3> D = degrees 
% <4> PTCH (ID indicating pitch of vessel) 
% <5> A = angular displacement 
% <6> Roll: oscillation of vessel about its longitudinal axis.  Roll to the starboard is positive.  Value reported to the nearest 0.1 degree. 
% <7> D = degrees 
% <8> ROLL (ID indicating roll of vessel) 
%
%
if iscell(sentences)
    sentences=char(sentences);
end
%Seperate sentences into A or B:
for i=1:length(sentences);yxxdr_commas(i,1)=length(strfind(sentences(i,:),','));end
%now that we know the number of commas per sentence, we can distinguish
%between types A and B.
sentencesA=sentences(yxxdr_commas==16,:);%this unit of 16 may need to be relaxed (i.e. >=10 etc.)
sentencesB=sentences(yxxdr_commas==8,:);
clear sentences;

for i=1:size(sentencesA,1)
    %XXa(i,:)=textscan(sentencesA(i,:),'%[$YXXDR] %c %f %c %s %c %f %c %s %c %f %c %s %c %f %c \r\n','delimiter',',');
    XXa(i,:)=textscan(sentencesA(i,:),'%s %s %f %s %s %s %f %s %s %s %f %s %s %s %f %s','delimiter',',');
end
YXXDR.A.WindChillTemperature=([XXa{:,3}]');
end






function [WIMDA] = wimdaread( sentences)
%wimdaread reads in NMEA $WIMDA-syntaxed sentences and returns a structured
%array of the results.
%
%%=======
% $WIMDA,<1>,<2>,<3>,<4>,<5>,<6>,<7>,<8>,<9>,<10>,<11>,<12>,<13>,<14>,<15>,<16>,<17>,<18>,<19>,<20>*hh<CR><LF>
% ======
% <1> Barometric pressure, inches of mercury, to the nearest 0.01 inch
% <2> I = inches of mercury
% <3> Barometric pressure, bars, to the nearest .001 bar
% <4> B = bars
% <5> Air temperature, degrees C, to the nearest 0.1 degree C
% <6> C = degrees C
% <7> Water temperature, degrees C (this field left blank by WeatherStation instrument)
% <8> C = degrees C (this field left blank by WeatherStation instrument)
% <9> Relative humidity, percent, to the nearest 0.1 percent
% <10> Absolute humidity, percent (this field left blank by WeatherStation instrument)
% <11> Dew point, degrees C, to the nearest 0.1 degree C
% <12> C = degrees C
% <13> Wind direction, degrees True, to the nearest 0.1 degree
% <14> T = true
% <15> Wind direction, degrees Magnetic, to the nearest 0.1 degree
% <16> M = magnetic
% <17> Wind speed, knots, to the nearest 0.1 knot
% <18> N = knots
% <19> Wind speed, meters per second, to the nearest 0.1 m/s
% <20> M = meters per second

if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end

for i=1:size(sentences,1)
    %XX(i,:)=textscan(sentences(i,:),'%[$WIMDA] %f %c %f %c %f %c %f %c %f %f %f %c %f %c %f %c %f %c %f %s \r\n','delimiter',',');
    XX(i,:)=textscan(sentences(i,:),'%s %f %s %f %s %f %s %f %s %f %f %f %s %f %s %f %s %f %s %f %s','delimiter',',');
end
WIMDA.pressueInches=[XX{:,2}]';
WIMDA.pressueBar=[XX{:,4}]';
WIMDA.AirTemperature=[XX{:,6}]';
WIMDA.TWD=[XX{:,14}]';
WIMDA.RWD=[XX{:,16}]';
WIMDA.TWS_knots=[XX{:,18}]';
WIMDA.TWS=[XX{:,20}]';
end








function [WIMWV] = wimwvread( sentences)
%wimwvread reads in NMEA $WIMWV-syntaxed sentences and returns a structured
%array of the results.
%
%=======
% $WIMWV,<1>,<2>,<3>,<4>,<5>*hh<CR><LF>
% ======
% <1> Wind angle, 0.0 to 359.9 degrees, in relation to the vessel�s bow/centerline, to the nearest 0.1 degree. If the data for this field is not valid, the field will be blank.
% <2> Reference:
% R = Relative (apparent wind, as felt when standing on the moving ship)
% T = Theoretical (calculated actual wind, as though the vessel were stationary)
% <3> Wind speed, to the nearest tenth of a unit. If the data for this field is not valid, the field will be blank.
% <4> Wind speed units:
% K = km/hr
% M = m/s
% N = knots
% S = statute miles/hr
% In the WeatherStation instrument, this field always contains "N" (knots).
% <5> Status: A = data valid; V = data invalid
if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end

for i=1:size(sentences,1)
%XX(i,:)=textscan(sentences(i,:),'%[$WIMWV] %f %c %f %c %s \r\n','delimiter',',');
XX(i,:)=textscan(sentences(i,:),'%s %f %s %f %s %s','delimiter',',');
end
if char([XX{1,3}]')=='R'   %relative wind values
   WIMWV.RWD=([XX{:,2}]');
   
    if char([XX{1,5}]')=='N'  
        WIMWV.RWS_knots=([XX{:,4}]');
    elseif char([XX{1,5}]')=='M'
        WIMWV.RWS_m_per_sec=([XX{:,4}]');
    elseif char([XX{1,5}]')=='K'
        WIMWV.RWS_kmh=([XX{:,4}]');
    end

elseif char([XX{1,3}]')=='T'   %True wind values (corrected for ship motion)
   WIMWV.TWD=([XX{:,2}]');
   
    if char([XX{1,5}]')=='N'  
        WIMWV.TWS_knots=([XX{:,4}]');
    elseif char([XX{1,5}]')=='M'
        WIMWV.TWS_m_per_sec=([XX{:,4}]');
    elseif char([XX{1,5}]')=='K'
        WIMWV.TWS_kmh=([XX{:,4}]');
    end
end
%we ignore the remaining datafields which would only be populated with more
%expensive AIRMAR units.
end







function [GPGGA] = gpggaread(sentences)
%gpggaread will read in the gpgga sentences and output data structure.
%=======
%  $GPGGA,<1>,<2>,<3>,<4>,<5>,<6>,<7>,<8>,<9>,<10>,<11>,<12>,<13>,<14>*hh<CR><LF>
%  ======
%  <1> UTC of position, in the form hhmmss
%  <2> Latitude, to the nearest .0001 minute
%  <3> N if field <2> is North Latitude
%      S if field <2> is South Latitude
%  <4> Longitude, to the nearest .0001 minute
%  <5> E if field <4> is East Longitude
%      W if field <4> is West Longitude
%  <6> GPS quality indicator:
%                            0 = Fix not available or invalid
%                            1 = GPS SPS Mode, fix valid
%                            2 = Differential GPS, SPS Mode, fix valid
%                            3 = GPS PPS Mode, fix valid
%                            4 = Real Time Kinematic (RTK)
%                            5 = Float RTK
%                            6 = Estimated (dead reckoning) Mode
%                            7 = Manual Input Mode
%                            8 = Simulator Mode
%      When providing data from the WX Series WeatherStation Sensor�s
%      internal GPS, the only valid values for the GPS quality indicator are 0,
%      1, and 2.
%  <7> Number of satellites in use, 0-12
%  <8> Horizontal dilution of precision (HDOP)
%  <9> Altitude relative to mean-sea-level (geoid), meters (to the nearest
%      whole meter)
% <10> M
% <11> Geoidal separation, meters (to the nearest whole meter). In the WX
%      Series WeatherStation Sensor, this field contains the separation data, if
%      available, otherwise, it is not provided, and appears as a null field.
% <12> M. In the WX Series WeatherStation Sensor, this field contains M, if
%      separation data is available, otherwise, it is not provided, and appears
%      as a null field.
% <13> Age of Differential GPS data, seconds. This field is not provided by the
%      WX Series WeatherStation Sensor, and appears as a null field.
% <14> Differential reference station ID, 0000-1023. This field is not provided
%      by the WX Series WeatherStation Sensor, and appears as a null field.
%
%error: strread: %q, %c, %[] or bit width format specifiers are not supported yet.
%error: called from
%    strread at line 329 column 7
%    textscan at line 321 column 8
%    airmar2mat>gpggaread at line 512 column 12
%    airmar2mat at line 76 column 18
% sentences=AIRMtext(airmar_raw.gpgga_index);

if iscell(sentences)       %convert to characters/strings if presented as cells
    sentences=char(sentences);
end


for i=1:size(sentences,1)
    %XX(i,:)=textscan(sentences(i,:),'%[$GPGGA] %s %f %c %f %c %f %f %f %f %c %f %c %s \r\n','delimiter',','); %Octave cannot deal with %c, %[]
    XX(i,:)=textscan(sentences(i,:),'%s %s %f %s %f %s %f %f %f %f %s %f %s %s %s','delimiter',',');
end
% <1>    <day fraction>
DayfractionUTC=char([XX{:,2}]');
GPGGA.DayfractionUTC=str2num(DayfractionUTC(:,1:2))./24  +  ...
    str2num(DayfractionUTC(:,3:4))./3600   +   ...
    str2num(DayfractionUTC(:,5:9))./86400;
% <2>   <latitude>
lat=([XX{:,3}]');
GPGGA.latitude=dm2degrees([floor(lat(:)./100) lat-100.*(floor(lat(:)./100))]);
%GPGGA.latitude(char(XX(:,4))=='S')=-abs(GPGGA.latitude(char(XX(:,4))=='S'));
GPGGA.latitude(char(vertcat(XX{:,4}))=='S')=-abs(GPGGA.latitude(char(vertcat(XX{:,4}))=='S'));
% <4>   <longitude>
lon=[XX{:,5}]';
GPGGA.longitude=dm2degrees([floor(lon(:)./100) lon-100.*(floor(lon(:)./100))]);
%GPGGA.longitude(char(XX(:,6))=='W')=-abs(GPGGA.longitude(char(XX(:,6))=='W'));
GPGGA.longitude(char(vertcat(XX{:,6}))=='W')=-abs(GPGGA.longitude(char(vertcat(XX{:,6}))=='W'));
% <6>   <GPS quality factor>
GPGGA.GPSquality=[XX{:,7}]';
% <7>
% <8>
% <9>   altitude
GPGGA.altitude=[XX{:,10}]';
% <10>
% <11>
end


