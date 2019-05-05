function [GPS] = gpsFix(GPS,PivotDate)
%gpsFix is a custom function written specifically for one purpose; to
%convert the latitude and longitude fields from the inaccurate numeric
%values transcribed from the NMEA sentences.
%
% If PivotDate is specified, the script is used to correct systime.
% PivotDate = DDMMYY;
%
% Brian Scanlon, 27 July 2018

if nargin<2
    %Interpret the lat and long data:
    GPS.lat=floor(GPS.lat./100)+mod(GPS.lat,100)/60;
    GPS.lon=floor(GPS.long./100)+mod(GPS.long,100)/60;
    GPS=rmfield(GPS, 'long');
    
    %Extract the correct times and convert them to matlab timestamps:
    
    %GPS time:
    dd=floor(GPS.gpsdate./10000);
    mm=floor(mod(GPS.gpsdate,10000)/100);
    yy=mod(GPS.gpsdate,100)+2000;
    H=floor((GPS.gpstime(end))/10000);
    M=floor(mod(GPS.gpstime,10000)/100);
    S=mod(GPS.gpstime,100);
    GPS.gpsdte=datenum(yy,mm,dd,H,M,S);
    GPS=rmfield(GPS,{'gpstime','gpsdate'});
    %Sys time:
    H=floor(GPS.systime/10000);
    M=floor(mod(GPS.systime,10000)/100);
    S=mod(GPS.systime,100);
    GPS.dte=datenum(yy,mm,dd,H,M,S);
    GPS=rmfield(GPS,{'systime'});
else  %Secondary utility, just correct systime (for use with IMU and RADAR data)
    dd=floor(PivotDate./10000);
    mm=floor(mod(PivotDate,10000)/100);
    yy=mod(PivotDate,100)+2000;
    H=floor(GPS.systime/10000);
    M=floor(mod(GPS.systime,10000)/100);
    S=mod(GPS.systime,100);
    GPS.dte=datenum(yy,mm,dd,H,M,S);
    GPS=rmfield(GPS,{'systime'});
    %test for +/-1day issues:
    if sum((abs(diff(GPS.dte))>8*mad(abs(diff(GPS.dte)))))>0
        badDTE=[((abs(diff(GPS.dte))>8*mad(abs(diff(GPS.dte))))) 1==0];
        badVals=[diff(GPS.dte) 0];
        corr=round(badVals(badDTE));
        GPS.dte(badDTE) = GPS.dte(badDTE) + corr;
    end
end

