function [gps] =gpsToMat(filepath)
%script that reads Moxa-logged GPS NMEA 0183-syntax sentences and 
%translates the data to matlab-readable format

%filepath='F:\FuelWave\Data\raw\171112_120000.00.GPS8';
%filepath='F:\FuelWave\Data\raw\171128_170000.00.GPS8';
start=now;

fid=fopen(filepath,'rt');
GPStext=textscan(fid,'%s',-1,'Delimiter','\n'); % scan each line as string
fclose(fid);
GPStext=char(GPStext{1});


gprmcIND=1;
%read through each sentence:
for i=1:size(GPStext,1)
    %split up the sentence into the relevant segments:
    sentence= strsplit(GPStext(i,:),{'\t',','});
    %distinguish between each $NMEA syntax and extract relevant data:   
   
    %-----------
    % &GPRMC
    %-----------
   if sentence{2}=='$GPRMC'
        %extract MOXA system time
        gps.gprmc.moxaDTE(gprmcIND,1)=datenum(['20' sentence{1}],'yyyymmddHHMMSS.FFF');
        %extract the GPS time (HHMMSS)
        gpsDate=sentence{3};
        gpsDay=sentence{11};
        %add '20' to year (convert yy to yyyy)
        gpsDay(7:8)=gpsDay(5:6);
        gpsDay(5:6)='20';
        gps.gprmc.gpsDTE(gprmcIND,1)=datenum([gpsDay gpsDate],'ddmmyyyyHHMMSS');
        %data validity
        if sentence{4}=='A'
            gps.gprmc.dataValidity(gprmcIND,1)='T';
        elseif sentence{4}=='V'
            gps.gprmc.dataValidity(gprmcIND,1)='F';
        end
        %Lat and Lon:
        lat=sentence{5};
        lon=sentence{7};
        gps.gprmc.Latitude(gprmcIND,1)=(str2double(lat(1:2)) + str2double(lat(3:end))/60);
        gps.gprmc.Longitude(gprmcIND,1)=(str2double(lon(1:3)) + str2double(lon(4:end))/60);
        if sentence{6}=='S'
            gps.gprmc.Latitude(gprmcIND,1)=gps.gprmc.Latitude(gprmcIND,1)*-1;
        end
        if sentence{8}=='W'
            gps.gprmc.Longitude(gprmcIND,1)=gps.gprmc.Longitude(gprmcIND,1)*-1;
        end
        %speeds
        gps.gprmc.speedOverGround_knots(gprmcIND,1)=str2double(sentence{9});
        gps.gprmc.CourseMadeGood_true(gprmcIND,1)=str2double(sentence{10});
        %magnetic variation
        gps.gprmc.MagneticVariation(gprmcIND,1)=str2double(sentence{12});
        checksum=sentence{13};
        %checksum (purge empty spaces suceeding sentence).
        checksum(checksum==' ')='';
        if checksum(1)=='W'
            gps.gprmc.MagneticVariation(gprmcIND,1)=gps.gprmc.MagneticVariation(gprmcIND,1)*-1;
        end
        gps.gprmc.checksum(gprmcIND,1)=str2double(checksum(end-1:end));
        %
        %Increase the counter index:
        gprmcIND=gprmcIND+1;
   end
    
    
    %-----------
    % &add additional NMEA syntaxes after this!
    %-----------
    
end



disp(['It took ' num2str((now-start)*1440) ' minutes to process GPS file']); 








