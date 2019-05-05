function [x,y]=latlon2merc(longitude,latitude,mapWidth,mapHeight)
%This function takes in spherical coordinates (latitude and longitude) and
%translates them into pixel value coordinates for a Mercator projection.
%
%brian Scanlon, 2017, NUI Galway
x = (longitude+180)*(mapWidth/360);

mercN = log(tan((pi/4)+( (-latitude*pi/180)/2) ));
y     = (mapHeight/2)-(mapHeight*mercN/(2*pi));

y=mapHeight-y;     %just to put origin in top right-hand corner
y(y<1)=1;
y(y>mapHeight)=mapHeight;

