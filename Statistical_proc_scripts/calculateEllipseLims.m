function [xLim,yLim] = calculateEllipseLims(ElipseRx,ElipseRy,steps)
%This function calculates points along the Elipse circumference in the
%first quadrant, spaced evenly using angle

theta = linspace(0, 90, steps)' .* (pi / 180);
sinalpha = sin(theta);
cosalpha = cos(theta);

xLim = (ElipseRx * cosalpha);
yLim = ( ElipseRy * sinalpha);
end
%========
