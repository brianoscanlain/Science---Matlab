% rad2deg.m
% 
% DEGREES (RADIANS)
%
% This function takes radians as input and converts to equivalent degrees.

% Copyright 01-30-2006 Brendan C. Wood -- Synchroverge

function degrees = rad2deg(radians)

degrees = radians * 360 / ( 2 * pi );