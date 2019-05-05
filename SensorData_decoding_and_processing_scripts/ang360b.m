function [a360] = ang360(ang)
% ang360 maps arbitrary angles [degrees] to [0, 360) degrees.

% 25.1.99, Oyvind.Breivik@nrsc.no.

a360 = rem(ang, 360) - 180*(sign(ang)-1);
