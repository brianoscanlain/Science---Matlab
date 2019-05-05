%TRANS.M
%
%subroutine OUT = trans(IN,ANGLES,IFLAG)
%
% Function from EDDYCORR toolbox
%
% June 27, 1996
%
% This subroutine rotates a vector from one cartesian basis to
% another, based upon the three Euler angles, defined as
% rotations around the reference axes, xyz. The axes in the
% rotated frame are x'y'z'.
%
%  IN     - 3XN matrix of input vector components
%
%  ANGLES - 3XN matrix of the angles phi, theta, psi, where:
%
%  phi   = angles(1,:) - rotation of x'y'z' about x axis (roll)
%  theta = angles(2,:) - rotation of x'y'z' about y axis (pitch)
%  psi   = angles(3,:) - rotation of x'y'z' about z axis (yaw)
%
%  OUT    - output vector components
%
% The transformation is given by
%
%    u' = [A]u      or     u = [transpose(A)]u'
%
% The argument IFLAG indicates in which direction the transform-
% ation is applied.
%
%   IFLAG = 0  "in" vector is transformed from x'y'z' to xyz
%   IFLAG = 1  "in" vector is transformed from xyz to x'y'z'
%
% For rotation of measurements from body coordinates to earth
% coordinates IFLAG is zero.
%
% This rotation is a 321 rotation. This means that the first rotation
% is around the 3 axis (z-axis, angle psi), the second rotation is then
% about the intermediate 2 axis (y-axis, theta), and the third rotation
% is about the intermediate 1 axis (x-axis, phi). 
%
% See Also: motcorr,  concat

function OUT = trans(IN,ANGLES,IFLAG)

  if nargin==2
   IFLAG=0;
  end

  p  = ANGLES(1,:);
  t  = ANGLES(2,:);
  ps = ANGLES(3,:);

  up = IN(1,:);
  vp = IN(2,:);
  wp = IN(3,:);

    if IFLAG          	% =1, from xyz to x'y'z'

      u = up.*cos(t).*cos(ps)                           + vp.*cos(t).*sin(ps)                           - wp.*sin(t);
      v = up.*(sin(p).*sin(t).*cos(ps)-cos(p).*sin(ps)) + vp.*(sin(p).*sin(t).*sin(ps)+cos(p).*cos(ps)) + wp.*(cos(t).*sin(p));
      w = up.*(cos(p).*sin(t).*cos(ps)+sin(p).*sin(ps)) + vp.*(cos(p).*sin(t).*sin(ps)-sin(p).*cos(ps)) + wp.*(cos(t).*cos(p));

    else		% =0, from x'y'z' to xyz

      u = up.*cos(t).*cos(ps) + vp.*(sin(p).*sin(t).*cos(ps)-cos(p).*sin(ps)) + wp.*(cos(p).*sin(t).*cos(ps)+sin(p).*sin(ps));
      v = up.*cos(t).*sin(ps) + vp.*(sin(p).*sin(t).*sin(ps)+cos(p).*cos(ps)) + wp.*(cos(p).*sin(t).*sin(ps)-sin(p).*cos(ps));
      w = up.*(-sin(t))       + vp.*(cos(t).*sin(p))                          + wp.*(cos(t).*cos(p));

    end;

  OUT = [u;v;w];


