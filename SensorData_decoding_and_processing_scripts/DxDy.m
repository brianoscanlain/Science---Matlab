function [dx,dy] = DxDy(vlength,roll,pitch)
%DXDY is a script to calculate the X and Y displacements of vector once
%it's roll and pitch values are known.
%   Input:
%   TangentData: vector length. 
%   roll: the roll of the vector (radians)
%   pitch: pitch of the vector   (radians)
%
%   Output:
%   dx: x-displacement from vector orgin to the end of the vector.
%   dy: y-displacement from vector origin to the end of the vector

%Negative roll required to calculate a positive dy:
dx=sin(pitch).*vlength';    
dy=-sin(roll).*vlength';
end

