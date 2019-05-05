function [pair] = RelativePolarCoords(theta, radius)
%RELATIVEPOLARCOORDS Calculates the relative distances and polar angles
%between all possible pairs, for a given set of polar locations.
%
% Input:
% ------
% theta   = (:,N) matrix of polar argument values for N given locations
% radius  = (:,N) matrix of polar radii for N given locations
%
% Output:
% -------
% pair.alpha_JK =  (:,N*(N-1)/2) matrix of relative angles between Jth-Kth pairs, 
%                               for all possible pairs. It represents the
%                               angle translation from J-th pair to Kth
%                               line [radians]
% pair.alpha_J  =  (:,N*(N-1)/2) matrix of angle for J-th pair [metres]
% pair.alpha_K  =  (:,N*(N-1)/2) matrix of angle for K-th pair [metres]
% pair.r_J      =  (:,N*(N-1)/2) matrix of length of J-th pair [metres]
% pair.r_K      =  (:,N*(N-1)/2) matrix of length of K-th pair [metres]
%
%
%--------------------------------------------------------------------------
%Brian Scanlon, 13 December 2017, NUIG
%-------------------------------------

%Debugging:
%theta=data.wdm.theta; radius= data.wdm.radius;

%Convert to Cartesian coordinates:
X = radius.*cos(theta);
Y = radius.*sin(theta);
%Calculate the distances of the joining lines:
l=0;x=[];y=[];
for j=1:size(X,2)-1
for k=(j+1):size(X,2)
   l=l+1;
   x(:,l)=X(:,k)-X(:,j); %x distance from j to k point
   y(:,l)=Y(:,k)-Y(:,j); %y distance from j to k point
end
end

i=sqrt(-1);
r=abs(x+i.*y); %distance between pairs!
a=atan2(y,x); %argument of the joining connecting the pair
%a(a<0)=a(a<0)+pi;
%find the 
l=0;pair.Ang_pair=[];
for j=1:size(X,2)-1
for k=(j+1):size(X,2)
   l=l+1;
   pair.alpha_JK(:,l)=a(:,k)-a(:,j); %angle going from J-th to K-th pair
   pair.alpha_J(:,l)=a(:,j);
   pair.alpha_K(:,l)=a(:,k);
   pair.r_J(:,l)=r(:,j); %length of Jth pair
   pair.r_K(:,l)=r(:,k); %length of Kth pair
end
end
pair.Ang_pair=pair.Ang_pair*180/pi;


% ii=find(pair.Ang_pair<0);
% pair.Ang_pair(ii)=pair.Ang_pair(ii)+360;
% ii=find(pair.Ang_pair>70 & pair.Ang_pair<110); %(find perpendicular~ish pairs)
% jj=find(pair.Ang_pair>250 & pair.Ang_pair<290);
% pair.isPerp=sort([ii jj]);


end

