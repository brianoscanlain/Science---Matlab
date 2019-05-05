function [theta,r,origin] = polarCoords(left,middle,right)
%POLARCOORDS Draw a circle through left, middle & right coordinates.
%calculate the distance vector between the points, it's slope, and then
%calculate the perpendicular line which passes through it's centre. Doing
%this for both lines, we solve where the two points meet. This intersection
%point will be the centre of the circle. Next, the left, middle and right
%coordinates will be translated into polar coordinates from this centre
%point.
%----------------------------------
%    We adopt Cartesian 2-d (x,y) coordinate system for sake of universality:
%
%   (y axis)
%    ^
%    |
%    |
%    |
%    O- - - - - > (x axis)
%
%   Input:
%   ------
%   leftCoord = the left-most coordinate [x y] 
%   middleCoord = middle coordinate [x y] 
%   rightCoord = right coordinate [x y] 
%
%   Output:
%   -------
%   theta = polar angle matrix [Nx3 array] [left middle right;...]
%   r = polar radius matrix [Nx1 array] [circle radius]
%   origin = origin of drawn circle matrix [Nx2 array] [x y]
%
%   Brian Scanlon, NUIG, 12 Dec 2017
%----------------------------------
% left=[-data.port.VN100.ship.dyF+CartLocMat(1,1), ...
%        data.port.VN100.ship.dxF+CartLocMat(1,2)];
% middle=[-data.bow.VN100.ship.dyF+CartLocMat(2,1),...
%        data.bow.VN100.ship.dxF+CartLocMat(2,2)];       
% right=[-data.stbd.VN100.ship.dyF+CartLocMat(3,1), ...%2nd  column is Bow (x,y)
%        data.stbd.VN100.ship.dxF+CartLocMat(3,2)];
% 
% left=[-data.port.VN100.ship.dy-3.22, data.port.VN100.ship.dx-.77];
% middle=[-data.bow.VN100.ship.dy+0,data.bow.VN100.ship.dx+0];
% right=[-data.stbd.VN100.ship.dy+3.22, data.stbd.VN100.ship.dx-.77];

%find the two midpoints of the joining lines:
midL=(left-middle)/2;
midL(:,2)=-midL(:,2);
midR=(right-middle)/2;
midR(:,2)=-midR(:,2);
%Find the slope perpendicular to the joining lines:
SlopeL=((middle(:,2)-left(:,2))./(middle(:,1)-left(:,1)));
SlopeR=((middle(:,2)-right(:,2))./(middle(:,1)-right(:,1)));
   
%Find the origin coordinates of the circle:
 origin(:,1) = (SlopeL.*SlopeR.*(right(:,2)-left(:,2))+...
     SlopeL.*(middle(:,1)+right(:,1))-...
     SlopeR.*(left(:,1)+middle(:,1)))./(2*(SlopeL-SlopeR));
    
origin(:,2)=-1./SlopeL.*(origin(:,1) - (left(:,1) + middle(:,1))./2 )...
    +(left(:,2) + middle(:,2))/2;

%Move the left, right and middle coordinates to the origin:
left=left-origin;
middle=middle-origin;
right=right-origin;

   
%Find the radius and theta angle:
[theta(:,1),r]=cart2pol(left(:,1),left(:,2));
[theta(:,2),~]=cart2pol(middle(:,1),middle(:,2));
[theta(:,3),~]=cart2pol(right(:,1),right(:,2));
end


% Plotting (for debugging
% R = mean(sqrt((origin(:,1)-left(:,1)).^2+(origin(:,2)-left(:,2)).^2));
% figure,
% plot(left(:,1),left(:,2),'.')
% hold on
% plot(middle(:,1),middle(:,2),'.')
% plot(right(:,1),right(:,2),'.')
% th = 0:pi/50:2*pi;
% xunit =mean(r) * cos(th);
% yunit = mean(r) * sin(th);
% h = plot(xunit, yunit);
% line([0 mean(left(:,1))],[0 mean(left(:,2))])
% line([0 mean(middle(:,1))],[0 mean(middle(:,2))])
% line([0 mean(right(:,1))],[0 mean(right(:,2))])
% line([mean(left(:,1)) mean(right(:,1))],[mean(left(:,2)) mean(right(:,2))])
% line([mean(left(:,1)) mean(middle(:,1))],[mean(left(:,2)) mean(middle(:,2))])
% line([mean(middle(:,1)) mean(right(:,1))],[mean(middle(:,2)) mean(right(:,2))])
% figure
% plot(midL(:,1),midL(:,2))
% hold on
% plot(midR(:,1),midR(:,2))
% plot(left(:,1),left(:,2),'.')
% plot(right(:,1),right(:,2),'.')
% plot(middle(:,1),middle(:,2),'.')
% end
% 
% figure
% plot(left(:,1),left(:,2),'.')
% hold on
% plot(right(:,1),right(:,2),'.')
% plot(-1.095,-1.87,'rs','markersize',30)
% plot(+1.095,-1.87,'rs','markersize',30)
% plot(+1.095,-1.87,'rs','markersize',30,'markerfacecolor',[1,0,0])
% plot(-1.095,-1.87,'rs','markersize',30,'markerfacecolor',[1,0,0])
% plot( middle(:,1),middle(:,2),'.')
% plot(0,0,'rs','markersize',30,'markerfacecolor',[1,0,0])