function [XFilt,YFilt,kx,ky] =ElypticFilt(X,Y,error)
%ElypticFilt distinguishes bad data based on an elyptical X-Y
%statistical distribution limit, and by the error codes (if any). The code first
%derives the centre of the elipse, using median statistics (the mean works
%good when the data is well behaved, but for contaminated data it can
%derive poor results).
%
% Input:
%-------
%      X  (unfiltered data [Nx1 array])
%      Y  (unfiltered data [Nx1 array])
%      error (errorcode data, specific for the symeo radar application)
%
% Output:
%--------
%     XFilt - filtered data
%     YFilt - Filtered data
%     kx - relation between std and mad of X, proxy for Normal distribution.
%     ky - relation between std and mad of Y, proxy for Normal distribution.
%
% Brian Scanlon, NUIG, 13th Feb 2018

%Check to see that the VargIn have the correct structure:
[~, dimX]=max(size(X));
[~, dimY]=max(size(Y));
if dimX>1
    X=X';
end
if dimY>1
    Y=Y';
end
%--------------------    


%Tuning parameters
Zscore=4;   %default was 4, however I set it to 6 to relax the eliptical tolerance
MAD_STD_scale_factor=1.4826;
MADzZcore=Zscore*MAD_STD_scale_factor;

% %Debugging and testing artifacts:
Debug=0;
if Debug==1
    X=port.symeo.distance;
    Y=port.symeo.velocity;
    error=port.symeo.error;
    close all;
    plot(X,Y,'.')
    set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
    xlabel('Distance [m]','fontsize',16,'fontname','times')
    ylabel('Velocity [m/s]','fontsize',16,'fontname','times')
end
% %----

if Debug~=1
if nargin>2
    X(error>0)=0./0; %Set the error codes above zero to Nan values (Error code 0 is good data)
    Y(error>0)=0./0; %Set the error codes above zero to Nan values (Error code 0 is good data)
end
end

%First we need to find the centre of the "good" data. To do this, I first
%start off with the median. Next we purge outliers, and keep doing this 
%until we converge toward a centre point. This will work if the 
%size(good data) >> size(bad data)


%=========================================================================================
%Find the centre of the main data, and also estimate the elipse radii using
%mean statistics (Normal Distribution):
%=========================================================================================
tolerance=1e-7;
x=X;
y=Y;
i=0;dX=1; dY=1;
while dX>tolerance || dY>tolerance
    avgX=mean(x(~isnan(x)));
    avgY= mean(y(~isnan(y)));
    stdX=std(x(~isnan(x)));
    stdY= std(y(~isnan(y)));
    
    %Remove the outermost layer of data:
    x=x(x>=avgX-3*stdX &x<avgX+3*stdY);
    y=y(y>=avgY-3*stdY &y<avgY+3*stdY);
    %recalculate the mean:
    avgXfew=nanmean(x);
    avgYfew=nanmean(y);
    %Calculate the resulting shift in the mean of the data:
    dX=abs(avgX-avgXfew);
    dY=abs(avgY-avgYfew);
    i=i+1;
    if Debug==1
        fprintf('\ni= %d, dDist = %08.7f, dVel = %08.7f, d=%03.1f, v=%03.1f',i,dX,dY,avgX,avgY)
    end
end
clear dis vel avgVfew avgDfew dVel dDist i
%========


%=========================================================================================
%Find the centre of the main data, and also estimate the elipse radii using
%median statistics:
%=========================================================================================
tolerance=1e-7;
xM=X; %for median calc
yM=Y; %for median calc
i=0;dXM=1; dYM=1;
while dXM>tolerance || dYM>tolerance
    %find medians:
    avgXM=median(xM(~isnan(xM)));
    avgYM=median(yM(~isnan(yM)));
    madX=median(abs(xM(~isnan(xM))-avgXM));
    madY=median(abs(yM(~isnan(yM))-avgYM));
    
    %Remove the outermost layer of data:
%     xM=xM(xM>=avgXM-3*madX &xM<avgXM+3*madX);
%     yM=yM(yM>=avgYM-3*madY &yM<avgYM+3*madY);
    in_are =( ( (abs(xM-avgXM)./(3*madX)).^2+(abs(yM-avgYM)./(3*madY)).^2 ) < 1 );
    xM=xM(in_are==1);
    yM=yM(in_are==1);
    
    %recalculate the mean:
    avgXMfew=median(xM(~isnan(xM)));
    avgYMfew=median(yM(~isnan(yM)));
    %Calculate the resulting shift in the mean of the data:
    dXM=abs(avgXM-avgXMfew);
    dYM=abs(avgYM-avgYMfew);
    i=i+1;
    if Debug==1
        fprintf('\ni= %d, dDistM = %08.7f, dVelM = %08.7f, d=%03.1f, v=%03.1f',i,dXM,dYM,avgXM,avgYM)
    end
end
clear disM velM avgDMfew avgVmfew dDistM dVelM i
%=========================================================================================
%Now let's assess the normal distribution of the data by comparing the
%median absolute deviation and standard deviation. std/mad =1.4826 if data
%is normal:
kx=stdX/madX;
ky=stdY/madY;
%It is not clear yet if this information can be actively used in the
%filter routine, it is an interesting evaluation parameter all the same.
ElipseCx=avgXM;
ElipseCy=avgYM;
ElipseRx=MADzZcore*madX;
ElipseRy=MADzZcore*madY;

if Debug==1
    [Ex,Ey]=calculateEllipse(ElipseCx,ElipseCy,ElipseRx,ElipseRy,0);
    hold on;
    plot(Ex,Ey,'--r','linewidth',1.4)
end
%========



%===================================================================
%  Distinguish data inside an elipsoid, at defined centre and radii.
%===================================================================
% To do this, I select bands of velocities, and assign distance limits for
% each of these bands. As the velocity bands come to the top of the Elipse,
% the distance cutoff limits tend toward one another (toward the avg
% Distance value)
XFilt=zeros(length(Y),1)./0; %create filtered copies of the original data
YFilt=zeros(length(Y),1)./0;
steps =10; % =1 will create a square block, =100 will create 100 velocity bands and thus 100 block elipse
[xLim,yLim]=calculateEllipseLims(ElipseRx,ElipseRy,steps); %calculate points on the elipse

for i=1:length(yLim)-1
    %We calculate the data indices for data falling between YBin(i) and YBin(i+1):
    Gd_indx=(abs(Y-ElipseCy)>=yLim(i) & abs(Y-ElipseCy)<yLim(i+1))...
        & (X>=ElipseCx-xLim(i) & X<ElipseCx+xLim(i));
    XFilt(Gd_indx)=X(Gd_indx); %set the bad data to Nan's
    YFilt(Gd_indx)=Y(Gd_indx); %set the bad data to Nan's
end
if Debug==1
    plot(XFilt,YFilt,'g.')
end
end
%========




%========Subfunction 1=====================================================
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

%========Subfunction 2=====================================================
function [X,Y] = calculateEllipse(x, y, a, b, angle, steps)
    % This functions returns points to draw an ellipse
    %  @param x         X coordinate
    %  @param y         Y coordinate
    %  @param a         Semimajor axis
    %  @param b         Semiminor axis
    %  @param angle     Angle of the ellipse (in degrees)
    narginchk(5, 6);
    if nargin<6, steps = 36; end

    beta = -angle * (pi / 180);
    sinbeta = sin(beta);
    cosbeta = cos(beta);
    alpha = linspace(0, 360, steps)' .* (pi / 180);
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
    X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
    Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);
    if nargout==1, X = [X Y]; end
end
%========