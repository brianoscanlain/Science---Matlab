function [XFilt,YFilt,kx,ky,PC,kx2nd,ky2nd,PC2nd,NansAboveThresh,Scenario] =ElypticFiltBiModal(X,Y,error,intrp,Zscore,SecondPass)
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
%      interp = 1: interpolates the output XFilt and YFilt signals
%             = 0: fills bad-data entries in output signals with Nan values
%      ZScore - The scale (in units of Zscore or sigma) upon which Ellyptical Filt
%               works. Set to 7 for Symeo data, and 12 for IMU data. As the
%               IMU is not normally distributed, as it is a function of
%               SHM. e.g. Y=sin(t); Y is not normally distributive... hence
%               the larger Zscore to allow for my code to handle such distributions.
%      SecondPass - 1 Perform Second Pass
%                   0 Just do the First pass.
%
% Output:
%--------
%     XFilt - filtered data
%     YFilt - Filtered data
%     kx - relation between std and mad of X, proxy for Normal distribution.
%     ky - relation between std and mad of Y, proxy for Normal distribution.
%     PC - percentage of data which fell outside filter constraints(bad
%          data found)
%     kx2nd - same as above, except for the Second pass
%     ky2nd - same as above, except for the Second pass
%     PC2nd - same as above, except for the Second pass
%     Scenario -    1 = well-behaved data (outliers <20%)
%                   2 = poorly-behaved data (outliers >20% & <10%)
%                   3 = Distinct bimodal distributive data
%                   4 = Manual Filtering (no modes found)
%
% Brian Scanlon, NUIG, 13th Feb 2018
% BS. Sept 2018: Added interp function, and added PC output.
% BS October 2018: full rework, added bimodal, manual and runmedian
% filtering functionalities.
%-------------------------------------------------------
% %Debugging and testing artifacts:
%-------------------------------------------------------
%Set the following coefficients in the case where we need to do manual
%filtering (this is a last resort!)
maxDist = 22;     %max expected distance
MinDist = 10;     %min expected distance
MinVel  = -5;     %min expected velocity (-5 m/s is a lot for Symeo data)
MaxVel  =  5;     %max expected velocity (+5 m/s is a lot for symeo data)
FracGoodContentThresh=0.80; %check that there is more than 80% of good data after the filter
FracGoodContentThreshRunMedian=0.50; %check that there is more than 50% of good data after the runmedian filter
runMedianWindowFast=70; %Moving window size used for the FAST runmedian in the First pass when more than 20% of bad data is present upon first pass.
runMedianWindowFastSP=30; %Moving window size used for the FAST runmedian in the Second pass
runMedianWindowSlow=1500; %Slow moving window size for Second Pass.
GapThreshold=500; %limit for interpolating series of Nan's longer than this limit
Debug=0;
VerboseEllipticOutput=0; %Outputs the iterative convergence
%---------------


if Debug==1
    X=data.radar.dist;
    Y=data.radar.vel;
    %X=data.IMU.acc(1,:);
    %Y=data.IMU.gyro(1,:);
    intrp=1;
    error=[]; %keep it empty
    Zscore=7;
    SecondPass=1;
    close all;
    %     plot(X,Y,'.')
    %     set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
    %     xlabel('Distance [m]','fontsize',16,'fontname','times')
    %     ylabel('Velocity [m/s]','fontsize',16,'fontname','times')
end
%Allocate empty Outputs incase interp=0, or secondpass=0, etc.
kx2nd=[]; ky2nd=[]; PC2nd=[]; NansAboveThresh=[];  Scenario=[];
%-------

%-------------------------------------------------------
%Check to see that the VargIn have the correct structure:
%-------------------------------------------------------
[~, dimX]=max(size(X));
[~, dimY]=max(size(Y));
if dimX>1
    X=X';
end
if dimY>1
    Y=Y';
end
%-------

%-------------------------------------------------------
%Check inputs for "error" and "intrp"
%-------------------------------------------------------
if Debug~=1
    if nargin>3
        if intrp~=0
            if intrp~=1
                error('input ''intrp'' must be a logical (== 1 | 0)');
            end
        end
    else
        intrp=0; %if not specified, just keep it turned off.
    end
    if nargin>2 && ~isempty(error)
        X(error>0)=0./0; %Set the error codes above zero to Nan values (Error code 0 is good data)
        Y(error>0)=0./0; %Set the error codes above zero to Nan values (Error code 0 is good data)
    end
    % Check "Zscore" input
    if nargin>4
        %we define Zscore!
        defineZ=0;
    else
        defineZ=1;
    end
else
    defineZ=1; %if debug=1
end
%---------


%-------------------------------------------------------
%Tuning parameters
%-------------------------------------------------------
if defineZ   %I have noticed that IMU data can have zscores up to 12 (i.e. a different distribution)
    Zscore=12;  %default was 4, however I set it to 6 to relax the eliptical tolerance.
    Zscore=7;
end         %For Symeo w/ Kalman, I relax it more as the Kalman can often
%go outside the natural distibution. While such data may not be
%real, it is itself a smoothened "fill", which we would not be
%able to match. Instead our fills can be compared with linear
%interps. Hence we try to favour Kalman extrapolation above our
%solutions.
%---------


%-------------------------------------------------------
%Filtering routine
%-------------------------------------------------------
%Firstly, all negative distances are set to Nan:
negInd=X<0;
X(negInd)=0./0;
Y(negInd)=0./0;
%--------------
%First pass is performed:
[XFilt,YFilt,kx,ky,PC] = ElypticFilter(X,Y,Zscore,Debug,VerboseEllipticOutput);
%If the First pass removes a significant amount of data, we need to check
%that the presence of bimodal distribution within the data.
fprintf('~~\nElliptic Filter: %05.2f%% good data returned on First pass.\n',sum(~isnan(XFilt))/length(X)*100)
if sum(~isnan(XFilt))/length(X)<FracGoodContentThresh
    fprintf('Checking for secondary distributive modes\n')
    %Here I recalculate the primary mode with a wider Zscore scaler so that
    %we can clean up the mode that much better.
    [XFilt1,YFilt1,kx,ky,PC] = ElypticFilter(X,Y,Zscore*2,Debug,VerboseEllipticOutput);
    %Rerun the ElypticFilter, but set the previously-found good data to
    %Nan:
    X2=X;  X2(~isnan(XFilt1))=0./0;
    Y2=Y;  Y2(~isnan(XFilt1))=0./0;
    [XFilt2,YFilt2,kx2,ky2,PC2] = ElypticFilter(X2,Y2,Zscore,Debug,VerboseEllipticOutput);
    %This is a very neat trick (If I don't say so myself) as it removes the
    %first mode of data and searches for the secondary mode. In cases where
    %kalman filtering is used at a hardware level, multiple modes can be
    %present in the data. This arises when the Kalman filter "latches" onto
    %a mode which it thinks is the appropriate range where good data is
    %found. With weak signals, the Kalman can jump onto different levels,
    %especially ocean wave elevation data!
    %
    %Now we must see if 2nd Bimode data is larger
    %
    if sum(~isnan(XFilt2))/length(X(isnan(XFilt)))>FracGoodContentThresh
        fprintf('Secondary bimode found, cleaning the time series now...\n')
        Scenario=3;
        if Debug
            fig1=figure;
            plot(X,'.') %plot raw data
            hold on
            plot(XFilt,'.','markersize',8) %Plot first guess filter (likely bimodal)
            plot(XFilt2,'.','markersize',8)%Plot filter results of second mode (with first mode removed)
        end
        if 1 %Apply a median filter to catch outliers!
            fprintf('Applying filter using median moving average on 2nd mode\n');
            XFiltruned=runmedian(XFilt2,runMedianWindowFast);
            XFiltnorm=XFilt2-XFiltruned;
            YFiltnorm=YFilt2-runmedian(YFilt2,runMedianWindowFast);
            %[XFilt3,YFilt3,kx2,ky2,PC2] = ElypticFilter([diff(XFilt2);0],[diff(YFilt2);0],12,Debug);
            [XFilt3,YFilt3,kxx,kyy,P] = ElypticFilter(XFiltnorm,YFiltnorm,Zscore*2,Debug,VerboseEllipticOutput);
            if sum(~isnan(XFilt3))/length(isnan(XFilt2))>FracGoodContentThreshRunMedian
                fprintf('runmedian successful\n');
                XFilt2(isnan(XFilt3)) = 0./0;
                YFilt2(isnan(YFilt3)) = 0./0;
                if Debug
                    figure(fig1);
                    plot(XFiltruned,'k--');
                    plot(XFilt2,'.','markersize',8);
                    legend('raw','1st guess','bimodal filtered','runmedian','runmedian filtered')
                end
                PC2=P;
                kx2=kxx;
                ky2=kyy;
                fprintf('runmedian not applied to second mode\n');
            else
            end
        end
        XFilt=XFilt2;
        YFilt=YFilt2;
        PC=PC2;
        kx=kx2;
        ky=ky2;
    else    %No sign of a second Mode!
        fprintf('No Secondary mode found.\n')
        %There are two scenarios, 1. the Ellyptical filter fails to notice any
        %modes, thus, we will need to Manually filter to see if we can recover something.
        %2. The Ellyptical filter works fine, but there is a lot of bad data. So as a second
        %measure, we can apply a RunMedian to clean up the data.
        if sum(~isnan(XFilt))/length(X)<0.1 && sum(~isnan(XFilt2))/length(X(isnan(XFilt)))<0.1
            fprintf('Filter cannot identify any modes.\n')
            fprintf('Applying manual filter on the data now...\n')
            rmInd= (X>maxDist) + (X<MinDist) + (Y<MinVel) + (Y>MaxVel)>0; %find data outside window
            fprintf('Manual Filtering: %05.2f%% good data found\n',sum(~rmInd)/length(rmInd)*100);
            %chop off outliers manually
            X2=X; X2(rmInd)=0./0;
            Y2=Y; Y2(rmInd)=0./0;
            [XFilt,YFilt,kx,ky,PC] = ElypticFilter(X2,Y2,Zscore,Debug,VerboseEllipticOutput);
            Scenario=4;
        else % Let's run a runnimg mean on the data and clean it up!
            % This is now commented out as I perform a compulsory
            % runmedian filter right after this nested loop series!
            fprintf('Applying filter using median moving average...\n');
            X2=XFilt-runmedian(XFilt,runMedianWindowFast);
            Y2=YFilt-runmedian(YFilt,runMedianWindowFast);
            [XFilt4,YFilt4,kx,ky,PC] = ElypticFilter(X2,Y2,Zscore,Debug,VerboseEllipticOutput);
            Scenario=2;
            if sum(~isnan(XFilt4))/length(XFilt)>FracGoodContentThreshRunMedian %We will only allow runmedian to filter half the data as bad.
                XFilt(isnan(XFilt4)) = 0./0;
                YFilt(isnan(YFilt4)) = 0./0;
            else %Here the last ellyptical filter found no signature of a distributive mode
                fprintf('runmedian not applied.\n');
            end
        end
        
    end
    
else %Data passes - Scenario 1 assigned!
    Scenario=1;
end
%---------


%----------------------------------------------------------------------
% Second pass (as Gandalf said (to the outliers), you shallll not pass!
%----------------------------------------------------------------------
%We require a second catch for bimodal data which are somewhat merged,
%and the elliptical filter cannot distinguish between them, but instead
%passes both modes as one. In such cases, I propose doing a running median
%on the X and Y data. We can then, for example, get the standard deviation
%of the runmedian, and check if it is similar (in order of magnitude) to
%that of the raw data. I suspect that the raw data will have a larger std
%in the cases where two modes are passed as one.
%   A simpler solution (which may be more successful) is to run an
%ellyptical filter on the "raw-runmedian(raw)" data.
if SecondPass
    X2=runmedian(XFilt,runMedianWindowSlow);
    Y2=runmedian(YFilt,runMedianWindowSlow);
    [XFilt5,YFilt5,kx2ndSlow,ky2ndSlow,PC2ndSlow] = ElypticFilter(XFilt-median(X2(~isnan(X2))),YFilt-median(Y2(~isnan(Y2))),Zscore/1.2,Debug,VerboseEllipticOutput);
    XFilt(isnan(XFilt5)) = 0./0; %This refines the selection of the first pass, enabling a good runmean to be found!
    YFilt(isnan(YFilt5)) = 0./0; %I haven't encountered this to ruin the data, so I directly modify the original first-pass data
    X4=runmedian(XFilt,runMedianWindowFastSP);
    Y4=runmedian(YFilt,runMedianWindowFastSP);
    [XFilt5,YFilt5,kx2nd,ky2nd,PC2nd] = ElypticFilter(X-median(X4(~isnan(X4))),Y-median(Y4(~isnan(Y4))),Zscore/1.5,Debug,VerboseEllipticOutput);
    if Debug
        figure
        plot(X,'.') %plot raw data
        hold on
        plot(XFilt,'.','markersize',8);
        hold on
        plot(X2,'k--','linewidth',1.3)
        plot(X4,'r--','linewidth',1.3)
        if VerboseEllipticOutput
            fprintf('distance: std(unfilt) / std(runmedian) = %05.2f%%\n',(nanstd(XFilt)/nanstd(X4))*100);
            fprintf('Velocity: std(unfilt) / std(runmedian) = %05.2f%%\n',(nanstd(YFilt)/nanstd(Y4))*100);
        end
    end
    %Check if the Second pass succeeded, if it removes more than half the
    %data, do nothing and maintain the old! For example, this returns false
    %if manual chop is performed
    if sum(~isnan(XFilt5))/length(XFilt)>FracGoodContentThreshRunMedian
        XFilt=X; XFilt(isnan(XFilt5))=0./0;
        YFilt=Y; YFilt(isnan(YFilt5))=0./0;
        fprintf('~~\nElliptic Filter: %05.2f%% good data returned on Second pass.\n',sum(~isnan(XFilt))/length(X)*100)
    else
        fprintf('~~\nSecond pass fails, reverting to last best: %05.2f%% good data returned.\n',sum(~isnan(XFilt))/length(X)*100)
        kx2nd=kx2ndSlow;
        ky2nd=ky2ndSlow;
        PC2nd=PC2ndSlow;
    end
end
%---------


%-------------------------------------------------------
%Interpolation and reconstruction of the time series:
%-------------------------------------------------------
%Here I have started some ideas for reconstructing the gappy time series.
%Some ideas would be to use resample(), interp1 and/or complex solutions
%such as kalman filtering or specialized wave-signal reconstruction
%techniques. Another idea is to only interpolate small gaps with interp1()
%and reconstruct medium gaps with a fancy technique, and then leave large
%gaps as Nan's (these gaps could then be passed through the code and
%handlered dynamically/intelligently - for example, the interpolation routine
%can shorten the time series, or the WDM could have a fancy decision process to
%find and weigh poor sections, and minimise their influence on the final mean
%(pair averaged) wave estimates.
if sum(isnan(XFilt))~=0
    if Debug
        [NansAboveThresh50,GapIdx,Gaplen] = GapFinder(XFilt,50);
        fprintf('%d gap(s) of Nan''s present...\n',length(GapIdx));
        fprintf('%05.2f%% Nan''s in gaps of length less than %d points\n',sum(~NansAboveThresh50)/length(isnan(XFilt))*100,50);
        [NansAboveThresh500,~,~] = GapFinder(XFilt,500);
        fprintf('%05.2f%% Nan''s in gaps of length greater than %d points\n',sum(NansAboveThresh500)/length(isnan(XFilt))*100,500);
        [NansAboveThresh5000,~,~] = GapFinder(XFilt,5000);
        fprintf('%05.2f%% Nan''s in gaps of length greater than %d points\n',sum(NansAboveThresh5000)/length(isnan(XFilt))*100,5000);
        if length(GapIdx)>30 %Only show a plot if there are more than (say) 30 gaps
            figure, plot(GapIdx,Gaplen); ylabel('Gap size')
        end
    end
    %GapThreshold=500;   %Now set in preamble
    [NansAboveThresh,GapIdx,Gaplen] = GapFinder(XFilt,GapThreshold);
    %Interpolate the Nans!
    if intrp && ~isempty(XFilt(~NansAboveThresh))
        x=linspace(1,length(XFilt),length(XFilt))';
        XFilt(~NansAboveThresh)=interp1(x(~isnan(XFilt)), XFilt(~isnan(XFilt)),x(~NansAboveThresh),'pchip','extrap');
        YFilt(~NansAboveThresh)=interp1(x(~isnan(YFilt)), YFilt(~isnan(YFilt)),x(~NansAboveThresh),'pchip','extrap');
        %fprintf('interpo finito!');
        %Prevent any erroneous extrapolation at the end of the data series:
        
        if GapIdx(end)+(Gaplen(end)-1)==length(XFilt)
            XFilt(GapIdx(end)-1:end)=0./0;
            YFilt(GapIdx(end)-1:end)=0./0;
        end
        %Also prevent erroneous extrapolation at the start of the data
        %series!
        if GapIdx(1)==1
            XFilt(GapIdx(1):Gaplen(1))=0./0;
            YFilt(GapIdx(1):Gaplen(1))=0./0;
        end
        
    end
else
    %No Nan's, do nothing!!
end
%Finally plot the difference!
if Debug
    figure,
    plot(X,'.')
    hold on
    plot(XFilt,'.')
    ylabel('Distance [m]')
    if exist('fileNm','var')
        title(['i= ' num2str(i) ' - ' fileNm(1:13) '   ' fileNm(15:end-4)])
    end
    if sum(isnan(XFilt))>0
        legend('raw',['filtered (' num2str(PC) '% bad data, ' ...
            num2str(round(sum(isnan(XFilt))/length(XFilt)*100*100)/100)  '% Nan''s)'] );
    end
end
NansAboveThresh=sum(NansAboveThresh);
%=========================================================================
end
%==========================================================================
%==========================================================================
%==========================================================================




















%==========================================================================
% Subfunction 1
%==========================================================================
function [XFilt,YFilt,kx,ky,PC] =ElypticFilter(X,Y,Zscore,Plotting,Verbose)
%
MAD_STD_scale_factor=1.4826;
MADzZcore=Zscore*MAD_STD_scale_factor;
%First we need to find the centre of the "good" data. To do this, I first
%start off with the median. Next we purge outliers, and keep doing this
%until we converge toward a centre point. This will work if the
%size(good data) >> size(bad data)
%========================================================================
%Find the centre of the main data, and also estimate the elipse radii using
%mean statistics (Normal Distribution):
%========================================================================
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
    if Verbose==1
        fprintf('i= %d, dDist = %08.7f, dVel = %08.7f, d=%03.1f, v=%03.1f\n',i,dX,dY,avgX,avgY)
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
    if Verbose==1
        fprintf('i= %d, dDistM = %08.7f, dVelM = %08.7f, d=%03.1f, v=%03.1f\n',i,dXM,dYM,avgXM,avgYM)
    end
end
clear disM velM avgDMfew avgVmfew dDistM dVelM i
%=========================================================================================
%Now let's assess the normal distribution of the data by comparing the
%median absolute deviation and standard deviation. std/mad =1.4826 if data
%is normal:
%kx=stdX/madX; %normaity scale for the elyptical data
%ky=stdY/madY;
%Let's output this for the original raw data:
kx=std(X(~isnan(X))/mad(X(~isnan(X))));
ky=std(Y(~isnan(X)))/mad(Y(~isnan(Y)));

%It is not clear yet if this information can be actively used in the
%filter routine, it is an interesting evaluation parameter all the same.
ElipseCx=avgXM;
ElipseCy=avgYM;
ElipseRx=MADzZcore*madX;
ElipseRy=MADzZcore*madY;

if Plotting==1
    [Ex,Ey]=calculateEllipse(ElipseCx,ElipseCy,ElipseRx,ElipseRy,0);
    figure
    plot(X,Y,'.')
    set(gca,'fontname','times','fontsize',12,'linewidth',1.20)
    xlabel('Distance [m]','fontsize',16,'fontname','times')
    ylabel('Velocity [m/s]','fontsize',16,'fontname','times')
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
if Plotting==1
    plot(XFilt,YFilt,'g.')
end

%Calculate the percentage content of bad data (which didn't make it through filter):
PC=100*((sum(isnan(XFilt)) + sum(isnan(YFilt)))/2)/length(XFilt); %take the mean, but they should both be the same!
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
% This function returns points to draw an ellipse
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




%========Subfunction 3=====================================================
function [NansAboveThresh,Gapstrt,Gaplen] = GapFinder(x,thres)
% This function finds Nan gaps in the X data
% Inputs:
% x = data to measure [Nx1] array
% thres = the threshold value
Nans=find(isnan(x));
prevIdx=0;
start=1;
gapIdx=1;
for i=1:length(Nans)
    if i==1
        Gapstrt(gapIdx)=Nans(i);
        counter=1;
    else
        if Nans(i)-prevIdx==1 %Nan is part of current gap series
            counter=counter+1;
        else
            Gaplen(gapIdx)=counter; %write count to GapLen
            counter=1; %reset counter
            gapIdx=gapIdx+1;
            Gapstrt(gapIdx)=Nans(i);
        end
        if i==length(Nans) %last nan:
            Gaplen(gapIdx)=counter; %write count to GapLen
        end
    end
    prevIdx=Nans(i);
end

%Find all the indices of x which are part of gaps greater than "thres"
NansAboveThresh=x.*0==1;
bigGaps=find(Gaplen>thres);
if ~isempty(bigGaps)
    for ii=1:length(bigGaps)
        NansAboveThresh(Gapstrt(bigGaps(ii)):Gapstrt(bigGaps(ii))+Gaplen(bigGaps(ii))-1)=1;
    end
end
end
%========


%========Subfunction 4=====================================================
function [Xout] = runmedian(x,window)
% This function acts as a running or moving median average window. Can
% handle Nan's. If the window has only Nan's, it returns a Nan for that
% particular datapoint!
%
%Brian Scanlon, October 2018

% Make sure the x input is the correct shape:
if size(x,1)>size(x,2)
    x=x';
    transPos=1;
else
    transPos=0;
end

%Prepad the start and end!
preWin=round(window/2);
postWin=round(window/2)-1;
Xout=[fliplr(x(1:preWin)) x fliplr(x(end-postWin:end))];


for i=1:length(x)
    buffer=(Xout(i:i+preWin+postWin));
    buffer=buffer(~isnan(buffer));
    if ~isempty(buffer)
        Xout(i+preWin)=median(buffer);
    else
        Xout(i+preWin)=0./0;
    end
end
%Chop off the pre and post padding.
Xout=Xout(preWin+1:end-postWin-1);

%Revert to initial shape (if reshaping was performed earlier)
if transPos
    Xout=Xout';
end
end
%========


%========Subfunction 4=====================================================
function [Xout] = mad(x)
% This function acts as a running or moving median average window. Can
% handle Nan's. If the window has only Nan's, it returns a Nan for that
% particular datapoint!
%
%Brian Scanlon, October 2018
Xout = median(abs(x-median(x)));
end