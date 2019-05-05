function doy = dnum2doy(dnum)

% DNUM2DOY - Returns the day-of-year for the datenumber
% 
% DNUM2DOY will return the day-of-year for the corresponding 
% datenumber.
%
% Use As: doy = dnum2doy(dnum)
%
% Inputs: dnum = datenumber 
%
% Output: doy  = day-of-year
%

% Brian Ward, WHOI, 24-04-2002

%dvec=datevec(dnum);
[yr,mn,dy,h,m,s]=datevec(dnum);
mn1=ones(size(dnum));
dy0=zeros(size(dnum));
doy=dnum-datenum(yr,mn1,dy0);
