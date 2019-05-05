function [dist,Lvl,w] = WellenzahlLocMax(Level,Distance,Phase)
%WellenzahlLocMax looks at the top three local maxima
%Brian Scanlon, NUIG Mar 2018
szD=size(Distance);
if szD(1)>szD(2)
    Distance=Distance';
end
szL=size(Level);
if szL(1)>szL(2)
    Level=Level';
end

szP=size(Phase);
if szP(1)>szP(2)
    Phase=Phase';
end

% %Interpolate The Level and Distance, and increase the resolution:
% DistanceR=interp1(1:Res:length(Distance)*Res,Distance,1:length(Distance)*Res,'pchip','extrap');
% LevelR=interp1(1:Res:length(Level)*Res,Level,1:length(Level)*Res,'pchip','extrap');

%Find the peaks and sort them out:
[L,I]=findpeaks(Level,'MinPeakDistance',10); %Have a minimum spacing between peaks
[Lsrt,Isrt]=sort(L);
Dsrt=I(Isrt);
%Find the three highest values:
Lvl=flip(Lsrt(end-2:end));
dist=flip(Distance(Dsrt(end-2:end)));
w=flip(Phase(Dsrt(end-2:end)));
end

