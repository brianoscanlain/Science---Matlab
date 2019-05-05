function [dist,Lvl,Phi] = WellenzahlLocMaxVel(Level,Distance,phi_angle)
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


%Find the peaks and sort them out:
[L,D]=findpeaks(Level,'MinPeakDistance',10); %Have a minimum spacing between peaks
[Lsrt,Isrt]=sort(L);
Dsrt=D(Isrt);
%Find the three highest values:
Lvl=flip(Lsrt(end-2:end));
dist=flip(Distance(Dsrt(end-2:end)));
Phi=flip(phi_angle(Dsrt(end-2:end)));

end

