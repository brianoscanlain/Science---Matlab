% WDM.m  Wavelet Directional Method 
% Driver for producing WDM directional wave spectra plots
%           from a specified array of 3 or more wave staffs.
%
%  
% M.Donelan / Bergen 1994.11.26; Modified by A.Babanin and M.Donelan
% 2012.05 and earlier
clf
clear

np = 6;% Set number of wave staffs.
npp = (np*(np-1))/2;

%  load data array with each wave staff time series as a column.
for run = [62]

if run > 99
   eval(['load s87',int2str(run)])
elseif run > 9
  eval(['load s870',int2str(run)])
else
  eval(['load s8700',int2str(run)])
end
if length(data) > 4095
    
%data=binavg(data,2);   % binavg;  Modify as necessary to subsample time series.

% Reduce data to a multiple of 4096 points: wavelet analysis done 4096
% point chunks; i.e. 2^12 points
tlen = fix(length(data)/4096);
   data=data(1:4096*tlen,:);

% Apply calibrations [to metres] to the wave staffs and set in desired order.

ws=[.58467*data(:,3) .57756*data(:,1) .59762*data(:,4) ...
      .58998*data(:,6) .57155*data(:,5) .57363*data(:,2) ];
  
      ws = detrend(ws);% Remove mean and trend.
%
%   a and R are the angle [rad] and position [m] of the wave staffs 
%       from the array centre; i.e. Polar coordinates of the staffs   
%   np is number of probes
%
A = [0 33:72:359]; dr=pi/180;  A = dr*A; % 28 Oct 87 DB  Corrected June 22 1994
R = [0 0.25 0.25 0.25 0.25 0.25] ;   % Radius of array (m). i.e. Polar coordinates of the staffs

X = R.*cos(A);
Y = R.*sin(A);

l=0;x=[];y=[];
for j=1:np-1
for k=(j+1):np
   l=l+1;
   x(l)=X(k)-X(j);
   y(l)=Y(k)-Y(j);
end
end

i=sqrt(-1);
r=abs(x+i*y);
a=atan2(y,x);

l=0;rr=[];
for j=1:npp-1
for k=(j+1):npp
   l=l+1;
   rr(l)=a(j)-a(k);
   csj(l)=cos(a(j));
   csk(l)=cos(a(k));
   snj(l)=sin(a(j));
   snk(l)=sin(a(k));
   rk(l)=r(k);
   rj(l)=r(j);
end
end
rr=rr*180/pi;
ii=find(rr<0);
rr(ii)=rr(ii)+360;
ii=find(rr>70 & rr<110);
jj=find(rr>250 & rr<290);
ij=sort([ii jj]);

rr
clear ii jj X Y x y rr a r
ns=4;n=4096;lf=.0625;hf=1;nv=4;% Set parameters required by Wavelet.m 

for i1=1:n:length(ws(:,1))-n+1
i1
  lp=log(lf)/log(2);lp=floor(lp);
  hp=log(hf)/log(2);hp=ceil(hp);

 
  AMP = [];
  for jj=1:np
     eval(['wx',int2str(jj),' = [];'])
     eval(['B',int2str(jj),' = [];'])
     eval(['AMP',int2str(jj),' = [];'])
     x=ws(i1:i1+n-1,jj);x=x(:);
     eval(['[wx',int2str(jj),',f]=wavelet(x,lp,hp,nv,ns);'])
     eval(['B',int2str(jj),'=angle(wx',int2str(jj),');'])
     eval(['AMP',int2str(jj),'=abs(wx',int2str(jj),');'])
  end
  AMP = AMP1;
  for jj = 2:np
      eval(['AMP = AMP + AMP',int2str(jj),';'])
      eval(['clear AMP',int2str(jj),' wx',int2str(jj)])
  end
  AMP = AMP/np;
  clear AMP1 wx1 x

mxmf=max(find(f < ns/2));
kkm=[];thm=[];kks=[];ths=[];
%********************************************
for mf=1:mxmf
%********************************************

l=0;b=[];
for j=1:np-1
for k=(j+1):np
   l=l+1;
   eval(['b(:,l)=B',int2str(k),'(:,mf)-B',int2str(j),'(:,mf);'])
end
end
ill=find(b(:) > pi); b(ill)=b(ill)-2*pi;
ill=find(b(:) < -pi); b(ill)=b(ill)+2*pi;
lb=length(b);
l=0;AA=[];
for j=1:npp-1
for k=(j+1):npp
   l=l+1;
   if length(find(ij == l)) ==0
      AA(:,l)=zeros(lb,1); bk(:,l)=zeros(lb,1);bj(:,l)=zeros(lb,1);
   else
      AA(:,l)=rk(l)/rj(l)*b(:,j)./b(:,k);
      bk(:,l)=b(:,k);bj(:,l)=b(:,j);
   end
end
end

th=[];ll=0;

for l=ij
   ll=ll+1;
   th(:,ll)=atan2((AA(:,l)*csk(l)-csj(l)),(snj(l)-AA(:,l)*snk(l)));
   kk(:,ll)=(snk(l)*bj(:,l)/rj(l)-snj(l)*bk(:,l)/rk(l)) ./ ...
           (csj(l)*snk(l) - csk(l)*snj(l))./cos(th(:,ll));
end
th = th+(kk<0)*pi;
kk = abs(kk);
thm = [thm meanang(th')'];
ths = [ths std(th')'];
kkm = [kkm mean(kk')'];
kks = [kks std(kk')'];

end        %%%  end of 'for mf=1:mxmf'  loop (frequency bin loop)
clear AA b bj bk ill th kk j k l

if i1 == 1 
   AAp=AMP;
   ddd=round(thm*180/pi);
   kkmp=kkm;
   eval(['save yw',int2str(run),' AAp ddd kkmp f'])
else
  eval(['load yw',int2str(run)])
  AAp=[AAp;AMP];
  kkmp=[kkmp;kkm];
  ddd=[ddd;round(thm*180/pi)];
  eval(['save yw',int2str(run),' AAp ddd kkmp f'])
end 
clear AAp ddd kkmp f
end

clear B1
 for jj = 2:np
      eval(['clear B',int2str(jj)])
  end

clear csj csk dr hf hp kkm kks lb ll lf mf mxmf mw A R
clear n rj rk snj snk thm ths wn t i1 ij jj lp data AMP

waveplots
wavenums
% wavelognums
else
end   % end of 'if length(data) > 4096' loop

end

