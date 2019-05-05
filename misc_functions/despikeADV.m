function [u1, Id] = despikeADV(u, hx, hy)
% This code removes spikes from the Accoustic Doppler Velocimeter (ADV)
% measurements. It uses bivariate kernel density function to separate the data
% cluster from the spike clusters. The bivariate kernel density function is
% generated using Botev et al. (2010)'s FFT based code. The elements
% identified as spikes are removed and replaced by the linearly interpolated
% values.
%
% Input variables: 
% u is time series of velocity [Nx1],
% hx and hy are the band widths. If no value is provided, the code will
% assume hx = hy = 0.01;
%
% Output variables:
% u1 is the despiked and reconstructed time series. 
% Id is the index of the elements that were identified as spikes and replaced. 
%
% References: 
% (1)Islam, M.R. and Zhu, D.Z. (2013), A Kernel Density Based algorithm to Despike ADV Data,
% accepted in Journal of Hydraulic Engineering, ASCE.
% (2)Botev,Z.I., Grotowski,J.F., and Kroese,D.P.(2010),Kernel Density Estimation Via Diffusion,
% Annals of Statistics, 38(5),2916-2957."
%

global N 

N = length(u); 
        
% Calculating du
for i= 2:N-1               
    db= u(i)-u(i-1);
    df = u(i+1)-u(i);
    if abs(db)>abs(df)
        du(i)=df;
    else
        du(i)=db;
    end
end
du=[ du  0]';
u1=u; w1=du;

% Axis rotation
th = atan2((N*sum(u1.*w1)-sum(u1)*sum(w1)),(N*sum(u1.*u1) -sum(u1)*sum(u1)));
     ut = (u1)*cos(th)+ (w1)*sin(th);
     wt = -(u1)*sin(th)+(w1)*cos(th);
data =[ut wt];
% Applying kernel density function using Botev et al.(2010)'s algorithm
if nargin <2
    hx = 0.01; hy = 0.01;
end
[d,X,Y] = kde2dBotev(data,N,hx,hy);

uf = X(1,:); wf = Y(:,1);
dp = max(max(d));
[wp,up] = find(d==dp);

% Calculating cut-off threshold
c1=0.4;c2=0.4;
[ul,uu]= cutoff1(dp,uf,c1,c2,d(wp,:),up);
[wl,wu]= cutoff1(dp,wf,c1,c2,d(:,up)',wp);

% Calculating axes of ellipse and identifying spikes
uu1= uu-0.5*(uu+ul); ul1= ul-0.5*(uu+ul);
wu1= wu-0.5*(wu+wl); wl1= wl-0.5*(wu+wl);
Ut1= ut-0.5*(uu+ul); Wt1= wt-0.5*(wu+wl);
F=zeros(N,1);

at= 0.5*(uu1-ul1); bt = 0.5*(wu1-wl1);
for i = 1:N
    if Ut1(i)>uu1| Ut1(i)<ul1
        F(i)=1;
    else
        we=sqrt((bt^2)*(1-(Ut1(i)^2)/(at^2)));
        if Wt1(i)>we | Wt1(i)<-we 
            F(i)=1;
        end
    end
end
Id= find(F>0);

%Replacing spikes by linearly interpolated values
if Id(1)==1                
   Id(1)=[]; u(1)= mean(u);
end
if Id(length(Id))==N
    Id(length(Id))=[]; u(N)=mean(u);
end
[u1]= velreplace(u,Id);
end

function [ul,uu]= cutoff1(dp,uf,c1,c2,f,Ip)

lf = length(f);

dk =[0 diff(f)]*256/dp;             
for i = Ip-1:-1:2
  if f(i)/f(Ip)<=c1 &&  abs(dk(i))<=c2
      i1=i;
      break;
  end
end

for i = Ip+1:1:lf-1
   if f(i)/f(Ip)<=c1 &&  abs(dk(i))<=c2 
        i2= i;
        break;
   else  % added by Brian Scanlon
       i2=lf-1;
   end
end


ul= uf(i1);
uu= uf(i2);
end

function [u]= velreplace(u,Id);
global N;
I1 = [Id;N];
  for i = 1: length(I1)-1
    for j = i: length(I1)-1
       if I1(j+1)-I1(j)>1
          break;
       end
    end
      ds=j-i+2; du1= (u(I1(j)+1)-u(I1(i)-1))/ds;
      u(I1(i))= u(I1(i)-1)+du1;
            
  end
end

function [density,X,Y]=kde2dBotev(data,N,hx,hy)

% This function is a modified version of the kernel density estimation code
% developed by Botev et al. (2010). The original source code is available
% at http://www.mathworks.com/matlabcentral/fileexchange/authors/27236
% The Reference is available at "Botev,Z.I., Grotowski,J.F., and Kroese,D.P.(2010),
% KERNEL DENSITY ESTIMATION VIA DIFFUSION,Annals of Statistics, 38(5),2916-2957."

global A2 I
n=256;

% Calculating scaled data
%N=size(data,1);
MAX_XY=max(data,[],1); MIN_XY=min(data,[],1); 
scaling =MAX_XY-MIN_XY;
 
transformed_data=(data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);

%bin the data uniformly using regular grid;
initial_data=ndhist(transformed_data,n);
% discrete cosine transform of initial data
a= dct2d(initial_data);
% now compute the optimal bandwidth^2
I=(0:n-1).^2; A2=a.^2;
t_star=fzero( @(t)(t-evolve(t)),[0,0.1]);
p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
t_y = hy^2; t_x = hx^2; 

% smooth the discrete cosine transform of initial data using t_star
a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
% now apply the inverse discrete cosine transform
if nargout>1
    density=idct2d(a_t)*(numel(a_t))/prod(scaling);
    
    [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
end
end
%#######################################
function  [out,time]=evolve(t)
global N
Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
time=(2*pi*N*Sum_func)^(-1/3);
out=(t-time)/time;
end
%#######################################
function out=func(s,t)
global N
if sum(s)<=4
    Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
    time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
    out=psi(s,time);
else
    out=psi(s,t);
end

end
%#######################################
function out=psi(s,Time)
global I A2
% s is a vector
w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
wx=w.*(I.^s(1));
wy=w.*(I.^s(2));
out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
end
%#######################################
function out=K(s)
out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
end
%#######################################
function data=dct2d(data)
% computes the 2 dimensional discrete cosine transform of data
% data is an nd cube
[nrows,ncols]= size(data);
if nrows~=ncols
    error('data is not a square array!')
end
% Compute weights to multiply DFT coefficients
w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
weight=w(:,ones(1,ncols));
data=dct1d(dct1d(data)')';
    function transform1d=dct1d(x)

        % Re-order the elements of the columns of x
        x = [ x(1:2:end,:); x(end:-2:2,:) ];

        % Multiply FFT by weights:
        transform1d = real(weight.* fft(x));
    end
end
%#######################################
function data = idct2d(data)
% computes the 2 dimensional inverse discrete cosine transform
[nrows,ncols]=size(data);
% Compute wieghts
w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
weights=w(:,ones(1,ncols));
data=idct1d(idct1d(data)');
    function out=idct1d(x)
        y = real(ifft(weights.*x));
        out = zeros(nrows,ncols);
        out(1:2:nrows,:) = y(1:nrows/2,:);
        out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
    end
end
%#######################################
function binned_data=ndhist(data,M)
% this function computes the histogram
% of an n-dimensional data set;
% 'data' is nrows by n columns
% M is the number of bins used in each dimension
% so that 'binned_data' is a hypercube with
% size length equal to M;
[nrows,ncols]=size(data);
bins=zeros(nrows,ncols);
for i=1:ncols
    [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
    bins(:,i) = min(bins(:,i),M);
end
% Combine the  vectors of 1D bin counts into a grid of nD bin
% counts.
binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));

end




















