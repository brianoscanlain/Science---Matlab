% Temperature related variables are evaulated and seen if any relationships 
% exist with the residuals:

cd('C:\Users\Brian\Desktop\Whitecap_Project\Images');

load data airt sst;
load Residuals.mat;

sst(1)=[]; airt(1)=[];  %remove the Nan data

% quick plot of the residuals:
[Ax,h1,h2]=plotyy(1:206,Residual_A,1:206,Residual_B);



set(Ax(:),'Fontsize',14)
set(Ax(:),'XLim',[1 206])

set(h1,'linewidth',2)
set(h2,'linewidth',2)

axes(Ax(1))
ylabel('Wa residual','fontsize', 14);

axes(Ax(2))
ylabel('Wb residual','fontsize', 14);

% exponential relationship between the residuals and Atmospheric stability?
% {\Delta}T = (sst-airt)

x1 = sst-airt;
x2 = sst;
y1= Residual_A;
y2= Residual_B;

% Y = Ae^{Bx}
% ln(Y) = ln(A) + Bx

%least squares fitting for exponentials::
% mathworld.wolfram.com/LeastSquaresFittingExponential.html
%

[a1 b1] = exp_fit_leastsq(x1,y1);
[a2 b2] = exp_fit_leastsq(x1,y2);

[a3 b3] = exp_fit_leastsq(x2,y2);
[a4 b4] = exp_fit_leastsq(x2,y2);

 


