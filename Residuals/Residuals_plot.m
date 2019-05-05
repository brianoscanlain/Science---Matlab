%Run relation_via_binned_u10.m first!!!
%
% tight subplot of remaining residuals
cd('C:\Users\Brian\Desktop\Whitecap_Project\Images');

load('data_wave_chl.mat', 'chl_flr')
load('data_wave_chl.mat', 'chl_sat')
load('WindHist.mat', 'WSslope')
ha = tight_subplot_gjs(3,1,[.05 .05],[.15 .1],[.15 .15]);
axes(ha(1));
[Ax,h1,h2]=plotyy(1:206,Resid_A,1:206,Resid_B);
set(h1,'linewidth',2)
set(h2,'linewidth',2)
set(Ax(:),'XLim',[1 206])
set(Ax(:),'Fontsize',12)
set(Ax(:),'XTickLabel',[])
set(get(Ax(1),'Ylabel'),'String','Wa(%) residual','FontSize',12)
set(get(Ax(2),'Ylabel'),'String','Wb(%) residual','FontSize',12)

axes(ha(1));
title('Residuals','fontsize',14)
axes(ha(2));
axes(ha(1));
axes(Ax(1));
axes(Ax(2));
axes(ha(2));

load('WindHist.mat', 'WDslope')
[Ax,h1,h2]=plotyy(1:206,WSslope(2:end),1:206,WDslope(2:end));
set(h1,'linewidth',2)
set(h2,'linewidth',2)
set(Ax(:),'XLim',[1 206])
set(Ax(:),'Fontsize',12)
set(Ax(:),'XTickLabel',[])
set(get(Ax(1),'Ylabel'),'String','Wind speed slope','FontSize',12)
set(get(Ax(2),'Ylabel'),'String','Wb(%) WInd dir slope','FontSize',12)
axes(ha(3));

[Ax,h1,h2]=plotyy(1:206,chl_flr(2:end),1:206,chl_sat(2:end));
set(h1,'linewidth',2)
set(h2,'linewidth',2)
set(Ax(:),'XLim',[1 206])
set(Ax(:),'Fontsize',12)

set(get(Ax(1),'Ylabel'),'String','Chl A (fluor)','FontSize',12)
set(get(Ax(2),'Ylabel'),'String','Chl A (sat)','FontSize',12)
xlabel('10 min image periods')