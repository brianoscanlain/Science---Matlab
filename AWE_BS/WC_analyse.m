%File to read the WC.txt file for the deployment 3, creating a matlab
%file of the whitecap values with a timestamp.
clear all;
cd('F:\Dropbox (AirSea Laboratory)\Whitecap_Data\SSWP_subjectivity_analysis\Proc1\raw_images3');

WC=load('WC.txt');

WC(:,2)=WC(:,2)+201100e+11;
times=num2str(WC(:,2),'%u');

dte_W=datenum(str2num(times(:,1:4)),str2num(times(:,5:6)),str2num(times(:,7:8)),...
    str2num(times(:,9:10)),str2num(times(:,11:12)),str2num(times(:,13:14)));

Wa=WC(:,3);
Wb=WC(:,4);
W=WC(:,5);

cd ..; cd ..;

save W_data3 Wa Wb 

% 
% dte_W1=dte_W;
% W1=W;
% Wa1=Wa;
% Wb1=Wb;
% 
% load('KnorrWhitecaps.mat')
% W=Knorr_W.raw.W;
% Wa=Knorr_W.raw.Wa;
% Wb=Knorr_W.raw.Wb;
% dte_W=Knorr_W.raw.dte;
% 
% 
% dte_total=[dte_W;dte_W1];
% W_total=[W;W1];
% Wa_total=[Wa;Wa1];
% Wb_total=[Wb;Wb1];
% % small error fix:
% dte_total(103943,1)=mean([dte_total(103943-1) dte_total(103943+1)]);
% 
% Knorr_W.raw.W=W_total;
% Knorr_W.raw.Wa=Wa_total;
% Knorr_W.raw.Wb=Wb_total;
% Knorr_W.raw.dte=dte_total;
% 
% 
% save KnorrWhitecaps Knorr_W