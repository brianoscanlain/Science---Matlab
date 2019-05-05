function [ax,hl1,hl2,hl3] = plotyyy(x1,y1,x2,y2,x3,y3,ylabels)
%PLOTYYY - extends plotyy to include a third y-axis
%
%Syntax: [ax,hl1,hl2,hl3] = plotyyy(x1,y1,x2,y2,x3,y3,ylabels)
%
%Input: x1,y1 are the xdata and ydata for the first axes' line
% x2,y2 are the xdata and ydata for the second axes' line
% x3,y3 are the xdata and ydata for the third axes' line
% ylabels is a 3x1 cell array containing the ylabel strings
%
%Output: ax is a 3x1 double array containing the axes' handles
% hl1 is the handle for line 1
% hl2 is the handle for line 2
% hl3 is the handle for line 3
%
%Example:
%x=0:10;
%y1=x; y2=x.^2; y3=x.^3;
%ylabels{1}='plain X';
%ylabels{2}='X-squared';
%ylabels{3}='X-cubed';
%[ax,hl1,hl2,hl3] = plotyyy(x,y1,x,y2,x,y3,ylabels);
%
%m-files required: plotyy (in Matlab 5.0 and above)

%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute
%Dept. of Fisheries and Oceans Canada
%email: nospam@thank.you
%Web: http://www.qc.dfo-mpo.gc.ca/iml/
%April 2000; Last revision: 07-Apr-2000

if nargin==6
   %Use empty strings for the ylabels
   ylabels{1}=' '; ylabels{2}=' '; ylabels{3}=' ';
elseif nargin > 7
   error('Too many input arguments')
elseif nargin < 6
   error('Not enough input arguments')
end

figure; clf;
set(gcf,'units','normalized')
[ax,hl1,hl2] = plotyy(x1,y1,x2,y2);
cfig=get(gcf,'color');
pos=[0.1 0.1 0.7 0.8];
offset = pos(3)/5.5;
%Reduce width of the two axes generated by plotyy
pos(3) = pos(3) - offset/2;
set(ax,'position',pos);

%Set the position of third axes
pos3=[pos(1) pos(2) pos(3)+offset pos(4)];

limx1=get(ax(1),'xlim');
xticks=get(ax(1),'xtick');
delta_xtick = xticks(end) - xticks(end-1);
limx3=[limx1(1) limx1(2) + delta_xtick];

ax(3)=axes('Position',pos3,...
   'Color','none','XColor','k','YColor','r',...
   'xtick',[],'xlim',limx3,'yaxislocation','right',...
   'XMinorTick','on','YMinorTick','on');

hl3=line(x3,y3,'Color','r','Parent',ax(3));
limy3=get(ax(3),'YLim');

%Hide unwanted portion of the x-axis line that lies
%between the end of the second and third axes
hl4=line([limx1(2) limx3(2)],[limy3(1) limy3(1)],...
   'Color',cfig,'Parent',ax(3),'Clipping','off');

%label all three y-axes
set(get(ax(1),'ylabel'),'string',ylabels{1})
set(get(ax(2),'ylabel'),'string',ylabels{2})
set(get(ax(3),'ylabel'),'string',ylabels{3}) 