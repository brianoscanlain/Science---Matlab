function [X]=DaFixer(Xseries,Yseries,type,limit,interval)

%% function [X]=DaFixer(Xseries,Yseries,type,limit,interval)
% The function looks at an interval of a timerseries of data and corrects
% any spikes depending on the selected 'type' and 'limit' using a linear
% fitting.
%
% 'interval' = [ min(Xdata) max(Xdata)] that you want to analyse. It's
%              optional to use. The two values are in units of Y rather than
%              indexes.
%
% 'type' is a string input, and is either:
%            'less' = if the Ydata is less than the specified 'limit'
%            'greater'= if the Ydata is greater than the specified 'limit'
%
% limit is a numerical threshold to distinguish between bad and good data.
% The User should plot(Xseries,Yseries) and observe any data that looks 
% suspicious. The user can then tune 'limit', 'interval' and 'type' to 
% correct the relevant data.
%
% The bad sectors are then filled using a linear fit. The function cannot 
% fix the first and last values of the Yseries data.
%
% Brian Scanlon, NUIGalway March,2013


%% Nargin

if nargin==4
    interval=[min(Xseries) max(Xseries)];

elseif nargin<4
    error('DaFixer::Not enough values inputted');
end

if ~(length(Xseries)==length(Yseries))
    error('DaFixer::X and Y series dimensions must agree!!');
end

%% Identifying bad sectors


if length(type)==4     %case type= 'less'
    BadInd=find(Yseries<limit);
else
    BadInd=find(Yseries>limit);
end
%apply the interval
X1=(find(Xseries>=interval(1)));
X2=(find(Xseries<interval(2)));
BadInd=BadInd(BadInd>=X1(1));
BadInd=BadInd(BadInd<X2(end));
clear X1 X2;

%% Failsafe incase Ind includes the 1st or last values of Yseries
if ~isempty((BadInd==1)) 
BadInd(BadInd==1)=[];
end
if  ~isempty(BadInd==length(Yseries))
    BadInd(BadInd==length(Yseries))=[];
end


%% Seperate data into individual groups if they are alongside eachother
count=1;
n=1;
Len=1;
%Allocate a maximum possible memory for BadSectors:
BadSectors=zeros(length(BadInd),2); % will be a 2 column matrix showing the
                                    %start and length of aligned bad sectors

while count<length(BadInd)
    
    if Len==1   %true for first Ind of each group
    BadSectors(n,1)=BadInd(count);
    end
    
    if BadInd(count+1)-BadInd(count)==1 %if the current and next Ind are 
       Len=Len+1;                       %together
    else                           %This marks the end of a group and so
        BadSectors(n,2)=Len;        %the data is written
        Len=1;
        n=n+1;
    end
    
    if count==length(BadInd)-1     %This is here for the final Ind value.
       if Len==1                   %We cannot check ind+1 so we look at 
                                   %Ind-1.
           BadSectors(n,1)=BadInd(count+1);    %Here last Ind is a single
       end
           BadSectors(n,2)=Len;                   
      
    end  
    
    count=count+1;
end

%shorten BadSectors (Shorten the memory allocation)

if ~isempty(BadSectors)
while (BadSectors(end,1)==0) && (BadSectors(end,2)==0)
    BadSectors(end,:)=[];
end
end

%% Correct the Yseries Data using a linear fit

for i=1:size(BadSectors,1)
    
    Y1=(BadSectors(i,1)-1);    %take the Ind before the 'ith' group
    Y2=Y1+BadSectors(i,2)+1;   %take the Ind after the 'ith' group
    
    Yseries(Y1:Y2)= linspace(Yseries(Y1),Yseries(Y2),BadSectors(i,2)+2);
    
   
end
    
 X=Yseries;   

