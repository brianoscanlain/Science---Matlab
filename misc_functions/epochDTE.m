function [ dNum,dStr ] = epochDTE( X )
%Converts Unix time [seconds] in to matlab datenum
%
% Brian Scanlon

if isstring(X)
   X=hex2dec(X); 
end

dNum=datenum(X/86400 + datenum(1970,1,1,0,0,0));
dStr=datestr(dNum);

end

