function [fcob,Fb,Fc] = cmorInfo(wname,precision)
%cmorInfo calculates the centre frequency of the complex Morlet wavelet.
%Precision = 10; %default

if strcmp(wname(1:4),'cmor')
    [psi,xval,Fb,Fc] = wavefunBS(wname,precision);
    Timerange=max(xval)-min(xval);
    sp=abs(fft(detrend(psi,0)));
    [~,ind]=max(sp);
    fcob= 1 / (Timerange/(ind-1) );
else
    error('This is not a morlet wavelet, please use matlab function if you want to try other wavelets');
    fcob=[];
    Fc=[];
    Fb=[];
end

end



%--------------------------------------------------------------------------
%----------------------------intwaveBS-------------------------------------
function [out1,xval] = intwaveBS(wname,in2)
precision = in2; %keep this between 1 and 20. matlab default=8
%Make sure it is cmor (complex morlet wavelet without scale function):
if strcmp(wname(1:4),'cmor')
    [psi,xval] = wavefunBS(wname,precision);
end

step = xval(2)-xval(1);
out1 = cumsum(psi)*step;    %function integrated value from n=0...ith value
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------





%--------------------------------------------------------------------------
%---------------------------wavefunBS--------------------------------------
function [out1,out2,Fb,Fc] = wavefunBS(wname,precision)
[~,~,~,bounds] =  ...
    wavemngr('fields',wname,'type','file','fn','bounds');
bounds=[-8 8];
iter=precision;  %=10     %2^precision are the number of points in mother wavelet
np = 2^iter;     %=1024 for iter=10
lb = bounds(1); %=-8
ub = bounds(2); %=+8
[out1,out2,Fb,Fc] = cmorwavfBS(lb,ub,np,wname);  %fname = cmorwavf
                                            %wname = cmor5-3
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end






%--------------------------------------------------------------------------
%-------------------------cmorwavfBS---------------------------------------
function [psi,X,Fb,Fc] = cmorwavfBS(LB,UB,N,Fb,Fc)
%LB = lower bound
%UB = upper bound
%N = Number of points (a.k.a 2^precision)
%Fb = function coefficient b, normalized frequency \omega_0 
%Fc = function coefficient c,                             
%                       ###(cmorFb-Fc)### 
%
%PSI(X) = ((pi*FB)^(-0.5))*exp(2*i*pi*FC*X)*exp(-(X^2)/FB)
nbIn = nargin;
switch nbIn
    case {0,1,2} ,
        error(message('Wavelet:FunctionInput:NotEnough_ArgNum'));
    case 3 ,
        Fc = 1; Fb = 1;
    case 4 ,
        if ischar(Fb)
            label = deblank(Fb);
            ind   = strncmpi('cmor',label,4);
            if isequal(ind,1)
                label(1:4) = [];
                len = length(label);
                if len>0
                    ind = strfind(label,'-');
                    if isempty(ind)
                        Fb = []; % error 
                    else
                        Fb = wstr2num(label(1:ind-1));
                        label(1:ind) = [];
                        Fc = wstr2num(label);    
                    end
                end
            else
                Fc = []; % error 
            end
        else
            Fb = []; Fc = []; % error 
        end
        
    case 5 ,
end       %at this stage we should be ready to calculate the mother wavelet
X = linspace(LB,UB,N);  % wavelet support.
psi = ((pi*Fb)^(-0.5))*exp(2*1i*pi*Fc*X).*exp(-(X.*X)/Fb);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



