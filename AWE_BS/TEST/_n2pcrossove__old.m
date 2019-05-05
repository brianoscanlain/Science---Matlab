function [N2P,P2N,N2P2, finalP2N] = n2pcrossover(LocalMax,threshold,PIP_sdx,PIP_s2dx,PIP_s3dx,PIP_s4dx,peak)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


          
%flip data and find index values up to the first derivative max:
if  length(LocalMax(:,1)) > 1
    indexNP = flipud(find(threshold < LocalMax(end,1)));
else
    indexNP = flipud(find(threshold < LocalMax(1,1)));
end
           
%Find the negative to positive crossing of the 1st derivative
N2P = [];

if ~isempty(indexNP)
   indexNP2 = indexNP((PIP_sdx(indexNP) < 0));     
   if ~isempty(indexNP2)
       if min(indexNP2)==min(indexNP)
           %its goin to be the 
           N2P=indexNP2(abs(diff(indexNP2))>1);
           N2P=max(threshold(N2P)>peak(1,1)); 
       else 
           %its going to be the max threshold(indexnp2)
           N2P=max(threshold(indexNP2));
       end

%       if res==1
%          if ~(threshold(N2P)>peak(2)-5 && threshold(N2P)<peak(2)+5)
%             N2P= max(indexNP2(  abs(diff(threshold(indexNP2)))>1  )>peak(2)-5);
%          end
%       end
    end
end
                

%Now find the first positive to negative crossing of the second
%derivative after the negative to positive crossing::
P2N=[];
if  ~isempty(N2P)
    indexPN = flipud(find(threshold > (N2P)));
    indexPN2= indexPN((PIP_s2dx(indexPN)>0)); 
    if ~isempty(indexPN2)
       P2N=indexPN2(1);
    end
end
                
                
%Now we find the negative to positive crossing of the 3rd derivative after
%P2N
N2P2=[];
if ~isempty(P2N)
          indexPN3= indexPN2(PIP_s3dx(indexPN2) > 0);
          if ~isempty(indexPN3)
             N2P2=indexPN3(1);
          end
%           processCode=5;      
end

finalP2N=[];
if ~isempty(N2P2)
          indexPN4=indexPN3(PIP_s4dx(indexPN3) > 0);
          if ~isempty(indexPN4)
          finalP2N=indexPN4(1);
          end
%           finalThresh(1,1)=threshold(finalP2N);
%           finalThresh(1,2) = sum(sum(editImage>finalThresh(1,1)))*100./(size(editImage,1)*size(editImage,2));
end

end

