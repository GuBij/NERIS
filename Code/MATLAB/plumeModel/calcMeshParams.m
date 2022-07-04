 function [cellSize,n,N]=calcMeshParams(l,l0,alpha)
   bool=1;
   N=length(alpha)-1;
   while bool
     r=(l/l0)^(1/N);
     n=ceil(alpha./((r.^(0:N))'*l0));
     if sum(n<2)>0
       N=N-1; bool=1;
     else
       bool=0;
       cellSize=zeros(N+1,1);
       cellSize(1)=alpha(1)/n(1);
       for i=2:(N+1)
        cellSize(i)=alpha(i)/n(i);
       end
%       alpha=n.*r.^(0:N)*l0/L;
     end
    end
 end

