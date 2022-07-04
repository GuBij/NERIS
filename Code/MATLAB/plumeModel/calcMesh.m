 function x=calcMesh(n,N,cellSize)
   noCells=sum(n)+1;
   x=zeros(noCells,1);
   dx=zeros(N+1,1);
   j1=1;
   for i=1:(N+1)
     j0=j1+1;
     j1=j1+n(i);
     x(j0:j1)=x(j0-1)+cellSize(i)*(1:n(i));
   end
   x=0.5*(x(2:noCells)+x(1:(noCells-1)));
 end

