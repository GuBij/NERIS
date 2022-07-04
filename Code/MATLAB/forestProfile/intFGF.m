function L=intFGF(z,option)

 load('treeParams.mat');

 n=length(z);
 L=zeros(n,1);

 if nargin==2 & strcmp(option,'derive')
  for j=1:n
   if z(j) < h & z(j) >= zc
     L(j)=duh/(eta*uh);
   elseif z(j) < h & z(j) >= zi
     y=(z(j)-h)/cdepth+1;
     L(j)=a1*y^a2*(a4-y)^a3;
   elseif z(j) < h & z(j) >= 10*z0g
     L(j)=TAD;
   elseif z(j) < h
     L(j)=TAD*10*z0g/z(j);
   end
  end
 else
  yc=(zc-h)/cdepth+1;
  yi=(zi-h)/cdepth+1;
  F21c=yc^(a2+1)*a4^a3*hypergeom([a2+1,-a3],a2+2,yc/a4);
  F21i=yi^(a2+1)*a4^a3*hypergeom([a2+1,-a3],a2+2,yi/a4);

  for j = 1:n
   if z(j) < h & z(j) >= zc
     L(j)=(h-z(j))*duh/(eta*uh);
   elseif z(j) < h & z(j) >= zi
     y=(z(j)-h)/cdepth+1;
     L(j)=(h-zc)*duh/(eta*uh)+a1*cdepth/(a2+1)*(F21c-y^(a2+1)*a4^a3*hypergeom([a2+1,-a3],a2+2,y/a4));
   elseif z(j) < h & z(j) >= 10*z0g
     L(j)=(h-zc)*duh/(eta*uh)+a1*cdepth/(a2+1)*(F21c-F21i);
     L(j)=L(j)+TAD*(zi-z(j));
   elseif z(j) < h
     L(j)=(h-zc)*duh/(eta*uh)+a1*cdepth/(a2+1)*(F21c-F21i)+TAD*(zi-10*z0g);
     L(j)=L(j)+10*z0g*TAD*log(10*z0g/z(j));
   end
  end
 end
