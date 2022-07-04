function tau=stressProfile(z,intFGF,option)

 load('treeParams.mat');
 load('treeParams.mat','mu');

 n=length(z);
 U=zeros(n,1);

 if nargin==3 & strcmp(option,'Wang')
  N=0.9; a1=0.7680; a2=0.9031; b0=0.4; b1=b0-kappa/log(h/z0g); b2=-4;
  a0=LAI/h; %mean(intFGF(z>=10*z0g & z < h)); %intFGF=intFGF(z,'derive');
  beta=uStar/uh;
  lc=2*beta^3/(Cd*a0);
  sh=lc/(lc^N+(kappa*h)^N)^(1/N);
  Lsbyh=kappa*sh/beta; Ls0byh=log(h/z0g);
  X=1/Lsbyh-1/Ls0byh*(Lsbyh/Ls0byh)^1.45;
  A=a2*X^2+a1*X;
%  FAI=a0*h;
%  fn=7.5*(0.25*FAI)^0.63*(1-0.25*FAI)^12;
%  Cu=1/beta; %(b0-b1*exp(b2*FAI));
%  A=Cd*Cu*FAI/((1-fn)^2*kappa*sh);
  gh=2*sqrt(A); g0=2*sqrt(A*z0g/h);
  Up=0; %(fn*uStar/h)/(Cd*Cu*a0);
  C1=uh/(besseli(0,gh)-besseli(0,g0)*besselk(0,gh)/besselk(0,g0));
%  C1=(2*uStar/(kappa*sh*gh))/(besseli(1,gh)+besseli(0,g0)*besselk(1,gh)/besselk(0,g0));
  C2=-(C1*besseli(0,g0)+Up)/besselk(0,g0);
 end

 for j=1:n
   if z(j) >= h
     tau(j)=uStar^2; %/kappa*(log((z(j)-d)/z0)+2/3*log(1+1.5*zStar/(mu*(z(j)-d)))*exp(-mu*(z(j)-d)/zStar)-2/3*log(1+1.5*zStar/(mu*z0))*exp(-mu*z0/zStar));
   elseif nargin==3 & strcmp(option,'Wang')
     g=2*sqrt(A*z(j)/h);
     tau(j)=0.5*g/z(j)*(C1*besseli(1,g)-C2*besselk(1,g))*kappa*z(j)*sh*uStar; %+fn*uStar^2*z(j)/h;
   else
     tau(j)=uStar^2*exp(-eta*intFGF(j));
   end
 end

end
