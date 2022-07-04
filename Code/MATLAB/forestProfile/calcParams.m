function params=calcParams(params)

 %initialize the required parameters
 LAI=params.LAI;
 Cd=params.Cd;
 extCoeff=params.extCoeff;
 h=params.h;
 canopyDepth=params.canopyDepth;
 standDensity=params.standDensity;
 DBH=params.DBH;
 zref=params.zref;
 uref=params.uref;

 kappa=0.4;
 lambda=1.5;

 %correction factor for the model of Mihailovic
 alpha=2*(Cd*LAI)^0.25; plantCover=1-exp(-LAI*extCoeff); corrFactor=plantCover*(alpha-1)+1;

 %roughness length
 if isfield(params,'z0')
   z0=params.z0;
   disp(sprintf('value z0 imposed; z0 = %3.2f\n', z0));
 else
   z0=0.071*h;
   params.z0=z0;
 end

 %thickness roughness sublayer
 if isfield(params,'zStar')
   zStar=params.zStar;
   disp(sprintf('value zStar imposed; zStar = %5.2f\n', zStar));
 else
   zStar=2*h;
   params.zStar=zStar;
 end

 %lad model
 a2=4.1920; a3=1.9264; a4=1; a1=LAI*(a2+1)/(hypergeom([a2+1,-a3],a2+2,1/a4)*canopyDepth);
 lad=@(y) a1*y.^a2*(a4-y).^a3;

 disp(sprintf('\nparameters LAD model: a1=%5.4f, a2=%5.4f, a3=%5.4f, a4=%5.4f\n',a1,a2,a3,a4));

 %trunk area density
 TAD=standDensity*0.5*pi*DBH*10^(-4);

 %determine the height zi at which lad equals tad 
 f=@(x) (lad(x)-TAD).^2;
 yi=fsolve(f,0.25);
 zi=(yi-1)*canopyDepth+h;
 params.zi=zi;

 %displacement height
 derivLad=@(y) a1*(a2*y.^(a2-1).*(a4-y).^a3-a3*y.^a2.*(a4-y).^(a3-1))/canopyDepth;
 d=fsolve(derivLad,0.5); d=(d-1)*canopyDepth+h;
 params.d=d;
 d=plantCover*alpha*d/(plantCover*(alpha-1)+1);

 %friction velocity and empirical parameter mu from De Ridder (2009)
 uh=1/(kappa*corrFactor)*log((corrFactor*h-plantCover*alpha*params.d)/(alpha^2*z0));
 f=@(x) [x(1)/kappa*(log((zref-d)/z0)+(log(1+lambda*(zStar-d)./(x(2)*(zref-d))).*exp(-x(2)*(zref-d)/(zStar-d))-log(1+lambda*(zStar-d)./(x(2)*z0)).*exp(-x(2)*z0/(zStar-d)))/lambda)-uref,x(1)/kappa*(log((h-d)/z0)+(log(1+lambda*(zStar-d)./(x(2)*(h-d))).*exp(-x(2)*(h-d)/(zStar-d))-log(1+lambda*(zStar-d)./(x(2)*z0)).*exp(-x(2)*z0/(zStar-d)))/lambda)-x(1)*uh];
 par=fsolve(f,[0.3,2.59]); %,optimoptions('fsolve','Algorithm','levenberg-marquardt'));
 uStar=par(1); mu=par(2);

 params.uStar=par(1); %gradU=uStar/(kappa*(800-d)), Utop=uStar/kappa*log((800-d)/z0), ustar=kappa*Utop/log(800/0.03)
 params.mu=par(2);  %TKE=ustar^2/sqrt(0.03), TDR=ustar^3/(kappa*69)

 %determine the height zc at which lad equals u'/(u*plantCover)
 duh=1/(kappa*(h-d))*(1-exp(-mu*(h-d)/(zStar-d)));
 ladTop=duh/(uh*plantCover);
 f=@(x) (lad(x)-ladTop).^2;
 yc=fsolve(f,0.75); %,optimoptions('fsolve','Algorithm','trust-region-reflective'));
 zc=(yc-1)*canopyDepth+h;

 params.zc=zc;

 %scale factor for the hypergeometric approximation
 hypergeomyc = hypergeom([a2+1,-a3],a2+2,yc/a4); 
 hypApprox=@(z) (1-z).^2+2/(a2+2).*(1-z).*z+2/((a2+3)*(a2+2)).*z.^2;
 maxIntegral = a1*a4^a3*(yc.^(a2+1).*hypergeomyc-yi.^(a2+1).*hypergeom([a2+1,-a3],a2+2,yi/a4))/(a2+1)*canopyDepth;
 maxIntegralApprox=a1*a4^a3*(yc.^(a2+1).*hypApprox(yc/a4)-yi.^(a2+1).*hypApprox(yi/a4))/(a2+1)*canopyDepth;
 scale = maxIntegral/maxIntegralApprox;

 params.scale=scale;

% ue=uStar/kappa*(log((ze-d)/z0)+(log(1+lambda*zStar./(mu*(ze-d))).*exp(-mu*(ze-d)/zStar)-log(1+lambda*zStar./(mu*z0)).*exp(-mu*z0/zStar))/lambda);
% duh=uStar/(kappa*(h-d))*(1-exp(-mu*(h-d)/zStar));
end
