function [L,uStar,class]=OFparams(u,T,T0,z0)

 eps=0.622;
 g=9.81;
 p=101300;
 RH=0.8;
 kappa=0.4;

 es=@(t) (1.0007+3.46*1e-06*p)*611.21*exp(17.502*t/(240.97+t)); %611.21*exp((18.678-t/234.5)*t/(257.14+t)); %Arden Buck equation
 rv=@(t) eps*RH*es(t)/(p-RH*es(t));
 Tv=@(t) (t+273.15)*(1+rv(t)/eps)/(1+rv(t));

 theta = Tv(T)*(p/(p-g*1.225*114))^(287/1004);
 theta0 = Tv(T0)*(p/(p-g*1.225*8))^(287/1004);

 uStar = u/log(69/z0);
 thetaStar = (theta-theta0)/log(114/8);
 L = uStar^2*theta0/(g*thetaStar);

 x0 = [uStar,L,thetaStar]'; %,z0]';
 x=picardOF(x0,u,theta,theta0,z0);

 uStar = x(1)*kappa;
 L = x(2);
 thetaStar = x(3)*kappa;

 if (abs(L) > 500)
   class = 3;
 elseif (L > 200)
   class = 2;
 elseif (L > 0)
   class = 1;
 elseif (L > -100)
   class = 6;
 elseif (L > -200)
   class = 5;
 elseif (L > -500)
   class = 4;
 end

end
