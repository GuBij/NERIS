function heff=plumeRise(T114,T8,uref,class,stackHeight)

 gradT=(T114-T8)/106;
 g=9.81;
 Ts=25+273.15;
 d=2.8;
 v=9.44/(pi*(0.5*d)^2);
 Tair=T8+gradT*(stackHeight-8)+273.15;
 F=0.25*g*v*d^2*(Ts - Tair)/Ts;
 s=g*(gradT+0.0098)/Tair;

 xf=0;
 if (F<55)
  xf=49*abs(F)^0.625;
 else
  xf=119*F^0.4;
 end

 if (class<3)
  if (xf > (1.84*uref/s^0.5))
    heff=2.4*(F/(uref*s))^(1/3);
  else
    heff=1.6*abs(F)/F*abs(F)^(1/3)*xf^(2/3)/uref;
  end
 else
   heff=1.6*abs(F)/F*abs(F)^(1/3)*xf^(2/3)/uref;
 end

 heff = heff + stackHeight;

end

