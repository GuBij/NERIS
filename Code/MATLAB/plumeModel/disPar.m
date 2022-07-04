function [sigy,sigz]=disPar(x,option,class)

 if strcmp(option,'BM')
%   A=0.418; a=0.796; B=0.520; b=0.711; %E3 (neutral)
%   A=0.586; a=0.796; B=0.700; b=0.711; %E4 (slightly unstable)
   switch class
	case 1
	  A=0.235; a=0.796; B=0.311; b=0.711;
	case 2
	  A=0.297; a=0.796; B=0.382; b=0.711;
	case 3
	  A=0.418; a=0.796; B=0.520; b=0.711;
	case 4
	  A=0.586; a=0.796; B=0.700; b=0.711;
	case 5
	  A=0.826; a=0.796; B=0.950; b=0.711;
	case 6
	  A=0.946; a=0.796; B=1.321; b=0.711;
	case 7
	  A=1.043; a=0.698; B=0.819; b=0.669;
   end
   sigy=A*x.^a; sigz=B*x.^b;
 elseif strcmp(option,'PG')
   switch class
	case 1
	  a=0.037; b=1170; c=0.134; d=0.022; e=0.7;
	case 2
	  a=0.0566; b=1070; c=0.137; d=0.0335; e=0.624;
	case 3
   	  a=0.0787; b=707; c=0.135; d=0.0475; e=0.465;
	case 4
	  a=0.134; b=283; c=0.134; d=0.0722; e=0.102;
	case 5
	  a=0.202; b=370; c=0.162; d=0.0962; e=-0.101;
	case 6
	  a=0.250; b=927; c=0.189; d=0.1020; e=-1.918;
   end
   sigy=a*x./(1+(x/b)).^c; sigz=d*x./(1+(x/b)).^e;
 elseif strcmp(option,'Brook')
   A=0.32; a=0.78; B=0.22; b=0.78;
   sigy=A*x.^a; sigz=B*x.^b;
 elseif strcmp(option,'Pasq')
   T=8.3333-0.72382*log(x/1000);
   sigy=1000*tan(T*pi/180)/2.15;
   n=length(x);
   sigz=zeros(n,1);
   for i=1:n
     if x(i)<300
	a=34.459; b=0.86974;
     elseif x(i) < 1000
	a=32.093; b=0.81066;
     elseif x(i) < 3000
	a=32.093; b=0.64403;
     end
   sigz(i)=a*(x(i)/1000)^b;
   end
 end

end
