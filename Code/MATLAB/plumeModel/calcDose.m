function dose=calcDose(modelInput)

  global grid Cfield;

  K=1.6*10^(-13+9); %conversion factor to nSv/h
  nuclide=modelInput.nuclide;
  locations=modelInput.monitoringPoints;
  L=modelInput.neighborhoodRadii;
  noMonitoringPoints=modelInput.noMonitoringPoints;
  N=grid.noRegions;

  dose=zeros(noMonitoringPoints,1);

  if strcmp(nuclide,'Ar-41')
    E=1.294;
    logE=log(E);
    A=59.4;
  else
    error('Nuclide %s is not implemented.',nuclide);
  end

  a1=sum([-8.46767*10^(-5),2.71806*10^(-4),2.93831*10^(-3),-6.30793*10^(-3),-1.84327*10^(-2),4.78301*10^(-2),-4.61222*10^(-2)].*(logE.^(6-(0:6))));
  a2=sum([-7.15631*10^(-4),-2.54597*10^(-3),9.31948*10^(-3),1.89129*10^(-2),-3.12917*10^(-2),1.49122*10^(-2),-1.02463*10^(-2)].*(logE.^(6-(0:6))));
  mu=exp(sum([-1.40465*10^(-4),3.07113*10^(-4),1.42298*10^(-2),-3.57795*10^(-3),-1.18921*10^(-1),-4.20821*10^(-1),-2.70365].*(logE.^(6-(0:6)))))/10*1.2041;
  A1=(a1+1)*mu;
  A2=(a2+1)*mu;
  if E<0.1
    sigair=exp(sum([-9.27823*10^(-3),2.62726*10^(-1),2.92391,1.58909*10,3.89774*10,3.06711*10].*(logE.^(5-(0:5)))))/10;
  else
    sigair=exp(sum([1.06134*10^(-5),-1.59201*10^(-3),7.51537*10^(-3),2.35368*10^(-2),-1.19158*10^(-1),-1.82464*10^(-1),-3.57954].*(logE.^(6-(0:6)))))/10;
  end

  for j=1:noMonitoringPoints
    CV=0; V=0;
    for id=1:N
      noCells=getfield(grid,strcat('region',num2str(id)),'size');
      cellCenters=getfield(grid,strcat('region',num2str(id)),'C');
      cellVolume=getfield(grid,strcat('region',num2str(id)),'V');
      C=getfield(Cfield,strcat('region',num2str(id)));
      for i=1:noCells
        r1=norm(cellCenters(i,:)-locations(j,:));
        r2=sqrt(sum((cellCenters(i,[1,3])-locations(j,[1,3])).^2)+(cellCenters(i,2)+locations(j,2))^2);
        CVi=C(i)*cellVolume;
        if r1<L(j) && r2<L(j)
	  CV=CV+2*CVi;
	  V=V+2*cellVolume;
        elseif r1<L(j) 
          CV=CV+CVi;
          V=V+cellVolume;
          Buildup=A*exp(-a1*mu*r2)+(1-A)*exp(-a2*mu*r2);
          dose(j)=dose(j)+K*sigair*E*CVi*Buildup*exp(-mu*r2)/(4*pi*r2^2);
        elseif r2<L(j)
          CV=CV+CVi;
          V=V+cellVolume;
          Buildup=A*exp(-a1*mu*r1)+(1-A)*exp(-a2*mu*r1);
          dose(j)=dose(j)+K*sigair*E*CVi*Buildup*exp(-mu*r1)/(4*pi*r1^2);
        else
	  Buildup1=A*exp(-a1*mu*r1)+(1-A)*exp(-a2*mu*r1);
	  Buildup2=A*exp(-a1*mu*r2)+(1-A)*exp(-a2*mu*r2);
	  dose(j)=dose(j)+K*sigair*E*CVi*(Buildup1*exp(-mu*r1)/(4*pi*r1^2)+Buildup2*exp(-mu*r2)/(4*pi*r2^2));
        end
      end
    end
    if V>0
     R=(3*V/(2*pi))^(1/3);
     dose(j)=dose(j)+0.5*K*sigair*E*CV/V*(-A*exp(-A1*R)/A1-(1-A)*exp(-A2*R)/A2+A/A1+(1-A)/A2);
    end
  end

end

