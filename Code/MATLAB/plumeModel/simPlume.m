function simPlume(fileName)

global axisGrid parameterization grid Cfield

disp(sprintf('\nInitialize model'));
modelInput=initialize({fileName});

Cfield.region1='';

disp('Calculate dose rate');
names=modelInput.dose.monitoringPointNames;
monitoringPoints=modelInput.dose.monitoringPoints;
noMonitoringPoints=modelInput.dose.noMonitoringPoints;
h=modelInput.conc.stackHeight;
href=modelInput.meteo.measurementHeight;
constWindSpeed=modelInput.meteo.constWindSpeed;
constWindDirection=modelInput.meteo.constWindDirection;
outputDir=modelInput.outputDir;
noCases=modelInput.meteo.noCases;
for c=1:noCases
 fileName=getfield(modelInput.meteo,strcat('case',num2str(c)),'file');
 disp(sprintf('-> process file %s',fileName));
 Uref=getfield(modelInput.meteo,strcat('case',num2str(c)),'refWindSpeed');
 T114=getfield(modelInput.meteo,strcat('case',num2str(c)),'refTemp114m');
 T8=getfield(modelInput.meteo,strcat('case',num2str(c)),'refTemp8m');
 theta=getfield(modelInput.meteo,strcat('case',num2str(c)),'windDirection');
 noTimes=length(theta);
 if constWindSpeed && constWindDirection
  noTimes=1;
  Uref=mean(Uref);
  theta=mean(theta);
  disp('   Constant wind speed is assumed');
  disp('   Constant wind direction is assumed');
 elseif constWindSpeed
  Uref=mean(Uref)*ones(noTimes,1);
  disp('   Constant wind speed is assumed');
 elseif constWindDirection
  theta=mean(theta)*ones(noTimes,1);
  disp('   Constant wind direction is assumed');
 end
 dose=zeros(noMonitoringPoints,noTimes);
 for t=1:noTimes
  disp(sprintf('   time %d',t));
  if strcmp(parameterization,'PG') || strcmp(parameterization,'Pasq')
    kappa=0.4; z0=0.03;
    [L,uStar,class]=OFparams(Uref(t),T114(t),T8(t),z0);
    U=uStar/kappa*(log(h/z0)-kernel(h/L,'M','integral')+kernel(z0/L,'M','integral'));
  else
    class = stab(T114(t),T8(t),Uref(t));
    U=windSpeed(h,parameterization,Uref(t),href,class);
  end
  heff=plumeRise(T114(t),T8(t),U,class,h);
  C=conc(modelInput.conc,class,heff);
  grid=mesh('full',axisGrid);
  N=grid.noRegions;
  Cfield=rmfield(Cfield,'region1');
  Cfield.region1='';
  for id=1:N
	if (isfield(Cfield,strcat('region',num2str(id))) &&  id>1)
	  Cfield=rmfield(Cfield,strcat('region',num2str(id)));
	end
  	ir=getfield(grid,strcat('region',num2str(id)),'range');
  	Ci=C(ir(1,1):ir(1,2),ir(2,1):ir(2,2),ir(3,1):ir(3,2));
  	Ci=permute(Ci,[2 1 3]);
  	Ci=reshape(Ci,[],1);
  	Cfield=setfield(Cfield,strcat('region',num2str(id)),Ci);
  end
  clear C Ci ir id;
  [modelInput.dose.monitoringPoints,modelInput.dose.neighborhoodRadii]=rotation(theta(t),monitoringPoints);
  dose(:,t)=calcDose(modelInput.dose);
  dose(:,t)=dose(:,t)/U;
 end
 writeToFile(dose,outputDir,fileName,names,'%6.3f');
end
end
%clear all;

