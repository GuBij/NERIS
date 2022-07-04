function modelInput=initialize(meteoFiles)

  global axisGrid parameterization;

  params;

  modelInput.meteo.measurementHeight=href;
  modelInput.meteo.constWindSpeed=constWindSpeed;
  modelInput.meteo.constWindDirection=constWindDirection;
  modelInput.conc.stackHeight=h;
  modelInput.conc.parameterization=parameterization;
  modelInput.dose.nuclide=nuclide;
  modelInput.conc.sourceStrength=Q;
  modelInput.outputDir=outputDir;

  noMeteoCases=length(meteoFiles);
  modelInput.meteo.noCases=noMeteoCases;
  for f=1:noMeteoCases
    fid=fopen(strcat(inputDir,'/meteo',meteoFiles{f},'.txt'),'r');
    if fid == -1
     error('Meteo file %s in directory %s not found.',meteoFiles{f},inputDir);
    end
    meteoFile=textscan(fid,'%f %f %f %f %f %f %f','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
    modelInput=setfield(modelInput,'meteo',strcat('case',num2str(f)),'file',meteoFiles{f});
    modelInput=setfield(modelInput,'meteo',strcat('case',num2str(f)),'refTemp114m',meteoFile{2});
    modelInput=setfield(modelInput,'meteo',strcat('case',num2str(f)),'refTemp8m',meteoFile{1}); 
    modelInput=setfield(modelInput,'meteo',strcat('case',num2str(f)),'refWindSpeed',meteoFile{3});
    modelInput=setfield(modelInput,'meteo',strcat('case',num2str(f)),'windDirection',meteoFile{4});
    fclose(fid);
  end

  if sum(alphax) ~= Lx
    error('Region division in the x-direction is not accurate; reconsider alphax');
  end
  if sum(alphay) ~= Ly
    error('Region division in the y-direction is not accurate; reconsider alphay');
  end

%{
  if sum(alphazg) ~= h
    error('Region division below the source is not accurate; reconsider alphazg');
  end
  if sum(alphazt) ~= Lz-h
    error('Region division above the source is not accurate; reconsider alphazt');
  end
%}

  gridInput.x.minSize=lx0;
  gridInput.x.maxSize=lx;
  gridInput.x.regionLengths=alphax;
  gridInput.y.minSize=ly0;
  gridInput.y.maxSize=ly;
  gridInput.y.regionLengths=alphay;
  gridInput.z.minSize=lz0;
  gridInput.z.belowSrc.maxSize=lzg;
  gridInput.z.belowSrc.regionLengths=alphazg;
  gridInput.z.aboveSrc.maxSize=lzt;
  gridInput.z.aboveSrc.regionLengths=alphazt;

  axisGrid=mesh('decompose',gridInput,h);

  if strcmp(monitoringPoints,'BR1')
    locationNames={'IMR/M11','IMR/M12','IMR/M07','IMR/M08','IMR/M13','IMR/M09','IMR/M10'};
    locations=zeros(7,3);
    locations(:,1)=[200.143503293714;200.155246861131;199.982903434886;199.855503938484;199.869570903673;199.904281147205;199.961095905096];
    locations(:,2)=[212.237160795659;211.962463226418;211.936320720149;211.93776306347;212.167760330885;212.263778117704;212.29068811741];

    srcLocation=[199.986,212.124]; %lambert coordinates of the BR1 chimney

    locations(:,1)=locations(:,1)-srcLocation(1);
    locations(:,2)=locations(:,2)-srcLocation(2);
    locations=locations*1000;
    locations(:,3)=1.5*ones(7,1);
  elseif strcmp(monitoringPoints,'fromFile')
    fid=fopen(strcat(inputDir,'/monitoringPoints'),'r');
    locations=textscan(fid,'%s %f %f %f %s','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
    fclose(fid);
    if fid == -1
     error('File monitoringPoints in directory %s not found.',inputDir);
    end
    locations=[locations{2},locations{3},locations{4}];
    locationNames={};
    for i=1:size(locations,1)
      locationNames{i}=strcat('S',num2str(i));
    end
  elseif strcmp(monitoringPoints,'test')
    locationNames={'IMR/M11','IMR/M10'};
    locations=zeros(2,3);
    locations(1,:)=[([200.143503293714,212.237160795659]-[199.986,212.124])*1000,1.5];
    locations(2,:)=[([199.961095905096,212.29068811741]-[199.986,212.124])*1000,1.5];
    normFactor=sqrt(sum(locations.^2,2));
    locations=100*[locations(1,:)/normFactor(1);locations(2,:)/normFactor(2)];
  else
    error('Monitoring point set %s is not implemented.',monitoringPoints);
  end

  if ~constWindDirection
    [minDist,i]=max(sqrt(sum(locations(:,1:2).^2,2)));
    if Lx < minDist
       error('Monitoring point %d falls potentially outside domain; Lx < %5.2f',i,minDist);
    end

    if Ly < minDist
%       error('Monitoring point %d falls outside domain; Ly < %5.2f',i,minDist);
    end

    if Lz < 1.5
       error('Monitoring points fall outside domain; Lz < 1.5');
    end
  end

  modelInput.dose.noMonitoringPoints=size(locations,1);
  modelInput.dose.monitoringPointNames=locationNames;
  modelInput.dose.monitoringPoints=locations;
end
