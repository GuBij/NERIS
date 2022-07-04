function maxC=pointFieldConc(U,Q,h,parameterization,option)

  monitoringPoints='meteo/monitoringPoints3';
  outputFile=strcat('pointFieldConc_src',num2str(h),'m_gaussPG');

  if nargin==5 && strcmp(option,'meteo')
   day=U;
   monitoringPoints=strcat('monitoringPoints/monitoringPoints',day);
   outputFile=strcat(day,'/pointFieldConc.out');
   fid=fopen(strcat('meteo/meteo',day),'r');
   meteoInCell=textscan(fid,'%f %*f','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
   fclose(fid);
   U=mean(meteoInCell{1});
  elseif isstr(U)
   U=str2num(U);
  end

  if isstr(Q)
    Q=str2num(Q);
  end

  if isstr(h)
    h=str2num(h);
  end

  fid=fopen(monitoringPoints,'r');
  xyz=textscan(fid,'%s %f %f %f %s','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
  fclose(fid);

  if fid == -1
   error('File %s not found.',monitoringPoints);
  end

  x=xyz{2}; 
  y=xyz{3};
  z=xyz{4};

  n=length(x);

  [sigy,sigz]=disPar(x,parameterization);
  Cfield=zeros(4,n);
  scale=Q/(2*pi*U*3600);
  for i=1:n
    if x(i)>0
     Cx=scale/(sigy(i)*sigz(i));
     Cy=exp(-0.5*(y(i)/sigy(i))^2);
     Cfield(4,i)=(exp(-0.5*((z(i)-h)/sigz(i))^2)+exp(-0.5*((z(i)+h)/sigz(i))^2))*Cy*Cx;
    end
    Cfield(1:3,i)=[x(i),y(i),z(i)];
  end

  [maxC,I] = max(Cfield(4,:)); maxC = ([Cfield(1:3,I);maxC])';
  columns={'x','y','z','C'};
  
  writeToFile(Cfield,'./output',outputFile,columns,'%e');

end

