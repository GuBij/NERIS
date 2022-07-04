function params=intialize(maxline)

 global bool dict meteoFile

 if bool == 1
  outputInterval_=dict{2}{strcmp(dict{1},'outputInterval')};
  stackHeight_=dict{2}{strcmp(dict{1},'stackHeight')};
  importStations_=dict{2}{strcmp(dict{1},'importStations')};

  if maxline == 0
    disp(sprintf('\nConstant wind direction is assumed.'));
    windDir=mean(meteoFile{2});
  else
    disp(sprintf('\nVariable wind direction is assumed.'));
    windDir=meteoFile{2}(maxline);
  end
  disp(sprintf('Wind direction: %4.2f\n',windDir));
  if strcmp(importStations_,'default')    
    stat=lambertToCart(windDir);
  else
    fileID=fopen(importStations_);
    stat=textscan(fileID,'%*s %f %f %*f %*s','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
    fclose(fileID);
    stat=[stat{1},stat{2}]/1000;
    stat=lambertToCart(windDir,stat,[0,0]); 
  end
  params.stations=stat;
  params.outputInterval=outputInterval_;
  params.stackHeight=stackHeight_;

 elseif bool == 0
  noi=size(dict{1},1);

  for i=1:noi
   if strcmp(dict{1}{i},'h')
	 params.h=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'LAI')
	params.LAI=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'extCoeff')
	params.extCoeff=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'standDensity')
	params.standDensity=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'DBH')
	params.DBH=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'zref')
        params.zref=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'Cd')
        params.Cd=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'z0')
	params.z0=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'z0g')
        params.z0g=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'canopyDepth')
        params.canopyDepth=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'zStar')
	params.zStar=str2num(dict{2}{i});
   elseif strcmp(dict{1}{i},'cstWindSpeed')
	if strcmp(dict{2}{i},'yes')
		cstWindSpeed=1;
	else
		cstWindSpeed=0;
	end
   end
  end

  if cstWindSpeed == 1
	disp(sprintf('\nConstant wind speed is assumed.'));
	windSpeed=meteoFile{1};
	params.uref=mean(windSpeed((end-24):end));
  else	
	disp(sprintf('\nVariable wind speed is assumed.'));
	params.uref=meteoFile{1}(maxline);
  end

 if nargin == 0
        disp('Constant wind direction is assumed.');
	meanWindDir=0; windDir=0;
 else
        disp('Variable wind direction is assumed.');
	windDir=meteoFile{2}(maxline);
	meanWindDir=mean(meteoFile{2});
 end
  params.windAngle=(meanWindDir-windDir)*pi/180;

  disp(sprintf('Wind speed: %4.2f',params.uref));
  disp(sprintf('Wind angle: %4.2f\n',params.windAngle));
 end

end
