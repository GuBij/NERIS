
 Q=1; %Table 3.6: Q = 1, Table 3.5: Q = 0.5712/1.5 (M), Q = 1 (GPG)
 threshold=10; %nSv/h
 measurements='GPG_10m/'; %'measurements/', 'GPG_10m/'
 modelName = 'GBM'; % F,LF,OF,LOF,GPG,GBM
 srcHeight = '10m'; %60m, 10m
 fac=2;
 files={'090117','100117','140217','130317','140317','280317','290317','110417','120417','130417','180417','190417','240417','030517','040517','200617'};
 
 workDir='./';
 N=0; Ntot=0; diff=0; middle=0; sqDiff=0; sumMeasured=0; sumSimulated=0; FAC2=0; 

 for i=1:length(files)
  measuredData=readFile(strcat(workDir,measurements,files{i},'.txt'));
  nocol=size(measuredData,2);
  measuredData(:,1)=[]; % remove the time column
  simulatedData=readFile(strcat(workDir,strcat(modelName,'_',srcHeight,'/'),files{i},'.txt'));
  simulatedData(:,1)=[]; % remove the time column
  simulatedData = Q*simulatedData;
 
  nanValues=find(isnan(measuredData));
  tooLow=find(measuredData < threshold);
 
  noel=numel(measuredData); Ntot=Ntot+noel; 
  nonRejectedValues=setdiff(1:noel,[nanValues;tooLow]); N=N+length(nonRejectedValues);
  measuredData=measuredData(nonRejectedValues);
  simulatedData=simulatedData(nonRejectedValues);

  ratio=measuredData./simulatedData;
  FAC2=FAC2+length(find(ratio > 1/fac & ratio < fac));
  diff=diff+sum(measuredData-simulatedData); 
  middle=middle+sum(0.5*(measuredData+simulatedData));
  sqDiff=sqDiff+sum((measuredData-simulatedData).^2);
  sumMeasured=sumMeasured+sum(measuredData);
  sumSimulated=sumSimulated+sum(simulatedData);
 end

 FB=diff/middle;
 NMSE=N*sqDiff/(sumMeasured*sumSimulated);
 FAC2=FAC2/N;

 disp(sprintf('\n%d values were not rejected',N));
 disp(sprintf('FB: %6.4f',FB));
 disp(sprintf('NMSE: %3.2f',NMSE));
 disp(sprintf('FAC2: %4.2f\n',FAC2));
