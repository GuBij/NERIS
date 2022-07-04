function plotStruct=SI()

 dataDir='SI/'; %data directory
 fileList='dirList.orig';

 fileID=fopen(fullfile(dataDir,fileList));
 fileNames_=textscan(fileID,'%s','Delimiter',{' ','\t'},'CommentStyle','#','MultipleDelimsAsOne',1);
 fclose(fileID);
 fileNames_=fileNames_{1};
 nof=length(fileNames_);
 nos=3;
 dose=zeros(3,1,nof);

 for i = 1:nof
   dosei=readFile(strcat(dataDir,fileNames_{i},'.txt')); dosei(1)=[]; dosei=dosei(1:nos)';
   dose(:,:,i)=dosei;
 end
 clear dosei;

 nopar=0.5*nof;
 if mod(nof,2)>0
  error('each parameter has to be changed twice; an odd number of parameters has been found.');
 end

 plotStruct.par1='';
 for i=1:nopar
   doseMax=max(dose(:,:,2*i),dose(:,:,2*i-1)); 
   doseMin=min(dose(:,:,2*i),dose(:,:,2*i-1));
   index=(doseMax-doseMin)./doseMax;
   fileName=fileNames_{2*i};
   par=''; k=1; kmax = length(fileName);
   while (k <= kmax && ~( fileName(k) > 47 && fileName(k) < 58 ))
     par=strcat(par,fileName(k));
     k=k+1;
   end
   if strcmp(par,'z') && strcmp(fileName(3),'g')
     par='z_{0,s}';
   elseif strcmp(par,'z') && strcmp(fileName(2),'0')
     par='z_{0}';
   elseif strcmp(par,'extCoeff')
     par='\gamma_{ext}';
   elseif strcmp(par,'canopyDepth')
     par='h_{c}';
   elseif strcmp(par,'Cd')
     par='C_{d}';
   elseif strcmp(par,'stackHeight')
     par='h_{eff}';
   elseif strcmp(par,'zStar')
     par='z_{\ast}';
   elseif strcmp(par,'U')
     par='u_{ref}';
   end
   plotStruct=setfield(plotStruct,strcat('par',num2str(i)),'name',par);
   plotStruct=setfield(plotStruct,strcat('par',num2str(i)),'value',index);
 end

end

