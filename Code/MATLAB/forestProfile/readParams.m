function params=readParams(maxline)

 global dict meteoFile

 fileID=fopen('./constant/inputDict');
 dict=textscan(fileID,'%s %s','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
 fclose(fileID);

 month_=dict{2}{strcmp(dict{1},'month')};
 day_=dict{2}{strcmp(dict{1},'day')};
 fileNames={strcat('./constant/meteo',month_,day_),strcat('./constant/meteo',day_,month_)};
 if sum(strcmp(dict{1},'year'))>0
   year_=dict{2}{strcmp(dict{1},'year')};
   fileNames{3}=strcat('./constant/meteo',day_,month_,year_);
 end
 fileID=-1;
 i=1;
 while (fileID < 0 & i<=length(fileNames))
   fileID=fopen(fileNames{i});
   i=i+1;
 end
 if fileID == -1
   error('Meteo file not found.');
 end
 meteoFile=textscan(fileID,'%f %f','Delimiter',{' ','\t'},'CommentStyle','%','MultipleDelimsAsOne',1);
 fclose(fileID);

 if nargin == 0
   params=initialize();
 else
   params=initialize(maxline);
 end
end

