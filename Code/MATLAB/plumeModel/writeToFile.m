function writeToFile(dose,outputDir,fileName,locationNames,colFormat)

 noo=length(locationNames);
 format='%d';
 for i=1:noo
  format=strcat(format,'\t',colFormat); %'%6.3f');
 end
 format=strcat(format,'\n');

 noTimes=size(dose,2);
 dose=[1:noTimes;dose];
 fileID=fopen(sprintf('%s/%s17.txt',outputDir,fileName),'w');
 fprintf(fileID,'period\t');
 for i=1:noo %(2*noo+1)
%   if i<=noo
     fprintf(fileID,'%s\t',locationNames{i});
%{
   elseif i==noo+1
     fprintf(fileID,'\n=======\t');
   else
     fprintf(fileID,'=======\t');

   end
%}
 end
 fprintf(fileID,'\n');
 fprintf(fileID,format,dose);
 fclose(fileID);

end
