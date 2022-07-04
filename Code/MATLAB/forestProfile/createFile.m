function createFile(params)

 global bool
 if bool == 0
	fileName='transportProperties';
 elseif bool == 1
	fileName='pollutantProperties';
 end 

 fileID=fopen(strcat('./constant/',fileName),'w');
 fprintf(fileID,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
 fprintf(fileID,'| =========                 |                                                 |\n');
 fprintf(fileID,'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
 fprintf(fileID,'|  \\\\    /   O peration     | Version:  4.1.x                                 |\n');
 fprintf(fileID,'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
 fprintf(fileID,'|    \\\\/     M anipulation  |                                                 |\n');
 fprintf(fileID,'\\*---------------------------------------------------------------------------*/\n');
 fprintf(fileID,'FoamFile\n');
 fprintf(fileID,'{\n');
 fprintf(fileID,'    version     4.1;\n');
 fprintf(fileID,'    format      ascii;\n');
 fprintf(fileID,'    class       dictionary;\n');
 fprintf(fileID,'    location    "constant";\n');
 fprintf(fileID,'    object      %s;\n',fileName);
 fprintf(fileID,'}\n');
 fprintf(fileID,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n');

 if bool == 1
  xy=params.stations;
  nos=length(xy(:,1));
  z=1.5;

  fprintf(fileID,'active              yes;\n\n');
  fprintf(fileID,'radiation           yes;\n\n');
  fprintf(fileID,'outputInterval      %s;\n\n',params.outputInterval);
  fprintf(fileID,'isotopes            "/vsc-hard-mounts/leuven-data/317/vsc31746/OpenFOAM/isotopes/";\n\n');
  fprintf(fileID,'monitoringPoints\n\n');
  fprintf(fileID,'(\n');
  for i=1:nos
   fprintf(fileID,'  ( %f\t%f\t%f )\n',xy(i,1),xy(i,2),z);
  end
  fprintf(fileID,');\n\n');
  fprintf(fileID,'sources\n');
  fprintf(fileID,'(\n');
  fprintf(fileID,'    source1\n');
  fprintf(fileID,'    {\n');
  fprintf(fileID,'        pollutant           Ar-41;\n\n');
  fprintf(fileID,'        volumeSourceModel   singleSource;\n');
  fprintf(fileID,'        source              on;\n');
  fprintf(fileID,'        dimensionless       off;\n');
  fprintf(fileID,'        timeDependence      off;\n\n');
  fprintf(fileID,'        singleSourceCoeffs\n');
  fprintf(fileID,'        {\n');
  fprintf(fileID,'            position    ( 0   0   %s);\n',params.stackHeight);
  fprintf(fileID,'            size        1;\n');
  fprintf(fileID,'            Lref        1;\n');
  fprintf(fileID,'            Vref        1;\n');
  fprintf(fileID,'        }\n');
  fprintf(fileID,'    }\n\n');
  fprintf(fileID,');\n\n');
 else
  fprintf(fileID,'uStar           %5.4f;\n',params.uStar);
  fprintf(fileID,'kappa           0.4;\n');
  fprintf(fileID,'d               %6.4f;\n',params.d);
  fprintf(fileID,'Cd              %4.2f;\n',params.Cd);
  fprintf(fileID,'LAI             %4.2f;\n',params.LAI);
  fprintf(fileID,'DBH             %4.3f;\n',params.DBH);
  fprintf(fileID,'treeHeight      %3.1f;\n',params.h);
  fprintf(fileID,'canopyDepth     %3.1f;\n',params.canopyDepth);
  fprintf(fileID,'extCoeff        %4.3f;\n',params.extCoeff);
  fprintf(fileID,'standDensity    %5.1f;\n',params.standDensity);
  fprintf(fileID,'mu              %5.4f;\n',params.mu);
  fprintf(fileID,'zc              %6.4f;\n',params.zc );
  fprintf(fileID,'zd              %6.4f;\n',params.zi);
  fprintf(fileID,'scale           %6.4f;\n',params.scale);
  if isfield(params,'zStar')
    fprintf(fileID,'zStar           %6.4f;\n',params.zStar);
  end
  if isfield(params,'z0')
    fprintf(fileID,'z0              %6.4f;\n',params.z0);
  end
  if isfield(params,'z0g')
    fprintf(fileID,'z0g             %6.4f;\n',params.z0g);
  else
    fprintf(fileID,'z0g             %6.4f;\n',0.25);
  end
  fprintf(fileID,'windDirection   (%5.4f %5.4f 0);\n\n',cos(params.windAngle),sin(params.windAngle));
  fprintf(fileID,'Sc              Sc [ 0 0  0 0 0 0 0 ] 0.9;\n');
  fprintf(fileID,'Dv              Dv [ 0 2 -1 0 0 0 0 ] 1e-20;\n\n');
  fprintf(fileID,'pollutant       yes;\n\n');
 end

  fprintf(fileID,'// ************************************************************************* //');
  fclose(fileID);
end

