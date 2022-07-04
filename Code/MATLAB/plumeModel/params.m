%DISPERSION
 h=60; %stack height
 href=69; %height of the wind speed measurement
 parameterization='BM'; % Three types: Bultynck-Malet ('BM'), Pasquill-Gifford ('PG') or Pasquill ('Pasq')
 Q=1.5*10^11; %source strength in Bq/h
 nuclide='Ar-41';
 monitoringPoints='BR1'; % employ pre-specified monitoring points ('BR1') or specified by file <inputDir>/monitoringPoints ('fromFile')

 inputDir='./meteo'; %input directory (meteo files)
 outputDir='./output'; %output directory

 constWindSpeed=0; % employ stationary wind speed? yes (1) or no (0)
 constWindDirection=0; % employ stationary wind direction? yes (1) or no (0)

%MESH
%x-direction
 Lx=2200; %2200; %480; %1100; % domain length along the wind direction
 lx=5; % cell size near the end of the domain along the wind direction
 lx0=5; % cell size at the source along the wind direction
 alphax=[1100;1100]; %[1100;1100]; %[240;240]; [550;550]; % the sum of the elements of alphax must equal Lx
 % if alphax has N elements, then the part of the x-axis corresponding with the i-th element of alphax has cell size lx0*(lx/lx0)^((i-1)/(N-1))

%y-direction
 Ly=300; % domain length along the cross wind direction
 ly=5; % cell size near the end of the domain along the cross wind direction
 ly0=5; % cell size at the source along the wind direction
 alphay=[150;150]; % the sum of the elements of alphay must equal Ly
 % if alphay has N elements, then the part of the y-axis corresponding with the i-th element of alphay has cell size ly0*(ly/ly0)^((i-1)/(N-1))

%z-direction
 Lz=200; % the height of the domain
 lzg=1; %3; %1; % cell size at the ground of the domain
 lzt=10; % cell at the top of the domain
 lz0=3; % cell size at stack height
 alphazg=[40;40]; % {30;30] % the sum of the elements of alphazg must equal the stack height (h)
 % alphazg=[5;5]; %[30;15;15]; %/h;
 % if alphazg has N elements, then the i-th element of alphazg has cell size lz0*(lzg/lz0)^((i-1)/(N-1))
 alphazt=[40;40;40]; %[30;50;60]; %[20;20;50;50;50]; %[30;50;60]; % the sum of the elements of alphazt must equal Lz minus the stack height (h)
 % if alphazt has N elements, then the part of the z-axis corresponding with the i-th element of alphazt has cell size lz0*(lzt/lz0)^((i-1)/(N-1))
