function u=windSpeed(z,option,uref,zref,class)

 if strcmp(option,'BM')
  switch class
    case 1
     p=0.53;
    case 2
     p=0.4;
    case 3
     p=0.33;
    case 4
     p=0.23;
    case 5
     p=0.16;
    case 6
     p=0.1;
    case 7
     p=0.33
  end
  u = uref*(z/zref)^p;

end
