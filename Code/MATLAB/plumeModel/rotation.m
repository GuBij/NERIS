function [points,L]=rotation(windDir,points)

 global axisGrid

 windDir=windDir-floor(windDir/360)*360;
 rotAngle=(-(90-abs(windDir)))*pi/180+pi;
 rotMat=[[cos(rotAngle),-sin(rotAngle)];[sin(rotAngle),cos(rotAngle)]];

 points(:,1:2)=points(:,1:2)*rotMat';

 noPoints=size(points,1);
 L=zeros(noPoints,1);
 axis={'x','y','z'};
 for i=1:noPoints
   coord=points(i,:);
   if coord(1)<0
     coord(1)=-coord(1);
   end
   if coord(2)<0
     coord(2)=-coord(2);
   end
   cellVolume=1;
   for j=1:3
     centers=getfield(axisGrid,axis{j},'C');
     sizes=getfield(axisGrid,axis{j},'regionSizes');
     sizings=getfield(axisGrid,axis{j},'regionSizings');
     cellid=max(sum(coord(j)>centers),1); %disp(sprintf('\ncellid=%d, coord(j)=%f, j=%d\n',cellid,coord(j),j));
     if abs(coord(j)-centers(cellid+1))<abs(coord(j)-centers(cellid))
       cellid=cellid+1;
     end
     coord(j)=centers(cellid);
     regionid=sum(cellid<=cumsum(sizes));
     cellVolume=cellVolume*sizings(regionid);
   end
   L(i)=4*(cellVolume)^(1/3); %points(i,:), disp(sprintf('\npoints(%d,1)=%f, points(%d,1)+L(i)=%f\n',i,points(i,1),i,points(i,1)+L(i))); %, points(%d,1)+L=%f, points(%d,2)=%d\n',i,points(i,1),i,points(i,1)+L,i,points(i,2)));
   if points(i,1)<0 && (points(i,1)+L(i))>0 && points(i,2)<0
     points(i,1)=-coord(1);
     points(i,2)=-coord(2);
   elseif points(i,1)<0 && (points(i,1)+L(i))>0
     points(i,1)=-coord(1);
     points(i,2)=coord(2);
   elseif points(i,1)>0 && points(i,2)<0
     points(i,1)=coord(1);
     points(i,2)=-coord(2);
   elseif points(i,1)>0
     points(i,1)=coord(1);
     points(i,2)=coord(2);
   end
   points(i,3)=coord(3);
  end

end
