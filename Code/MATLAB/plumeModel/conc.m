function C=conc(modelInput,class,heff)

  global axisGrid;

  Q=modelInput.sourceStrength;
%  h=modelInput.stackHeight;
  parameterization=modelInput.parameterization;
  x=axisGrid.x.C;
  y=axisGrid.y.C;
  z=axisGrid.z.C;
  nx=sum(axisGrid.x.regionSizes);
  ny=sum(axisGrid.y.regionSizes);
  nz=sum(axisGrid.z.regionSizes);

  [sigy,sigz]=disPar(x,parameterization,class); 
  C=zeros(nx,ny,nz);
  scale=Q/(2*pi);
  for i=1:nx
    Cx=scale/(sigy(i)*sigz(i));
    for j=1:ny
      Cy=exp(-0.5*(y(j)/sigy(i))^2);
      for k=1:nz
        C(i,j,k)=(exp(-0.5*((z(k)-heff)/sigz(i))^2)+exp(-0.5*((z(k)+heff)/sigz(i))^2))*Cy*Cx;
      end
    end
  end

end
