function grid=mesh(type,gridInput,refinementHeight)

 if strcmp(type,'decompose')

   lx0=gridInput.x.minSize; lx=gridInput.x.maxSize; alphax=gridInput.x.regionLengths;
   ly0=gridInput.y.minSize; ly=gridInput.y.maxSize; alphay=gridInput.y.regionLengths;
   lz0=gridInput.z.minSize; lzg=gridInput.z.belowSrc.maxSize; alphazg=gridInput.z.belowSrc.regionLengths;
   lzt=gridInput.z.aboveSrc.maxSize; alphazt=gridInput.z.aboveSrc.regionLengths;

   [dy,ny,Ny]=calcMeshParams(ly,ly0,alphay);
   [dx,nx,Nx]=calcMeshParams(lx,lx0,alphax);
   [dzg,nzg,Nzg]=calcMeshParams(lzg,lz0,alphazg);
   [dzt,nzt,Nzt]=calcMeshParams(lzt,lz0,alphazt);

   y=calcMesh(ny,Ny,dy);
   x=calcMesh(nx,Nx,dx);
   zg=calcMesh(nzg,Nzg,dzg);
   zt=calcMesh(nzt,Nzt,dzt);
   zg=zg(end-(0:(end-1))); dzg=dzg(end-(0:(end-1))); nzg=nzg(end-(0:(end-1)));
%   zt=zt(end-(0:(end-1))); dzt=dzt(end-(0:(end-1))); nzt=nzt(end-(0:(end-1)));
   z=[refinementHeight-zg;zt+refinementHeight]; dz=[dzg;dzt]; nz=[nzg;nzt];
%   z=[zt+refinementHeight;refinementHeight-zg]; dz=[dzt;dzg]; nz=[nzt;nzg];

   grid.x.C=x; grid.x.noRegions=Nx+1; grid.x.regionSizes=nx; grid.x.regionSizings=dx;
   grid.y.C=y; grid.y.noRegions=Ny+1; grid.y.regionSizes=ny; grid.y.regionSizings=dy;
   grid.z.C=z; grid.z.noRegions=Nzg+Nzt+2; grid.z.regionSizes=nz; grid.z.regionSizings=dz;

   disp(sprintf('-> mesh has %d cells\n',2*sum(nx)*sum(ny)*sum(nz)));
%{
   figure(1);
   [X,Y]=meshgrid(x,y);
   plot(X,Y,'ok');
   figure(2);
   [X,Z]=meshgrid(x,z);
   plot(X,Z,'ok');
   figure(3);
   [Y,Z]=meshgrid(y,z);
   plot(Y,Z,'ok');
%}
 elseif strcmp(type,'full')
   x=gridInput.x.C; Nx=gridInput.x.noRegions; nx=gridInput.x.regionSizes; dx=gridInput.x.regionSizings;
   y=gridInput.y.C; Ny=gridInput.y.noRegions; ny=gridInput.y.regionSizes; dy=gridInput.y.regionSizings;
   z=gridInput.z.C; Nz=gridInput.z.noRegions; nz=gridInput.z.regionSizes; dz=gridInput.z.regionSizings;

   grid.noRegions=0;
   regionid=0;
   i1=0;
   for i=1:Nx
    i0=i1+1;
    i1=i1+nx(i);
    j1=0;
    for j=1:Ny
     j0=j1+1;
     j1=j1+ny(j);
     k1=0;
     for k=1:Nz
	k0=k1+1;
	k1=k1+nz(k);
        regionid=regionid+1;
	[X,Y,Z]=meshgrid(x(i0:i1),y(j0:j1),z(k0:k1));
	X=reshape(X,[],1);
	Y=reshape(Y,[],1);
	Z=reshape(Z,[],1);

        grid=setfield(grid,strcat('region',num2str(regionid)),'size',nx(i)*ny(j)*nz(k));
	grid=setfield(grid,strcat('region',num2str(regionid)),'C',[X,Y,Z]);
        grid=setfield(grid,strcat('region',num2str(regionid)),'V',dx(i)*dy(j)*dz(k));
	grid=setfield(grid,strcat('region',num2str(regionid)),'range',[[i0,i1];[j0,j1];[k0,k1]]);
     end
    end
   end
   grid.noRegions=regionid;
 end

end
