function Mesh = createMeshCylindrical(L,N)
Nx=N(1);
Ny=N(2);
d=1e-6; % distance of the first face from center
dx=L(1)/Nx;                         %cell witdh (uniform grid) 
dy=L(2)/Ny;                         %cell witdh (uniform grid)
Mesh.Csizex=ones(Nx+2,1)*dx;      %Cell size (distance between cell centers (can be manualy inputed)
Mesh.Csizey=ones(Ny+2,1)*dy;      %Cell size (distance between cell centers (can be manualy inputed)
Mesh.Clocx=[1:Nx]'*dx-dx/2+d;       %Cell center location
Mesh.Clocy=[1:Ny]'*dy-dy/2;       %Cell center location
Mesh.Flocx=[0:Nx]'*dx+d;            %Face locations
Mesh.Flocy=[0:Ny]'*dy;            %Face locations
Mesh.Nx=Nx;
Mesh.Ny=Ny;

Mesh.numbering=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
Mesh.c=Mesh.numbering([1,end], [1,end]);
%% no need to calc here
Dr=diff(Mesh.Flocx);
Dt=diff(Mesh.Flocy);
R=0.5*(Mesh.Flocx(1:end-1)+Mesh.Flocx(2:end)); %average of rn and rs
Mesh.V=ones(Nx+2,Ny+2);
Mesh.V(2:Nx+1,2:Ny+1)=(R.*Dr)*Dt';
% Mesh.V=[Mesh.V(1,:); Mesh.V; Mesh.V(end,:)]; %adding ghoast cell volume ,x(dummy)
% Mesh.V=[Mesh.V(:,1) Mesh.V Mesh.V(:,end)]; %adding ghoast cell volume ,y(dummy)

