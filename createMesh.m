function Mesh = createMesh(L,Nx)
dx=L/Nx;                         %cell witdh (uniform grid) 
Mesh.Csize=ones(Nx+2,1)*dx;      %Cell size (distance between cell centers (can be manualy inputed)
Mesh.Cloc=[0:Nx+1]'*dx-dx/2;       %Cell center location
Mesh.Floc=[0:Nx]'*dx;            %Face locations
Mesh.Nx=Nx;
Mesh.V=


