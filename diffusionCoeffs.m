 function [A, B] = diffusionCoeffs(m,prop,f)
% matrix A is implicit terms
% matrix B is explicit terms
Dx=m.Csize;
Nx=m.Nx;
dx=0.5*(Dx(1:end-1)+Dx(2:end));
G=1:length(Dx);
AW=prop.k./dx(1:end-1);
AE=prop.k./dx(2:end);
AP=f*(AE+AW);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1),Nx,1); % main diagonal x
iix(1:3*Nx) = repmat(rowx_index,3,1);
jjx(1:3*Nx) = [reshape(G(1:Nx),Nx,1); reshape(G(2:Nx+1),Nx,1); reshape(G(3:Nx+2),Nx,1)];
sA(1:3*Nx) = [-f*AW; f*(AE+AW); -f*AE];
sB(1:3*Nx) = [(1-f)*AW; -(1-f)*(AE+AW); (1-f)*AE];

% build the sparse matrix 
kx = 3*Nx;
A = sparse(iix(1:kx), jjx(1:kx), sA(1:kx), Nx+2, Nx+2);
B = sparse(iix(1:kx), jjx(1:kx), sB(1:kx), Nx+2, Nx+2);