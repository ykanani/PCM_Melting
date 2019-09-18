function [A, B] = sourceCoeffsCylindrical(m,dt,Ste)
% matrix A is implicit terms
% matrix B is explicit terms (RHS)
Dr=diff(m.Flocx);
Dt=diff(m.Flocy);
R=0.5*(m.Flocx(1:end-1)+m.Flocx(2:end)); %average of RW and RE
Nx=m.Nx;
Ny=m.Ny;
IntNodes = reshape(m.numbering(2:Nx+1,2:Ny+1),Nx*Ny,1); % main diagonal (only internal cells)
V=(R.*Dr)*Dt';
b=V/dt/Ste ;
% b=Dr.*R/dt;

A=[];
B=sparse(IntNodes,IntNodes,b,(Nx+2)*(Ny+2),(Nx+2)*(Ny+2));