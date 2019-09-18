function [A, B] = transientCoeffs(m,dt,prop)
% matrix A is implicit terms
% matrix B is explicit terms
Dx=diff(m.Floc);
Nx=m.Nx;
a=prop.c*prop.rho*Dx/dt;
b=prop.c*prop.rho*Dx/dt;
IntNodes=2:Nx+1;
A=sparse(IntNodes,IntNodes,a,Nx+2,Nx+2);
B=sparse(IntNodes,IntNodes,b,Nx+2,Nx+2);