function [A, RHS] = BoundaryCondition(m,BC)

Nr=m.Nx;
Nt=m.Ny;
nb=2*(Nr+Nt); %number of boundary nodes
% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);
% define the RHS column vector
RHS = zeros((Nr+2)*(Nt+2), 1);

q=0;
G=reshape(1:(Nr+2)*(Nt+2), Nr+2, Nt+2);

% assign value to the corner nodes (useless cells)
q = 1:4;
ii(q) = m.c(:); jj(q) = m.c(:);
s(q) = 1; BCRHS(m.c(:)) = 0;

%Left boundary (High Theta) 
i=2:Nr+1;
j=Nt+2;
q=q(end)+(1:Nr);
ii(q) = G(i,j);  jj(q) = G(i,j);   s(q) = BC.left.A/2 + BC.left.B/m.Csizey(end); 
q=q(end)+(1:Nr);
ii(q) = G(i,j);  jj(q) = G(i,j-1); s(q) = BC.left.A/2 - BC.left.B/m.Csizey(end);
RHS(G(i,j)) = BC.left.C; 
% Right boundary (low Theta) 
i=2:Nr+1;
j = 1;
q=q(end)+(1:Nr);
ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = BC.right.A/2 + BC.right.B/m.Csizey(1);
q=q(end)+(1:Nr);
ii(q) = G(i,j);  jj(q) = G(i,j+1);    s(q) = BC.right.A/2 - BC.right.B/m.Csizey(1);
RHS(G(i,j)) = BC.right.C;


% Top boundary (High R ) 
i=Nr+2;
j =2:Nt+1;
q=q(end)+(1:Nt);
ii(q) = G(i,j);  jj(q) = G(i,j);  s(q) = BC.top.A/2 + BC.top.B/m.Csizex(end);
q=q(end)+(1:Nt);
ii(q) = G(i,j);  jj(q) = G(i-1,j);s(q) = BC.top.A/2 - BC.top.B/m.Csizex(end);
RHS(G(i,j)) = BC.top.C;

%Bottom boundary (low r) 
i=1;
j=2:Nt+1;
q=q(end)+(1:Nt);
ii(q) = G(i,j);  jj(q) = G(i,j);   s(q) = BC.bottom.A/2 + BC.bottom.B/m.Csizex(1); 
q=q(end)+(1:Nt);
ii(q) = G(i,j);  jj(q) = G(i+1,j); s(q) = BC.bottom.A/2 - BC.bottom.B/m.Csizex(1);
RHS(G(i,j)) = BC.bottom.C; 


A=sparse(ii(1:q(end)), jj(1:q(end)), s(1:q(end)), ...
    (Nr+2)*(Nt+2), (Nr+2)*(Nt+2));