 function [A, B] = diffusionCoeffsCylindrical(m,f)
% matrix A is implicit terms (LHS)
% matrix B is explicit terms (RHS)
Nr=m.Nx;
Nt=m.Ny;
G=reshape(1:(Nr+2)*(Nt+2), Nr+2, Nt+2);

Dr=repmat(m.Csizex,1,Nt+2);
Dt=repmat(m.Csizey',Nr+2,1);
rf=repmat(m.Flocx,1,Nt);
rp=repmat(m.Clocx,1,Nt);

dr=0.5*(Dr(1:end-1,:)+Dr(2:end,:));
dt=0.5*(Dt(:,1:end-1)+Dt(:,2:end));

% build the sparse matrix based on the numbering system
iix = zeros(3*(Nr+2)*(Nt+2),1);	iiy = zeros(3*(Nr+2)*(Nt+2),1);
jjx = zeros(3*(Nr+2)*(Nt+2),1);	jjy = zeros(3*(Nr+2)*(Nt+2),1);
sx = zeros(3*(Nr+2)*(Nt+2),1);	sy = zeros(3*(Nr+2)*(Nt+2),1);
NrNt = Nr*Nt;	mny = Nr*Nt;

DN=rf(2:Nr+1,:).*Dt(2:Nr+1,2:Nt+1)./dr(2:Nr+1,2:Nt+1);
DS=rf(1:Nr,:).*Dt(2:Nr+1,2:Nt+1)./dr(1:Nr,2:Nt+1);

DW=Dr(2:Nr+1,2:Nt+1)./rp./dt(2:Nr+1,2:Nt+1);
DE=Dr(2:Nr+1,2:Nt+1)./rp./dt(2:Nr+1,1:Nt);

AN=reshape(DN,NrNt,1);
AS=reshape(DS,NrNt,1);
AW=reshape(DW,NrNt,1);
AE=reshape(DE,NrNt,1);

AP=reshape((AE+AW+AN+AS),NrNt,1);


%% node numbering scheme
%   --> Nt
%|  1        Nr+3  . . .
%|  2        Nr+4
%Nr .        .
%   .        .
%   .        .     . . . Nr*Nt
%   Nr+2

% Matrix A is NrNt * NrNt, each row blongs to one node, the row number is
% the node number, the column numbers are the node numbers. The coeeficient
% is applied to the node numbers (col number) coresponding to N,W,E,S which
% can be easily determined from the numering scheme matrix
%%
% build the sparse matrix based on the numbering system
row_index = reshape(G(2:Nr+1,2:Nt+1),NrNt,1); % index of internal nodes 
ii(1:5*NrNt) = repmat(row_index,5,1);         % row number is identical for all 5 coeff

%%% jj= node number associated to [AE AW AP AS AN]
jj(1:5*NrNt) = [reshape(G(2:Nr+1,1:Nt)  ,NrNt,1); ... AE
                reshape(G(2:Nr+1,3:Nt+2),NrNt,1); ... AW
                reshape(G(2:Nr+1,2:Nt+1),NrNt,1); ... AP
                reshape(G(1:Nr  ,2:Nt+1),NrNt,1); ... AS
                reshape(G(3:Nr+2,2:Nt+1),NrNt,1)  ... AN
                ];

sA(1:5*NrNt) = f*[-AE; -AW; AP; -AS; -AN];
sB(1:5*NrNt) = (1-f)*[AE; AW; -AP; AS; AN];




% build the sparse matrix 
k = 5*NrNt;
A = sparse(ii(1:k), jj(1:k), sA(1:k), (Nr+2)*(Nt+2), (Nr+2)*(Nt+2));
B = sparse(ii(1:k), jj(1:k), sB(1:k), (Nr+2)*(Nt+2), (Nr+2)*(Nt+2));