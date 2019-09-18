function [Tmelt Tavg]=solve(options,Ste,Bi,PCM)
%%% Simulation paramters
%All t parameters are Fo number
%All T values are theta = (T-Tm)/(Tinf-Tm)
t0=options.t0;
t1=options.t1;
dt=options.dt;                               %deltaT
Ly=options.Ly;%pi;                              %Lengh of domain (radius)
Lx=options.Lx;%pi;                              %Lengh of domain (radius)
Nx=options.Nx;                                  %number of cells (radial)
Ny=options.Ny;                                  %number of cells (theta)
Ncel=(Nx+2)*(Ny+2);
thresh=options.thresh;



%%% Creating Uniform Mesh
%%% *-----|-----o-----|-----o-----|-----o-----|-----*
%%% ?     ?      ?                       ?     ?     ?
%%% T1          T2          T3          T4          T5
%%%      x=0                                  x=L    
%%% Cx1        Cx2                     CxN+1       CxN+2
%%%      Fx1        xe,Re       xw,Rw              FxN+1 
%%%       ??Csize1????Csize2????Csize3??            (Dx) or (Dr)
%%% ? dx1,dr1 ?? dx2,dr2 ?? dx3,dr3 ?? dx4,dr4 ?     (dx) or (dr) 

Mesh=createMeshCylindrical([Lx Ly],[Nx Ny]);                  %Creating Uniform Mesh
%%% Setting Boundary Conditions
% B.C. equation:
%(right)          A*(Tb+Tg)/2+B*(Tb-Tg)/dx = C,
%(left)           -A*(Tb+Tg)/2-B*(Tb-Tg)/dx = -C,
% Tb=T2 or TN+1 (Boundary point, cell center)
% Tg=ghoast cell


%Note: The area in already included in the U calculations, U only needs to 
%be multiplied by the temperature difference to provide Q 






%Initialize
Tinit=options.Tinit;
T=Tinit*ones(Nx+2,Ny+2); % two zeros corresponds to the ghoasted cells and will be corrected according to the B.C.


%%
% %%%% Left Boundary - High Theta (adiabatic, symmetry)
BC.left.A=zeros(1,Nx);
BC.left.B=ones(1,Nx);
BC.left.C=zeros(1,Nx);

%%%% Right Boundary - Low Theta (adiabatic, symmetry)
BC.right.A=zeros(1,Nx);
BC.right.B=ones(1,Nx);               %k*A , because U is already multiplied by area , k(end) k of end cell (boundary
BC.right.C=zeros(1,Nx);

% %%%% Top Boundary - High R, (adiabatic, symmetry)
BC.top.A=Bi;
BC.top.B=1.*ones(1,Ny);               %k*A , because U is already multiplied by area Lx/R =1
BC.top.C=Bi*Tinit.*ones(1,Ny);

% %%%% Left Boundary - Low R, (adiabatic, symmetry)
BC.bottom.A=zeros(1,Ny);
BC.bottom.B=ones(1,Ny);
BC.bottom.C=zeros(1,Ny);

[B, RHS] = BoundaryCondition(Mesh,BC);  % B is multiplied by T , change sign to negative for RHS

[Itran, Etran] = transientCoeffsCylindrical(Mesh,dt);
[IdiffImplicit, EdiffImplicit] = diffusionCoeffsCylindrical(Mesh,1);
[Is, Es] = sourceCoeffsCylindrical(Mesh,dt,Ste);

%Apply boundary condition
Intnodes=Mesh.numbering(2:Nx+1,2:Ny+1);
T(:)=(B+sparse(Intnodes(:),Intnodes(:),1,Ncel,Ncel))\(RHS+sparse(Intnodes(:),1,T(Intnodes(:)),Ncel,1));




%% creating mesh grid and adding boundary points
% Tplot=figure;

Tfull=calcBoundaryPoints(T,'average');
% % % [~,hcT]=contourf(X,Y,Tfull,10);
% % % axis equal
% % % set(hcT,'LineColor','none')
% % % colormap('jet')
% % % colorbar
% % % caxis([-1 1])
% % % 
% % % Gammaplot=figure;
% % % [~, ~, results.gamma(:,:,1)]= calcState(Tfull,PCM,0);
% % % 
% % % [~,hcG]=contourf(X,Y,results.gamma(:,:,1),10);
% % % 
% % % axis equal
% % % set(hcG,'LineColor','none')
% % % % colormap('jet')
% % % colorbar
% % % caxis([0 1])



Tmelt=0;
%%
% pb = CmdLineProgressBar(['Running Ste=' num2str(Ste) 'Bi=' num2str(Bi)]);
Qtot=0;
time=t0;
for t=t0+dt:dt:t1
%     disp(t)
    %%% update Boundary Conditions

    %%%% TOP Boundary (wall+convective)
    Air.T=0.5*(options.Tlow-options.Thigh)*(sign(sin(2*pi*t/options.period))+1)+options.Thigh;
   
    BC.top.A=Bi;
    BC.top.B=1.*ones(1,Ny);               %k*A , because U is already multiplied by area
    BC.top.C=Bi*Air.T.*ones(1,Ny);
    [B, RHS] = BoundaryCondition(Mesh,BC);  % B is multiplied by T , change sign to negative for RHS
    %% gauss-sidel loop
    
    S=[];
    Q=[];
    T2=[];
    % Initial Estimate of Q (fully explicit)
    Q=-Ste*dt*IdiffImplicit*T(:)./(Mesh.V(:));
    % Initial Estimate of S
    [S, T2, dHdT0, dHdT1]= calcS(T,reshape(Q,Nx+2,Ny+2),PCM,0);
    e=1;

    %% 1 in m, 2 is m+1
    %% 0 is old (previous time step)

    while (e>thresh)
      T1=T2;
      %update B.C.
      %%% my method
      Nodes = reshape(Mesh.numbering(1:Nx+2,1:Ny+2),(Nx+2)*(Ny+2),1); % main diagonal (only internal cells)
      Sp=sparse(Nodes,Nodes,-Es*dHdT1(:),(Nx+2)*(Ny+2),(Nx+2)*(Ny+2));
      Sc=Es*S(:)-Sp*T1(:);
      T2(:)=(Itran+IdiffImplicit+B)\((Etran+EdiffImplicit)*T(:)+RHS+Es*S(:)); %implicit solver
      
      e=max(max(abs(T2-T1)));
      % Estimate Q (fully explicit)
      Q=-Ste*dt*IdiffImplicit*T2(:)./(Mesh.V(:));
      
      % update S
      [S, ~, dHdT0, dHdT1]= calcS(T,reshape(Q,Nx+2,Ny+2),PCM,0);
      
    end
%     Tint=I\(E*T+RHS+Es*S); %implicit solver
    T=T2;
    [~, ~, gamma]= calcState(calcBoundaryPoints(T,'average'),PCM,0);
    gammaAvg=volAvg(Mesh,gamma);
    q2d=calcBoundaryPoints(reshape(-Ste*IdiffImplicit*T(:),Nx+2,Ny+2),'zerograd'); %Watts ;%reshape(Q,Nx+2,Ny+2);
    qtot=squeeze(sum(sum(abs(q2d(2:Nx+1,2:Ny+1)),1),2));  % Watts
    Qtot(end+1)= qtot/(pi*options.Lx^2*(options.Ly/(2*pi)));
    time(end+1)= t;
    trapz(time,Qtot);
    Tavg = volAvg(Mesh,T);
%     pb.print(gammaAvg,0.999)
    if (Tavg>0.499)%(trapz(time,Qtot)>1+0.999*Ste)%(Qtot(end)/dt < 1e-3) %(sum(Qtot)>1+0.999*Ste)   %(gammaAvg>0.999 && Tmelt==0)
        Tmelt=t;
        pb = CmdLineProgressBar(['Finished Running Ste=' num2str(Ste) 'Bi=' num2str(Bi)]);
        pb.delete;
%         fprintf('Finished Running Ste=%2.2f, Bi=%2.2f',Ste, Bi)
        break
    end

end
fprintf('Finished Running Ste=%2.2f, Bi=%2.2f',Ste, Bi)


