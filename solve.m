function [stringOut, results]=solve(options,Ste,Bi,meltTime,fraction,Nperiods)
%%% Simulation paramters
%All t parameters are Fo number
%All T values are theta = (T-Tm)/(Tinf-Tm)
options.period=fraction*meltTime*2;
t0=options.t0;
t1=options.period*Nperiods;
nperiodConv=Nperiods;
%options.dt=meltTime/round(meltTime/.01);
dt = options.period/round(options.period/options.dt);                               %deltaT
Ly=options.Ly;%pi;                              %Lengh of domain (radius)
Lx=options.Lx;%pi;                              %Lengh of domain (radius)
Nx=options.Nx;                                  %number of cells (radial)
Ny=options.Ny;                                  %number of cells (theta)
Ncel=(Nx+2)*(Ny+2);
thresh=options.thresh;

e=options.e;%(tl-ts)/2
PCM = calcPCMcurve(e,Ste);

options.Ste = Ste;
options.Bi = Bi;
options.fraction = fraction;
options.meltTime = meltTime;
options.Nperiods = Nperiods;


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
[R,Theta]=ndgrid([Mesh.Flocx(1); Mesh.Clocx; Mesh.Flocx(end)],[Mesh.Flocy(1); Mesh.Clocy; Mesh.Flocy(end)]);
X=R.*cos(Theta);
Y=R.*sin(Theta);
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



results.T(:,:,1)=Tfull;
results.S(:,:,1)=Tfull*0;
results.Q(:,:,1)=Tfull*0;
results.pcmk(:,:,1)=Tfull*0;
results.gamma(:,:,1)=Tfull*0;
results.gammaAvg(1)=volAvg(Mesh,results.gamma(:,:,1));
results.t=0;
results.Mesh=Mesh;
results.AirT=Tinit;
results.options=options;
results.Tmelt=0;
results.time=t0;
fprintf('\nStarted Running Ste=%2.2f, Bi=%2.2f\n',Ste, Bi)
pause(0.5)
%%
% fprintf('Started Running Ste=%2.2f, Bi=%2.2f \n',Ste, Bi)
% pb = CmdLineProgressBar(['Running Ste=' num2str(Ste) 'Bi=' num2str(Bi) '...t/T=']);
for t=t0+dt:dt:t1
%     disp(t)
    %%% update Boundary Conditions
%     pb.print(t/options.period,options.t1/options.period)
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
    results.T(:,:,end+1)=calcBoundaryPoints(T,'average');
    results.S(:,:,end+1)=calcBoundaryPoints(reshape(Ste*Es*S(:),Nx+2,Ny+2),'zerograd');
    results.Q(:,:,end+1)=calcBoundaryPoints(reshape(-Ste*IdiffImplicit*T(:),Nx+2,Ny+2),'zerograd') ;
    results.AirT(end+1)=Air.T;
    [~, ~, results.gamma(:,:,end+1)]= calcState(calcBoundaryPoints(T,'average'),PCM,0);
    results.gammaAvg(end+1)=volAvg(Mesh,results.gamma(:,:,end));
    results.time(end+1)=t;
    if (results.gammaAvg(end)>0.999 && results.Tmelt==0)
        results.Tmelt=t;
    end
    %%% check priodicity in time
    
    if (rem(t,options.period) < dt/2)
        nperiod = round(t/options.period);
        if nperiod>1
            oldPriod = find(results.time >= (nperiod-2)*options.period & results.time < (nperiod-1)*options.period);
            newPriod = find(results.time >= (nperiod-1)*options.period & results.time < (nperiod)*options.period);
            gmeanOld=mean(results.gammaAvg(oldPriod));
            gmeanNew=mean(results.gammaAvg(newPriod));
            qNew=trapz(results.time(newPriod),abs(squeeze(sum(sum((results.Q(2:Nx+1,2:Ny+1,newPriod)),1),2))))/(pi*options.Lx^2*(options.Ly/(2*pi)))/(max(results.time(newPriod))-min(results.time(newPriod)));
            qOld=trapz(results.time(oldPriod),abs(squeeze(sum(sum((results.Q(2:Nx+1,2:Ny+1,oldPriod)),1),2))))/(pi*options.Lx^2*(options.Ly/(2*pi)))/(max(results.time(oldPriod))-min(results.time(oldPriod)));
% % %             fprintf('Bi=%2.2f,Ste=%2.2f,fraction=%2.2f,Nperiods=%d \t e = %2.5f,gammaAvg=%2.5f,dQ=%2.5f  \n',Bi,Ste,fraction,nperiod,abs(gmeanOld-gmeanNew),gmeanNew,(qNew-qOld))
% % %             %fprintf('Bi=%2.2f,Ste=%2.2f,fraction=%2.2f,Nperiods=%d \t e = %2.5f,gammaAvg=%2.5f,Q=%2.5f  \n',Bi,Ste,fraction,nperiod,abs(gmeanOld-gmeanNew),gmeanNew,qNew)
% % %             if (abs(gmeanNew-0.5)<0.005 )%(abs(gmeanOld-gmeanNew) <0.00001 && nperiod>=10)
% % %                 nperiodConv = nperiod;
% % %             end
% % %             if (nperiod > nperiodConv +5)
% % %                 fprintf('Bi=%2.2f,Ste=%2.2f,fraction=%2.2f,Nperiods=%d \t e = %2.5f,gammaAvg=%2.5f,dQ=%2.5f  \n',Bi,Ste,fraction,nperiod,abs(gmeanOld-gmeanNew),gmeanNew,abs(qNew-qOld))
% % %                 break;
% % %             end
        end
    end
%     hcT.ZData=calcBoundaryPoints(T');
% %     figure(Gammaplot)
% %     [~,h]=contourf(X,Y,results.gamma(:,:,end),50);
% %     axis equal
% %     set(h,'LineColor','none')
% %     % colormap('jet')
% %     colorbar
% %     caxis([0 1])
% %     figure(Tplot)
% %     [~,h]=contourf(X,Y,calcBoundaryPoints(T,'average'),50);
% %     axis equal
% %     set(h,'LineColor','none')
% %     colormap('jet')
% %     colorbar
% %     caxis([Air.Tlow Air.Thigh])
% %     drawnow
end
nperiodFinal = round(t/options.period);
% results.QavePerVolume = trapz(results.time,squeeze(sum(sum(abs(results.Q(2:Nx+1,2:Ny+1,:)),1),2)))/(pi*Lx^2*(Ly/(2*pi)))/(options.t1/options.period*2);
%% calc averages
q.tot=abs(squeeze(sum(sum((results.Q(2:Nx+1,2:Ny+1,:)),1),2)));  % Watts
q.lat=abs(squeeze(sum(sum((results.S(2:Nx+1,2:Ny+1,:)),1),2)));  % J/s , Watts
% Tsurf = squeeze(results.T(end,2:Ny+1,:));
q.sen=q.tot-q.lat;
 %% averaging range
Tmin=nperiodFinal -5;
Tmax=nperiodFinal;
indices = find(results.time>=Tmin*options.period & results.time<Tmax*options.period);
Qtot = trapz(results.time(indices),q.tot(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)))/((Tmax-Tmin)*2);
Qlat = trapz(results.time(indices),q.lat(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)))/((Tmax-Tmin)*2);
qtot = Qtot /(fraction*meltTime);
qlat = Qlat /(fraction*meltTime);
% Qtot,Qlat,Ste,Bi,Period,MeltingTime,Fraction,dt
stringOut=sprintf('%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f\n',...
                        Qtot,Qlat,qtot,qlat,Ste,Bi,meltTime,fraction,Nperiods,Tmin,Tmax,dt);
    
fprintf(stringOut)
fprintf('Finished Running Ste=%2.2f, Bi=%2.2f',Ste, Bi)


