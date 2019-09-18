clear all
close all
%%% Simulation paramters
t0=0;
dt=1;                               %deltaT
Nx=20;                                  %number of cells
thresh=1e-4;
f=1;                                  %time discritization scheme, f=1 fully implicit, f=0 explicit


%%% PCM PCMerties
PCM.name='RT35HC';
PCM.c=2e3;                % specific heat of the meterial
PCM.rho=(7700+880)/2;      % Density of Meterial\
PCM.k=0.2;                % thermal conductivity
PCM.hfg=240e3;            % Heat of fusion
PCM.D=2.78e-3;            %PCM diameter  
PCM.dhh=[3 4 3 5 5 7 11 29 108 60 11 3 2 3 2 3]*1e3; % Partial Enthanlpy J/kg
PCM.dhc=[3 2 3 2 2 2 13 94 72 27 7 6 4 4 3 3]*1e3;
PCM.dHTref=27:42;
[PCM.T, PCM.H, PCM.Htot]=calcPCMcurve(PCM);
%%% External flow (air) parameters
Air.Tlow=25.3;
Air.Thigh=42.9;

Air.nu=1.5571e-5;                       % Kinematic Viscosity, m^2/s
Air.rho=1.1845;                         % Density, kg/m^3
Air.cp=1.0063E+3;                       % Specific Heat J/kg.K
Air.k=0.025969;                         % Conductivity W/m.K
Air.alfa=Air.k/(Air.cp*Air.rho);        % Thermal Diffusivity m^2/s
Air.Pr=Air.nu/Air.alfa;                 %Prandl number

                                          
%%% HDPE (tube) parameters
HDPE.D=3.18e-3;                            % HDPE Tube Diameter
HDPE.cp=2250;                           % Specific Heat J/kg.K
HDPE.rho=970;                           % Density, kg/m^3
HDPE.k=0.49;                            % Conductivity W/m.K

%%% reading kapilow's data
tmelt=10;
tfreeze=250;
kapilow=readKapilow([tmelt tfreeze]); 

%%% calc heat tranfer coeff at t=0;
Air.V=9;                                 %Air velocity m/s
Air.h=cylinderH(HDPE.D,Air);             % Air Heat Transfer Coefficient
%%% Overal heat transfer coefficient
Air.R =1/(HDPE.D/2*Air.h);              % Thermal Resistance of Air (K/W)
HDPE.R=log(HDPE.D/PCM.D)/HDPE.k;        % Thermal Resistance of HDPE (K/W)
Air.U=1/(Air.R+HDPE.R);                     % overall heat transfer coefficient, W/K

%Note: The area in already included in the U calculations, U only needs to 
%be multiplied by the temperature difference to provide Q 

%%% Creating Uniform Mesh
%%% *-----|-----o-----|-----o-----|-----o-----|-----*
%%% ?     ?      ?                       ?     ?     ?
%%% T1          T2          T3          T4          T5
%%%      x=0                                  x=L    
%%% Cx1        Cx2                     CxN+1       CxN+2
%%%      Fx1        xe,Re       xw,Rw              FxN+1 
%%%       ??Csize1????Csize2????Csize3??            (Dx) or (Dr)
%%% ? dx1,dr1 ?? dx2,dr2 ?? dx3,dr3 ?? dx4,dr4 ?     (dx) or (dr) 
L=PCM.D/2;                              %Lengh of domain (radius)
Mesh=createMeshCylindrical(L,Nx);                  %Creating Uniform Mesh
%%% Setting Boundary Conditions
% B.C. equation:
%(right)          A*(Tb+Tg)/2+B*(Tb-Tg)/dx = C,
%(left)           -A*(Tb+Tg)/2-B*(Tb-Tg)/dx = -C,
% Tb=T2 or TN+1 (Boundary point, cell center)
% Tg=ghoast cell

% %%%% Left Boundary (adiabatic, symmetry)
BC.left.A=0;
BC.left.B=1;
BC.left.C=0;

%%%% Right Boundary (wall+convective)
BC.right.A=Air.U;
BC.right.B=PCM.k*(PCM.D/2);               %k*A , because U is already multiplied by area
BC.right.C=Air.U*Air.Thigh;

[B, RHS] = BoundaryCondition(Mesh,BC);  % B is multiplied by T , change sign to negative for RHS

[Itran, Etran] = transientCoeffsCylindrical(Mesh,dt,PCM);
[IdiffImplicit, EdiffImplicit] = diffusionCoeffsCylindrical(Mesh,PCM,1);
[IdiffExplicit, EdiffExplicit] = diffusionCoeffsCylindrical(Mesh,PCM,0);
[Is, Es] = sourceCoeffsCylindrical(Mesh,dt,PCM);


I=Itran+IdiffImplicit+B;
E=Etran+EdiffImplicit;

%Initialize
Tinit=Air.Tlow;
T=[0; Tinit*ones(Nx,1); 0]; % two zeros corresponds to the ghoasted cells and will be corrected according to the B.C.
%Apply boundary condition
T=(B+sparse(2:Nx+1,2:Nx+1,1,Nx+2,Nx+2))\(RHS+T);

%%solve
figure(1)

%Temperature plot
subplot(2,2,1)
hold on
Tline=plot(Mesh.Cloc,T,'ok','MarkerFaceColor','b');
hold on
Tanaline=plot(Mesh.Cloc,T,'--k');
Tbc=plot(Mesh.Floc(end),(T(end-1)+T(end))/2,'ok','MarkerFaceColor','b');
xlabel('Radius')
ylabel('Temperature')
axis([Mesh.Floc(1) Mesh.Floc(end) Air.Tlow-1 Air.Thigh+1])
%phase fraction plot
subplot(2,2,3)
Gammaline=plot(Mesh.Cloc,ones(1,Nx+2),'ok','MarkerFaceColor','r');
xlabel('Radius')
ylabel('Phase Fraction')
text1=text(0.1,0.7,'$$\overline{a}$$', 'Units', 'normalized','FontSize',18,'Interpreter','latex');
axis([Mesh.Floc(1) Mesh.Floc(end) -0.05 1.05])
%Average phase fraction plot
subplot(2,2,2)
AvgGline=plot(0,0,'k','LineWidth',2);
hold on
plot(kapilow.gammat,kapilow.gamma,'^b');
xlabel('time (s)')
ylabel('Average Phase Fraction')
text2=text(0.1,0.7,'$$T_{air}$$', 'Units', 'normalized','FontSize',18,'Interpreter','latex');
axis([0 600 -0.05 1.05]);


%Overall heat transfer Coef
subplot(2,2,4)
Uline=plot(0,0,'k','LineWidth',2);
hold on
plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');
xlabel('time (s)')
ylabel('Average Phase Fraction')
axis([0 600 0 250])

gammaAvg=[];
Uplot=[];
% figure(2)
% hold on   
for t=t0+dt:dt:600
%     disp(t)
    %%% update Boundary Conditions
    % %%%% Left Boundary (adiabatic, symmetry)
    BC.left.A=0;
    BC.left.B=1;
    BC.left.C=0;
    %%%% Right Boundary (wall+convective)

    Air.V=8.8;
    if (t<tmelt)
        Air.T=Air.Tlow;
    elseif (t<tfreeze)
        Air.T=Air.Thigh;
    else
        Air.T=Air.Tlow;
    end
    Air.h=cylinderH(PCM.D,Air);             % Air Heat Transfer Coefficient
% % %     Air.h=100;
    BC.right.A=Air.U;
    BC.right.B=PCM.k*(PCM.D/2);               %k*A , because U is already multiplied by area
    BC.right.C=Air.U*Air.T;
    [B, RHS] = BoundaryCondition(Mesh,BC);  % B is multiplied by T , change sign to negative for RHS
    %%% update coef matrix
    I=Itran+IdiffImplicit+B;
    E=Etran+EdiffImplicit;
    %%% gauss-sidel loop
    
    S=[];
    Q=[];
    Tint=[];
    
%     % Initial Estimate of T
%     Tint=I\(E*T+RHS);
    % Initial Estimate of Q (fully explicit)
    Q=dt*EdiffExplicit*T./(Mesh.V*PCM.rho);
    % Initial Estimate of S
    [S, Tint]= calcS(T,Q,PCM,0);
    e=1;

%     Tint=T';

    while (e>thresh)
      TintOld=Tint;
      %%% my method
      Tint=I\(E*T+RHS+Es*S); %implicit solver
      %%% voller's method
%       G=IdiffImplicit-eye(Nx+2,Nx+2).*IdiffImplicit;
%       for k=2:Nx+1
%           Tint(k)=(Etran(k,:)*T-G(k,:)*Tint+Es(k,:)*S)/(IdiffImplicit(k,k)+Etran(k,k));
%       end
%       %Apply boundary condition
%       Tint=(B+sparse(2:Nx+1,2:Nx+1,1,Nx+2,Nx+2))\(RHS+[0; Tint(2:Nx+1); 0]);
      e=mean(abs(Tint-TintOld));
      % Estimate Q (fully explicit)
      Q=dt*EdiffExplicit*Tint./(Mesh.V*PCM.rho);
      % update S
      S= calcS(T,Q,PCM,0);
    end
%     Tint=I\(E*T+RHS+Es*S); %implicit solver
    T=Tint;
%     [Qtot Qlat]=CalcQ(Mesh,T,PCM);
    
    q.tot=sum(EdiffExplicit*T);  % Watts
    q.lat=sum(Es*S);  % J/s , Watts
    q.sen=q.tot+q.lat;
    q.surf=PCM.D/2*PCM.k*(T(end)-T(end-1))/Mesh.Csize(end);
       
    disp(struct2table(q))
    
    U.tot=q.tot/abs(35-Air.T)/(HDPE.D/2);
    U.lat=q.lat/abs(35-Air.T)/(HDPE.D/2);
    U.sen=q.sen/abs(35-Air.T)/(HDPE.D/2);
    Uplot(end+1)=abs(U.lat);
    disp(struct2table(U))
%     figure(2)
    [H, Htot, gamma]= calcState(T,PCM,0);
    gammaAvg(end+1)=sum(gamma(2:Nx+1).*Mesh.V(2:Nx+1))/sum(Mesh.V(2:Nx+1));
    if (true)%rem(t,10)==0)
        set(Tline,'YData',T)
        set(Tbc,'YData',(T(end-1)+T(end))/2);
        set(Gammaline,'YData',gamma)
        set(AvgGline,'YData',gammaAvg,'XData',t0+dt:dt:t)
        set(text1,'String',['$$\overline{\gamma}=$$ ' num2str(round(gammaAvg(end),2))  ]);%  '%']);
        set(text2,'String',['$$t=$$' num2str(t,2) char(10) '$$T_{air}=$$'   num2str(Air.T,2) char(10) '$$V_{air}=$$'   num2str(Air.V,2)  ]);%  '%']);
        set(Uline,'YData',Uplot,'XData',t0+dt:dt:t);
    
    end
    
    
    
% %     %exact solution (only sensible)
% %     Thetaana=zeros(Nx+2,1);
% %     Bi=(Air.U/(PCM.D/2))*L/PCM.k;
% %     lambda=lambdan([-1 2000],Bi);
% %     for n=1:length(lambda)
% %         %An=(-2*(-1)^n*Bi*sqrt(lambda(n)^2+Bi^2))/(lambda(n)*(lambda(n)^2+Bi*(1+Bi)));
% %         kn=lambda(n);
% %         An=(2/kn)*(besselj(1,kn))/(besselj(0,kn)^2+besselj(1,kn)^2);
% %         Fo=(t-t0)*PCM.k/(PCM.rho*PCM.c)/L^2;
% %         Thetaana=Thetaana+An*besselj(0,kn*Mesh.Cloc/L)*exp(-kn^2*Fo);
% %     end
% %     Tana=Thetaana*(Tinit-Air.T)+Air.T;
% %     set(Tanaline,'YData',Tana)
    
    
    refreshdata
    drawnow

end


