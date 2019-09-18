clear all
close all
%%% PCM PCMerties
PCM.name='RT35HC';
% PCM.c=2.1e3;                % specific heat of the meterial
% PCM.cs=4e3;                % specific heat of the meterial

PCM.rhol=770;     % Density of Meterial\
PCM.rhos=880;      % Density of Meterial\
PCM.kl=0.14;                % thermal conductivity
PCM.ks=0.36;                % thermal conductivity
PCM.hfg=220e3;            % Heat of fusion
PCM.D=2.78e-3;            % PCM diameter  
PCM.dhh=[3 4 3 5 5 7 11 29 108 60 11 3 2 3 2 3]*1e3; % Partial Enthanlpy J/kg
PCM.dhc=[3 2 3 2 2 2 13 94 72 27 7 6 4 4 3 3]*1e3;
PCM.dHTref=27:42;
[PCM.T, PCM.Hlat, PCM.Htot, PCM.Hsen, PCM.cref, PCM.cs, PCM.cl, PCM.Hfit]=calcPCMcurve(PCM);
%%% External flow (air) parameters
Air.Tlow=25.3;
Air.Thigh=42.9;

Air.nu=1.5571e-5;                       % Kinematic Viscosity, m^2/s
Air.rho=1.1845;                         % Density, kg/m^3
Air.cp=1.0063E+3;                       % Specific Heat J/kg.K
Air.k=0.025969;                         % Conductivity W/m.K
Air.alfa=Air.k/(Air.cp*Air.rho);        % Thermal Diffusivity m^2/s
Air.Pr=Air.nu/Air.alfa;                 %Prandl number
Air.V=8.8;                                 %Air velocity m/s
                                          
%%% HDPE (tube) parameters
HDPE.D=3.18e-3;                            % HDPE Tube Diameter
HDPE.cp=2250;                           % Specific Heat J/kg.K
HDPE.rho=970;                           % Density, kg/m^3
HDPE.k=0.49;                            % Conductivity W/m.K


options.t0=0;
options.t1=800;
options.dt=1;
options.Nx=20;  %radial
options.Ny=30;  %theta
options.Lx=PCM.D/2; 
options.Ly=pi;
options.period=400;
options.thresh=1e-6;
datafile=['PhaseChange-Nx' num2str(options.Nx) '-dt' num2str(options.dt) '-k' num2str(PCM.kl) '-' num2str(PCM.ks) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];
results=solve(options,PCM,Air,HDPE);
save(datafile)




