clear all
% close all
%%% PCM PCMerties
PCM.name='none';
% PCM.c=2.1e3;                % specific heat of the meterial
% PCM.cs=4e3;                % specific heat of the meterial
Ste=0.0682;
Bi=1.428943793553002;


options.t0=0;

options.Ny=1;  %theta
options.Ly=pi;
options.Lx=1;
options.thresh=1e-6;
options.Tlow=0.5; 
options.Thigh=-.5; % it is always 1
options.Tinit=-0.5;
options.period=40; % this is Fo number
options.t1=0.5*options.period;

e=0.01;%(tl-ts)/2
[PCM.T, PCM.Hlat, PCM.Htot, PCM.Hsen, PCM.gamma]=calcPCMcurve(e,Ste);
options.Nx=40;  %radial
options.dt=0.05;
datafile=['PhaseChangeNonDim-Nx' num2str(options.Nx) '-dt' num2str(options.dt) '-Bi' num2str(Bi) '-Ste' num2str(Ste) 'T' num2str(options.period) 'e' num2str(e) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];
results=solve(options,Ste,Bi,PCM);
save(datafile)
results.Tmelt
datafiles{1}=datafile

options.dt=0.01;
datafile=['PhaseChangeNonDim-Nx' num2str(options.Nx) '-dt' num2str(options.dt) '-Bi' num2str(Bi) '-Ste' num2str(Ste) 'T' num2str(options.period) 'e' num2str(e) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];
results=solve(options,Ste,Bi,PCM);
save(datafile)
results.Tmelt
datafiles{2}=datafile





postprocess(datafiles,options,0)
