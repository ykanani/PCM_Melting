clear all
% close all
%%% PCM PCMerties
SteArray=fliplr(logspace(-2,0,25));
BiArray=fliplr(logspace(-1,1,50));
options.SteArray=[0.1];
options.BiArray=[1];
% options.SteArray=logspace(-2,0,50);%[0.1 0.2];%
% options.BiArray=logspace(-1,1,100);%[0.5 1 ];%
[Ste, Bi]=meshgrid(options.SteArray,options.BiArray);

options.t0=0;

options.Ny=1;  %theta
options.Ly=pi;
options.Lx=1;
options.thresh=1e-6;
options.Tlow=0.5; 
options.Thigh=-.5; % it is always 1
options.Tinit=0;

Tmelt = @(Bi,Ste) (Ste.*Bi+1+3.5*Ste+0.5*Bi)./Ste./Bi;

options.period=[];%40; % this is Fo number
options.t1=[];



options.Nx=40;  %radial
options.dt=0.01;
options.e=0.05;

datafile=['Ste-Bi-Effect-PhaseChangeNonDimQ-Nx' num2str(options.Nx)  'e' num2str(options.e) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];

fractions=[0.1 0.5 1];
config=repmat([Ste(:),Bi(:)],length(fractions),1);
fraction=repelem(fractions,length(Ste(:)))';
% fraction=repelem(3,length(Ste(:)))';

meltingTime = Tmelt(Bi,Ste);
config=[config meltingTime*ones(length(fraction),1) fraction]; %% for fully melting period = 2 * meltTime
N=size(config,1);

for i=1:N
    [string{i}, results{i}]=solve(options,config(i,1),config(i,2),config(i,3),config(i,4),200);
end

save(datafile)
% datafiles{1}=datafile
% postprocess(datafiles,options,0)
