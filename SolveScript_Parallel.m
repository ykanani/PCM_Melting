clear all
% close all
%%% PCM PCMerties
PCM.name='none';
% PCM.c=2.1e3;                % specific heat of the meterial
% PCM.cs=4e3;                % specific heat of the meterial
SteArray=fliplr(logspace(-2,0,25));
BiArray=fliplr(logspace(-1,1,50));
options.SteArray=SteArray(1:3:end);
options.BiArray=[BiArray(1) BiArray(2:2:end)];
[Ste, Bi]=meshgrid(options.SteArray,options.BiArray);

options.t0=0;

options.Ny=1;  %theta
options.Ly=pi;
options.Lx=1;
options.thresh=1e-6;
options.Tlow=0.5; 
options.Thigh=-.5; % it is always 1
options.Tinit=0;

a1 =      0.9578;
a2 =       0.485;
b =      -1.007;
b1 =      -1.006;
b2 =      -1.003;
c1 =       3.417;
c2 =      0.9678;


options.period=[];%40; % this is Fo number
options.t1=[];



options.Nx=40;  %radial
options.dt=0.01;
options.e=0.01;

datafile=['PhaseChangeNonDimQ-Nx' num2str(options.Nx)  'e' num2str(options.e) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];

fractions=[1 2 3 4 5 6 7 8 9 10 12 15 30]/10;
config=repmat([Ste(:),Bi(:)],length(fractions),1);
fraction=repelem(fractions,length(Ste(:)))';
% fraction=repelem(3,length(Ste(:)))';

meltingTime = (a1*config(:,1).^b1+c1).*config(:,2).^b+(a2*config(:,1).^b2+c2);
config=[config meltingTime fraction]; %% for fully melting period = 2 * meltTime
N=size(config,1);


Nworkers=15;
configWorker=cell(Nworkers,1);
c=1;
for i=1:N
    configWorker{c}(end+1,:)=config(i,:);
    c=c+1;
    if c>Nworkers
        c=1;
    end
end

WaitMessage = parfor_wait(N, 'Waitbar', true); 
spmd
    fileID = fopen(['CPU_' num2str(labindex) '_Date' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.txt'],'w');
    results=cell(size(configWorker{labindex},1),1);
    for i=1:size(configWorker{labindex},1)
        config=configWorker{labindex}(i,:);
        [string,~]=solve(options,config(1),config(2),config(3),config(4),50);
        fprintf(fileID,string);
        WaitMessage.Send
    end
    fclose(fileID);
end
WaitMessage.Destroy
% dsave(datafile)



