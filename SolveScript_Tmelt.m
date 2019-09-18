clear all
% close all
%%% PCM PCMerties
PCM.name='none';
% PCM.c=2.1e3;                % specific heat of the meterial
% PCM.cs=4e3;                % specific heat of the meterial
SteArray=fliplr(logspace(-2,0,25));
BiArray=fliplr(logspace(-1,1,50));
[Ste, Bi]=meshgrid(SteArray,BiArray);

options.t0=0;

options.Ny=1;  %theta
options.Ly=pi;
options.Lx=1;
options.thresh=1e-7;
options.Tlow=0.5; 
options.Thigh=-.5; % it is always 1
options.Tinit=-0.5;
options.period=10000; % this is Fo number
options.t1=0.5*options.period;

e=0.01;%(tl-ts)/2

options.Nx=40;  %radial
options.dt=0.01;

datafile=['PhaseChangeNonDim-Nx' num2str(options.Nx) '-dt' num2str(options.dt) 'e' num2str(e) '-date-' datestr(now,'mm-dd-yyyy_HH-MM-SS') '.mat'];
tmelt=zeros(1,length(Ste(:)));
Tavg=zeros(1,length(Ste(:)));
tic
parfor i=1:length(Ste(:))
    [tmelt(i), Tavg(i)]=solveMinimal(options,Ste(i),Bi(i),calcPCMcurve(e,Ste(i)));
end
meltingTime=reshape(tmelt,size(Ste));
TmeltContour=figure('DefaultAxesFontSize',15);
[~,h]=contourf(log10(Ste),log10(Bi),meltingTime);

axis equal
set(h,'LineColor','none')
colormap(jet)
c = colorbar;
c.Label.String = 'Fo_{melt}';
c.Label.FontSize=14;

Fo=meltingTime(:);
Ste1=Ste(:);
Bi1=Bi(:);
[fitresult, gof] = createFitBiSte(Bi1, Ste1, Fo)

fig1=figure('DefaultAxesFontSize',12);
% orginal % % FoPredicted = (fitresult.a1*Ste.^fitresult.b1+fitresult.c1).*Bi.^fitresult.b+(fitresult.a2*Ste.^fitresult.b2+fitresult.c2); %  Ste.^fitresult.a.*(fitresult.b+fitresult.c*Bi.^fitresult.d);
FoPredicted =1+(1+3.5*Ste+0.5*Bi)./(Ste.*Bi); %(1*Ste.^(-1)+3.5).*Bi.^(-1)+(0.5*Ste.^(-1)+1); %  Ste.^fitresult.a.*(fitresult.b+fitresult.c*Bi.^fitresult.d);

%h=surface(Ste,Bi,fitresult.a*Ste.^fitresult.b.*Bi.^fitresult.c+Ste.^fitresult.d);
h=surface(Ste,Bi,meltingTime,'DisplayName','Predicted');

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'ZScale', 'log')
axis([min(Ste1) max(Ste1) min(Bi1) max(Bi1) min(Fo) max(Fo)])
set(h,'EdgeColor','none')
view( 25.2, 47.6 );
hold on
xlabel('$Ste$','Interpreter','latex')
ylabel('$Bi$','Interpreter','latex')
zlabel('$\tilde{\mathcal{T}}_\mathrm{ref}$','Interpreter','latex')
plot3(Ste1,Bi1,FoPredicted(:),'.r','DisplayName','Eq. \ref{41}')
error=(abs(FoPredicted(:)-Fo)./Fo)*100;
max(error)
mean(error)
min(error)
c = colorbar;
colormap(jet)
grid('on')
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;
c.Label.String = '$\tilde{\mathcal{T}}_\mathrm{ref}$';
view([64 20])
legend('show','location','best')
ax=gca;
ax.MinorGridAlpha=0.8;
ax.MinorGridColor=[0 0 0];
ax.GridAlpha=0.8;
ax.GridColor=[0 0 0];

save(datafile)

