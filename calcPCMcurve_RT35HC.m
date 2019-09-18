function [T, Hlat, Htot, Hsen, PCMc, cs, cl, fitresult]=calcPCMcurve(PCM)

% PCM.cs=3000;
% PCM.cl=2100;
Tref=[PCM.dHTref-0.5 PCM.dHTref(end)+0.5];
HtotrefHot=[0 cumsum(PCM.dhh)];%+PCM.cs*Tref(1);
HtotrefCool=[0 cumsum(PCM.dhc)];
Htotref=(HtotrefHot+HtotrefCool)/2;
ps=polyfit(Tref(1:7),Htotref(1:7),1);
cs=ps(1);
pl=polyfit(Tref(end-5:end),Htotref(end-5:end),1);
cl=pl(1);
p3=polyfit([Tref(8) Tref(end-6)], [cs cl],1);
Htotref=Htotref+cs*Tref(1);
HtotrefHot=HtotrefHot+cs*Tref(1);
HtotrefCool=HtotrefCool+cs*Tref(1);
cref=zeros(1,length(Tref));
cref(1:7)=cs;
cref(8:end-6)=p3(1)*Tref(8:end-6)+p3(2);
cref(end-7:end)=cl;
Hsenref=cumsum(cref.*[Tref(1) diff(Tref)]);
Hpartial=Htotref-Hsenref;

% PCM.hfg=267510;
% Hpartial=(Htotref-cs*Tref)./(1+(cl-cs)/PCM.hfg*Tref);



[xData, yData] = prepareCurveData( Tref, Hpartial/PCM.hfg );

% Set up fittype and options.
ft = fittype( 'a*erf((x-b)/c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.5 35 0.4505 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% a =0.4694;  
% b =35.09;  
% c =1.286;  
% d =0.4754;
a =fitresult.a;  
b =fitresult.b;  
c =fitresult.c;  
d =fitresult.d;
T=0:0.001:60;
Hlat=(a*erf((T-b)/c)+d)*PCM.hfg;
rho=(PCM.rhol+PCM.rhos)/2;
gamma=Hlat./max(Hlat);
PCMc=gamma*cl+(1-gamma)*cs;
Hsen=cumsum(PCMc.*[0 diff(T)]);
% Htot=H+cumsum(PCMc.*[0 diff(T)]);
Htot=Hlat+Hsen; 
figure('DefaultAxesFontSize',15)
plot(Tref,HtotrefHot,'ok');%,'MarkerFaceColor','k')
hold on
plot(Tref,HtotrefCool,'^k');%,'MarkerFaceColor','k')
plot(Tref,Htotref,'*k');%,'MarkerFaceColor','k')
hold on
plot(Tref,Hpartial,'o')
plot(T,Htot,'--k')
plot(T,Hlat,'--r')
%%
figure('DefaultAxesFontSize',15)
% plot(Tref,HtotrefHot,'or','MarkerFaceColor','r','MarkerSize',3)
% hold on
% plot(Tref,HtotrefCool,'ob','MarkerFaceColor','b','MarkerSize',3)
plot(Tref,Htotref,'ok','LineWidth',1.5);%,'MarkerFaceColor',[.2 .3 1])
hold on
xlabel('Temperature ({\circ}C)')
ylabel('h-h_{ref} (J/kg)')
plot(Tref,Hpartial,'xk','LineWidth',1.5);%,'MarkerFaceColor','k')
plot(Tref,Hsenref,'>k','LineWidth',1.5);%,'MarkerFaceColor','k')
plot(T,Htot,'k')
plot(T,Hlat,'--k')
plot(T,Hsen,'-.k')
axis([25 43 -inf inf])
%'h_{tot,hot}' 'h_{tot,cold}'
legend({ 'h_{tot}' 'h_{lat}' 'h_{sen}' 'h_{tot}' 'h_{lat}' 'h_{sen}'},'location','best')