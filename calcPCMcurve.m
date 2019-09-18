function PCM=calcPCMcurve(e,Ste)
%% ALL H values are H/L enthalpy/latent heat of fusion
% PCM.cs=3000;
% PCM.cl=2100;
T=-2:e/10:2;

Tm=0;
Ts=-e;
Tl=e;
gamma=zeros(size(T));


gamma(T<Ts)=0;
gamma(T>Tl)=1;
gamma((T>=Ts) & (T<=Tl))=(T((T>=Ts) & (T<=Tl))-Ts)/(Tl-Ts);

Hlat=gamma;



Hsen=cumtrapz(T,Ste*ones(1,length(T)));

Htot=Hsen+Hlat;

% 
% figure
% plot(T,Htot,'--k')
% hold on
% plot(T,Hlat,'--r')

PCM.T=T;
PCM.Hlat=Hlat;
PCM.Htot=Htot;
PCM.Hsen=Hsen;
PCM.gamma=gamma;