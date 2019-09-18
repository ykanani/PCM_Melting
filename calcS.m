function [S, T, dHdT0, dHdT1]=calcS(T0,Q,PCM,ifplot)
%%% Htot    Total Enthalpy/L (Ste*T+H)/L
%%% H    Latent Enthalpy/L or gamma
%%% find the Total Enthalpy associated with T0

Nx=size(T0,1);
Ny=size(T0,2);
dhdt=[0 diff(PCM.Hlat)./diff(PCM.T)];
% dhdt=zeros(size(PCM.T));
% dhdt(PCM.T==0)=1e8;
H0=interp1(PCM.T,PCM.Hlat,T0,'linear','extrap'); %Hlatent
dHdT0=interp1(PCM.T,dhdt,T0,'linear','extrap');

Hsen0=interp1(PCM.T,PCM.Hsen,T0,'linear','extrap');
Htot0=H0+Hsen0;

Htot1=Htot0+Q;

% myfun=@(x,h) h-(a*erf((x-b)/c)+d)*PCM.hfg-interp1(PCM.T,PCM.Hsen,x);
% T1=zeros(size(T0));
% for i=1:length(T0(:))
%     T1(i)=fzero(@(x) myfun(x,Htot1(i)),T0(i));
% end

T1=interp1(PCM.Htot,PCM.T,Htot1,'linear','extrap');


% T1=(erfinv(((Htot1/240000-d)/a))*c+b;
H1=interp1(PCM.T,PCM.Hlat,T1,'linear','extrap'); %Hlatent
dHdT1=interp1(PCM.T,dhdt,T1,'linear','extrap'); 

gamma0=H0;
gamma1=H1;



S=gamma0-gamma1;   % the coef matrix already has 1/Ste , this is gamma0-gamma1


T=T0;
T(abs(S)>0.001)=T1(abs(S)>0.001);



if (ifplot==1)
    
    
    plot(PCM.T,PCM.Htot,'k')
    hold on
    plot(T0(2:Nx-1,2:Ny-1),Htot0(2:Nx-1,2:Ny-1),'*r')
    plot(T1(2:Nx-1,2:Ny-1),Htot1(2:Nx-1,2:Ny-1),'or')

    plot(PCM.T,PCM.Hlat,'-.k')
    plot(T0(2:Nx-1,2:Ny-1),H0(2:Nx-1,2:Ny-1),'*b')
    plot(T1(2:Nx-1,2:Ny-1),H1(2:Nx-1,2:Ny-1),'ob')
end
%%% fitting curve (not used)
% ft=fittype('(a*erf(b*(x-c)/(d-e))+f)', 'independent', 'x', 'dependent', 'y' );
% f = fit(Tref',H'/1e6,ft,'StartPoint',[1 1 mean(Tref) max(Tref) min(Tref) 1])