function [H, Htot, gamma]=calcState(T,PCM,ifplot)
Nx=length(T)-2;

Hsen=interp1(PCM.T,PCM.Hsen,T,'linear','extrap');
H=interp1(PCM.T,PCM.Hlat,T,'linear','extrap'); %Hlatent
Htot=H+Hsen;

gamma=H;

Nx=size(T,1);
Ny=size(T,2);


if (ifplot==1)  
    plot(PCM.T,PCM.Htot,'k')
    plot(T(2:Nx+1),Htot(2:Nx+1),'or')
    plot(PCM.T,PCM.H,'-.k')
    plot(T(2:Nx+1),H(2:Nx+1),'ob')
end