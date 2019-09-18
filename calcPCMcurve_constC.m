function [T, H, Htot]=calcPCMcurve(PCM)
Tref=[PCM.dHTref-0.5 PCM.dHTref(end)+0.5];
Hpartial=PCM.dhh-PCM.c*diff(Tref);            % partial latent heat J/kg.K
Href=[0 cumsum(Hpartial)];                       % latent heat content J/kg.K
%%% adding lower and higher bound temperatures
n=50;
Tref=[fliplr(Tref(1)-1:-1:Tref(1)-n) Tref Tref(end)+0.5:Tref(end)+n];
Href=[zeros(1,n) Href ones(1,n)*Href(end)];

Htotref=PCM.c*Tref+Href;

T=0:0.001:60;
Htot=interp1(Tref,Htotref,T,'spline');
H=interp1(Tref,Href,T,'spline');