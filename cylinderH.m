function h=cylinderH(D,Air,angles)

ReD=D*Air.V/Air.nu;
Pr=Air.Pr;
Nu=NuDistribute(ReD,angles/pi*180);
%Nu=ones(1,length(Nu))*mean(Nu); %uniform distribution
%Churchill and Bernstein
NuAvg=0.3+(0.62*ReD^0.5*Pr^(1/3))/(1+(0.4*Pr)^(2/3))^(1/4)*(1+(ReD/282000)^(5/8))^(4/5);
%Hilpert
if (ReD < 4)
    C=0.989; m=0.330;
elseif (ReD < 40)
    C=0.911; m=0.385;
elseif (ReD < 4000)
    C=0.683; m=0.466;
elseif (ReD < 40000)
    C=0.193; m=0.618;
else 
    C=0.027; m=0.805;
end
NuAvg=C*ReD^m*Pr^(1/3);


h=Nu*Air.k/D;