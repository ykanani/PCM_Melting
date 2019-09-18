function Nu=NuDistribute(Re,angle)
Nu500=csvread('500.csv');
Nu3840=csvread('3840.csv');
Nu8290=csvread('8290.csv');
Nu10000=csvread('10000.csv');
N500=sortrows(Nu500,1);
N3840=sortrows(Nu3840,1);
N8290=sortrows(Nu8290,1);
N10000=sortrows(Nu10000,1);


NN500=interp1(N500(:,1),N500(:,2),0:0.1:180,'linear','extrap');
NN3840=interp1(N3840(:,1),N3840(:,2),0:0.1:180,'linear','extrap');
NN8290=interp1(N8290(:,1),N8290(:,2),0:0.1:180,'linear','extrap');
NN10000=interp1(N10000(:,1),N10000(:,2),0:0.1:180,'linear','extrap');
figure
plot(0:0.1:180,NN500,'r')
hold on
plot(0:0.1:180,NN3840,'b')
plot(0:0.1:180,NN8290,'b')
plot(0:0.1:180,NN10000,'b')

% plot(N8290(:,1),N8290(:,2))
% plot(N10000(:,1),N10000(:,2))
legend show
Nu=[];

for i=1:length(angle)
    Nu0=interp1(0:0.1:180,NN500,angle(i));
    Nu1=interp1(0:0.1:180,NN3840,angle(i));
    Nu(i)=interp1([500 3840],[Nu0 Nu1],Re,'linear','extrap');
end

plot(angle,Nu,'k')

