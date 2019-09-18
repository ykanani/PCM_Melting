 function kapilow=readKapilow(t)

path='kapilow/';
tmelt=t(1)+8;
tsolid=t(2);
temp=csvread([path 'F4-Melting-Gamma-8.8.csv']);
temptm=temp(:,1);
tempGm=temp(:,2);

temp=csvread([path 'F4-Freezing-Gamma-8.8.csv']);
temptf=temp(:,1);
tempGf=temp(:,2);

kapilow.gamma=[tempGm; tempGf];
kapilow.gammat=[temptm+tmelt ;temptf+tsolid];
kapilow.gammaError=0.8338623634939601-0.7872870210282072;

temp=csvread([path 'kapilow-v9-UvsTime-Melting.csv']);
temptm=temp(:,1);
tempUm=temp(:,2);

temp=csvread([path 'kapilow-v9-UvsTime-Solid.csv']);
tempts=temp(:,1);
tempUs=temp(:,2);

kapilow.U=[tempUm; tempUs];
kapilow.Ut=[temptm+tmelt; tempts+tsolid];
kapilow.UError=214.9207145444918-239.9557822181349;


temp=csvread([path 'AirT.csv']);
kapilow.T=temp(:,2);
kapilow.Tt=temp(:,1);
