% function data=readData(folder)
clear all
close all
folder = '4thRun'
list=dir([folder '\*.txt']);
data=[];
for i=1:length(list)
    data=[data ; table2array(readtable([folder '\' list(i).name])) ];
end
% removing Bi<1
% indices=find(data(:,3)<0.1 | data(:,4)<2);
% data(indices,:)=[];
Qtot=data(:,1);
Qlat=data(:,2);
qtot=data(:,3);
qlat=data(:,4);
Ste=data(:,5);
Bi=data(:,6);
meltTime=data(:,7);
fraction=data(:,8);
QtotByQmax=Qtot./(1+Ste);
qtotByQmax=qtot./(1+Ste);


% % % % % % % % %%fractions > 1
% % % % % % % % ifone=find(fraction>0.95);
% % % % % % % % Qtoti=Qtot(ifone);
% % % % % % % % Qmaxi=1+Ste(ifone);
% % % % % % % % error=abs(Qtoti-Qmaxi)./Qmaxi*100;
% % % % % % % % max(error)
% % % % % % % % min(error)
% % % % % % % % 
% % % % % % % % iminBiSte=find(Ste==min(Ste) & Bi==min(Bi));
% % % % % % % % figure
% % % % % % % % plot(fraction(iminBiSte),QtotByQmax(iminBiSte));
% % % % % % % % hold on
% % % % % % % % imaxBiSte=find(Ste==max(Ste) & Bi==max(Bi));
% % % % % % % % plot(fraction(imaxBiSte),QtotByQmax(imaxBiSte),'--');
% % % % % % % % axis([0 1 0 1.5])
% % % % % % % % 
% % % % % % % % figure()
% % % % % % % % plot(fraction,QtotByQmax,'.')
% % % % % % % % 
% % % % % % % % figure()
% % % % % % % % ifrac=find(fraction==0.5 & Bi>1.676 & Bi<1.677);
% % % % % % % % plot3(Bi(ifrac),Ste(ifrac),QtotByQmax(ifrac),'r*')
% % % % % % % % xlabel('Bi')
% % % % % % % % ylabel('Ste')
% % % % % % % % Bi1=Bi(ifrac);
% % % % % % % % Ste1=Ste(ifrac);
% % % % % % % % QtotbyQmax1=QtotByQmax(ifrac);
% % % % % % % % hold on
% % % % % % % % ifrac=find(fraction==0.5);
% % % % % % % % plot3(Bi(ifrac),Ste(ifrac),QtotByQmax(ifrac),'b*')
% % % % % % % % xlabel('Bi')
% % % % % % % % ylabel('Ste')
% % % % % % % % 
% % % % % % % % 
a=[];
b=[];
c=[];
Stearray=[];
Biarray=[];
SteValues = unique(Ste);
fracValues = unique(fraction);
BiValues = unique(Bi);
% % % % % for i=1:length(SteValues)
% % % % %     for j=1:length(BiValues)
% % % % %         indices = find(Ste==SteValues(i) & Bi==BiValues(j));
% % % % %         if isempty(indices)
% % % % %             break
% % % % %         end
% % % % %         fracRef=fraction(indices);
% % % % %         QtotbyQmaxRef=QtotByQmax(indices);
% % % % %         
% % % % %         
% % % % %         [fitresult, gof] = createFitFraction(fracRef, QtotbyQmaxRef);
% % % % %         a(end+1)=fitresult.a;
% % % % %         b(end+1)=fitresult.b;
% % % % %         c(end+1)=fitresult.c;
% % % % %         Biarray(end+1)=BiValues(j);
% % % % %         Stearray(end+1)=SteValues(i);
% % % % %     end
% % % % % end

% (1-b(1)*exp((b(2)*Bi^b(3)*Ste^b(4))*x^(b(5)*Bi^b(6)*Ste^b(7))))


% fun = @(b,X) (1-b(1).*exp(-(b(2).*X(:,1).^b(3).*X(:,2).^b(4)).*X(:,3).^(b(5).*X(:,1).^b(6).*X(:,2).^b(7))));%                  b(1).*X(:,1).^b(2).*X(:,2).^b(3).*X(:,3).^b(4)
% fun = @(b,X) (1-b(1).*exp(b(2).*X(:,3).^b(3)));% 
b0=ones(1,7);
options.MaxIter=10000;
options.RobustWgtFun='bisquare';
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit([Bi Ste fraction],QtotByQmax,@fitmodel,b0,options);

QtotFitted=fitmodel(beta,[Bi Ste fraction]);
SSE=sum(R.^2);
SSTO=sum((QtotByQmax-mean(QtotByQmax)).^2);
R2=1-SSE/SSTO;
max((QtotFitted-QtotByQmax))
figure()
plot3(Bi,fraction,QtotByQmax,'bo')
hold on
plot3(Bi,fraction,QtotFitted,'r*')
xlabel('Bi')
ylabel('fraction')
error=abs(QtotFitted-QtotByQmax)./QtotByQmax*100;
max(error)
mean(error)
figure
plot(error,'*')

[beta2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@fitmodel,b0,[Bi Ste fraction],QtotByQmax);
QtotFitted2=fitmodel(beta2,[Bi Ste fraction]);
SSE=sum(R.^2);
SSTO=sum((QtotByQmax-mean(QtotByQmax)).^2);
R2=1-SSE/SSTO;
max((QtotFitted2-QtotByQmax))
figure()
plot3(Bi,fraction,QtotByQmax,'bo')
hold on
plot3(Bi,fraction,QtotFitted2,'r*')
xlabel('Bi')
ylabel('fraction')
error=abs(QtotFitted2-QtotByQmax)./QtotByQmax*100;
max(error)
mean(error)
% % % 
% % % 
%%
% % % % % param = QtotByQmax;
% % % % % % X=[ones(size(Bi)) Bi Ste fraction Bi.*Ste Bi.*fraction Ste.*fraction Bi.^2 Ste.^2 fraction.^2 (Bi.*Ste).^2 (Bi.*fraction).^2 (Ste.*fraction).^2];
% % % % % X=[ones(size(Bi)) Bi Bi.^2 Bi.^3 fraction fraction.^2 fraction.^3 Ste Ste.^2 Ste.^3 fraction.*Ste fraction.^2.*Ste fraction.*Ste.^2 fraction.^2.*Ste.^2 ];
% % % % % b=regress(param,X)
% % % % % % res=b(1)+ b(2)*Bi+ b(3)*Ste +b(4)*fraction +b(5)*Bi.*Ste +b(6)*Bi.*fraction +b(7)*Ste.*fraction+...
% % % % % %     b(8)*Bi.^2 +b(9)*Ste.^2 +b(10)*fraction.^2+...
% % % % % %     b(11)*(Bi.*Ste).^2 +b(12)*(Bi.*fraction).^2 +b(13)*(Ste.*fraction).^2;
% % % % % 
% % % % % %%this is reference
% % % % % % res=b(1)+ ...
% % % % % %     b(2)*Bi+ b(3)*Bi.^2 +b(4)*Bi.^3+ ...
% % % % % %     b(5)*fraction+ b(6)*fraction.^2 +b(7)*fraction.^3+ ...
% % % % % %     b(8)*Ste+ b(9)*Ste.^2 +b(10)*Ste.^3 +...
% % % % % %     b(11)*fraction.*Ste + b(12)*fraction.^2.*Ste + b(13)*fraction.*Ste.^2 +b(14)*fraction.^2.*Ste.^2; 
% % % % % res=dot(repmat(b',size(X,1),1),X,2);
% % % % % 
% % % % % error=abs(res-param)./param*100;
% % % % % max(abs(error))
% % % % % 
% % % % % mean(abs(error))
% % % % % figure
% % % % % plot(error,'*')
% % % 
% % % param = Qtot;
% % % X=[ones(size(Bi)) Bi Ste fraction Bi.*Ste Bi.*fraction Ste.*fraction  Ste.^2 fraction.^2 (Ste.*fraction).^2];
% % % b=regress(param,X)
% % % res=b(1)+ b(2)*Bi+ b(3)*Ste +b(4)*fraction +b(5)*Bi.*Ste +b(6)*Bi.*fraction +b(7)*Ste.*fraction+...
% % %     +b(8)*Ste.^2 +b(9)*fraction.^2+...
% % %     b(10)*(Ste.*fraction).^2;
% % % 
% % % error=(res-param)./param*100;
% % % max(abs(error))
% % % 
% % % mean(abs(error))
% % % plot(error)

    