function postprocess(saveddatafiles,options,ifvideo)
lstyles={'-' '--' '-.' ':' '-' '--' '-.' ':' };

%%% reading kapilow's data

kapilow=readKapilow([options.period/2 options.period]); 
 

%% Average phase fraction plot
Gplot=figure('DefaultAxesFontSize',15);
plot(kapilow.gammat,kapilow.gamma,'^b','DisplayName','Kapilow et al. \cite{Kapilow2018}');
hold on
% plot(kapilow.gammat,kapilow.gamma+kapilow.gammaError,'-');
% errorbar(kapilow.gammat,kapilow.gamma,kapilow.gammaError*ones(length(kapilow.gamma),1),'^b','MarkerFaceColor','b')
xlabel('time (s)')
ylabel('\gamma')
axis([0 options.t1 -0.05 1.05]);
%axis([0 600 -0.05 1.05]);
yyaxis right
ylabel('T_{\infty} (^\circ C)')

% text2=text(0.7,0.9,'$$T_{air}$$', 'Units', 'normalized','FontSize',18,'Interpreter','latex');

for k=1:length(saveddatafiles)
    
end



%% Overall heat transfer Coef (KAPILOW)
UplotKipilow=figure('DefaultAxesFontSize',15);

% shadedErrorBar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1))
errorbar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1),'ob','MarkerFaceColor','b','DisplayName','U_{lat} (Kapilow)')
title('Kapilo''s Approach U=Q/(T_{melt}-T_{\infty})','Interpreter','tex')
hold on
% plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');

xlabel('time (s)')
ylabel('Overal Heat Transfer Coef. (W/m^2K)')
axis([0 options.t1 0 250])
yyaxis right
ylabel('Air Temperature (^\circ C)')


hold on


%% Overall heat transfer Coef (CURRENT)
UplotCurrent=figure('DefaultAxesFontSize',15);

% shadedErrorBar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1))
% errorbar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1),'ob','MarkerFaceColor','b','DisplayName','U_{lat} (Kapilow)')
plot(0,0)
title('Current Approach U=Q/(T_{core}-T_{\infty})','Interpreter','tex')
hold on
% plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');

xlabel('time (s)')
ylabel('Overal Heat Transfer Coef. (W/m^2K)')
axis([0 options.t1 0 250])
yyaxis right
ylabel('Air Temperature (^\circ C)')

%% Total Q
Qplot=figure('DefaultAxesFontSize',15);
plot(0,0)
title('Q (J/s)','Interpreter','tex')
hold on
% plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');

xlabel('time (s)')
ylabel('Overal Heat Transfer Coef. (W/m^2K)')
axis([0 options.t1 -inf inf])
yyaxis right
ylabel('Air Temperature (^\circ C)')


hold on

%% T contour

Tcontour=figure;


%% Gamma contour

Gcontour=figure;
Kcontour=figure;
q.totMean=[];
q.latMean=[];
q.senMean=[];
period=[];
%%
for k=1:length(saveddatafiles)
    load(['' saveddatafiles{k}])
    Air.T=results.AirT;
    Mesh=results.Mesh;
    % creating mesh grid and adding boundary points
    [R,Theta]=ndgrid([Mesh.Flocx(1); Mesh.Clocx; Mesh.Flocx(end)],[Mesh.Flocy(1); Mesh.Clocy; Mesh.Flocy(end)]);
    X=R.*cos(Theta);
    Y=R.*sin(Theta);
    %%%
    Nx=options.Nx;
    Ny=options.Ny;
    Mesh=results.Mesh;                  %Creating Uniform Mesh
    t=options.t0:options.dt:options.t1;
    
    
    q.tot=squeeze(sum(sum(abs(results.Q(2:Nx+1,2:Ny+1,:)),1),2));  % Watts
    q.lat=squeeze(sum(sum(abs(results.S(2:Nx+1,2:Ny+1,:)),1),2));  % J/s , Watts
    q.sen=q.tot-q.lat;
    q.totMean(k) = mean(q.tot);
    q.latMean(k) = mean(q.lat);
    q.senMean(k) = mean(q.sen);
    period(k) = options.period;
%     q.surf=PCM.D/2*PCM.k*(results.T(end,end)-results.T(end,end-1))/Mesh.Csize(end);
    
    %%% Uplolt (KAPILOW)
    U.tot=q.tot./abs(35-Air.T')/(Mesh.Flocy(end)*HDPE.D/2);
    U.lat=q.lat./abs(35-Air.T')/(Mesh.Flocy(end)*HDPE.D/2);
    U.sen=q.sen./abs(35-Air.T')/(Mesh.Flocy(end)*HDPE.D/2);
    figure(UplotKipilow)
    yyaxis left
    plot(t,U.lat,'k','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','U_{lat}');
    plot(t,U.tot,'r','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','U_{tot}');
    yyaxis right
    plot(t,Air.T,':','LineWidth',2,'DisplayName','T_{air}');
    axis([0 options.t1 Air.Tlow-5 Air.Thigh+5])
    
    %%% Uplolt (CURRENT) 
    %%- using exterem value
    Un.tot=q.tot./squeeze(max(max(abs(results.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2))/(Mesh.Flocy(end)*HDPE.D/2);
    Un.lat=q.lat./squeeze(max(max(abs(results.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2))/(Mesh.Flocy(end)*HDPE.D/2);
    Un.sen=q.sen./squeeze(max(max(abs(results.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2))/(Mesh.Flocy(end)*HDPE.D/2);
    %% using core temp
%     Un.tot=q.tot./abs(squeeze(results.T(1,1,:))-Air.T')/(HDPE.D/2);
%     Un.lat=q.lat./abs(squeeze(results.T(1,1,:))-Air.T')/(HDPE.D/2);
%     Un.sen=q.sen./abs(squeeze(results.T(1,1,:))-Air.T')/(HDPE.D/2);
%     
    
    figure(UplotCurrent)
    yyaxis left
    plot(t,Un.lat,'k','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','U_{lat}');
    plot(t,Un.tot,'r','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','U_{tot}');
    yyaxis right
    plot(t,Air.T,'--','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','T_{air}');
    axis([0 options.t1 Air.Tlow-5 Air.Thigh+5])
    PCM.Rl=log(PCM.D/2/Mesh.Flocx(1))/PCM.kl;        % Thermal Resistance of PCM (K/W)
    PCM.Rs=log(PCM.D/2/Mesh.Flocx(1))/PCM.ks;        % Thermal Resistance of PCM (K/W)
    Uanal=1/(mean(results.Air.R)+results.HDPE.R+PCM.Rl)/(HDPE.D/2);
    Uanas=1/(mean(results.Air.R)+results.HDPE.R+PCM.Rs)/(HDPE.D/2);
    plot([t(1) t(end)],[Uanal Uanal],'--k');
    plot([t(1) t(end)],[Uanas Uanas],'--k');
    %%% Gamma plolt
    figure(Gplot)
    
    yyaxis left
    plot(t,results.gammaAvg,'k','LineStyle',lstyles{k},'LineWidth',2,'DisplayName','Current');
    yyaxis right
    plot(t,Air.T,':','LineWidth',2,'DisplayName','T_{\infty}');
    axis([0 options.t1 Air.Tlow-0.05*(Air.Thigh-Air.Tlow) Air.Thigh+0.05*(Air.Thigh-Air.Tlow)])
%    axis([0 600 Air.Tlow-0.05*(Air.Thigh-Air.Tlow) Air.Thigh+0.05*(Air.Thigh-Air.Tlow)])
    
    %%% Q plot
    figure(Qplot)
    yyaxis left
    plot(t,q.tot,'LineWidth',2);
    yyaxis right
    plot(t,Air.T,'--','LineWidth',2);
    axis([0 options.t1 Air.Tlow-0.05*(Air.Thigh-Air.Tlow) Air.Thigh+0.05*(Air.Thigh-Air.Tlow)])
    
    
    if(ifvideo)
        % video
        % create the video writer with 5 fps
        writerObjT = VideoWriter('T.avi');
        writerObjT.FrameRate = 5;
        open(writerObjT);
        writerObjG = VideoWriter('G.avi');
        writerObjG.FrameRate = 5;
        open(writerObjG);
        for kk=1:size(results.T,3)
            t(kk)
            figure(Tcontour)
            [~,h]=contourf(X,Y,results.T(:,:,kk),256);
            text(0,max(max(Y))*1.2,'Temperature Distribution','FontSize',14,'HorizontalAlignment','center');

            set(gca,'visible','off')
            axis equal
            set(h,'LineColor','none')
            colormap(jet)
            c = colorbar;
            c.Label.String = 'T_{PCM}';
            c.Label.FontSize=14;
            caxis([Air.Tlow Air.Thigh])
            F1(kk) = getframe(gcf);
%             drawnow
%             writeVideo(writerObjT, F1(kk));
            figure(Gcontour)
            [~,h]=contourf(X,Y,results.gamma(:,:,kk),256);
            text(0,max(max(Y))*1.2,'Phase fraction Distribution','FontSize',14,'HorizontalAlignment','center');

            axis equal
            set(gca,'visible','off')
            set(h,'LineColor','none')
            colormap(cool)
            colorbar
            c = colorbar;
            c.Label.String = '\gamma';
            c.Label.FontSize=14;
            caxis([0 .2])
            F2(kk) = getframe(gcf);
%             drawnow
%             writeVideo(writerObjG, F2(kk));



            figure(Kcontour)
            [~,h]=contourf(X,Y,results.pcmk(:,:,kk),256);
            text(0,max(max(Y))*1.2,'Conductivity Distribution','FontSize',14,'HorizontalAlignment','center');

            set(gca,'visible','off')
            axis equal
            set(h,'LineColor','none')
            colormap(jet)
            c = colorbar;
            c.Label.String = 'k_{PCM}';
            c.Label.FontSize=14;
            caxis([PCM.kl PCM.ks])     
            drawnow

        end
        close(writerObjT);
        close(writerObjG);
    end
end

figure(UplotKipilow)
legend('show')


figure(UplotCurrent)
legend('show')

figure(Gplot)
legend('show')
print('results/Gamma','-dsvg','-r300')

figure(UplotKipilow)
print('results/U','-dsvg','-r300')

figure()
plot(period,q.totMean,'o')
hold on
plot(period,q.latMean,'*')
plot(period,q.senMean,'^')

% set the seconds per image
figure('DefaultAxesFontSize',18);
plot(t,squeeze(results.T(2,2,:)),'LineWidth',2,'DisplayNAme','T_{r=0}')
hold on
plot(t,Air.T,'--r','LineWidth',2,'DisplayNAme','T_{infty}')
axis([0 2000 20 45])
xlabel('time (s)')
ylabel('Temeprature')
yyaxis right
ylabel('Q')
plot(t,q.tot,'-.k','LineWidth',2,'DisplayNAme','Q_{tot}')
plot(t,q.lat,'-.g','LineWidth',2,'DisplayNAme','Q_{tot}')
legend('show')
figure('DefaultAxesFontSize',18);
plot(t,squeeze(max(max(abs(results.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2)),'LineWidth',2,'DisplayNAme','\Delta T')
xlabel('time (s)')
ylabel('Temeprature')
yyaxis right
ylabel('Q_{tot}')
plot(t,q.tot,'-.k','LineWidth',2,'DisplayNAme','Q_{tot}')
legend('show')
axis([0 2000 -inf inf])

%% analytical solution
L=PCM.D/2;                     %Lengh of domain (radius)
Ueair=1/(mean(results.Air.R)+results.HDPE.R)/(HDPE.D/2);
kpcm=results.pcmk(2,2,444);
rho=(PCM.rhol+PCM.rhos)/2;
Nx=options.Nx;
Ny=options.Ny;
Bi=Ueair*L/kpcm;


Tinf=42.9;
Tinit=38;
figure('DefaultAxesFontSize',18);
Tanaline=plot(results.Mesh.Clocx,results.Mesh.Clocx,'--k');
axis([-inf inf 30 38])

Q=[];
deltaT=[];
U=[];
figure('DefaultAxesFontSize',18);
Qline=plot(0:100,0:100);
ylabel('Q')
xlabel('Time')
yyaxis right

DTline=plot(0:100,0:100)
ylabel('\Delta T')
figure()
Uline=plot(0:100,0:100)
xlabel('Time')
ylabel('U')
for t=0:100
    lambda=lambdan([-1 2000],Bi);
    Thetaana=zeros(Nx,1);
    derivative=zeros(Nx,1);
    for n=1:length(lambda)
        %An=(-2*(-1)^n*Bi*sqrt(lambda(n)^2+Bi^2))/(lambda(n)*(lambda(n)^2+Bi*(1+Bi)));
        kn=lambda(n);
        An=(2/kn)*(besselj(1,kn))/(besselj(0,kn)^2+besselj(1,kn)^2);
        Fo=t*kpcm/(rho*PCM.cl)/L^2;
        Thetaana=Thetaana+An*besselj(0,kn*results.Mesh.Clocx/L)*exp(-kn^2*Fo);
        derivative = derivative + An *exp(-kn^2*Fo) * (-besselj(1,kn*results.Mesh.Clocx/L)*kn/L);
    end
    deriv = derivative * (Tinit-Tinf);
    Q(end+1) = -deriv(end)*kpcm;
    
    Tana=Thetaana*(Tinit-Tinf)+Tinf;
    deltaT(end+1)=Tana(1) - Tinf;
    U(end+1)= Q/deltaT;
    set(Tanaline,'YData',Tana)
    set(Uline,'YData',U)
    set(Qline,'YData',Q)
    set(DTline,'YData',deltaT)
    set(Uline,'XData',0:t)
    set(Qline,'XData',0:t)
    set(DTline,'XData',0:t)
    refreshdata
    drawnow
end
%%%% animation

% load([datafiles{1}])
% %Temperature plot
% anim=figure('DefaultAxesFontSize',15);
% hold on
% Tline=plot(Mesh.Cloc,results.T(1,:),'-k','LineWidth',2,'MarkerFaceColor','b');
% hold on
% % Tbc=plot(Mesh.Floc(end),(results.T(1,end-1)+results.T(1,end))/2,'ok','MarkerFaceColor','b');
% xlabel('Radius')
% ylabel('Temperature')
% axis([Mesh.Floc(1) Mesh.Floc(end) Air.Tlow-1 Air.Thigh+1])
% %phase fraction plot
% yyaxis right
% Gammaline=plot(Mesh.Cloc,ones(1,Nx+2),'--','LineWidth',2,'MarkerFaceColor','r');
% xlabel('Radius')
% ylabel('Phase Fraction')
% text1=text(0.1,0.7,'$$\overline{a}$$', 'Units', 'normalized','FontSize',18,'Interpreter','latex');
% axis([Mesh.Floc(1) Mesh.Floc(end) -0.05 1.05])
% 
% gammaAvg=[];
% for t=options.t0+options.dt:options.dt:options.t1
%     i=round((t-options.t0)/options.dt);
%     [~, ~, gamma]= calcState(results.T(i,:),PCM,0);
%     gammaAvg(end+1)=sum(gamma(2:Nx+1)'.*Mesh.V(2:Nx+1))/sum(Mesh.V(2:Nx+1));
%     
%     if (true)%rem(t,10)==0)
%         set(Tline,'YData',results.T(i,:))
% %         set(Tbc,'YData',(results.T(i,end-1)+results.T(i,end))/2);
%         set(Gammaline,'YData',gamma)
% 
%         set(text1,'String',['$$\overline{\gamma}=$$ ' num2str(round(gammaAvg(end),2)) char(10)  '$$t=$$ ' num2str(t,3) char(10) '$$T_{air}=$$'   num2str(results.AirT(i),3)  ]);%  '%']);
% 
% 
%     
%     end
%     
%     refreshdata
%     drawnow
% 
% end
% 
% 



