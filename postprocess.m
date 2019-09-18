function postprocess(saveddatafiles,options,ifvideo)
lstyles={'-' '--' '-.' ':' '-' '--' '-.' ':' };
lwidth={2 2 2 2 1 1 1 1};
lcolor={[0 0 0] [0 0 0] [0 0 0] [0 0 0] [1 1 1] [1 1 1] [1 1 1] [1 1 1]}; 
%%% reading - just to get Tlow and Thigh
load(['' saveddatafiles{1}])

kapilow=readKapilow([options.period/2 options.period]); 
 
Tinf=0.5*(results.options.Tlow-results.options.Thigh)*(sign(sin(2*pi*(0:0.01:options.t1/options.period)))+1)+results.options.Thigh;
TaxisLow=min([results.options.Tlow results.options.Thigh])-0.1;
TaxisHigh=max([results.options.Tlow results.options.Thigh])+0.1;
%% Average phase fraction plot
Gplot=figure('DefaultAxesFontSize',15);
set(gca,'YColor',[0 0 0])
yyaxis right

plot((0:0.01:options.t1/options.period),Tinf,'b','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
ylabel('T_{\infty} (^\circ C)')
set(gca,'YColor',[0 0 0])
axis([0 options.t1/options.period TaxisLow TaxisHigh])
% plot(kapilow.gammat,kapilow.gamma+kapilow.gammaError,'-');
% errorbar(kapilow.gammat,kapilow.gamma,kapilow.gammaError*ones(length(kapilow.gamma),1),'^b','MarkerFaceColor','b')


% text2=text(0.7,0.9,'$$T_{air}$$', 'Units', 'normalized','FontSize',18,'Interpreter','latex');

for k=1:length(saveddatafiles)
    
end



% % % % %% Overall heat transfer Coef (KAPILOW)
% % % % UplotKipilow=figure('DefaultAxesFontSize',15);
% % % % 
% % % % % shadedErrorBar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1))
% % % % errorbar(kapilow.Ut,kapilow.U,kapilow.UError*ones(length(kapilow.U),1),'ob','MarkerFaceColor','b','DisplayName','U_{lat} (Kapilow)')
% % % % title('Kapilo''s Approach U=Q/(T_{melt}-T_{\infty})','Interpreter','tex')
% % % % hold on
% % % % % plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');
% % % % 
% % % % xlabel('time (s)')
% % % % ylabel('Overal Heat Transfer Coef. (W/m^2K)')
% % % % axis([0 options.t1 0 250])
% % % % set(gca,'YColor',[0 0 0])
% % % % yyaxis right
% % % % set(gca,'YColor',[0 0 0])
% % % % ylabel('Air Temperature (^\circ C)')
% % % % 
% % % % 
% % % % hold on
% % % % 
% % % % 
%% Total Q
Qplot=figure('DefaultAxesFontSize',15);
set(gca,'YColor',[0 0 0])
yyaxis right
% % % % 
% % % % 
% % % % % plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');
plot((0:0.01:options.t1/options.period),Tinf,'-','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
set(gca,'YColor',[0 0 0])
ylabel('T_{\infty} (^\circ C)')
axis([0 options.t1/options.period results.options.Tlow-0.05*(1-results.options.Tlow) 1+0.05*(1-results.options.Tlow)])
% % % % 
% % % % 
% % % % 
% % % % hold on
% % % % 
% % % % %% latent Q
% % % % Qplotlat=figure('DefaultAxesFontSize',15);
% % % % set(gca,'YColor',[0 0 0])
% % % % yyaxis right
% % % % 
% % % % 
% % % % % plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');
% % % % plot((0:0.01:options.t1/options.period),Tinf,'-','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
% % % % set(gca,'YColor',[0 0 0])
% % % % ylabel('T_{\infty} (^\circ C)')
% % % % axis([0 options.t1/options.period results.options.Tlow-0.05*(1-results.options.Tlow) 1+0.05*(1-results.options.Tlow)])
% % % % 
% % % % 
% % % % 
% % % % hold on
% % % % %% T contour
% % % % 
% % % % Tcontour=figure;
% % % % 
% % % % 
% % % % %% Gamma contour
% % % % 
% % % % Gcontour=figure;
% % % % Kcontour=figure;
% % % % q.totMean=[];
% % % % q.latMean=[];
% % % % q.senMean=[];
% % % % period=[];
%%
for k=1:length(saveddatafiles)
    load(['' saveddatafiles{k}])
    results.options.dt
    results.Tmelt
    
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
     %% averaging range
    indices = find(t>=20*options.period & t<30*options.period);
    q.totMean(k) = mean(q.tot(indices));
    q.latMean(k) = mean(q.lat(indices));
    q.senMean(k) = mean(q.sen(indices));
    period(k) = options.period;
%     q.surf=PCM.D/2*PCM.k*(results.T(end,end)-results.T(end,end-1))/Mesh.Csize(end);
    
    
    
    
    %%% Gamma plolt
    figure(Gplot)
    
    yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(results.time/options.period,results.gammaAvg,'k','LineStyle',lstyles{k},'LineWidth',lwidth{k},'DisplayName',['$\gamma$, $\mathcal{T}$=' num2str(options.period) 's']);
    hold on
    xlabel('$t/\mathcal{T}$ ','Interpreter','latex')
    ylabel('\gamma','FontSize',24)
    axis([0 options.t1/options.period -0.05 1.05]);
    
    
%    axis([0 600 Air.Tlow-0.05*(Air.Thigh-Air.Tlow) Air.Thigh+0.05*(Air.Thigh-Air.Tlow)])
    
    %%% Q plot
    figure(Qplot)
    yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(t/options.period,q.tot,'k','LineStyle',lstyles{k},'LineWidth',lwidth{k},'DisplayName',['$Q_{tot}$, $\mathcal{T}$=' num2str(options.period) 's']);
    hold on
    %plot(t/options.period,q.lat,'LineStyle',lstyles{k},'Color',[0.8 0.8 0.8],'LineWidth',lwidth{k},'DisplayName',['$Q_{lat}$, $\mathcal{T}$=' num2str(options.period) 's']);
    
    xlabel('$t/\mathcal{T}$ ','Interpreter','latex')
    ylabel('Heat Transfer Rate (W)')
    axis([0 options.t1/options.period 0 22])
% % % %     
% % % %     
% % % %     %%% Q plot
% % % %     figure(Qplotlat)
% % % %     yyaxis left
% % % %     set(gca,'YColor',[0 0 0])
% % % % 
% % % %     plot(t/options.period,q.lat,'LineStyle',lstyles{k},'Color',[0 0 0],'LineWidth',lwidth{k},'DisplayName',['$Q_{lat}$, $\mathcal{T}$=' num2str(options.period) 's']);
% % % %     hold on
% % % %     xlabel('$t/\mathcal{T}$ ','Interpreter','latex')
% % % %     ylabel('Heat Transfer Rate (W)')
% % % %     axis([0 options.t1/options.period 0 22])
% % % %     
% % % % 
% % % %     
    %% contour plots
% %     figure('DefaultAxesFontSize',15);
% %     set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 2])
% %     for jj=1:8
% %         %subplot(2,4,jj)
% %         subplot_tight(2, 4, jj, [0.01]);
% %         e=abs(t/options.period-(29+(jj-1)*0.125));
% %         cindex = find(e==min(e));%+[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875]+0.01 & t/options.period>29+[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875]-0.01);
% %         [~,h]=contourf(X,Y,results.gamma(:,:,cindex(1)),50);
% %         t(cindex(1))/options.period
% %         axis equal
% %         set(gca,'visible','off')
% %         set(h,'LineColor','none')
% %         colormap((summer))
% %     %     colorbar
% %     %     c = colorbar;
% %     %     c.Label.String = '\gamma';
% %     %     c.Label.FontSize=14;
% %         caxis([0 1])
% %     end
% %     
% %     %%gamma colorbar 
% %     figure('DefaultAxesFontSize',15);
% %     set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 6])
% %     e=abs(t/options.period-(29+(1-1)*0.125));
% %     cindex = find(e==min(e));%+[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875]+0.01 & t/options.period>29+[0 0.125 0.25 0.375 0.5 0.625 0.75 0.875]-0.01);
% %     [~,h]=contourf(X,Y,results.gamma(:,:,cindex(1)),50);
% %     t(cindex(1))/options.period
% %     axis equal
% %     set(gca,'visible','off')
% %     set(h,'LineColor','none')
% %     colormap((summer))
% %     c = colorbar('southoutside');
% %     c.Label.String = '\gamma';
% %     c.Label.FontSize=20;
% %     c.FontSize=14;
% %     caxis([0 1])


    
    text(0,max(max(Y))*1.2,'Phase fraction Distribution','FontSize',14,'HorizontalAlignment','center');
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

% % % % figure(UplotKipilow)
% % % % legend('show')


figure(Gplot)
lGplot=legend('show')
lGplot.Interpreter='latex';
lGplot.Location='best';
% print('results/Gamma','-dsvg','-r300')

% % % % figure(UplotKipilow)
% % % % % print('results/U','-dsvg','-r300')
% % % % 
% % % % figure(Qplot)
% % % % lQplot=legend('show');
% % % % lQplot.Interpreter='latex';
% % % % lQplot.Location='best';
% % % % 
% % % % figure(Qplotlat)
% % % % lQplotlat=legend('show');
% % % % lQplotlat.Interpreter='latex';
% % % % lQplotlat.Location='best';
% % % % 
% % % % figure('DefaultAxesFontSize',15);
% % % % plot(period,q.totMean,'-ok','LineWidth',2,'DisplayNAme','Q_{tot}')
% % % % hold on
% % % % plot(period,q.latMean,'--ok','LineWidth',2,'DisplayNAme','Q_{lat}')
% % % % plot(period,q.senMean,'-.ok','LineWidth',2,'DisplayNAme','Q_{sen}')
% % % % 
% % % % xlabel('$\mathcal{T}$ (s)','Interpreter','latex')
% % % % ylabel('Average Heat Transfer rate (W)')
% % % % legend('show')
% set the seconds per image
% figure('DefaultAxesFontSize',18);
% plot(t,squeeze(results.T(2,2,:)),'LineWidth',2,'DisplayNAme','T_{r=0}')
% hold on
% plot(t,Air.T,'--r','LineWidth',2,'DisplayNAme','T_{infty}')
% axis([0 2000 20 45])
% xlabel('time (s)')
% ylabel('Temeprature')
% yyaxis right
% ylabel('Q')
% plot(t,q.tot,'-.k','LineWidth',2,'DisplayNAme','Q_{tot}')
% plot(t,q.lat,'-.g','LineWidth',2,'DisplayNAme','Q_{tot}')
% legend('show')
% figure('DefaultAxesFontSize',18);
% plot(t,squeeze(max(max(abs(results.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2)),'LineWidth',2,'DisplayNAme','\Delta T')
% xlabel('time (s)')
% ylabel('Temeprature')
% yyaxis right
% ylabel('Q_{tot}')
% plot(t,q.tot,'-.k','LineWidth',2,'DisplayNAme','Q_{tot}')
% legend('show')
% axis([0 2000 -inf inf])

%% analytical solution
% % % % L=PCM.D/2;                     %Lengh of domain (radius)
% % % % Ueair=1/(mean(results.Air.R)+results.HDPE.R)/(HDPE.D/2);
% % % % kpcm=results.pcmk(2,2,444);
% % % % rho=(PCM.rhol+PCM.rhos)/2;
% % % % Nx=options.Nx;
% % % % Ny=options.Ny;
% % % % Bi=Ueair*L/kpcm;
% % % % 
% % % % 
% % % % Tinf=42.9;
% % % % Tinit=38;
% % % % figure('DefaultAxesFontSize',18);
% % % % Tanaline=plot(results.Mesh.Clocx,results.Mesh.Clocx,'--k');
% % % % axis([-inf inf 30 38])
% % % % 
% % % % Q=[];
% % % % deltaT=[];
% % % % U=[];
% % % % figure('DefaultAxesFontSize',18);
% % % % Qline=plot(0:100,0:100);
% % % % ylabel('Q')
% % % % xlabel('Time')
% % % % yyaxis right
% % % % 
% % % % DTline=plot(0:100,0:100)
% % % % ylabel('\Delta T')
% % % % figure()
% % % % Uline=plot(0:100,0:100)
% % % % xlabel('Time')
% % % % ylabel('U')
% % % % for t=0:100
% % % %     lambda=lambdan([-1 2000],Bi);
% % % %     Thetaana=zeros(Nx,1);
% % % %     derivative=zeros(Nx,1);
% % % %     for n=1:length(lambda)
% % % %         %An=(-2*(-1)^n*Bi*sqrt(lambda(n)^2+Bi^2))/(lambda(n)*(lambda(n)^2+Bi*(1+Bi)));
% % % %         kn=lambda(n);
% % % %         An=(2/kn)*(besselj(1,kn))/(besselj(0,kn)^2+besselj(1,kn)^2);
% % % %         Fo=t*kpcm/(rho*PCM.cl)/L^2;
% % % %         Thetaana=Thetaana+An*besselj(0,kn*results.Mesh.Clocx/L)*exp(-kn^2*Fo);
% % % %         derivative = derivative + An *exp(-kn^2*Fo) * (-besselj(1,kn*results.Mesh.Clocx/L)*kn/L);
% % % %     end
% % % %     deriv = derivative * (Tinit-Tinf);
% % % %     Q(end+1) = -deriv(end)*kpcm;
% % % %     
% % % %     Tana=Thetaana*(Tinit-Tinf)+Tinf;
% % % %     deltaT(end+1)=Tana(1) - Tinf;
% % % %     U(end+1)= Q/deltaT;
% % % %     set(Tanaline,'YData',Tana)
% % % %     set(Uline,'YData',U)
% % % %     set(Qline,'YData',Q)
% % % %     set(DTline,'YData',deltaT)
% % % %     set(Uline,'XData',0:t)
% % % %     set(Qline,'XData',0:t)
% % % %     set(DTline,'XData',0:t)
% % % %     refreshdata
% % % %     drawnow
% % % % end
% % % % %%%% animation

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



