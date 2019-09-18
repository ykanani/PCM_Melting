function postprocessResults(results,NperiodsAvg,ifvideo)
lstyles=repmat({'-' '--' '-.' ':'},1,10);
lwidth=1.5*ones(1,100);
cTemp=linspecer(10);
colors=zeros(10,3);%[repelem(cTemp(:,1),4) repelem(cTemp(:,2),4)  repelem(cTemp(:,3),4) ];%
colors2=['r' 'b' 'g'];
r1=results{1};
options=r1.options;
options.t1=options.period/4;
Tinf=0.5*(options.Tlow-options.Thigh)*(sign(sin(2*pi*(0:0.01:options.t1/options.period)))+1)+options.Thigh;
TaxisLow=min([options.Tlow options.Thigh]);
TaxisHigh=max([options.Tlow options.Thigh]);

%% Average phase fraction plot
Gplot=figure('DefaultAxesFontSize',15);
Gplot2=figure('DefaultAxesFontSize',15);
%set(gca,'YColor',[0 0 0])
% yyaxis right
% 
% plot((0:0.01:options.t1/options.period)*options.period,Tinf,'b','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
% ylabel('T_{\infty} (^\circ C)')
% set(gca,'YColor',[0 0 0])
axis([-inf inf TaxisLow TaxisHigh])



%% Centerline Temp. plot
Tplot=figure('DefaultAxesFontSize',15);
set(gca,'YColor',[0 0 0])
% yyaxis right
% 
% plot((0:0.01:options.t1/options.period)*options.period,Tinf,'b','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
% ylabel('T_{\infty} (^\circ C)')
% set(gca,'YColor',[0 0 0])
axis([-inf inf TaxisLow TaxisHigh])

Tplot2=figure('DefaultAxesFontSize',15);
set(gca,'YColor',[0 0 0])
axis([-inf inf TaxisLow TaxisHigh])
%% Average phase fraction plot
GTplot=figure('DefaultAxesFontSize',15);
set(gca,'YColor',[0 0 0])
yyaxis right
axis([-inf inf TaxisLow TaxisHigh])
set(gca,'YColor',[0 0 0])
% 
% plot((0:0.01:options.t1/options.period)*options.period,Tinf,'b','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
% ylabel('T_{\infty} (^\circ C)')
% set(gca,'YColor',[0 0 0])

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

%% Total Q2
Q2plot=figure('DefaultAxesFontSize',15);
%set(gca,'YColor',[0 0 0])
% yyaxis right
% % % % 
% % % % 
% % % % % plot(kapilow.Ut,kapilow.U,'ob','MarkerFaceColor','b');
%plot((0:0.01:options.t1/options.period),Tinf,'-','LineWidth',1,'Color',[0.3 0.7 0.9],'DisplayName','$T_{\infty}$');
%set(gca,'YColor',[0 0 0])
%ylabel('T_{\infty} (^\circ C)')
%axis([options.Nperiods-NperiodsAvg options.t1/options.period TaxisLow TaxisHigh])

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
% % % % axis([0 options.t1/options.period results{k}.options.Tlow-0.05*(1-results{k}.options.Tlow) 1+0.05*(1-results{k}.options.Tlow)])
% % % % 
% % % % 
% % % % 
% % % % hold on
%% T contour


                

%% Gamma contour

Ganim=figure;
Kcontour=figure;
q.totMean=zeros(1,length(results));
q.latMean=zeros(1,length(results));
q.senMean=zeros(1,length(results));
% % % % period=[];
%%
for k=1:length(results)
    fprintf('###############################\n#############################\n')
    fprintf('Ste = %2.2f\n',results{k}.options.Ste)
    fprintf('Bi = %2.2f\n',results{k}.options.Bi)
%     fprintf('dt/T = %2.5f\n',results{k}.options.dt/results{k}.options.period)
    fprintf('Period = %2.2f\n',results{k}.options.period)
    fprintf('Melt Time = %2.2f\n',results{k}.options.meltTime)
    fprintf('Fraction = %2.2f\n',results{k}.options.fraction)

%     results{k}.Tmelt
    
    options = results{k}.options;
    options.t1=options.Nperiods*options.period;
    Air.T=results{k}.AirT;
    Mesh=results{k}.Mesh;
    t = results{k}.time - results{k}.time(end) + options.meltTime*2 ;
    % creating mesh grid and adding boundary points
    [R,Theta]=ndgrid([Mesh.Flocx(1); Mesh.Clocx; Mesh.Flocx(end)],[Mesh.Flocy(1); Mesh.Clocy; Mesh.Flocy(end)]);
    X=R.*cos(Theta);
    Y=R.*sin(Theta);
    %%%
    Nx=options.Nx;
    Ny=options.Ny;
    Mesh=results{k}.Mesh;                  %Creating Uniform Mesh
%     t=options.t0:options.dt:options.t1;
    
   
    q.tot=abs(squeeze(sum(sum((results{k}.Q(2:Nx+1,2:Ny+1,:)),1),2)));  % Watts
    q.lat=abs(squeeze(sum(sum((results{k}.S(2:Nx+1,2:Ny+1,:)),1),2)));  % J/s , Watts
    q.sen=q.tot-q.lat;
     %% averaging range
     
    indices = find(t>0);%find(t>=Tmin*options.period & t<Tmax*options.period);
%     indices=indices(1:20:end);
    tmin=0;%min(results{k}.time(indices));
    tmax=t(end);
    q.totMean(k) = mean(q.tot(indices));
    q.latMean(k) = mean(q.lat(indices));
    q.senMean(k) = mean(q.sen(indices));
    period(k) = options.period;
    QavePerVolume = trapz(t(indices),q.tot(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)))/((tmax-tmin)*2);
    QavePerVolumeLat = options.Ste*trapz(t(indices),q.lat(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)))/((tmax-tmin)*2);
    QtotPerVolume = trapz(t(indices),q.tot(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)));
    QtotPerVolumeLat = options.Ste*trapz(t(indices),q.lat(indices))/(pi*options.Lx^2*(options.Ly/(2*pi)));
    1+options.Ste;
%     fprintf('Qtot = %2.3f,Qlat = %2.3f,  1+Ste=%2.3f, ratio =%2.3f\n',QtotPerVolume,QtotPerVolumeLat,1+options.Ste,QavePerVolume/(1+options.Ste));
    
    
    %trapz(results{1}.time,q.tot)/4/(pi/2)*.1
%     q.surf=PCM.D/2*PCM.k*(results{k}.T(end,end)-results{k}.T(end,end-1))/Mesh.Csize(end);
    
    %% melt front detection
     %%% Gamma plolt
    figure(Gplot)
    
    %yyaxis left
%     set(gca,'YColor',[0 0 0])
    plot(t(indices)/options.meltTime,results{k}.gammaAvg(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$\zeta=$' num2str(options.fraction)]);
    hold on
    xlabel('$\tilde{t}/\tilde{\mathcal{T}}_\mathrm{ref}$ ','Interpreter','latex','FontSize',18)
    ylabel('\gamma','FontSize',18)
    axis([0 2 0 1]);
    figure(Gplot2)
% %     subplot(5,1,5)

    plot(t(indices)/options.meltTime,results{k}.gammaAvg(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',0.5,...
        'color',[0.5 0.5 0.5],'DisplayName',['$\gamma,\zeta=$' num2str(options.fraction)]);
    hold on
% %     subplot(5,1,1:4)
    hold on
    Xi = zeroCrossing( Y(:,2), squeeze(results{k}.T(:,2,indices)) );
    times=repmat(t(indices)/options.meltTime,[size(Xi,1) 1]);
    XT=[Xi(:),times(:)];
    XT(any(isnan(XT), 2),:)=[];
    dist = pdist2(XT,XT);

    N = size(XT,1);
    %res = NaN(1,N);
    res=[];
    res(1) = 1; % first point is first row in data matrix
    dist(:,1) = Inf;
    oldpoint = 1;
    isfirst=1;
    for ii=2:N
        
        [mindist, closest_idx] = min(dist(oldpoint,:));
        if mindist>0.07
            if isfirst
                plot(XT(res,2),XT(res,1),'-','color',colors(k,:),'LineStyle',lstyles{k},'LineWidth',lwidth(k),'DisplayName',['$\tilde{r}_{interface},\zeta=$' num2str(options.fraction)]);
                isfirst = 0;
            else
                plot(XT(res,2),XT(res,1),'-','color',colors(k,:),'LineStyle',lstyles{k},'LineWidth',lwidth(k),'DisplayName','','HandleVisibility','off');
            end
            res=[];
        end
        res(end+1) = closest_idx;
        oldpoint=closest_idx;
        dist(:,closest_idx) = Inf;
    end
    plot(XT(res,2),XT(res,1),'-','color',colors(k,:),'LineStyle',lstyles{k},'LineWidth',lwidth(k),'DisplayName','','HandleVisibility','off');
    
    xlabel('$\tilde{t}/\tilde{\mathcal{T}}_\mathrm{ref}$ ','Interpreter','latex','FontSize',18)
    ylabel('$\gamma,\tilde{r}$','Interpreter','latex','FontSize',18)
    axis([0 2 0 1]);
    
    %%% Center  plolt
    figure(Tplot)
    
    %yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(t(indices)/options.meltTime,squeeze(results{k}.T(end-1,2,indices)),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$\zeta=$' num2str(options.fraction)]);
    hold on
    xlabel('$\tilde{t}/\tilde{\mathcal{T}}_\mathrm{ref}$ ','Interpreter','latex','FontSize',18)
    ylabel('$\theta_{\tilde{r}=1}$','Interpreter','latex','FontSize',18)
    axis([0 2 TaxisLow TaxisHigh]);
    
    figure(Tplot2)
    
    %yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(t(indices)/options.meltTime,squeeze(results{k}.T(1,2,indices)),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$\zeta=$' num2str(options.fraction)]);
    hold on
    xlabel('$\tilde{t}/\tilde{\mathcal{T}}_\mathrm{ref}$ ','Interpreter','latex','FontSize',18)
    ylabel('$\theta_{\tilde{r}=0}$','Interpreter','latex','FontSize',18)
    axis([0 2 TaxisLow TaxisHigh]);
    
    %%% Gamma plolt
    figure(GTplot)
    
    yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(t(indices)/options.meltTime,results{k}.gammaAvg(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$Ste=$' num2str(options.Ste) '$,Bi=$' num2str(options.Bi)]);
    hold on
    xlabel('$\tilde{t}Ste$ ','Interpreter','latex','FontSize',24)
    ylabel('\gamma','FontSize',24)
    axis([0 2 0 1]);
    
   
    yyaxis right
    set(gca,'YColor',[0 0 0])
    plot(t(indices),squeeze(results{k}.T(end-1,2,indices)),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color','r','DisplayName','','HandleVisibility','off');
    hold on
    %xlabel('$\tilde{t}$ ','Interpreter','latex','FontSize',24)
    ylabel('\theta','FontSize',24)
    axis([-inf inf TaxisLow TaxisHigh]);
%    axis([0 600 Air.Tlow-0.05*(Air.Thigh-Air.Tlow) Air.Thigh+0.05*(Air.Thigh-Air.Tlow)])
    
    %%% Q plot
    figure(Qplot)
%     yyaxis left
%     set(gca,'YColor',[0 0 0])
    Qdot = q.tot(indices)/((pi*options.Lx^2*(options.Ly/(2*pi))));
    plot(t(indices)/options.meltTime,Qdot,'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$\zeta=$' num2str(options.fraction)]);
    hold on
    %plot(t/options.period,q.lat,'LineStyle',lstyles{k},'Color',[0.8 0.8 0.8],'LineWidth',lwidth{k},'DisplayName',['$Q_{lat}$, $\mathcal{T}$=' num2str(options.period) 's']);
    
    xlabel('$\tilde{t}/\tilde{\mathcal{T}}_\mathrm{ref}$','Interpreter','latex');
    ylabel('$|\tilde{q}|$','Interpreter','latex')
%     yyaxis right
%     set(gca,'YColor',[0 0 0])
%     plot(t(indices)/options.meltTime,results{k}.gammaAvg(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
%         'color',[0.5 0.5 0.5],'DisplayName','','HandleVisibility','off');
%     ylabel('\gamma','FontSize',18)
    axis([0 2 -inf inf ]);
    
    
    %%total transfered Q
    Qtot = mean(Qdot) * (tmax-tmin)/2;
    fprintf('qdot AVG = %2.3f, Dt=%2.3f, Qtot=%2.3f\n',mean(Qdot),tmax-tmin,Qtot);
    %%% Q2 plot
    figure(Q2plot)
    yyaxis left
    set(gca,'YColor',[0 0 0])
    plot(t(indices),q.tot(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName',['$Ste=$' num2str(options.Ste) '$,Bi=$' num2str(options.Bi)]);
    hold on
    
    %plot(t/options.period,q.lat,'LineStyle',lstyles{k},'Color',[0.8 0.8 0.8],'LineWidth',lwidth{k},'DisplayName',['$Q_{lat}$, $\mathcal{T}$=' num2str(options.period) 's']);
    %set(gca, 'YScale', 'log')
    %set(gca, 'XScale', 'log')
    xlabel('$\tilde{t}$','Interpreter','latex','FontSize',18) %/\tilde{\mathcal{T}}_\mathrm{ref}$ '
    ylabel('$\tilde{q}$','Interpreter','latex','FontSize',18)
    yyaxis right
    set(gca,'YColor',[0 0 0])
    plot(results{k}.time(indices),results{k}.gammaAvg(indices),'MarkerEdgeColor','none','LineStyle',lstyles{k},'LineWidth',lwidth(k),...
        'color',colors(k,:),'DisplayName','','HandleVisibility','off');
    ylabel('\gamma','FontSize',18)
    axis([-inf inf 0 1]);

    
    if(ifvideo)
        % video
        % create the video writer with 5 fps
        anim=figure('DefaultAxesFontSize',15);
        yyaxis left;
        set(gca,'YColor',[0 0 0])
        Tline = plot(Y(:,2),results{k}.T(:,2,1),'MarkerEdgeColor','none','LineStyle',lstyles{1},'LineWidth',lwidth(1),...
                    'color',colors(1,:),'DisplayName','$\theta$');
        axis([0 1 TaxisLow TaxisHigh])
        hold on
        Mline = plot(0,0,'pr');
        yyaxis right
        set(gca,'YColor',[0 0 0])
        Gline = plot(Y(:,2),results{k}.gamma(:,2,1),'MarkerEdgeColor','none','LineStyle',lstyles{1},'LineWidth',lwidth(1),...
                    'color','r','DisplayName','$\gamma$');
        
        axis([0 1 0 1])
        animplot=legend('show');
        animplot.Interpreter='latex';
        animplot.Location='best';
        writerObjTG = VideoWriter('TG.avi');
        writerObjTG.FrameRate = 5;
        open(writerObjTG);
        
        for kk=indices
            
            set(Tline,'YData',results{k}.T(:,2,kk))
            set(Gline,'YData',results{k}.gamma(:,2,kk))
            Xi= zeroCrossing( Y(:,2), results{k}.T(:,2,kk) );
            set(Mline,'YData',zeros(size(Xi)));
            set(Mline,'XData',Xi);
            F1(kk) = getframe(gcf);
            refreshdata
            drawnow

        end
        close(writerObjTG);
    end
end

% % % % figure(UplotKipilow)
% % % % legend('show')


figure(Gplot)
lGplot=legend('show');
lGplot.Interpreter='latex';
lGplot.Location='best';

figure(Gplot2)
lGplot2=legend('show');
lGplot2.Interpreter='latex';
lGplot2.Location='best';


figure(Tplot)
lTplot=legend('show');
lTplot.Interpreter='latex';
lTplot.Location='best';

figure(Tplot2)
lTplot2=legend('show');
lTplot2.Interpreter='latex';
lTplot2.Location='best';

figure(Qplot)
lQplot=legend('show');
lQplot.Interpreter='latex';
lQplot.Location='best';

figure(Q2plot)
lQ2plot=legend('show');
lQ2plot.Interpreter='latex';
lQ2plot.Location='best';

figure(GTplot)
lGTplot=legend('show');
lGTplot.Interpreter='latex';
lGTplot.Location='best';
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
% plot(t,squeeze(results{k}.T(2,2,:)),'LineWidth',2,'DisplayNAme','T_{r=0}')
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
% plot(t,squeeze(max(max(abs(results{k}.T-repmat(shiftdim(Air.T',-2),Nx+2,Ny+2,1)),[],1),[],2)),'LineWidth',2,'DisplayNAme','\Delta T')
% xlabel('time (s)')
% ylabel('Temeprature')
% yyaxis right
% ylabel('Q_{tot}')
% plot(t,q.tot,'-.k','LineWidth',2,'DisplayNAme','Q_{tot}')
% legend('show')
% axis([0 2000 -inf inf])

%% analytical solution
% % % % L=PCM.D/2;                     %Lengh of domain (radius)
% % % % Ueair=1/(mean(results{k}.Air.R)+results{k}.HDPE.R)/(HDPE.D/2);
% % % % kpcm=results{k}.pcmk(2,2,444);
% % % % rho=(PCM.rhol+PCM.rhos)/2;
% % % % Nx=options.Nx;
% % % % Ny=options.Ny;
% % % % Bi=Ueair*L/kpcm;
% % % % 
% % % % 
% % % % Tinf=42.9;
% % % % Tinit=38;
% % % % figure('DefaultAxesFontSize',18);
% % % % Tanaline=plot(results{k}.Mesh.Clocx,results{k}.Mesh.Clocx,'--k');
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
% % % %         Thetaana=Thetaana+An*besselj(0,kn*results{k}.Mesh.Clocx/L)*exp(-kn^2*Fo);
% % % %         derivative = derivative + An *exp(-kn^2*Fo) * (-besselj(1,kn*results{k}.Mesh.Clocx/L)*kn/L);
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
% Tline=plot(Mesh.Cloc,results{k}.T(1,:),'-k','LineWidth',2,'MarkerFaceColor','b');
% hold on
% % Tbc=plot(Mesh.Floc(end),(results{k}.T(1,end-1)+results{k}.T(1,end))/2,'ok','MarkerFaceColor','b');
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
%     [~, ~, gamma]= calcState(results{k}.T(i,:),PCM,0);
%     gammaAvg(end+1)=sum(gamma(2:Nx+1)'.*Mesh.V(2:Nx+1))/sum(Mesh.V(2:Nx+1));
%     
%     if (true)%rem(t,10)==0)
%         set(Tline,'YData',results{k}.T(i,:))
% %         set(Tbc,'YData',(results{k}.T(i,end-1)+results{k}.T(i,end))/2);
%         set(Gammaline,'YData',gamma)
% 
%         set(text1,'String',['$$\overline{\gamma}=$$ ' num2str(round(gammaAvg(end),2)) char(10)  '$$t=$$ ' num2str(t,3) char(10) '$$T_{air}=$$'   num2str(results{k}.AirT(i),3)  ]);%  '%']);
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



