params=updateConstDir();
treeParams;
uStar=params.uStar;
plantCover=1-exp(-params.LAI*params.extCoeff);
clear params
save('treeParams.mat');
%clearvars -except plantCover uh duh z0 h uStar

heights = [0.26,0.62,0.84,0.96,1,2.07]*h;
winds = [0.6,0.7,1.8,2.5,2.6,5.8];
sdu=[0.02,0.5,1.7,3.0,3.0,3.5];
stress=[0.3,0.8,0.9,1.0,1.1^2,1.1];
z=linspace(0,3*h,1000);
L=intFGF(z);
U=UProfile(z,L); 
gamma_1 = eta*(alpha-1)+1;
Z0 = alpha^2*z0/gamma_1;
z_M = linspace(Z0+d,2*h,100);
z_PGR = linspace(z0+d,3*h,1000);
U_M = uStar/(kappa*gamma_1)*log((z_M-d)/Z0);
U_PGR = uStar/kappa*(log((z_PGR-d)/z0)+2/3*log(1+1.5*zStar./(mu*(z_PGR-d))).*exp(-mu*(z_PGR-d)/zStar)-2/3*log(1+1.5*zStar/(mu*z0))*exp(-mu*z0/zStar)); 
z_Yi = z(z<=h); U_Yi = U(z<=h);

 [X,Z]=meshgrid(U,z);
 y=min(max(Z,cdepth),h);
 Y=(y-cdepth)/(h-cdepth); clear y;
 lad=a1*Y.^a2.*(a4-Y).^a3;
 nmax=sum(z<=h);
 lad(1:nmax,:)=max(lad(1:nmax,:),TAD);

 subplot(1,3,[1,2]);
 norm=winds(end-1);
 pcolor(X/norm,Z/h,lad); hold on;
 ZeroToOne=linspace(0,1,64)'; OneToZero=ZeroToOne(64-(0:63));
 cmap=[OneToZero,ones(64,1),OneToZero];
 colormap(cmap);
 caxis([0,max(max(lad))]);
 shading interp;
 axis tight;

 handle = zeros(1,3);
 plot(U_PGR/norm,z_PGR/h,'--k'); hold on; 
 plot(U_Yi/norm,z_Yi/h,'-k'); hold on; 
 ax=gca; xax = ax.XLim; xax=linspace(xax(1),xax(2),100);
 plot(xax,z0g/h*ones(1,100),'.k','MarkerSize',3); hold on;
 plot(xax,(d+z0)/h*ones(1,100),'.k','MarkerSize',3); hold on;
 plot(xax,ones(1,100),'.k','MarkerSize',3); hold on;
 plot(xax,2*ones(1,100),'.k','MarkerSize',3); hold on;
 hold off;
 yticks([z0g/h,(d+z0)/h,1,2]); 
 yticklabels({'$z_{0,s}$','$d_{\rm tr}+z_{0}$','$h_{\rm tr}$','$z_{\ast}$'}); 
 xlabel('wind speed','FontSize',18,'interpreter','latex');
 set(gca,'FontSize',16,'xtick',[],'TickLabelInterpreter','latex');
 ylh = ylabel('height','FontSize',18,'interpreter','latex');
 ylh.Position(1) = ylh.Position(1) + abs(ylh.Position(1) * 0.5);
 posax = ax.Position;
 annotation('doublearrow','Position',[posax(1)+0.8*posax(3),posax(2),0,1/3*posax(4)]);
 annotation('textbox','String','inside canopy','Interpreter','latex','FontSize',16,'EdgeColor','none','BackgroundColor','w','Position',[posax(1)+0.45*posax(3),posax(2)+2*z0/h*posax(4),0.1*posax(3),0.1*posax(4)],'FitBoxToText','on');
 annotation('doublearrow','Position',[posax(1)+0.12*posax(3),posax(2)+1/3*posax(4),0,2/3*posax(4)]);
 annotation('textbox','String','above canopy','Interpreter','latex','FontSize',16,'EdgeColor','none','BackgroundColor','w','Position',[posax(1)+0.02*posax(3),posax(2)+0.7*posax(4),0.1*posax(3),0.1*posax(4)],'FitBoxToText','on')
 annotation('textarrow','String','class $C^{1}$ condition','FontSize',16,'interpreter','latex','TextBackgroundColor','w','Position',[posax(1)+U_Yi(end)/U_PGR(end)*posax(3),posax(2)+1.6/3*posax(4),0,-0.6/3*posax(4)]);
 annotation('textbox','String','$\bar{u}(z)=\bar{u}_{h}e^{-\eta \int _{z}^{h_{\rm tr}} \mathcal{F}(s) {\rm d}s}$','interpreter','latex','FontSize',15,'EdgeColor','none','BackgroundColor','none','Position',[posax(1)+0.09*posax(3),posax(2)+0.02*posax(4),0.1*posax(3),0.1*posax(4)],'FitBoxToText','on');
 annotation('textbox','String','$\frac{{\rm d}\bar{u}}{{\rm d}Z} = \frac{u_{\ast}}{\kappa Z} \Phi _{M} \varphi _{M}$','interpreter','latex','FontSize',16,'EdgeColor','none','BackgroundColor','none','Position',[posax(1)+0.18*posax(3),posax(2)+0.9*posax(4),0.1*posax(3),0.1*posax(4)],'FitBoxToText','on');
 annotation('textbox','String','$\frac{{\rm d}\theta _{v}}{{\rm d}Z} = \frac{\theta _{\ast}}{\kappa Z} \Phi _{H} \varphi _{H}$','interpreter','latex','FontSize',16,'EdgeColor','none','BackgroundColor','none','Position',[posax(1)+0.18*posax(3),posax(2)+0.8*posax(4),0.1*posax(3),0.1*posax(4)],'FitBoxToText','on');

 subplot(1,3,3);
 z_Yi = z(z<h); z_Yi(1)=[];
 L=intFGF(z_Yi,'derive'); 
 s=pcolor([[0,max(L(2:end))/L(1)];[0,max(L(2:end))/L(1)]],[[z(1)/h,z(1)/h];[z(end)/h,z(end)/h]],[[0 0];[0 0]]); hold on;
 s.FaceColor=[1 1 1]; 
 s.EdgeColor=[1 1 1];
 plot(L(2:end)/max(L),z_Yi(2:end)/h,'-k'); hold on; 
 ax=gca; xax = ax.XLim; nopts = 30; xax=linspace(xax(1),xax(2),nopts);
 plot(xax,10*z0/h*ones(1,nopts),'.k','MarkerSize',3); hold on;
 plot(xax,zi/h*ones(1,nopts),'.k','MarkerSize',3); hold on;
 plot(xax,zc/h*ones(1,nopts),'.k','MarkerSize',3); hold on;
 plot(xax,ones(1,nopts),'.k','MarkerSize',3); hold on;
 hold off;
 ylim([z(1)/h,z(end)/h]);
 yticks([10*z0/h,zi/h,zc/h*0.95,1.05]);
 yticklabels({'$z_{1}$','$z_{2}$','$z_{3}$','$h_{\rm tr}$'});
 xlabel('$\mathcal{F}(z)$','FontSize',18,'interpreter','latex');
 set(gca,'FontSize',16,'xtick',[],'TickLabelInterpreter','latex'); 
 ylabel('');
 box off;

 set(gcf,'Units','Inches');
 pos = get(gcf,'Position');
 set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,'forestSketch','-dpdf'); %,'-r300');

% set(gcf,'renderer','painters');
% saveas(gcf,'forestSketch.pdf'); %,'epsc');

