noModels = 6;
noSources = 2;
modelNames={'LB-F','LB-OF','LT-OF','LT-F','G-F','G-OF'};
markers={'d','p','o','s','^','v'};
FB = zeros(noModels,noSources);
NMSE = zeros(noModels,noSources);
writeNames = 0;

dy = 0.2;
dx = 0.02;
FB(:,1)=[1.21;1.18;1.10;0.85;1.05;1.10];
FB(:,2)=[0.43;0.38;0.27;-0.06;0.20;0.26];
NMSE(:,1)=[4.0;3.67;3.02;1.68;2.67;3.02];
NMSE(:,2)=[0.59;0.51;0.37;0.28;0.34;0.38];
delta_y = zeros(noModels,1); delta_y(end) = delta_y(end)+1.5*dy;
delta_x = zeros(noModels,1); delta_x(end) = 10*dx; delta_x([2,3]) = -4*dx; 

for i=1:noModels
 plot(FB(i,1),NMSE(i,1),strcat('k',markers{i}),'MarkerFaceColor','k'); hold on;
end
if (writeNames)
  text(FB(:,1)-20*dx+delta_x,NMSE(:,1)+delta_y,modelNames,'Color','k','interpreter','latex','FontSize',16);
end
for i=1:noModels
 plot(FB(i,2),NMSE(i,2),strcat('k',markers{i}),'MarkerEdgeColor','k','LineWidth',1.2); hold on;
end
if (writeNames)
 delta_x = [2;0;-18;-20;-14;-14]*dx; 
 delta_y(6) = 0; delta_y([2,3,5,6]) = [2;2.5;1;4.5]*dy; 
 text(FB(:,2)+delta_x,NMSE(:,2)+delta_y,modelNames,'Color','k','interpreter','latex','FontSize',16);
end
x=linspace(-1.5,1.5);
NMSE_min = 4*x.^2./(4-x.^2);
plot(x,NMSE_min,'-k','LineWidth',1.2); hold off;
xlabel('FB [-]','FontSize',18,'interpreter','latex');
ylabel('NMSE [-]','FontSize',18,'interpreter','latex'); 
set(gca,'FontSize',14);
set(gcf,'renderer','painters','Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(gcf,'NMSE_FB','-dpdf');
%saveas(gcf,'NMSE_FB.eps','epsc');
