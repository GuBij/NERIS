params=updateConstDir();
treeParams;
uStar=params.uStar;
plantCover=1-exp(-params.LAI*params.extCoeff);
clear params
save('treeParams.mat');
clearvars -except plantCover uh duh z0 h uStar

quantity='windSpeed'; %'stress'; %'windSpeed';
heights = [0.26,0.62,0.84,0.96,1,2.07]*h;
winds = [0.6,0.7,1.8,2.5,2.6,5.8];
sdu=[0.02,0.5,1.7,3.0,3.0,3.5];
stress=[0.3,0.8,0.9,1.0,1.1^2,1.1];
z=linspace(0,2.3*h,1000);
L=intFGF(z);
U=UProfile(z,L); %Ubar=mean(U);

%L=intFGF(z,'derive');
U_cstLAD=UProfile(z,L,'Wang');
tau=stressProfile(z,L);
tau_cstLAD=stressProfile(z,L,'Wang');

if strcmp(quantity,'windSpeed')
 norm=winds(end-1);
 plot(U/norm,z/h,'-k','LineWidth',1.1); hold on; %z/h %+4
 plot(U_cstLAD/norm,z/h,'--k','LineWidth',1.1); hold on;
 plot(winds/norm,heights/h,'ok','MarkerFaceColor','k'); hold off;
 xlim([0,max(2.3,winds(end)/norm)]);
 ylim([0,max(2.3,winds(end)/norm)]);
 xlabel('$u/u_h$ [-]','FontSize',18,'interpreter','latex');
elseif strcmp(quantity,'stress')
 norm=stress(end-1);
 plot(tau/norm,z/h,'-k','LineWidth',1.1); hold on;
 plot(tau_cstLAD/norm,z/h,'--k','LineWidth',1.1); hold on;
 plot(stress/norm,heights/h,'ok','MarkerFaceColor','k'); hold off;
 xlim([0,max(2.3,stress(end)/norm)]);
 ylim([0,max(2.3,stress(end)/norm)]);
 xlabel('$\tau/\tau_h$ [-]','FontSize',18,'interpreter','latex');
end
set(gca,'FontSize',18);
ylabel('$z/h$ [-]','FontSize',18,'interpreter','latex');
set(gcf,'renderer','painters');
%saveas(gcf,'canopyProfiletau_Sellier.eps','epsc'); %'LIDARwind.png'); %'NERIS_windProfiles.pdf'); %canopyProfiletau_Sellier.png');
