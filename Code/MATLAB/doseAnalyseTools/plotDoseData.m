function plotDoseData()

 Q=0.5712/1.5;
 threshold = 10;
 workDir='./'; 
 mloc = 'measurements/'; 
 sloc = {'LF_60m/','F_60m/'}; 
 files={'0405'}; 
 C = {'k','w'};
 alpha = [0.5,0.5];

 nor = 0;

 for j=1:length(sloc)
  for i=1:1 %length(files)
    mData=readFile(strcat(workDir,mloc,files{i},'17.txt'));
    mData(:,1) = [];
    mData = mData(:,4:5); mData = mData(:); [mData,order] = sort(mData);

    sData=readFile(strcat(workDir,sloc{j},files{i},'17.txt'));
    sData(:,1)=[]; % remove the time column
    sData = sData(:,4:5); sData = sData(:); sData = sData(order);
  end

  mData = [mData;flip(mData)];
  sData = [sData;Q*flip(sData)];
  patch(mData,sData,C{j},'FaceAlpha',alpha(j)); hold on;
 end
 h=refline(1,0); hold off;
 h.LineWidth = 1.2;
 h.Color = 'k';
 box on;
 xlabel('$\dot{d}_{\gamma ,o}$ $[nSv/h]$','interpreter','latex','FontSize',18);
 ylabel('$\dot{d}_{\gamma ,p}$ $[nSv/h]$','interpreter','latex','FontSize',18);

 axis([0 140 0 140]);
 set(gca,'FontSize',14);
 set(gcf,'Units','Inches');
 pos = get(gcf,'Position');
 set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,'NERIS0405_PMF','-dpdf','-r300');

% set(gcf,'renderer','painters');
% saveas(gcf,'NERIS0405_PMF.pdf'); %,'epsc');

end
