plotClass = 'BM'; % 'BM', 'OF'
BMclass = 'BM.txt';
OFclass = 'MOST_Forest.txt'; % 'MOST_OF.txt', 'MOST_Forest.txt'
labels = {'S','SS','N','SU','U','VU','A'};
rgbmatrix = [[0 0 0];[1 1 1]*64/255;[1 1 1]*128/255;[1 1 1]*160/255;[1 1 1]*192/255;[1 1 1]*224/255;[1 1 1]];
nocol=1;

dataBM=readFile(BMclass,0,1);
dataOF=readFile(OFclass,0,nocol);

stabDistr = zeros(7,2); %8 ipv 7

for i=1:length(dataBM)

 if dataBM(i) == 7
	stabDistr(3,1) = stabDistr(3,1)+1;
 else
	stabDistr(dataBM(i),1) = stabDistr(dataBM(i),1)+1;
 end

 L = dataOF(i,nocol);
 if (abs(L) > 500)
	stabDistr(3,2) = stabDistr(3,2)+1;
 elseif (L >= 200)
	stabDistr(2,2) = stabDistr(2,2)+1;
 elseif (L >= 50)
	stabDistr(1,2) = stabDistr(1,2)+1;
 elseif (L < -50 && L >= -100)
	stabDistr(6,2) = stabDistr(6,2)+1;
 elseif (L >= -200 && L < -50)
	stabDistr(5,2) = stabDistr(5,2)+1;
 elseif (L >= -500 && L < -50)
	stabDistr(4,2) = stabDistr(4,2)+1;
 else
	stabDistr(7,2) = stabDistr(7,2)+1; %8 ipv 7
 end
end

zero1 = {};
zero2 = {};

for i=1:length(stabDistr(:,1))
  if (stabDistr(i,1) == 0)
	stabDistr(i,1)=2^(-64);
	zero1{length(zero1)+1} = i;
  end

  if (stabDistr(i,2) == 0)
	stabDistr(i,2)=2^(-64);
	zero2{length(zero2)+1} = i;
  end
end

figure('units','normalized','outerposition',[0 0 1 1]);

if strcmp(plotClass,'BM')
 p=pie((stabDistr(:,1))');
 prct = findobj(p,'Type','text');
 for i = 1:length(stabDistr(:,1))
   prct(i).String = strcat(labels{i},' (',prct(i).String,')');
 end
 if (length(zero1)>0)
  for i=1:length(zero1)
    prct(zero1{i}).String = '';
  end
 end
elseif strcmp(plotClass,'OF')
 p=pie((stabDistr(:,2))');
 prct = findobj(p,'Type','text');
 for i = 1:length(stabDistr(:,2))
   prct(i).String = strcat(labels{i},' (',prct(i).String,')');
 end
 if (length(zero2)>0)
  for i=1:length(zero2)
    prct(zero2{i}).String = '';
  end
 end
end
set(prct,'FontSize',20);
for k = 1 : length(stabDistr(:,1))
   pieColorMap = rgbmatrix(k,:);  % Color for this segment.
   set(p(k*2-1),'FaceColor', pieColorMap);
end

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition',[ 0 0 1 1 ]);

set(gcf,'renderer','painters');
%saveas(gcf,'pieChart_Forest.eps','epsc');

