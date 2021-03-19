function fig = mkScatterPlot(Eavg,Iavg,Ethres,Ithres)
      
figure('Position',[654 850 679 590]);set(gcf,'color','w')
subplot(5,5,[6:9,11:14,16:19,21:24])          
failures =intersect(find(Iavg<=Ithres),find(Eavg>=Ethres));   % neither epsc nor ipsc
successes = setdiff(1:length(Iavg),failures); % everything else
% failures are gray  
plot(-Eavg(failures),Iavg(failures),'o','color',[0.8 0.8 0.8],'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','none')
hold on;
% successes are black
plot(-Eavg(successes),Iavg(successes),'ko','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','none')
axis([-1 -3.*mean(Eavg(successes)) -1 3.*mean(Iavg(successes))])
tempAx=axis();
vline(-Ethres)
hline(Ithres)
xlabel('E amplitude')
ylabel('I amplitude')
xticks([0,-mean(Eavg(successes)),-2.*mean(Eavg(successes))])
xticklabels({'0','1','2'})
yticks([0,mean(Iavg(successes)),2.*mean(Iavg(successes))])
yticklabels({'0','1','2'})
set(gca,'box','off','FontSize',20)

subplot(5,5,[1:4])
histogram(-Eavg(failures),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor',[0.8 0.8 0.8],'LineWidth',2)
hold on;
histogram(-Eavg(successes),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor','r','LineWidth',2)
xlim([tempAx(1), tempAx(2)])
yax1=ylim;
xticks([])
set(gca,'box','off','FontSize',20)

subplot(5,5,[10,15,20,25])
histogram(Iavg(failures),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor',[0.8 0.8 0.8],'LineWidth',2)
hold on;
histogram(Iavg(successes),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor','b','LineWidth',2)
xlim([tempAx(3), tempAx(4)])
yax2=ylim;
view([90 -90])
yticks([])
xticks([])
set(gca,'box','off','FontSize',20)
finalY=max([yax1,yax2]);
finalY=80;
subplot(5,5,[1:4])
ylim([0 finalY])
yticks([50])
subplot(5,5,[10,15,20,25])
ylim([0 finalY])
yticks([50])
fig=gcf;
