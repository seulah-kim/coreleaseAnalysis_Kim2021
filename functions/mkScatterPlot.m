function fig = mkScatterPlot(Eavg,Iavg,Ethres,Ithres)
% mScatterPlot function accepts amplitudes of E, I and threshold of each to make scatter plots as shown in Fig 2.c in Kim et al. 2021
      
figure('Position',[654 850 679 590]);set(gcf,'color','w')
subplot(5,5,[6:9,11:14,16:19,21:24])          
failures =intersect(find(Iavg<=Ithres),find(Eavg>=Ethres));   % neither epsc nor ipsc
successes = setdiff(1:length(Iavg),failures); % everything else
% failure trials are gray circles
plot(-Eavg(failures),Iavg(failures),'o','color',[0.8 0.8 0.8],'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','none')
hold on;
% successes trials are black circles
plot(-Eavg(successes),Iavg(successes),'ko','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','none')
axis([-1 -3.*mean(Eavg(successes)) -1 3.*mean(Iavg(successes))]) 
tempAx=axis();
vline(-Ethres)
hline(Ithres)
xlabel('E amplitude')
ylabel('I amplitude')
xticks([0,-mean(Eavg(successes)),-2.*mean(Eavg(successes))]) % ticks are labeled as normalized amplitudes of success trials
xticklabels({'0','1','2'})
yticks([0,mean(Iavg(successes)),2.*mean(Iavg(successes))]) % ticks are labeled as normalized amplitudes of success trials
yticklabels({'0','1','2'})
set(gca,'box','off','FontSize',20)

% histograms at the top
subplot(5,5,[1:4])
histogram(-Eavg(failures),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor',[0.8 0.8 0.8],'LineWidth',2) % failure trials are in gray
hold on;
histogram(-Eavg(successes),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor','r','LineWidth',2) % success trials are in red
xlim([tempAx(1), tempAx(2)])
yax1=ylim;
xticks([])
set(gca,'box','off','FontSize',20)

% histograms on the right
subplot(5,5,[10,15,20,25])
histogram(Iavg(failures),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor',[0.8 0.8 0.8],'LineWidth',2) % failure trials are in gray
hold on;
histogram(Iavg(successes),'Normalization','count','DisplayStyle','stairs','BinWidth',1,'EdgeColor','b','LineWidth',2) % success trials are in blue
xlim([tempAx(1), tempAx(2)])
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
