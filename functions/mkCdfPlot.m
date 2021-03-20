function fig = mkCdfPlot(Eavg,Iavg,Ethres,Ithres)
% mkCdfPlot function accepts amplitudes of E, I and threshold of each to make cdf plots as shown in Fig 2.e in Kim et al. 2021

Isuccess = find(Iavg>Ithres); % I
Ifail = setdiff(1:length(Iavg),Isuccess); % no I
Esuccess = find(Eavg<Ethres); % E
Efail = setdiff(1:length(Eavg),Esuccess); % no E
failures =intersect(find(Iavg<=Ithres),find(Eavg>=Ethres)); % neither epsc nor ipsc
successes = setdiff(1:length(Iavg),failures); % everything else

f1=cdfplot(Eavg(Isuccess)./mean(Eavg(successes))); % normalized to amplitudes of success trials
f1.Color = 'r';
f1.LineStyle = '-';

hold on;
f2=cdfplot(Eavg(Ifail)./mean(Eavg(successes))); % normalized to amplitudes of success trials
f2.Color = 'r';
f2.LineStyle = ':';
[H1,P1] = kstest2(-Eavg(Isuccess),-Eavg(Ifail));

f3=cdfplot(Iavg(Esuccess)./mean(Iavg(successes))); % normalized to amplitudes of success trials
f3.Color = 'b';
f3.LineStyle = '-';
f4=cdfplot(Iavg(Efail)./mean(Iavg(successes))); % normalized to amplitudes of success trials
f4.Color = 'b';
f4.LineStyle = ':';

xlim([0 3])
xticks([1,2])
xticklabels({'1','2'})
ylabel('cdf')
xlabel('i')
set(gca,'box','off','FontSize',20)
title ''
fig = gcf;
