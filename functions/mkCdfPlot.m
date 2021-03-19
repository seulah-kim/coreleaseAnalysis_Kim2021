function fig = mkCdfPlot(Eavg,Iavg,Ethres,Ithres)

Isuccess = find(Iavg>Ithres);
Ifail = setdiff(1:length(Iavg),Isuccess);
Esuccess = find(Eavg<Ethres);
Efail = setdiff(1:length(Eavg),Esuccess);
failures =intersect(find(Iavg<=Ithres),find(Eavg>=Ethres));   % neither epsc nor ipsc
successes = setdiff(1:length(Iavg),failures); % everything else

f1=cdfplot(Eavg(Isuccess)./mean(Eavg(successes)));
f1.Color = 'r';
f1.LineStyle = '-';

hold on;
f2=cdfplot(Eavg(Ifail)./mean(Eavg(successes)));
f2.Color = 'r';
f2.LineStyle = ':';
[H1,P1] = kstest2(-Eavg(Isuccess),-Eavg(Ifail));

f3=cdfplot(Iavg(Esuccess)./mean(Iavg(successes)));
f3.Color = 'b';
f3.LineStyle = '-';
f4=cdfplot(Iavg(Efail)./mean(Iavg(successes)));
f4.Color = 'b';
f4.LineStyle = ':';
[H2,P2] = kstest2(Iavg(Esuccess),Iavg(Efail));
xlim([0 3])
xticks([1,2])
xticklabels({'1','2'})
ylabel('cdf')
xlabel('i')
set(gca,'box','off','FontSize',20)
title ''
fig = gcf;