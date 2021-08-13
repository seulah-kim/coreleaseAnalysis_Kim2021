function [outliers,trialsMout,OutTrialsM] = filterSweepsQC(epoch,trialsM,figureOn)
% This function will remove sweeps that are 3 MAD away from median. It will sort
% sweeps that pass quality check and those that do not into separate
% groups.
% It will then save all the sorted parameters into separate files

% load the sweeps that belong to the specified epoch
load(strcat(['physParamsEpoch',num2str(epoch),'.mat']))

% estimate initialRaccess
initRaccess = median(Raccess(1:length(Raccess)/3));

% apply moving median filter to detect Raccess that exceeds 25% of original
% Raccess
RaOutliers =find(abs(movmedian(Raccess,20)-initRaccess)>(initRaccess.*0.25));

% filter anything that are baseline outliers
persweep=find(abs(mean(trialsM(:,1:end/2),2)-mean(trialsM(:,(end/2+1):end),2))>30); % 30pA difference from 1st half and last half
acrosssweep= find(abs(mean(trialsM,2)-median(mean(trialsM,2)))>30); % 30pA difference from the median baseline value
baselineOutlier = union(persweep,acrosssweep); 

% combine the two criteria
outliers = union(RaOutliers,baselineOutlier);

% update the sweepNums that are included
sweepNum = [1:size(trialsM,1)];
sweepNum(outliers)=[];

% save sweeps that failed QC
OutTrialsM = trialsM(outliers,:);
% and keep the sweeps that passed QC
trialsM(outliers,:)=[];
trialsMout = trialsM;

if figureOn
    figure
    subplot(2,2,1);
    plot(Raccess,'ko');hold on;plot(outliers,Raccess(outliers),'ro');legend('Raccess (MOhm)')
    xlabel('Sweep #')
    ylim([0 50])
    subplot(2,2,2);plot(Rm,'ko');hold on;plot(outliers,Rm(outliers),'ro')
    ylim([0 2000]);
    legend('Rm (MOhm)')
    xlabel('Sweep #')
    subplot(2,2,3);plot(tau,'ko');hold on;plot(outliers,tau(outliers),'ro')
    ylim([0 10]);legend('{\tau}_m (ms)')
    xlabel('Sweep #')
    subplot(2,2,4);plot(Cm,'ko');hold on;plot(outliers,Cm(outliers),'ro')
    ylim([0 50]);legend('Cm (pF)')
    xlabel('Sweep #')

    set(gcf, 'Color', 'w');
end

% params of sweeps that fail QC get saved
failRaccess = Raccess(outliers);
failRm = Rm(outliers);
failtau = tau(outliers);
failCm = Cm(outliers);

save(sprintf('physParamsEpoch%d-failQC.mat',epoch),'failRm','failCm','failtau','failRaccess','outliers');

% params of sweeps that pass QC get saved
Raccess(outliers)=[];
Rm(outliers)=[];
tau(outliers)=[];
Cm(outliers)=[];

save(sprintf('physParamsEpoch%d-passQC.mat',epoch),'Rm','Cm','tau','Raccess','sweepNum');

