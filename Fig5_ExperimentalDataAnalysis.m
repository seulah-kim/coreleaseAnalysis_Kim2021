% Fig5_ExperimentalDataAnalysis.m
% 
% This script demonstrates analysis workflow shown in Supplemental Figure 
% 4A from Kim et al. 2021. It reproduces Figures 5 in Kim et al. 2021 using
% two example datasets.

%% load raw data
clear all 
clc
% select main folder where we can find subfolder names & excel file
[~,mainpath] = uiputfile('*.*','Select main data folder', 'mainpath.mat');
% move to folder
cd(mainpath)
% define epoch number
% epoch = 6; % for co-packaging example
epoch = 9; % for independent example

% loads the average trace for the given epoch and saves trial numbers
load(strcat(['AD0_e',num2str(epoch),'p2avg']))
traceNums = evalin('base',['AD0_e',num2str(epoch),'p2avg.UserData.Components']);

% saves pattern sequences for the given epoch
patternSeq=str2num(valueFromHeaderString('state.DMD.patternsString',...
  evalin('base',['AD0_e',num2str(epoch),'p2avg.UserData.headerString'])));

% concatenate matrix of individual traces for the given epoch
trialsM=[];
for i=1:length(traceNums)
    load([traceNums{i},'.mat'])
    trialsM(i,:) = eval([traceNums{i},'.data']);
end

clearvars AD0_*     
%% filter out sweeps that are outside 25% range
clearvars sweepNum

% requires physParamsEpoch*.mat file
if isempty(dir(sprintf('physParamsEpoch%d-passQC.mat',epoch)))
    [outliers,trialsMout] = filterSweepsQC(epoch,trialsM,0);
else
    load(sprintf('physParamsEpoch%d-passQC.mat',epoch),'sweepNum')
    load(sprintf('physParamsEpoch%d-failQC.mat',epoch),'outliers')
    trialsMout = trialsM(sweepNum,:);   
end

%% pre-process
% baseline current subtraction (first 300ms) 
base = mean(trialsMout(:,1:2999),2);
baseM = repmat(base,1,size(trialsMout,2));
trialsM_subtracted = trialsMout-baseM;

% pre-process
y=preprocess(trialsM_subtracted); 

%% concatenate matrix for all putatitve spots
clearvars tlist v MAD correctedSweeps
numPattern= length(patternSeq);

tlist = 300:100:(300+(numPattern-1)*100);   % stimulation times, each targeting different spatial location
intv = [-30 70];  % time window with respect to stim onset
len = diff(intv)*10;

% loop through each stimulation time points
for ii = 1:length(tlist)
    clearvars rawSweeps preprocessedSweeps baseLoc offset
    
    t = tlist(ii); % in ms  
    baseLoc = mean(y(:,((t-30)*10+1):(t)*10),2);   % baseline calculated from 30ms before photostim.
    % outlier traces are removed
    rawSweeps = trialsM(:,[(t+intv(1))*10+1:(t+intv(2))*10]);
    rawSweeps(outliers,:)=[];

    % baseline subtract again to make sure the trace is centered on y axis
    preprocessedSweeps = y(:,[(t+intv(1))*10+1:(t+intv(2))*10])-repmat(baseLoc,1,len);
    % 13ms before photostim as baseline trace
    preprocessedSweeps_base = y(:,[(t+intv(1)-13)*10+1:(t+intv(2))*10])-repmat(baseLoc,1,len+130);
    % measure offset around stimulation period
    offset = mean(preprocessedSweeps(:,abs(intv(1))*10:(abs(intv(1))*10+50)),2);

    % correct the raw traces with offset calculted over stim artifact
    correctedSweeps{ii} = preprocessedSweeps-offset;
    correctedSweeps_base{ii}= preprocessedSweeps_base-offset;
    MAD(ii,:) = mad(correctedSweeps{ii},1);  % calculate mad acros sweeps
end
%% sort data into hotspots vs. null spots
clearvars ipt pkWin_MAD
thresholdFactor=3;
hotspots=find(sum(isoutlier(MAD,'median',1,'ThresholdFactor',thresholdFactor),2)>50); % find spots that exceed thresholdFactor scaled MAD from the median for at least 5ms

%% determine analysis time window per hotspot based on changepoint analysis
clearvars pkWin ipt
cpaWind= 0; % option to widen the analysis window (unit in 0.1 ms), default is 0

for u=1:length(hotspots)
   clearvars timeRange 
    
   % changepoint analysis
   ipt(u,:)=findchangepts([correctedSweeps{hotspots(u)}],'Statistic','rms','MaxNumChanges',2); % rms works better than std or mean
   
   figure
   plot([correctedSweeps{hotspots(u)}]')
   hold on;
   vline(ipt(u,1),'k');vline(ipt(u,2),'k');
   vline(ipt(u,1)-10,':k');vline(ipt(u,2)+10,':k');
   vline(ipt(u,1)-30,':r');vline(ipt(u,2)+30,':r');
   vline(ipt(u,1)-50,':b');vline(ipt(u,2)+50,':b');
   set(gcf,'color','w')
   title(sprintf('spot #%d, cpa window increase =%.2fms',hotspots(u),cpaWind/10))
%    savefig(sprintf('changePoint_spot%d',hotspots(u))) % uncomment to save
%    figures
   close all
   
   ipt(u,1)=ipt(u,1)-cpaWind;   % final analysis time window
   ipt(u,2)=ipt(u,2)+cpaWind;
end

%% build a noise model for a given cell using parametric fitting
% pd1 is a gaussian model fit to null spots
% 
% pd2 is a gaussian model fit to -30ms to -5ms period prior to
% photostimulation onset
% 
% pd3 is a gaussian model fit to pooled datapoints of null spots and -30ms 
% to -5ms period prior to photostimulation onset

% below for pd1
notHotspots = 1:size(MAD);
notHotspots(hotspots)=[];
notHotspotData=[];
for ui = 1:length(notHotspots)  
    notHotspotData=[notHotspotData;correctedSweeps{notHotspots(ui)}]; 
end
notHotspotData=notHotspotData(:,355:end);   % concatenate portion not containing stim artifact
notHotspotData=reshape(notHotspotData,1,[]);    % all datapoints for null spots

figure;subplot(1,3,1);set(gcf,'color','w')
histogram(notHotspotData,'Normalization','pdf')
TF=isoutlier(notHotspotData);
MADthres(1)=3.*1.4826.*median(abs(notHotspotData-median(notHotspotData))); % based on the isoutlier function definition

hold on; histogram(notHotspotData(~TF),'Normalization','pdf');
pd1 = fitdist(notHotspotData(~TF)','Normal');

xrange=[min(notHotspotData),max(notHotspotData)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd1,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('null spots')

% below for pd2
preStimData=[];
for uj = 1:size(MAD)   
    preStimData=[preStimData;correctedSweeps{uj}];    
end

preStimData=preStimData(:,1:250);
preStimData=reshape(preStimData,1,[]);
subplot(1,3,2);histogram(preStimData,'Normalization','pdf')
TF2=isoutlier(preStimData);
MADthres(2)=3.*1.4826.*median(abs(preStimData-median(preStimData))); % based on the isoutlier function definition

hold on; histogram(preStimData(~TF2),'Normalization','pdf');
pd2 = fitdist(preStimData(~TF2)','Normal');
xrange=[min(preStimData),max(preStimData)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd2,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('pre-stim period')

% below for pd3
subplot(1,3,3);
allNull=[notHotspotData(~TF),preStimData(~TF2)];
histogram(allNull,'Normalization','pdf')
pd3 = fitdist(allNull','Normal');
xrange=[min(allNull),max(allNull)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd3,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('pre-stim period + null spots')

% save these parameters for a given cell
save('HotspotTimeWin','ipt','hotspots','cpaWind')
save('NoiseModel','pd1','pd2','pd3')

close all
%% EPSC vs. IPSC peak comparison within the time window defined by changepoint analysis 
clearvars avgEonly avgItoo pEPSC pEgivenI pI pIgivenE pboth

for u=1:length(hotspots)
    % make a separate folder per hotspot
    mkdir(sprintf('coreleaseAnalysis-spot%d',hotspots(u))) 
    
    clearvars E I Eloc Iloc adjInd_E Eavg adjInd_I Iavg tempMat tempMat_base E_null I_null baseline_u elim_u
    tempMat= correctedSweeps{hotspots(u)};
    tempMat_base=correctedSweeps_base{hotspots(u)};
    baseline_u = [mean(tempMat(:,1:100),2),mean(tempMat(:,(end-100):end),2)];
    baseDev=3;  % deviation threshold from subtracted baseline 
    elim_u = find(sum(abs(baseline_u)>baseDev,2)==2);
    tempMat(elim_u,:)=[];
    tempMat_base(elim_u,:)=[];
    
    % calculate baseline peak sizes(-30ms to photostim onset)
    [E_null,E_nullLoc]=min(tempMat(:,1:(abs(intv(1))*10-5)),[],2); 
    [I_null,I_nullLoc]=max(tempMat(:,1:(abs(intv(1))*10-5)),[],2);
    % calculate EPSC and IPSC peak sizes within time window
    [E,Eloc]=min(tempMat(:,ipt(u,1):ipt(u,2)),[],2);
    [I,Iloc]=max(tempMat(:,ipt(u,1):ipt(u,2)),[],2);
    
    % adjusted peak location indices
    adjInd_E=Eloc+ipt(u,1)-1;
    adjInd_I=Iloc+ipt(u,1)-1;
   
    % avg around peak point
    [Eavg,Iavg,E_nullavg,I_nullavg] = extractPk(tempMat,adjInd_E,adjInd_I,tempMat_base,E_nullLoc,I_nullLoc); 
    
    % set threshold based on the noise model (2*standard devation of the
    % symmetric noise)
    load('NoiseModel.mat','pd3')
    Ethres=-2*pd3.sigma;
    Ithres=2*pd3.sigma;

    % calculate the probability of each quadrant
    pEPSC(u)=length(find(Eavg<Ethres))./length(Eavg);  % EPSC peak that is larger or equal to threshold 
    pEgivenI(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(find(Iavg>Ithres)); % conditional probability, given IPSC
    pEgivenNoI(u)=length(intersect(find(Eavg<Ethres),find(Iavg<=Ithres)))./length(find(Iavg<=Ithres)); % conditional probability, given no IPSC
    pIPSC(u)=length(find(Iavg>Ithres))./length(Iavg);  % EPSC peak that is larger or equal to threshold -- for a given data
    pIgivenE(u)=length(intersect(find(Iavg>Ithres),find(Eavg<Ethres)))./length(find(Eavg<Ethres)); % conditional probability, given EPSC
    pIgivenNoE(u)=length(intersect(find(Iavg>Ithres),find(Eavg>=Ethres)))./length(find(Eavg>=Ethres)); % conditional probability, given no EPSC
    pboth(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(Eavg);
    % calculate correlation
    EICorr(u)=corr2(-Eavg,Iavg);
    EICorr_null(u)=corr2(-E_nullavg,I_nullavg);
    
    % same for null distribution
    pEPSCn(u)=length(find(E_nullavg<Ethres))./length(E_nullavg);  % EPSC peak that is larger or equal to threshold
    pEgivenIn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg>Ithres)))./length(find(I_nullavg>Ithres)); % conditional probability, given IPSC
    pEgivenNoIn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg<=Ithres)))./length(find(I_nullavg<=Ithres)); % conditional probability, given no IPSC
    pIPSCn(u)=length(find(I_nullavg>Ithres))./length(I_nullavg);  % EPSC peak that is larger or equal to threshold
    pIgivenEn(u)=length(intersect(find(I_nullavg>Ithres),find(E_nullavg<Ethres)))./length(find(E_nullavg<Ethres)); % conditional probability, given EPSC
    pIgivenNoEn(u)=length(intersect(find(I_nullavg>Ithres),find(E_nullavg>=Ethres)))./length(find(E_nullavg>=Ethres)); % conditional probability, given no EPSC
    pbothn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg>Ithres)))./length(E_nullavg);
  
    % bootstrap probabilities to classify hotspots
    nBoot=10000;
    bootpEPSC=zeros(nBoot,1);
    bootpIPSC=zeros(nBoot,1);
    bootpEPSCn=zeros(nBoot,1);
    bootpIPSCn=zeros(nBoot,1);
    for boot_i=1:nBoot
        clearvars tempSweep tempEavg tempIavg tempE_nullavg tempI_nullavg
        % resample the sweeps, with replacement
        tempSweep=randi(length(Eavg),size(Eavg));
        % resample the Eavg and Iavg based on resampled sweep numbers
        tempEavg=Eavg(tempSweep);
        tempIavg=Iavg(tempSweep);
        tempE_nullavg=E_nullavg(tempSweep);
        tempI_nullavg=I_nullavg(tempSweep);
        
        % compute bootstrapped p(E) and p(I)
        bootpEPSC(boot_i)=length(find(tempEavg<Ethres))./length(Eavg);  % EPSC peak that is larger or equal to threshold 
        bootpIPSC(boot_i)=length(find(tempIavg>Ithres))./length(Iavg);  % EPSC peak that is larger or equal to threshold -- for a given data
        bootpEPSCn(boot_i)=length(find(tempE_nullavg<Ethres))./length(E_nullavg);  % EPSC peak that is larger or equal to threshold
        bootpIPSCn(boot_i)=length(find(tempI_nullavg>Ithres))./length(I_nullavg);  % EPSC peak that is larger or equal to threshold
    end
        
    % calculate a 95% confidence interval for the simulation:
    myAlpha = 0.05;     % MATLAB convention for determining 95% CI
    idxHi = ceil(nBoot * (1 - myAlpha/2));
    idxLo = floor(nBoot * (myAlpha/2));
    sortedbootpEPSCn = sort(bootpEPSCn);
    sortedbootpIPSCn = sort(bootpIPSCn);
    bootCI_pEPSCn(u,:) = [sortedbootpEPSCn(idxLo),sortedbootpEPSCn(idxHi)];
    bootCI_pIPSCn(u,:) = [sortedbootpIPSCn(idxLo),sortedbootpIPSCn(idxHi)];
    
    % classify hotspots and their subtypes based on the confidence interval
    evokedEsite(u)=pEPSC(u)>bootCI_pEPSCn(u,2);
    evokedIsite(u)=pIPSC(u)>bootCI_pIPSCn(u,2);
          
end 

% save these parameters of all hotspots for a given cell
save('statistics','pEPSCn','pEPSC','pEgivenIn','pEgivenI','pEgivenNoIn','pEgivenNoI','pIPSCn','pIPSC','pIgivenEn','pIgivenE','pIgivenNoEn','pIgivenNoE','EICorr','EICorr_null',...
    'bootCI_pEPSCn','bootCI_pIPSCn','evokedEsite','evokedIsite')

%% testing statistical independence - shuffled bootstraping and model simulations

% load the active hotspots, already sorted out failure spots
load('HotspotTimeWin.mat')

allspots=[];sitenum=[];allpEPSC=[];allpIPSC=[];

% load the params
load('HotspotTimeWin.mat','hotspots')
load('statistics.mat','evokedEsite','evokedIsite','pEPSC','pIPSC')
temp=[evokedEsite;evokedIsite];
allspots=[allspots,temp];
allpEPSC=[allpEPSC,pEPSC];
allpIPSC=[allpIPSC,pIPSC];
sitenum=[sitenum,hotspots'];

Eonly=allspots(1,:)==1 & allspots(2,:)==0 & allpEPSC<0.95;  % updated this on 6/4/20 but should not impact Fig 4
Ionly=allspots(2,:)==1 & allspots(1,:)==0 & allpIPSC<0.95;
CoRelease=allspots(2,:)==1 & allspots(1,:)==1 & allpEPSC<0.95 & allpIPSC<0.95 & ipt(:,1)'>301;    % added ipt lower limit 
spotNum=sitenum(find(CoRelease));
iptLen = diff(ipt,[],2);
iptLen = iptLen(find(CoRelease));
save('minimalCoReleaseSites','spotNum','iptLen');
counter_spoti =1;
delEcdf50=[];
delIcdf50=[];

% do the following analysis for co-release hotspots 
for spot_i = find(CoRelease)

    % find u for given hotspot
    u=spot_i;

    mkdir(sprintf('coreleaseAnalysis-spot%d',hotspots(u)))
    % compute Eavg and Iavg
    clearvars E I Eloc Iloc adjInd_E Eavg adjInd_I Iavg tempMat tempMat_base E_null I_null baseline_u elim_u
    tempMat= correctedSweeps{hotspots(u)};
    tempMat_base=correctedSweeps_base{hotspots(u)};
    baseline_u = [median(tempMat(:,1:100),2),median(tempMat(:,(end-100):end),2)];   % changed to median from mean to avoid case where spontaneous activity drives the mean value
    baseDev=3;  % deviation threshold from subtracted baseline -- previously 1.5, 1.30.20
    elim_u = find(sum(abs(baseline_u)>baseDev,2)==2);
    tempMat(elim_u,:)=[];
    tempMat_base(elim_u,:)=[];

    % calculate baseline peak sizes(-30ms to photostim onset)
    [E_null,E_nullLoc]=min(tempMat(:,1:(abs(intv(1))*10-5)),[],2);  
    [I_null,I_nullLoc]=max(tempMat(:,1:(abs(intv(1))*10-5)),[],2);

    % calculate EPSC and IPSC peak sizes within time window
    [E,Eloc]=min(tempMat(:,ipt(u,1):ipt(u,2)),[],2);
    [I,Iloc]=max(tempMat(:,ipt(u,1):ipt(u,2)),[],2);

    adjInd_E=Eloc+ipt(u,1)-1;
    adjInd_I=Iloc+ipt(u,1)-1;

    clearvars *base_beforepk
    % avg around peak point
    [Eavg,Iavg,~,~] = extractPk(tempMat,adjInd_E,adjInd_I,tempMat_base,E_nullLoc,I_nullLoc); 

    % set threshold based on the noise model (2*standard devation of the
    % symmetric noise)
    load('NoiseModel.mat','pd3')
    Ethres=-2*pd3.sigma;
    Ithres=2*pd3.sigma;

    % calculate the probability of each quadrant
    pEPSC(u)=length(find(Eavg<Ethres))./length(Eavg);  % EPSC peak that is larger or equal to threshold 
    pEgivenI(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(find(Iavg>Ithres)); % conditional probability, given IPSC
    pEgivenNoI(u)=length(intersect(find(Eavg<Ethres),find(Iavg<=Ithres)))./length(find(Iavg<=Ithres)); % conditional probability, given no IPSC
    pIPSC(u)=length(find(Iavg>Ithres))./length(Iavg);  % EPSC peak that is larger or equal to threshold -- for a given data
    pIgivenE(u)=length(intersect(find(Iavg>Ithres),find(Eavg<Ethres)))./length(find(Eavg<Ethres)); % conditional probability, given EPSC
    pIgivenNoE(u)=length(intersect(find(Iavg>Ithres),find(Eavg>=Ethres)))./length(find(Eavg>=Ethres)); % conditional probability, given no EPSC
    pboth(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(Eavg);
    
    % trial classification
    failures =intersect(find(Iavg<=Ithres),find(Eavg>=Ethres));   % neither epsc nor ipsc
    successes = setdiff(1:length(Iavg),failures); % everything else
    nSuccess = length(successes);
    nTrials(counter_spoti) = length(Eavg);
    pSuccess = nSuccess./nTrials(counter_spoti);

    if nSuccess <10
        fprintf('SKIPPED spot %d : did not meet the minimum nunber of success trials\n',hotspots(u))
        break
    end    

    rng(0)  % for repeatability

    if length(successes)<9
        % show less traces
        sweepn=[successes,randsample(failures,3)];    % to show enough of success trials in low pr
        sweepn=sweepn(randperm(length(sweepn)));
    elseif length(failures)<3
        % show less traces
        sweepn=[randsample(successes,9),failures];    % to show enough of success trials in low pr
        sweepn=sweepn(randperm(length(sweepn)));
    else
        sweepn=[randsample(successes,9),randsample(failures,3)];    % to show enough of success trials in low pr
        sweepn=sweepn(randperm(length(sweepn)));
    end
      
    %% figure 5a: example of peak detection across sweeps
    xx =linspace(intv(1),intv(2),length(tempMat));
    x0=0;
    y0=20;
    
    figure('Position',[547 586 256 518]);set(gcf,'color','w')
    subplot(6,1,[1:5])
    for (counter=1:length(sweepn))
        plot(xx+x0*(counter-1),tempMat(sweepn(counter),:)+y0*(counter-1),'k');
        hold on;
        plot(adjInd_E(sweepn(counter))/10+intv(1)+x0*(counter-1),Eavg(sweepn(counter))+y0*(counter-1),'ro');
        plot(adjInd_I(sweepn(counter))/10+intv(1)+x0*(counter-1),Iavg(sweepn(counter))+y0*(counter-1),'bo');
    end
    Yrange=ylim;

    fill([ipt(u,1)/10+intv(1),ipt(u,1)/10+intv(1),ipt(u,2)/10+intv(1),ipt(u,2)/10+intv(1)],...
    [Yrange(1),Yrange(2),Yrange(2),Yrange(1)],'k','FaceAlpha',0.1,'EdgeAlpha',0);   % gray shaded area over analysis window

    sweepNum=traceNums{1};
    sweepNum(3)='3';    % AD3 is the channel for phodiode
    laser=load(sweepNum);
    laserOn=find(laser.(sweepNum).data(((hotspots(u)+2)*1000):(hotspots(u)+2)*1000+1000)>1);

    fill([min(laserOn)*0.1,min(laserOn)*0.1,max(laserOn)*0.1,max(laserOn)*0.1],...
    [Yrange(1),Yrange(2),Yrange(2),Yrange(1)],'b','FaceAlpha',0.1,'EdgeAlpha',0);   % blue shaded area over photostimulation

    sb=scalebar;    % put a scalebar
    sb.YLen=10;
    sb.XLen=10;

    set(gca,'box','off','FontSize',20)
    axis off
    Xrange=xlim;

    subplot(6,1,6)
    histogram(adjInd_E/10+intv(1),'Normalization','probability','DisplayStyle','stairs','BinWidth',3,'EdgeColor','r','LineWidth',1)
    hold on;
    histogram(adjInd_I/10+intv(1),'Normalization','probability','DisplayStyle','stairs','BinWidth',3,'EdgeColor','b','LineWidth',1)
    xlim([Xrange])
    axis off

    saveas(gcf,sprintf('coreleaseAnalysis-spot%d/Fig5a.pdf',hotspots(u)));
    %% figure 5b: scatterplots
    figure('Position',[654 850 679 590]);set(gcf,'color','w')
    subplot(5,5,[6:9,11:14,16:19,21:24])          

    fig1=mkScatterPlot(Eavg,Iavg,Ethres,Ithres);
    saveas(fig1,sprintf('coreleaseAnalysis-spot%d/Fig5b_scatter.pdf',hotspots(u)));
    
    close all
    %% figure 5d: distribution of imin and imax (w and w/o EPSC and IPSC)
    figure('Position',[654 850 300 300]); set(gcf,'color','w')         
    fig2=mkCdfPlot(Eavg,Iavg,Ethres,Ithres);
    saveas(gcf,sprintf('coreleaseAnalysis-spot%d/Fig5d_cdf.pdf',hotspots(u)));
    close all
    %% figure 5c: testing statistical independence and run simulation
    clearvars pboth_temp pEPSC_temp pIPSC_temp pboth_temp_null pIPSC_temp_null EICorr_s EICorr_all EICorrpval_s EICorrpval_all EICorrnull
    nboot=10000;
    
    % bootstrapping, with replacement
    rng default
    TrialLen = length(Eavg);

    parfor boot_i=1:nboot % changed to parfor
        rng(boot_i) % for repeatability

        Eavg_temp=Eavg;
        Iavg_temp=Iavg;

        bootidx=randi(length(Eavg),size(Eavg)); % with replacement

        % bootstrapped sample, with matching n, capturing variability in
        % the data
        Eavg_temp=Eavg_temp(bootidx);
        Iavg_temp=Iavg_temp(bootidx);
        % shuffle Iavg_temp with respect to Eavg_temp
        Iavg_temp_shuffled=Iavg_temp(randperm(length(Iavg_temp)));

        % calculate the probability of each quadrant
        pEPSC_temp(boot_i)=length(find(Eavg_temp<Ethres))./length(Eavg_temp);  % EPSC peak that is larger or equal to threshold 
        pIPSC_temp(boot_i)=length(find(Iavg_temp>Ithres))./length(Iavg_temp);  % EPSC peak that is larger or equal to threshold -- for a given data
        pIPSC_temp_null(boot_i)=length(find(Iavg_temp_shuffled>Ithres))./length(Iavg_temp_shuffled); 
        pboth_temp(boot_i)=length(intersect(find(Eavg_temp<Ethres),find(Iavg_temp>Ithres)))./length(Eavg_temp);
        pboth_temp_null(boot_i)=length(intersect(find(Eavg_temp<Ethres),find(Iavg_temp_shuffled>Ithres)))./length(Eavg_temp);

        % calculate correlation
        failures =intersect(find(Iavg_temp<=Ithres),find(Eavg_temp>=Ethres));   % neither epsc nor ipsc
        successes = setdiff(1:length(Iavg_temp),failures); % everything else 

        [r,p]=corrcoef(-Eavg_temp(successes),Iavg_temp(successes));   % for successes
        EICorr_s(boot_i)=r(2,1);
        EICorrpval_s(boot_i)=p(2,1);
        [r2,p2]=corrcoef(-Eavg_temp,Iavg_temp);   % for all 
        EICorr_all(boot_i)=r2(2,1);
        EICorrpval_all(boot_i)=p2(2,1);

        subset = Iavg_temp(successes);
        EICorrnull(boot_i)=corr2(-Eavg_temp(successes),subset(randperm(length(successes))));  % shuffling of success

    end

    % run simulation using p(E), p(I), and nSweep  
    rng default

    % binomial model
    [pEPSC_cop,pIPSC_cop,pboth_cop,EICorr_cop_all,EICorr_cop_s,EICorrpval_cop_all,EICorrpval_cop_s,EICorrnull_cop,...
     pEPSC_ind,pIPSC_ind,pboth_ind,EICorr_ind_all,EICorr_ind_s,EICorrpval_ind_all,EICorrpval_ind_s,EICorrnull_ind]=repeatSimAnalysis(500,[],pEPSC(u),pIPSC(u),TrialLen);    

    pProd_cop=pEPSC_cop.*pIPSC_cop; % pE * pI distribution for co-packaging
    pProd_ind=pEPSC_ind.*pIPSC_ind; % pE * pI distribution for independent
    
    save(sprintf('coreleaseAnalysis-spot%d/simResults_confProb.mat',hotspots(u)),'pEPSC_cop','pIPSC_cop','pboth_cop','EICorr_cop_all','EICorr_cop_s','EICorrpval_cop_all','EICorrpval_cop_s','EICorrnull_cop',...
       'pEPSC_ind','pIPSC_ind','pboth_ind','EICorr_ind_all','EICorr_ind_s','EICorrpval_ind_all','EICorrpval_ind_s','EICorrnull_ind')
    save(sprintf('coreleaseAnalysis-spot%d/bootstrapResults.mat',hotspots(u)),'pEPSC_temp','pIPSC_temp','pboth_temp','pboth_temp_null','pIPSC_temp_null','EICorr_all','EICorr_s','EICorrnull')

    % plot results    
    figure('Position',[79 826 1287 460]); set(gcf,'color','w')
    subplot(2,2,1)
    histogram(pboth_temp,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(pEPSC_temp.*pIPSC_temp,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',3,'LineStyle','-');xlim([0 1])
    ylabel('prob')
    legend('p(E,I)','p(E)p(I)')
    legend boxoff
    yticks([0 ,0.4])
    xticks([-0.5, 0, 0.5])
    ylim([0 ,0.4])
    set(gca,'box','off','FontSize',20)
    title('data')

    subplot(2,2,3)
    histogram(pboth_temp_null,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(pEPSC_temp.*pIPSC_temp_null,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-');xlim([0 1]);
    ylabel('prob')
    yticks([0 ,0.4])
    xticks([-0.5, 0, 0.5])
    ylim([0 ,0.4])
    set(gca,'box','off','FontSize',20)
    title('shuffled data')

    subplot(2,2,4)
    histogram(pboth_cop,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(pEPSC_cop.*pIPSC_cop,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',3,'LineStyle','-');xlim([0 1])
    ylabel('prob')
    yticks([0 ,0.4])
    xticks([-0.5, 0, 0.5])
    ylim([0 ,0.4])
    set(gca,'box','off','FontSize',20)
    title('copackage model')

    subplot(2,2,2)
    histogram(pboth_ind,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(pEPSC_ind.*pIPSC_ind,'Normalization','probability','BinWidth',0.03,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-');xlim([0 1]);
    ylabel('prob')
    yticks([0 ,0.4])
    xticks([-0.5, 0, 0.5])
    ylim([0 ,0.4])
    set(gca,'box','off','FontSize',20)
    title('independent model')
    saveas(gcf,sprintf('coreleaseAnalysis-spot%d/Fig5c_prob.pdf',hotspots(u)));

    %% figure 5e: correlation plot

    figure('Position',[654 850 1287 460]);set(gcf,'color','w');

    subplot(1,3,1)
    histogram(EICorr_all,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(EICorr_s,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2);
    histogram(EICorrnull,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.85 0.85 0.85],'LineWidth',2,'LineStyle','-');xlim([-1 1]);

    ylabel('prob')
    xlabel('correlation')
    yticks([0 ,0.6])
    xticks([-1, 0, 1])
    ylim([0 ,0.6])
    xlim([-1,1])
    set(gca,'box','off','FontSize',20)
    title('data')

    subplot(1,3,3)
    histogram(EICorr_cop_all,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(EICorr_cop_s,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2);
    histogram(EICorrnull_cop,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.85 0.85 0.85],'LineWidth',2,'LineStyle','-');xlim([-1 1]);

    ylabel('prob')
    xlabel('correlation')
    yticks([0 ,0.6])
    xticks([-1, 0, 1])
    ylim([0 ,0.6])
    xlim([-1,1])
    set(gca,'box','off','FontSize',20)
    title('co-packaging')

    subplot(1,3,2)
    histogram(EICorr_ind_all,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor','k','LineStyle','-','LineWidth',2);hold on;
    histogram(EICorr_ind_s,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2)
    histogram(EICorrnull_ind,'Normalization','probability','BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',[0.85 0.85 0.85],'LineWidth',2); xlim([-1 1]);
    ylabel('prob')
    legend('all','successes','shuffled')
    legend boxoff
    xlabel('correlation')
    yticks([0 ,0.6])
    xticks([-1, 0, 1])
    ylim([0 ,0.6])
    xlim([-1,1])
    set(gca,'box','off','FontSize',20)
    title('independent')
    saveas(gcf,sprintf('coreleaseAnalysis-spot%d/Fig5e_corr.pdf',hotspots(u)));

    counter_spoti = counter_spoti+1;   
end
 