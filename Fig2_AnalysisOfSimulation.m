% Fig2_AnalysisOfSimulation.m
% 
% This script reproduces figures 2.c-f published in Kim et al. 2021.
%

clear all
% parameters
nTrials=200; % number of trials
match_pr=0; % success rate is matched for two models for visualization of scatter plot

% simulate two release models
rng default % for reproducibility
rng(1)
[iNet_cop,iNet_ind,t]=runSimCorelease(nTrials,match_pr);
     
% pre-process
iNet_ind=preprocess(iNet_ind); % input matrix (observations x time)
iNet_cop=preprocess(iNet_cop);

% concatenate
len= 150;
preT=50; % stim Onset (matches param in runSimCorelease.m)
t=t(:,(preT-10):(preT-10+len));
iNet_ind=iNet_ind(:,(preT-10):(preT-10+len));
iNet_cop=iNet_cop(:,(preT-10):(preT-10+len));

% Trial-by-trial peak extraction 
% for co-packaging case
[~,ElocCop]=min(iNet_cop,[],2);
[~,IlocCop]=max(iNet_cop,[],2);      
% avg around peak point
for sweep_i = 1:size(iNet_cop,1)
    % to better estimate peak size post-stim period, average 0.5ms before and after peak and subtract from 10ms before 3ms of the peak
    Eavg_cop(sweep_i) = mean(iNet_cop(sweep_i,(max(ElocCop(sweep_i)-5,1)):(min(ElocCop(sweep_i)+5,size(iNet_cop,2)))),2);
    Iavg_cop(sweep_i) = mean(iNet_cop(sweep_i,(max(IlocCop(sweep_i)-5,1)):(min(IlocCop(sweep_i)+5,size(iNet_cop,2)))),2);

end

% for independent case
[~,ElocInd]=min(iNet_ind,[],2);
[~,IlocInd]=max(iNet_ind,[],2);
% avg around peak point
for sweep_i = 1:size(iNet_ind,1)
    % to better estimate peak size post-stim period, average 0.5ms before and after peak and subtract from 10ms before 3ms of the peak
    Eavg_ind(sweep_i) = mean(iNet_ind(sweep_i,(max(ElocInd(sweep_i)-5,1)):(min(ElocInd(sweep_i)+5,size(iNet_ind,2)))),2);
    Iavg_ind(sweep_i) = mean(iNet_ind(sweep_i,(max(IlocInd(sweep_i)-5,1)):(min(IlocInd(sweep_i)+5,size(iNet_ind,2)))),2);
end

%% Figure 2.c: scatter plot

sigma = 0.05;  % noise level - 3% of signal
signalSize= 10; % 20pA for both EPSC and IPSC

% threshold
Ethres=-2*sigma*signalSize;
Ithres=2*sigma*signalSize;

% For Independent model
fig1=mkScatterPlot(Eavg_ind,Iavg_ind,Ethres,Ithres);
saveas(fig1,'simulatedFig2c_scatterInd.pdf')

% For Co-package model
fig2=mkScatterPlot(Eavg_cop,Iavg_cop,Ethres,Ithres);
saveas(fig2,'simulatedFig2c_scatterCo.pdf')

%% Figure 2.e: cdf of amplitude distribution

figure('Position',[85 653 700 241]); set(gcf,'color','w')
% For Independent model
subplot(1,2,1)
fig2=mkCdfPlot(Eavg_ind,Iavg_ind,Ethres,Ithres);
title('independent')
% For Co-package model
subplot(1,2,2)
mkCdfPlot(Eavg_cop,Iavg_cop,Ethres,Ithres);
title('co-packaging')

saveas(fig2,'simulatedFig2e.pdf')

%% Figure 2.d & f: repeating simulations (nTrial=200) for 1000 runs
% this section could be skipped if simResults.mat file already exists

clear all
rng default % for reproducibility
[pEPSC_cop,pIPSC_cop,pboth_cop,EICorr_cop_all,EICorr_cop_s,EICorrpval_cop_all,EICorrpval_cop_s,EICorrnull_cop,diffVar_cop,...
   pEPSC_ind,pIPSC_ind,pboth_ind,EICorr_ind_all,EICorr_ind_s,EICorrpval_ind_all,EICorrpval_ind_s,EICorrnull_ind,diffVar_ind]=repeatSimAnalysis(1000);    % 200 nTrials per simulation

% distributions predicted by statistical independence
pProd_cop=pEPSC_cop.*pIPSC_cop;
pProd_ind=pEPSC_ind.*pIPSC_ind;

save('simResults.mat')  % save all the variables from simulations
%% Figure 2.d: distribution of joint probability and that predicted by statistical independence
clear all
load('simResults.mat')  % load the .mat file if this already exists in the directory

% We plot below
figure('Position',[654 850 1015 346]);set(gcf,'color','w');
subplot(1,2,1)
set(gcf,'color','w');histogram(pProd_ind,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.8 0.8 0.8],'LineWidth',2,'BinWidth',0.02);hold on; 
histogram(pboth_ind,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2,'BinWidth',0.02);xlim([0 0.65]);
legend('p(E)p(I)','p(E,I)')
ylabel('prob')
ylim([0 0.3])
yticks([0 ,0.3])
xticks([0, 0.5])
set(gca,'box','off','FontSize',20)
title('independent')

subplot(1,2,2)
histogram(pProd_cop,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.8 0.8 0.8],'LineWidth',2,'BinWidth',0.02);hold on;
histogram(pboth_cop,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2,'BinWidth',0.02);xlim([0 0.65]);
legend('p(E)p(I)','p(E,I)')
ylabel('prob')
ylim([0 0.3])
yticks([0 ,0.3])
xticks([0, 0.5])
set(gca,'box','off','FontSize',20)
title('co-packaging')

% saveas(gcf,'simulatedFig1d_prob.pdf')
%% Figure 2.f: correlation distributions

% We plot below
figure('Position',[654 850 1015 346]);set(gcf,'color','w');
subplot(1,2,1)
set(gcf,'color','w');histogram(EICorr_ind_all,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineStyle','-','BinWidth',0.05,'LineWidth',2);hold on;
histogram(EICorr_ind_s,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineStyle','-','BinWidth',0.05,'LineWidth',2)
histogram(EICorrnull_ind,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.85 0.85 0.85],'BinWidth',0.05,'LineWidth',2); xlim([-1 1]);
ylabel('prob')
legend('all','successes','shuffled')
xlabel('correlation')
yticks([0 ,1])
xticks([-0.5, 0, 0.5])
ylim([0 ,1])
set(gca,'box','off','FontSize',20)
title('independent')

subplot(1,2,2)
histogram(EICorr_cop_all,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineStyle','-','BinWidth',0.05,'LineWidth',2);hold on;
histogram(EICorr_cop_s,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineStyle','-','BinWidth',0.05,'LineWidth',2);
histogram(EICorrnull_cop,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.85 0.85 0.85],'LineWidth',2,'BinWidth',0.05,'LineStyle','-');xlim([-1 1]);
ylabel('prob')
legend('all','successes','shuffled')
xlabel('correlation')
yticks([0 ,1])
xticks([-0.5, 0, 0.5])
ylim([0 ,1])
set(gca,'box','off','FontSize',20)
title('co-packaging')

% saveas(gcf,'simulatedFig1f_corr.pdf')