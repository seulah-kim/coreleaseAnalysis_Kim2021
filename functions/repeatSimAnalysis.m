function [pEPSC_cop,pIPSC_cop,pboth_cop,EICorr_cop_all,EICorr_cop_s,EICorrpval_cop_all,EICorrpval_cop_s,EICorrnull_cop,...
   pEPSC_ind,pIPSC_ind,pboth_ind,EICorr_ind_all,EICorr_ind_s,EICorrpval_ind_all,EICorrpval_ind_s,EICorrnull_ind]=repeatSimAnalysis(varargin)
% This function performs nSim runs of simulations of glutamate/GABA corelease under independent and co-packaging release models.
% Outputs variables include probability and correlation distributions resulting from nSim runs of simulation, using the two models.

nSim=varargin{1};
if length(varargin)>1
    match_pr=varargin{2}; 
    pE = varargin{3};
    pI = varargin{4};
    nTrials = varargin{5};
else
    nTrials=200; % number of trials
    match_pr =1; % we match 
    pE = NaN;
    pI = NaN;
end
nInputs = nargin;   % number of input variables

parfor sim_i = 1:nSim
    rng(sim_i)  % for reproducibility   

    % simulate two release models
    switch nInputs
        case 1  % for simulation with only run number indicated 
            [iNet_cop,iNet_ind,t]=runSimCorelease(nTrials,match_pr);       
        otherwise % for simulation with other parameters specified
            [iNet_cop,iNet_ind,t]=runSimCorelease(nTrials,[],pE,pI);
    end
    % pre-process
    iNet_ind=preprocess(iNet_ind); % input matrix (observations x time)
    iNet_cop=preprocess(iNet_cop);

    % concatenate
    len= 150;
    preT=50; % stim Onset (matches param in runSimCorelease.m)
    t=t(:,(preT-10):(preT-10+len));
    iNet_ind=iNet_ind(:,(preT-10):(preT-10+len));
    iNet_cop=iNet_cop(:,(preT-10):(preT-10+len));

    %---------------- Trial-by-trial peak extraction ---------------------
    % for co-packaging case
    [~,ElocCop]=min(iNet_cop,[],2);
    [~,IlocCop]=max(iNet_cop,[],2);      
    
    % pre-allocating allows parfor
    Eavg_cop = NaN(1,size(iNet_cop,1));  
    Iavg_cop = NaN(1,size(iNet_cop,1)); 
    % avg around peak point
    for sweep_i = 1:size(iNet_cop,1)
        % to better estimate peak size post-stim period, average 0.5ms before and after peak and subtract from 10ms before 3ms of the peak
        Eavg_cop(sweep_i) = mean(iNet_cop(sweep_i,(max(ElocCop(sweep_i)-5,1)):(min(ElocCop(sweep_i)+5,size(iNet_cop,2)))),2);
        Iavg_cop(sweep_i) = mean(iNet_cop(sweep_i,(max(IlocCop(sweep_i)-5,1)):(min(IlocCop(sweep_i)+5,size(iNet_cop,2)))),2);

    end
    
    % pre-allocating allows parfor
    Eavg_ind = NaN(1,size(iNet_ind,1)); 
    Iavg_ind = NaN(1,size(iNet_ind,1));
    % for independent case
    [~,ElocInd]=min(iNet_ind,[],2);
    [~,IlocInd]=max(iNet_ind,[],2);
    % avg around peak point
    for sweep_i = 1:size(iNet_ind,1)
        % to better estimate peak size post-stim period, average 0.5ms before and after peak and subtract from 10ms before 3ms of the peak
        Eavg_ind(sweep_i) = mean(iNet_ind(sweep_i,(max(ElocInd(sweep_i)-5,1)):(min(ElocInd(sweep_i)+5,size(iNet_ind,2)))),2);
        Iavg_ind(sweep_i) = mean(iNet_ind(sweep_i,(max(IlocInd(sweep_i)-5,1)):(min(IlocInd(sweep_i)+5,size(iNet_ind,2)))),2);
    end
    
    %--------- Calculation of probabilities and correlations  --------------
    sigma = 0.05;  % noise level - 3% of signal
    signalSize= 10; % 20pA for both EPSC and IPSC

    % threshold
    Ethres=-2*sigma*signalSize;
    Ithres=2*sigma*signalSize;
    
    % For Co-package model
    failures_cop =intersect(find(Iavg_cop<=Ithres),find(Eavg_cop>=Ethres));   % neither epsc nor ipsc
    successes_cop = setdiff(1:length(Iavg_cop),failures_cop); % everything else 
    
    % calculate the probability of each quadrant
    pEPSC_cop(sim_i)=length(find(Eavg_cop<Ethres))./length(Eavg_cop);  % EPSC peak that is larger or equal to threshold 
    pIPSC_cop(sim_i)=length(find(Iavg_cop>Ithres))./length(Iavg_cop);  % EPSC peak that is larger or equal to threshold -- for a given data
    pboth_cop(sim_i)=length(intersect(find(Eavg_cop<Ethres),find(Iavg_cop>Ithres)))./length(Eavg_cop);
    
    % calculate correlation
    [r,p]=corrcoef(-Eavg_cop(successes_cop),Iavg_cop(successes_cop));   % for successes
    EICorr_cop_s(sim_i)=r(2,1);
    EICorrpval_cop_s(sim_i)=p(2,1);
    [r2,p2]=corrcoef(-Eavg_cop,Iavg_cop);   % for all 
    EICorr_cop_all(sim_i)=r2(2,1);
    EICorrpval_cop_all(sim_i)=p2(2,1);

    subset = Iavg_cop(successes_cop);
    EICorrnull_cop(sim_i)=corr2(-Eavg_cop(successes_cop),subset(randperm(length(successes_cop))));  % shuffling of success

%     % calculate variance difference
%     Isuccess = find(Iavg_cop>Ithres);
%     Ifail = setdiff(1:length(Iavg_cop),Isuccess);
%     Esuccess = find(Eavg_cop<Ethres);
%     Efail = setdiff(1:length(Eavg_cop),Esuccess);
%     bothsuccess = intersect(Isuccess,Esuccess); 
%     Eonly=intersect(Esuccess,Ifail);
%     Ionly=intersect(Efail,Isuccess);
%     diffVar_cop(sim_i)=var(Eavg_cop(bothsuccess),Iavg_cop(bothsuccess))-var(Eavg_cop(Eonly))-var(Iavg_cop(Ionly));
%     
    % for Independent model
    failures_ind =intersect(find(Iavg_ind<=Ithres),find(Eavg_ind>=Ethres));   % neither epsc nor ipsc
    successes_ind = setdiff(1:length(Iavg_ind),failures_ind); % everything else
    
    % calculate the probability of each quadrant
    pEPSC_ind(sim_i)=length(find(Eavg_ind<Ethres))./length(Eavg_ind);  % EPSC peak that is larger or equal to threshold 
    pIPSC_ind(sim_i)=length(find(Iavg_ind>Ithres))./length(Iavg_ind);  % EPSC peak that is larger or equal to threshold -- for a given data
    pboth_ind(sim_i)=length(intersect(find(Eavg_ind<Ethres),find(Iavg_ind>Ithres)))./length(Eavg_ind);
    
    % calculate correlation
    [r,p]=corrcoef(-Eavg_ind(successes_ind),Iavg_ind(successes_ind));   % for successes
    EICorr_ind_s(sim_i)=r(2,1);
    EICorrpval_ind_s(sim_i)=p(2,1);
    [r2,p2]=corrcoef(-Eavg_ind,Iavg_ind);   % for all
    EICorr_ind_all(sim_i)=r2(2,1);
    EICorrpval_ind_all(sim_i)=p2(2,1);
    
    subset = Iavg_ind(successes_ind);
    EICorrnull_ind(sim_i)=corr2(-Eavg_ind(successes_ind),subset(randperm(length(successes_ind))));    % correlation of shuffled data 
    
%     % calculate variance difference
%     Isuccess = find(Iavg_ind>Ithres);
%     Ifail = setdiff(1:length(Iavg_ind),Isuccess);
%     Esuccess = find(Eavg_ind<Ethres);
%     Efail = setdiff(1:length(Eavg_ind),Esuccess);
%     bothsuccess = intersect(Isuccess,Esuccess);
%     Eonly=intersect(Esuccess,Ifail);
%     Ionly=intersect(Efail,Isuccess);
%     diffVar_ind(sim_i)=var(Eavg_ind(bothsuccess),Iavg_ind(bothsuccess))-var(Eavg_ind(Eonly))-var(Iavg_ind(Ionly));

end
