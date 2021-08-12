function [iNet_cop,iNet_ind,t]=runSimCorelease(nTrials,match_pr,pE,pI)
% This function accepts two input parameters and simulates independent and co-packaging co-release of glutamate and GABA.
% Details of this model is described under Methods section in Kim et al. 2021.
%
% INPUTS:
%   nTrial determines number of trials (observations) per run. 
%   match_pr is a flag used to determine release probabilities of two models
%      match_pr = 1, release probabilites are set to 0.5 (default) for both models
%      match_pr = 0, release probability is set to 0.5 for independent and 0.75 for co-packaging. This is to match the failure rate in both models.
% OUTPUTS:
%   iNet_cop is a matrix of simulated net post-synaptic current based on the co-packaging model. Dimension is (nTrials x time).
%   iNet_ind is a matrix of simulated net post-synaptic current based on the independent model. Dimension is (nTrials x time).
%   t is a time array that matches the column dimension of iNet_cop and iNet_ind.

% parameters 
sigma = 0.05;  % noise level - 5% of signal
simT=500; % how long we stimulate after the NT release
preT=50; % stim Onset 

% make alpha function
TR=0:simT;
nPoints=simT+preT+1;
t = linspace(0,nPoints/10,nPoints);  % timescale, each point is 0.1ms
ratioKinetics = 3; % number of times slower in IPSC
signalSize= 10; % 10pA for both EPSC and IPSC

% simplifed EPSC current (this is equivalent to the normalized version of Eqn. 2) 
Ie=zeros(1,nPoints);
Ie(preT+TR) = TR.*(exp(-(TR)/10));  % perform a unit step function 
m=max(Ie(:)); % this is tau/e 
Ie=-signalSize./m.*Ie;  % we scale by the signal size
% simplifed IPSC current (this is equivalent to the normalized version of Eqn. 2) 
Ii=zeros(1,nPoints);
Ii(preT+TR) = TR.*(exp(-TR/(10*ratioKinetics)));  % perform a unit step function 
m=max(Ii(:)); % this is tau/e 
Ii=signalSize./m.*Ii; % we scale by the signal size
 
EIscale=1.5; % E to I ratio

% initialize
EPSC_cop = zeros(nTrials,length(Ie));
IPSC_cop = zeros(nTrials,length(Ii));
EPSC_ind = zeros(nTrials,length(Ie));
IPSC_ind = zeros(nTrials,length(Ii));

% two alternate simulation types, depending on specific parameters that we
% define
switch nargin
    case 2  % if pE and pI are not specified
        pReleaseInd=0.5;    % release probability of independent release model
        pReleaseCop=1-(1-pReleaseInd).^2;   % release probability of co-packaging release model

        if match_pr % flag for matching pr for two models
            pReleaseCop = pReleaseInd;
        end
        
        for n=1:nTrials
            % co-packing model
            if rand<=pReleaseCop % did release occur?
                sigmaVesicle=1+0.3*randn; % noise in vesicle content, centered around 1, 30% change
                EPSC_cop(n,:)=sigmaVesicle.*EIscale.*Ie;
                IPSC_cop(n,:)=sigmaVesicle.*Ii;
            end
            % independent model
            if rand<=pReleaseInd % did release occur?
                sigmaVesicle=1+0.3.*randn; % noise in vesicle content, centered around 1, 30% change
                EPSC_ind(n,:)=sigmaVesicle.*EIscale.*Ie;
                IPSC_ind(n,:)=sigmaVesicle.*Ii;
            end
        end
    otherwise % if we specify pE and pI
        % assume same for Glu and GABA so that I can run one simulation and compare
        % for co-packaging versus independent release using the same trials
        pReleaseInd_E=pE;
        pReleaseInd_I=pI;
        pReleaseCop=sqrt(pE.*pI);   % to match the p(E) and p(I)
        
        for n=1:nTrials
                % use two dies
                die1 = rand;
                die2 = rand;
              % co-packaging model
              if die1<=pReleaseCop % did release occur?
                  sigmaVesicle=1+0.3*randn; % noise in vesicle content, centered around 1, 30% change
                 
                  EPSC_cop(n,:)=sigmaVesicle.*EIscale.*Ie;
                  IPSC_cop(n,:)=sigmaVesicle.*Ii;
              end

              % independent model - two separate dies are used
              if die1<=pReleaseInd_E
                  sigmaVesicle=1+0.3.*randn; % noise in vesicle content, centered around 1, 30% change

                  EPSC_ind(n,:)=sigmaVesicle.*EIscale.*Ie;
              end

              if die2<=pReleaseInd_I
                  sigmaVesicle=1+0.3.*randn; % noise in vesicle content, centered around 1, 30% change
                  % independent model
                  IPSC_ind(n,:)=sigmaVesicle.*Ii;
              end
         end 
end

% shuffle the order in the IPSC to generate independent case
IPSC_rand=IPSC_ind(randperm(nTrials),:);

% sort trials based on the minimum amplitude size
[~,sortind_cop]=sort(min(EPSC_cop,[],2));
[~,sortind_ind]=sort(min(EPSC_ind,[],2));

% % Add uncorrelated gaussian noise to individual observations
% iNet_cop = EPSC_cop(sortind_cop,:)+IPSC_cop(sortind_cop,:)+sigma.*signalSize.*randn(size(EPSC_cop));
% iNet_ind = EPSC_ind(sortind_ind,:)+IPSC_rand(sortind_ind,:)+sigma.*signalSize.*randn(size(EPSC_ind));

% Add 3kHz corner frequency
fc = 3000;  % cutoff freq
fs = 10000;  % sampling rate
% [b,a] = butter(6,fc/(fs/2));    % 6th order lowpass Butterworth filter
[b,a] = besself(4,fc/(fs/2));   % 4th order (pole) lowpass Bessel filter
noiseIn_cop=randn(size(EPSC_cop));
noiseIn_ind=randn(size(EPSC_ind));
noiseOut_cop=filter(b,a,noiseIn_cop);
noiseOut_ind=filter(b,a,noiseIn_ind);
iNet_cop = EPSC_cop(sortind_cop,:)+IPSC_cop(sortind_cop,:)+sigma.*signalSize.*noiseOut_cop;
iNet_ind = EPSC_ind(sortind_ind,:)+IPSC_rand(sortind_ind,:)+sigma.*signalSize.*noiseOut_ind;



