function w = PBS( HX,Y,R )
%% Efficient implementation of the Particle Batch Smoother
% presented in Margulis et al. (2015; JHM).
% N.B. The observation errors are assumed to be uncorrelated (diagonal R)
% and Gaussian.
%
% Dimensions: No = Number of observations in the batch to assimilate.
%             Np = Number of parameters to update.
%             Ne = Number of ensemble members (particles).
%
% -----------------------------------------------------------------------
% Inputs:
%
%
% HX   => No x Ne matrix containing an ensemble of Ne predicted
%         observation column vectors each with No entries.
%
% Y     => No x 1 vector containing the batch of (unperturbed) observations.
%
% R     => No x No observation error variance matrix; this may also be
%         specified as a scalar corresponding to the constant variance of
%         all the observations in the case that these are all from the same
%         instrument.
%
% -----------------------------------------------------------------------
% Outputs:
%
% w     => 1 x Ne vector containing the ensemble of posterior weights,
%         the prior weights are implicitly 1/N_e.
%
% -----------------------------------------------------------------------
% See e.g. https://jblevins.org/log/log-sum-exp for underflow issue.
%
% Code by Kristoffer Aalstad (Feb. 2019)



% Calculate the diagonal of the inverse obs. error covariance.
No=size(Y,1);
if numel(R)==No
    if size(R,2)==No
        Rinv=R.^(-1);
    else
        Rinv=(R').^(-1);
    end
elseif numel(R)==1
    Rinv=(1/R).*ones(1,No);
else
    error('Expected numel(R)=No or scalar R')
end


% Calculate the likelihood.
Inn=repmat(Y,1,size(HX,2))-HX;   % Innovation.
EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
LLH=-0.5.*EObj; % log-likelihoods.
normc=logsumexp(LLH,2);

% NB! The likelihood coefficient (1/sqrt(2*pi...)) is
% omitted because it drops out in the normalization
% of the likelihood. Including it (very small term) would lead
% to problems with FP division.


% Calculate the posterior weights as the normalized likelihood.
logw=LLH-normc;
w=exp(logw); % Posterior weights.

% Need "log-sum-exp" trick to overcome numerical issues for small R/large
% number of obs.







    function s = logsumexp(a, dim)
        % Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
        % Default is dim = 1 (columns).
        % logsumexp(a, 2) will sum across rows instead of columns.
        % Unlike matlab's "sum", it will not switch the summing direction
        % if you provide a row vector.
        
        % Written by Tom Minka
        % (c) Microsoft Corporation. All rights reserved.
        
        if nargin < 2
            dim = 1;
        end
        
        % subtract the largest in each column
        [y, i] = max(a,[],dim);
        dims = ones(1,ndims(a));
        dims(dim) = size(a,dim);
        a = a - repmat(y, dims);
        s = y + log(sum(exp(a),dim));
        i = find(~isfinite(y));
        if ~isempty(i)
            s(i) = y(i);
        end
    end






end


