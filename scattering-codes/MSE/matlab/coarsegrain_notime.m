function [cgu] = coarsegrain_notime(data,tau)

    % ___________________________________________________________________________
    % Description: 
    % Coarse grain time series using scale tau 
    % ___________________________________________________________________________
    % Input parameters: 
    % data [1D array]     Data array holding time series
    % tau [int]           Scale of decimation
    % ___________________________________________________________________________
    
% Get length of time series
N = length(data) ;
N_tau = fix(N/tau);

% Initialise output array 
cgu = zeros(N_tau-1, 1);

% Coarse grain 
for j=1:N_tau;
    cgu(j) = (1/tau) * sum(data((j-1)*tau +1 : (j)*tau));
  
end 
