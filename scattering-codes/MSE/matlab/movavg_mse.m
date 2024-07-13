function [E] = movavg_mse(data,rval,scale, verbose)
    
    % ___________________________________________________________________________
    % Description: 
    % Function calculates MSE for data array using moving average method to
    % compare averaged segments of time series 
    % ___________________________________________________________________________
    % Input parameters: 
    % data [1D array]     Data array holding time series
    % rval [float]        Threshold value
    % scale [int]         Max scale for MSE
    % verbose [str]       Verbose 'yes/no'
    % ___________________________________________________________________________
    % Other functions required for this function: 
    %   moving_mean
    %   SampEn
    % __________________________________________________________________________


    % Standardise the data: mean zero and standard deviation 1
    data = zscore(data); 
    r = rval*std(data); % Note in this case this line doesn't really do anything because the std is 1

    % Define empty arrays
    mmean = [];
    E = [];

    % Calculate moving mean for each scale and then sample entropy 
    for tau = 1:scale;
    mm = movingmean(data, tau);
    E(tau) = SampEn(mm,r,tau);

    % If verbose, print update: 
    if verbose == 'yes'
        disp(strcat('Completed:  ', num2str(tau), '/', num2str(scale)))
    end
    
end
