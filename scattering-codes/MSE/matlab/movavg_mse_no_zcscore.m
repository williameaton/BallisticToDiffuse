%This function calls moving average and multi-scale sample entropy
%inputs: data: data vector
%        scale: scale numbers
%output: E: entropy vector for all the scale factors
%
function [E] = movavg_mse_no_zcscore(data,rval,scale)

r = rval;  %define threshold -- this is actually meaningless as the data has standard deviation 1
		     %therefore, rval is going to be the threshold value
             
mmean = [];
E = [];
for tau = 1:scale;
 cg = movingmean(data, tau);
E(tau) = SampEn(cg,r,tau);
end
