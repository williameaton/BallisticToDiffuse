function [] = moving_mse(r, tau, scale, input_fn, output_dir, output_fn)



% Create output file file:
mkdir(output_dir);
% Output path and name
output_name = strcat(output_dir, '/', output_fn);




% Loop through channels
for ch = 0:2


temp_input_fn = strcat('SLICES_', input_fn, '_ch', num2str(ch), '.txt');

% Verbose load data 
disp('Loading data from')
disp(temp_input_fn)
data = textread(temp_input_fn); 
disp('Load completed.') 


size_data = size(data)
no_slices = size_data(1)

out = zeros(no_slices, scale);

for i=1:no_slices;
    
    % Get data slice
    data_slice = data(i,:);
    cg = coarsegrain_notime(data_slice, tau);

    % Calculating moving MSE with no zscore
    out(i,:) = movavg_mse_no_zcscore(cg, r, scale);

    
    % optional update on slices completed 
     if rem(i, 10)== 0 
        disp(strcat('Completed:', num2str(i),'/', num2str(no_slices)))
    end
end

% Output MMSE 
writematrix(out,strcat(output_name, '_ch', num2str(ch), '.txt' ),'Delimiter',',')
disp('Completed:')
disp(output_name)
disp('_________________________________________')

end
disp('**************************************************')
disp('****************  FINISHED RUN   *****************')
disp('**************************************************')



end 
