
%function [] = mars_MSE(r, tau, scale, input_path, output_folder, output_file)

    % ___________________________________________________________________________
    % Description: 
    % Function calculating the multiscale entropy for Mars data 
    % ___________________________________________________________________________
    % Input parameters: 
    % r [float]           MSE Threshold value e.g. 0.1
    % tau [int]           Coarse graining decimation value e.g. 30 
    % scale [int]         Max. scale for MSE calculations e.g. 50
    % input_path [str]    Input file path and name e.g. 'example_folder/file.ascii'
    % output_folder [str] Output folder path  e.g. 'example_folder'
    % output_file [str]   Output file name  e.g. 'example_file.ascii'
    % ___________________________________________________________________________
    % Other functions required for this function: 
    %   coarsegrain_notime
    %   movavg_mse
    % ___________________________________________________________________________

    r = 0.1
    tau = 30
    scale = 50
    input_path = "./example_mars_data" 
    output_folder = "./TEST/"
    output_file = "test.txt"
	

	% If variables are accidently inputted as strings: 
	if isa(tau, 'char')
	tau = str2num(tau)
	end 
	
	if isa(scale, 'char')
	scale = str2num(scale)
	end 
    
  
  % ______________________________________________________________
  
    % Create output file file:
    mkdir(output_folder)

    % Define default output path and name
    output_name = strcat(output_folder, '/', output_file);


    % Append filetype: 
    % This needs more flexibility in input names
    input_path = strcat(input_path, '.ascii');
    
    % Verbose updates
    disp('Loading data from')
    disp(input_path)

    data = textread(input_path); 
    disp('Data load completed.') 


	% Check number of channels and whether row-major or column-major data: 
    % Get data dimensions: 
    d_dim = size(data);
    
    
    no_channels = min(d_dim); % Assuming that there are more data points in time series than channels
    % Transpose data if is [time series, channels] instead of [channels, time series]
    if no_channels == d_dim(2);
        data = transpose(data);
    end 
    
    % Initialise entropy matrix (MSE output matrix)
    entropy = zeros(no_channels, scale);
    
    for ch = 1:no_channels;
            temp_data = data(ch, :);

	    % NORMALISE DATA 
	    norm_data = temp_data/(max(abs(temp_data)));

	    % Coarse grain the time series - this is before any entropy calculation simply to decimate the timeseries for computational efficiency: 
	    disp('Apply Coarse Grain')

	    cg = coarsegrain_notime(norm_data, tau);
	    disp('Coarse Grain completed')

	    % Calculate entropy over different scales using moving-average         
	    entropy(ch, :) = movavg_mse(cg, r, scale, 'yes');

	    disp(strcat('Completed channel ',' ', num2str(ch), '/', num2str(no_channels)))

    end
    
    
    % Output calculated entropy data: 
    writematrix(entropy, output_name,'Delimiter',',')
       
    disp('Written to:')
    disp(output_name)
    disp('_________________________________________') 
    disp('**************************************************')
    disp('****************  FINISHED RUN   *****************')
    disp('**************************************************')
