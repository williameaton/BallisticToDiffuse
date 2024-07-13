clear all
close all
clc
format long

%stns = ["b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y"];
%chls = [ "R" , "T", "Z"];   % Define channels 

stns = ["b"];
chls = [ "R" , "T", "Z"];  

r = 0.2;                    % Threshold  
r_str = num2str(r);         
delta_v = 'p-0.2';
no_stns = length(stns);

m_min = 1 ;
m_max = 200;

num_m = m_max - m_min + 1; % Number of scales to calculates 

rawdata_dir = strcat('../../data/','/MSE_Pslices_250/',delta_v, '/');

%m = ["1.5", "8"; "2", "4"; "2", "5"; "2", "6"; "2", "7"; "2", "8"; "2.5", "6"; "2.5", "7"; "2.5", "8"];
sim = ["0.5", "2"];

loop_len = size(sim);

for mfp_i = 1:loop_len(1);
    
    rad      = sim(mfp_i,1);
    mfp      = sim(mfp_i,2);
    sim_name = strcat(delta_v, '_2hz_', mfp, '_mfp_', rad,'_rad');

    
    
    % Create output file file:
    out_folder_path = strcat('./entropy/', delta_v ,'/', mfp,'_', rad);
    mkdir(out_folder_path);
    
    for stn = 1:1;
        for ch = 1:3;
      
            channel = chls(ch);
            input_path = strcat(rawdata_dir, '/', sim_name, '/ST', num2str(stns(stn)), '.', channel);
            output_name = strcat(out_folder_path, '/', 'ST', num2str(stns(stn)), '.', channel, '_r', num2str(r));
    
            
            disp('Loading data from')
            disp(input_path);
            data = importdata(input_path); 
            time = data(:,1);
            d = data(:,2);
            disp('Load completed.') 
   
            % Initialise output entropy array
            e = zeros(num_m, 1);
            ctr = 1;

            
            norm_data = d/(max(abs(d)));
            for m=m_min:m_max;
                disp(strcat('m: ', num2str(m), '/', num2str(m_max)))
                

                %mm = movingmean(norm_data, m);
                %e(ctr) = SampEn(mm,r,m);
                
                cg_data = coarsegrain_notime(d, m);
                e(ctr) = sampen(cg_data, m, r);
                disp(e(ctr))

                %cg_t = coarsegrain_notime(time, m);
                %cg_data = coarsegrain_notime(norm_data, m)
                % Calculate sample entropy 
                %e(ctr) = SampEn(cg_data, r, m);
                
                % Update index counter 
                ctr = ctr + 1;
            end

        output_matrix = [m_min:m_max'; e']';

        writematrix(output_matrix, output_name,'Delimiter',' ', 'FileType','text');

        end 
        disp('Completed:')
        disp(output_name)
        disp('_______________________________________')
        
    end
end 

disp('************************************************')
disp('**************  FINISHED RUN   *****************')
disp('************************************************')
