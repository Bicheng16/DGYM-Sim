function Save_Result(prop, num_time_steps, current_sim_t, rst_fd)

    % Save variables into a folder
    current = 1;
    
    % Check if the folder exist
    ex_flag = exist(rst_fd);    
    if ex_flag ~= 7
        return
    end
    file_path= strcat(pwd, '\', rst_fd);
    
    % Write the time steps and report data in a text file
    t_dat = 'time.dat';
    if num_time_steps == 0
        % Delete all files of *.mat and *.PLT
        delete(fullfile(file_path, '\*.PLT'));
        delete(fullfile(file_path, '\*.mat'));
        
        t_fid = fopen(fullfile(file_path, t_dat), 'w');
        % Copy tecplot file
        plt_file_path = strcat(pwd, '\input\include\');
        copyfile(fullfile(plt_file_path, '\MESH_1.PLT'), fullfile(file_path, 'MESH_1.PLT'));        
    else
        t_fid = fopen(fullfile(file_path, t_dat), 'a+'); 
    end
    
    fprintf(t_fid, '%13.7e\n', current_sim_t);    
    fclose(t_fid);
    
    % Save P, Sw, Sg as example
    P = prop(current).cell_P;
    Sw = prop(current).cell_Sw;
    Sg = prop(current).cell_Sg;
    
    % Save variable in terms of 'timestep_variableName.mat'
    current_P = strcat(num2str(num_time_steps), '_P.mat');
    current_Sw = strcat(num2str(num_time_steps), '_Sw.mat');
    current_Sg = strcat(num2str(num_time_steps), '_Sg.mat');
    save(fullfile(file_path, current_P), 'P');    
    save(fullfile(file_path, current_Sw), 'Sw');    
    save(fullfile(file_path, current_Sg), 'Sg'); 
        
end