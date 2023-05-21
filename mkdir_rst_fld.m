function [rst_fd, rsm_fd] = mkdir_rst_fld()

    % Creat a folder to save mat file for post-processing

    rst_fd = 'rst_mat';
    ex_flag = exist(rst_fd);
    if ex_flag ~= 7
        mkdir(rst_fd);
    end
    
    rsm_fd = 'rsm';
    ex_flag = exist(rsm_fd);
    if ex_flag ~= 7
        mkdir(rsm_fd);
    end    
 
end