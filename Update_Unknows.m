function [Ni, Nw, P, Pwf] = Update_Unknows(x_vec, delta_x_vec, nCells, nComps, nWells)
    
    % Number of reservoir variables
    nRVar = nCells * (nComps + 2);
    
    % Add increment to the original 
    x_vec_temp_0 = x_vec + delta_x_vec;
    
    % Eliminate round off error 
    delta_x_vec = x_vec_temp_0 - x_vec;
    x_vec_temp_0 = x_vec + delta_x_vec;
    
    % Get Pwf 
    Pwf = zeros(nWells, 1);
    Pwf(1 : nWells, 1) = x_vec_temp_0(nRVar+1 : end);    
    
    % Assign back to primary variables
    x_vec_temp_1 = x_vec_temp_0(1 : nRVar);
    x_vec_temp_1 = reshape(x_vec_temp_1, nComps+2, nCells);
    
    Ni_dummy = x_vec_temp_1(1:nComps, :);
    Ni = reshape(Ni_dummy.', nCells, nComps);
    Nw = reshape(x_vec_temp_1(nComps+1, :), nCells, 1);   
    P = reshape(x_vec_temp_1(nComps+2, :), nCells, 1);

    % Get changing trend of each primary variable
    % var_search_dir_1 = ones(size(delta_x_vec));
    % var_search_dir_1(delta_x_vec < 0) = -1;
    
end