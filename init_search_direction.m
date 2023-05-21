function [var_search_dir] = init_search_direction(nComps, nCells, nWells, well)

    % Define searching direction of each primary variable
    dir_Ni  =  1.0 .* ones(nCells, nComps);
    dir_Nw  =  1.0 .* ones(nCells, 1);
    dir_P   = -1.0 .* ones(nCells, 1);
    dir_Pwf = -1.0 .* ones(nWells, 1);       % Producer for timestep 1
    
    % Update for injector
    dir_Pwf(well.Prod_Flag ~= 1) = 1.0;
    
    dir_Ni_dummy = reshape(dir_Ni', nComps, nCells);
    dir_Nw_dummy = reshape(dir_Nw, 1, nCells);   
    dir_P_dummy = reshape(dir_P, 1, nCells);
    
    % Reservoir variable
    var_search_dir = [dir_Ni_dummy; dir_Nw_dummy; dir_P_dummy];
    var_search_dir = reshape(var_search_dir, nCells*(nComps+2), 1);  

    % Well variable at the end
    var_search_dir = [var_search_dir; dir_Pwf];
    
end
