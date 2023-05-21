
function [var_search_dir_1] = update_search_direction(nComps, nCells, prop)

    % Define temporary variable
    current = 1;
    next    = 2;
    
    % Get change of each primary variable
    dir_Ni  = prop(next).cell_Ni - prop(current).cell_Ni;
    dir_Nw  = prop(next).cell_Nw - prop(current).cell_Nw;
    dir_P   = prop(next).cell_P  - prop(current).cell_P;
    dir_Pwf = prop(next).well_P  - prop(current).well_P;
    
    dir_Ni(dir_Ni < 0) = -1; dir_Ni(dir_Ni >= 0) = 1;
    dir_Nw(dir_Nw < 0) = -1; dir_Nw(dir_Nw >= 0) = 1;
    dir_P(dir_P < 0)   = -1; dir_P(dir_P >= 0)   = 1;
    dir_Pwf(dir_Pwf < 0) = -1; dir_Pwf(dir_Pwf >= 0) = 1;
   
    dir_Ni_dummy = reshape(dir_Ni', nComps, nCells);
    dir_Nw_dummy = reshape(dir_Nw, 1, nCells);   
    dir_P_dummy = reshape(dir_P, 1, nCells);
    
    % Assemble variable
    % Reservoir
    var_search_dir_1 = [dir_Ni_dummy; dir_Nw_dummy; dir_P_dummy];
    var_search_dir_1 = reshape(var_search_dir_1, nCells*(nComps+2), 1);  

    % Well 
    var_search_dir_1 = [var_search_dir_1; dir_Pwf];    
    
end