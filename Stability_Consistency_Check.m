function [prop, flag_Chop_dt_k, well_ctr, delta_x_vec] = Stability_Consistency_Check(init, ...
          prop, fluid, well, sim_ctr, x_vec, delta_x_vec, Ni_hold, Nw_hold, P_hold, Pwf_hold, newton_iter, well_ctr)

    % Define temporary variables
    current = 1;
    next = 2;
    nComps = fluid.nComps;
    nCells = init.nCells;
    nWells = well.nWells;
    
    flag_Chop_dt_k = zeros(10, 1);
    
    flag_Chop_dt_k(:) = false;
       
    %% Convergence check 
    % Increment check
    nRVar = (nComps+2)*nCells;
    delta_x_vec_tmp = reshape(delta_x_vec(1 : nRVar), nComps+2, nCells);
    delta_P_vec_tmp = reshape(delta_x_vec_tmp(nComps+2, :), nCells, 1); 
    delta_Ni_vec_tmp = reshape(delta_x_vec_tmp(1 : nComps, :), nCells, nComps); 
    delta_Nw_vec_tmp = reshape(delta_x_vec_tmp(nComps+1, :), nCells, 1);
    delta_Pwf_vec_tmp = delta_x_vec(nRVar+1 : end, 1); 
       
    dP_max = max(abs(delta_P_vec_tmp(:)));  
    
    % Check Pwf
    % Initialize newton damping factor
    damping_coef = 1;
    for iWell = 1 : nWells
        % Producer damping checking
        if prop(next).well_P(iWell) < well.Target_min_BHP(iWell) && well.Prod_Flag(iWell) == 1
            damping_coef_1 = (well.Target_min_BHP(iWell) - Pwf_hold(iWell)) / delta_Pwf_vec_tmp(iWell);
            damping_coef = min([damping_coef, damping_coef_1]);
            well_ctr(current).BHP_Ctr(iWell) = 1;
        end        
        % Injector damping checking
        if prop(next).well_P(iWell) > well.Target_max_BHP(iWell) && ( well.Inje_Wat_Flag(iWell) == 1 || well.Inje_Gas_Flag(iWell) == 1)
            damping_coef_1 = (well.Target_max_BHP(iWell) - Pwf_hold(iWell)) / delta_Pwf_vec_tmp(iWell);
            damping_coef = min([damping_coef, damping_coef_1]);
            well_ctr(current).BHP_Ctr(iWell) = 1;
        end        
    end
       
    % If damping is required
    if damping_coef ~= 1
        delta_x_vec = delta_x_vec .* damping_coef; 
        [prop(next).cell_Ni, prop(next).cell_Nw, prop(next).cell_P, prop(next).well_P] = Update_Unknows(x_vec, delta_x_vec, nCells, nComps, nWells);
    end
                 
    % Max and min ranges for key properties
    % Ensure saturations and compositions are between 0 and 1, and pressure>14.7
    p_min_constraint  = min(prop(next).cell_P); 
    p_NaN_Inf_constraint = sum( isnan(prop(next).cell_P) + isinf(prop(next).cell_P) ) ...
                         + sum( isnan(prop(next).well_P) + isinf(prop(next).well_P) );
      
    % Sometimes negative mass, 
    prop(next).cell_Ni(prop(next).cell_Ni < 0) = 1.0d-7;
    prop(next).cell_Nw(prop(next).cell_Nw < 0) = 1.0d-7;
                         
    Ni_min_constraint    = min( prop(next).cell_Ni(:) );
    Nw_min_constraint    = min( prop(next).cell_Nw(:) );
    
    %% If constraints, material balance, physical ranges, or max iterations are violated, chop dt
    Constraints_ok_flag = (Ni_min_constraint < 0) ...
                        + (Nw_min_constraint < 0) ...
                        + (p_min_constraint < 14.7d0) ...  
                        + (p_NaN_Inf_constraint > 0) ...
                        + (newton_iter > sim_ctr.iter_max); 

    if Constraints_ok_flag ~= 0
        flag_Chop_dt_k(1) = true;
        if (Ni_min_constraint < 0)
            disp('Minimum Ni is less than zero')
        end
        if (Nw_min_constraint < 0)
            disp('Minimum Nw is less than zero')
        end
        if (p_min_constraint < 0)
            disp('Minimum P is less than 14.7')
        end        
        if (p_NaN_Inf_constraint > 0)
            disp('Nonphysical P exist')
        end        
        if (newton_iter > sim_ctr.iter_max)
            disp( strcat('Newton iteration is ', num2str(newton_iter), ', going beyond max iteration' ) )
        end        
    end
    
    % If any of well shut off all perforations, cut timestep
    for iWell = 1 : nWells
        if sum( well_ctr(current).wconne_well( well_ctr(current).wconne_well == iWell ) ) == 0
            flag_Chop_dt_k(2) = true;
            disp( strcat( 'All perforations in WELL ', well.name(iWell), ' is closed due to crossflow!') )
        end
    end

end