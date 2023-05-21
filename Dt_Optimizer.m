function [dt_1, Reduce_dt_counter_1] = Dt_Optimizer(dt_0, fluid, init, well, prop, ...
    sim_ctr, current_sim_t, var_search_dir, Reduce_dt_counter)

    % Define temporary variables
    current = 1;
    next = 2;
    nWells  = well.nWells;
    [nCells, nComps] = size(prop(current).cell_Ni);
    nRVar = nCells * (nComps + 2);
      
    % Memory allocation
    dir_Pwf = zeros(nWells, 1);
    dir_P   = zeros(nCells, 1);
    flag_dt_back = zeros(nWells+5, 1);
    Reduce_dt_counter_1 = zeros(nWells+5, 1);
    
    flag_dt_back(:) = false;
    
    %% Decide timestep size multipler
    %     omega = 0.5;
    %     max_dP = 200;   % Maximum change of pressure per timestep
    %     max_dPwf = 200; % Maximum change of BHP per timestep
    %     
    %     dNi = prop(next).cell_Ni - Ni_hold;
    %     dNw = prop(next).cell_Nw - Nw_hold;
    %     dP  = prop(next).cell_P  - P_hold;
    %     dPwf = prop(next).well_P - Pwf_hold;
    %     
    %     dt_mul_1 = (1 + omega) .* max_dP ./ (dP + omega *  max_dP);  
    %     dt_mul_2 = (1 + omega) .* max_dPwf ./ (dPwf + omega *  max_dPwf);
    %     dt_mul = min([min(dt_mul_1(:)), min(dt_mul_2(:))]);
    %     dt_mul = min([dt_mul, sim_ctr.dt_mul]);

    %% Normal timestep size increase
    % If the current timestep converge good, predict next timestep size
    % Converged tstep, increase dt
    if ( dt_0 * sim_ctr.dt_mul ) < sim_ctr.dt_max
        dt_1 = dt_0 * sim_ctr.dt_mul;
    else
        dt_1 = sim_ctr.dt_max;
    end
           
    %% Timestep size decrease based on BHP drop and phase appearance prediction
    % Get variable searching direction 
    dir_Pwf(1 : nWells, 1) = var_search_dir(nRVar+1 : end);
    var_search_dir_1 = var_search_dir(1 : nRVar);
    var_search_dir_1 = reshape(var_search_dir_1, nComps+2, nCells);
    dir_P(:) = reshape(var_search_dir_1(nComps+2, :), nCells, 1);    
        
    %% Check well bottom hole pressure
    % Get variable searching direction
    if sim_ctr.dbhp_thresh > 0.0d0
        for iWell = 1 : well.nWells
            % distance
            distance = prop(current).well_P(iWell) - well.Target_min_BHP(iWell); 

            %(distance < sim_ctr.dbhp_thresh * relax_r) && ...
            %(distance > sim_ctr.dbhp_thresh * relax_l) && ...

            if  (well.Prod_Flag(iWell)) && ...     % Producer
                (distance > 0) && ...              % Distance check 
                (distance < sim_ctr.dbhp_thresh ) && ...
                (dir_Pwf(iWell) == -1) && ...      % Changing trend check
                (Reduce_dt_counter(iWell) < 1)     % Previous cut check

                % Producer BHP drop and close to minimum BHP
                flag_dt_back(iWell) = true;
                Reduce_dt_counter_1(iWell) = 1; 

            end
            if  (~well.Prod_Flag(iWell)) && ...   % Injector
                (-distance > 0) && ...            % Distance check
                (-distance < sim_ctr.dbhp_thresh) && ...
                (dir_Pwf(iWell) == +1) && ...     % Changing trend check
                (Reduce_dt_counter(iWell) < 1)    % Previous cut check

                % Injector BHP increase and close maximum BHP
                flag_dt_back(iWell) = true;
                Reduce_dt_counter_1(iWell) = 1;

            end
        end
    end

    %% Check if any cell pressure is close to bubble point 
    if sim_ctr.dpp_thresh > 0.0d0
        [P_min, iCell]  = min(prop(current).cell_P); 
        Pb = init.res_vle_chk(iCell, nComps+1);

        if  isequal(fluid.fluidtype, 'VLE') && ...  
            ( P_min > Pb ) && ...
            ( P_min - Pb < sim_ctr.dpp_thresh) && ...
            ( dir_P(iCell) == -1 ) && ...
            ( Reduce_dt_counter(nWells+1) < 1) 
        
            % One cell pressure is close to Pb and will decrease soon to Pb
            flag_dt_back(nWells+1) = true;
            Reduce_dt_counter_1(nWells+1) = 1;

        end

    end
    
    if any(flag_dt_back) == true
        dt_1 = dt_0 / (sim_ctr.dt_mul)^sim_ctr.dt_back;
    end
    
    %% Set final dt to match final_Sim_t
    if( current_sim_t + dt_1 ) > sim_ctr.t_total
        dt_1 = sim_ctr.t_total - current_sim_t;
    end 
       
end