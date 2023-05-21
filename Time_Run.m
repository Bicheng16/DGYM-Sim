function [flag_run_finish_ok, Chop_time_counter] = Time_Run(init, fluid, res, rockfluid, sim_ctr, well, well_ctr, prop, rst, CPU_time_start, input)

    disp('...Start of Time Stepping')
    disp('|===========================================================')

    format longE;    
    Jacobian_cond = [];
    Newtons = [];
    
    % Creat a folder to save files (mat and result) for post-processing
    [rst_fd, rsm_fd] = mkdir_rst_fld();
    % Creat a report file to store summary
    output = createReport(input);
    % Save data at the initial status: timestep = 0;
    Save_Result(prop, 0, 0, rst_fd);
    % Write result report 
    writeReport(output, rsm_fd, rst, fluid, init, 0);
    
    current_sim_t = 0;
    current = 1;
    next = 2;
    prev = 3;
    nCells = init.nCells;
    nConnes = init.nConnes;
    nComps = fluid.nComps;
    nWells = well.nWells;
    nPerfs = sum(well.nperf_cell);
    vol_cof = 5.61458144;                   % Conversion from bbl to cubic feet, original: 5.614581442311112 
    
    % Initialize time step size
    dt_k = sim_ctr.dt_init;
    % Estimated maximum number of Jacobian nonzero values
    Jacobian_nnz_est = (2 * nConnes + nCells) * (nComps + 2) ^ 2 ...
                     +  2 * nPerfs * (nComps + 2) + nWells;
    % Initialize time, newton iteration records
    Chop_time_counter = 0;
    num_time_steps = 0;
    num_newton_steps = 0;

    % Dt reducer counter
    Reduce_dt_counter = zeros(nWells+5, 1);
    
    % Nodes and connection range to calculate residual or RHS
    Res_Nodes  = 1 : nCells;
    Res_Connes = 1 : nConnes; 
        
    % Initialize search direction of each primary variable
    [var_search_dir_0] = init_search_direction(nComps, nCells, nWells, well);
    
    % For Slippage flow Multiplier
    GKapp_data = [];
        
    %% Time step loop
    while current_sim_t < sim_ctr.t_total
                
        % Initialize newton-step loop
        newton_iter = 0;
        error = 1000;
        error_0 = error; % To check the monotonousity of the convergence.
        
        % Initialize time step chop flag 
        flag_Chop_dt = zeros(10, 1);  

        % Initialize newton-step variables
        prop(next).cell_Ni = prop(current).cell_Ni;
        prop(next).cell_Nw = prop(current).cell_Nw;        
        prop(next).cell_P = prop(current).cell_P;
        prop(next).cell_T = prop(current).cell_T; 
        prop(next).well_P = prop(current).well_P;
        
        % Initialize well control flags
        % Perforation open/shutin: continous within different newtons or in a single timestep 
        % Rate/BHP control: continuou between different newtons and different timesteps 

        if num_time_steps == 65 
           a = 0; 
        end        
        
        if num_time_steps > 0
            % well_ctr(current).wconne_well = well.wconne_well;
            % Update check point for VLE 
            % [init] = Check_PSat(init, well, prop, fluid, res);                       
            % Calculate wellbore pressure loss due to gravity 
            [prop] = Cal_Well_Gradient(well, prop, fluid, well_ctr, init);
        end
        
%         disp('Kapp')
%         prop(current).cell_GKapp
%         temp = zeros(nCells, nComps+1);
%         temp(:, 1) = prop(current).cell_P;
%         temp(:, 2:end) = prop(current).cell_GKapp(:, :);
%         GKapp_data = [GKapp_data; temp];

%         disp('Pressure');
%         prop(current).cell_P;
                       
        well_ctr(prev) = well_ctr(current);   % Backup a copy
                                         
        %% Newton loop
        while (error >= sim_ctr.tol) && (newton_iter <= sim_ctr.iter_max) % Recommend error: 0.1
            
            % Iteration counter
            newton_iter = newton_iter + 1;
            
            if newton_iter == 1 && num_time_steps == 7 
               a = 1; 
            end
                                       
            % Initialize dummy struct to import or export variables 
            dummy = alloc_inherit_var;

            if newton_iter == 1
                % Initial guess from previous timestep
                prop(next) = prop(current);
                dummy = update_inherit_var(true, dummy, next, prop);
                % variable searching direction for timestep cut
                var_search_dir_current = var_search_dir_0; 
            else
                % Update cell properties based on primary variables from
                % previous newton
                [dummy, prop] = Cell_Update_RHS(Res_Nodes, dummy, init, fluid, res, rockfluid, prop);
            end

            % Assemble unknows
            [x_vec, nVar] = Assemble_Unknows(prop(next).cell_Ni, prop(next).cell_Nw, prop(next).cell_P, prop(next).well_P);

            % Backup primary variables
            Ni_hold = prop(next).cell_Ni;
            Nw_hold = prop(next).cell_Nw;
            P_hold = prop(next).cell_P;
            T_hold = prop(next).cell_T;    % Flag for non-isothermal case in the future
            Pwf_hold = prop(next).well_P;
                                                                                                     
            % Calculate Right Hand Side (RHS) for each primary functions
            % for reservoir part            
            [Node_Residual_0, Conne_Residual_0, R_well_0, ~, ~, well_ctr] = ...
            Residual_Pieces(Res_Nodes, Res_Connes, dt_k, init, fluid, res, well, prop, well_ctr);
        
            % Make well control consistent
            well_ctr(current) = well_ctr(next);  
                                                                               
            % Assemble Residual: RHS = [R_c; R_wat; R_vol]' * nCells; 
            [RHS_0, ~, ~, ~, ~, ~, ~] = Residual_Assemble( Conne_Residual_0, Node_Residual_0, R_well_0,    ...
             nCells, nComps, init.conne_cell1, init.conne_cell2);

            % Error of mass balance
            error = norm(RHS_0, 2);
                              
            % Clear newton-level secondary variables
            [prop] = zero_dynamic_properties(prop, nCells, nConnes, nComps, nPerfs);         
         
            % Evaluate and Assemble Jacobian
            [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Assemble(Jacobian_nnz_est, RHS_0, Node_Residual_0, Conne_Residual_0, R_well_0, ...
            well_ctr, var_search_dir_0, dt_k, dummy, sim_ctr, init, fluid, res, rockfluid, well, prop, Ni_hold, Nw_hold, P_hold, T_hold, Pwf_hold);

            % Evaluate and Assemble Jacobian
%             if (newton_iter > 1 && error_0 > error && error < 10)
%                 Triplet_I = Triplet_I0 ;
%                 Triplet_J = Triplet_J0; 
%                 Triplet_Val = Triplet_Val0;
%                 nonzeros_JVal = nonzeros_JVal0;                
%             else        
%                 [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Assemble(Jacobian_nnz_est, RHS_0, Node_Residual_0, Conne_Residual_0, R_well_0, ...
%                 well_ctr, var_search_dir_0, dt_k, dummy, sim_ctr, init, fluid, res, rockfluid, well, prop, Ni_hold, Nw_hold, P_hold, T_hold, Pwf_hold);
%                 Triplet_I0 = Triplet_I;
%                 Triplet_J0 = Triplet_J; 
%                 Triplet_Val0 = Triplet_Val;
%                 nonzeros_JVal0 = nonzeros_JVal;        
%             end
            
            % Solve the linear system: 4 vector feed into the solver
            [delta_x_vec, flag_Chop_dt(1), Jacobian_cond] = LinearSolver(Triplet_I, Triplet_J, Triplet_Val, ...
            nonzeros_JVal, nVar, -RHS_0, sim_ctr.solver, Jacobian_cond);
        
            % Update primary variables in the current newton iteration
            [prop(next).cell_Ni, prop(next).cell_Nw, prop(next).cell_P, prop(next).well_P] = Update_Unknows(x_vec, delta_x_vec, nCells, nComps, nWells);
            
            % Clear newton-level secondary variables after solve linear system AX = b
            [prop] = zero_dynamic_properties(prop, nCells, nConnes, nComps, nPerfs);
                       
            % Stability and consistency check 
            [prop, flag_Chop_dt_k, well_ctr, delta_x_vec] = Stability_Consistency_Check(init, ...
            prop, fluid, well, sim_ctr, x_vec, delta_x_vec, Ni_hold, Nw_hold, P_hold, Pwf_hold, newton_iter, well_ctr);        
                                                         
            % Error convergence
            if newton_iter == 1
                error_0 = error;
            else
                % Only quadratic convergence is treated as correct searching direction
                flag_correct_search = ( error_0 > error ); %(error_0/error > 10);  
                if flag_correct_search
                    var_search_dir_0(:) = 1;
                    var_search_dir_0(delta_x_vec < 0) = -1;                    
                end
                error_0 = error;
            end            
            
            % If necessary to chop time steps, exit from newton loop
            if any(flag_Chop_dt_k) == 1
                break;
            end
                                    
        end
        
        %%
        %% Result report
        %%
                
        if any(flag_Chop_dt_k) == 1
            
            % Solution is not convergent, cut timestep and go back again
            prop(next).cell_Ni = prop(current).cell_Ni; 
            prop(next).cell_Nw = prop(current).cell_Nw; 
            prop(next).cell_P = prop(current).cell_P;
            prop(next).cell_T = prop(current).cell_T; 
            
            % Get well control at previous timestep
            well_ctr(current) = well_ctr(prev); 
            
            % Get variable searching direction at previous timestep
            var_search_dir_0 = var_search_dir_current;            
                       
            if ( dt_k / 2 ) < sim_ctr.dt_min
                disp('...ERROR: delta_t has been reduced to the minimum allowable specified value');
                disp('...ERROR: simulation will be terminated');
                flag_run_finish_ok = false;
                return;
            else
                disp(  strcat('......INFO: delta_t will be reduced from: ', num2str(dt_k), ' days, to: ', num2str(dt_k/2), ' days')    );
                dt_k = dt_k / 2;
                Chop_time_counter = Chop_time_counter + 1;
            end
                        
        else
                        
            % Successful run: record number of time steps and number of newton iterations 
            num_time_steps = num_time_steps + 1;
            num_newton_steps = num_newton_steps + newton_iter;
            Newtons = [Newtons newton_iter]; 
              
            current_sim_t = current_sim_t + dt_k;
            CPU_time = etime( clock, CPU_time_start );
            
            % Initialize dummy struct to import or export variables 
            dummy = alloc_inherit_var;
            
            % Update cell properties first
            [~, prop] = Cell_Update_RHS(Res_Nodes, dummy, init, fluid, res, rockfluid, prop);
            
            % Update the changing trend of each primary
            % varialbe between the current and previous timestep, for more
            % accurate secant slope in newton
            [var_search_dir_1] = update_search_direction(nComps, nCells, prop);
            var_search_dir_0 = var_search_dir_1;
            
            %[vol_total_error, Balance_vol] = Volumetric_Check(prop, init, 2);            
            
            % Calculate Right Hand Side (RHS) for each primary functions
            [Node_Residual, Conne_Residual, R_well, W_SS_Residual, perf_P, well_ctr] = ...
             Residual_Pieces(Res_Nodes, Res_Connes, dt_k, init, fluid, res, well, prop, well_ctr);
         
            % Make well control consistent
            well_ctr(current) = well_ctr(next);         
                                         
            % Calculate material balance error
            [~, ~, ~, ~, Balance_hc, Balance_wat, Balance_vol] = Residual_Assemble( Conne_Residual, ...
            Node_Residual, R_well, nCells, nComps, init.conne_cell1, init.conne_cell2);            
                                
            % Assign prop(next) back to prop(current)
            prop(next).perf_P = perf_P;
            prop(current) = prop(next);
                                    
            % Update result variables            
            % Book keeping cumulative fluid production
            FOPT_prev = rst.FOPT; FGPT_prev = rst.FGPT; FWPT_prev = rst.FWPT; 
            FOIT_prev = rst.FOIT; FGIT_prev = rst.FGIT; FWIT_prev = rst.FWIT;
            
            rst = alloc_rst_var(nComps, nWells);
            
            % (1): Mainly fluid in place
            [rst] = Result_Var_Cals(next, rst, fluid, prop, init);

            % Other rst variables
            rst.dt = dt_k;
            rst.newtons = newton_iter;

            % (2): Mainly production parameters
            for iWell = 1 : nWells

                % Get bottom hole pressure 
                rst.well_P(iWell, 1) = prop(current).well_P(iWell, 1); 

                % W_SS_i, W_SS_wat are from Reservoir_Residual_Piece   

                % W_SS_Residual(:, 1:nComps) represents lb-mol of hydrocarbon produced from each perforated cell
                n_HC_moles_well = sum(sum( W_SS_Residual(well.wconne_well == iWell, 1:nComps), 2 )); 
                % W_SS_Residual(:, 1+nComps) represent water rate: lb-mol/day
                n_wat_moles_well = sum(W_SS_Residual(well.wconne_well == iWell, 1+nComps));
                                
                % Producer
                if well.Prod_Flag(iWell) == 1
                    
                    % Calculate well stream composition for the current producer
                    ni = abs( sum( W_SS_Residual(well.wconne_well == iWell, 1:nComps), 1 ) );
                    Zi_wellstream = ni ./ sum(ni);
                    Zi_wellstream(isnan(Zi_wellstream) | isinf(Zi_wellstream)) = 0; 
                    
                    % Update well stream composition again
                    well.Gas_Comp(iWell, :) = Zi_wellstream;
                                        
                    % Flash well stream composition to surface standard condition
                    % Surface flash calculation
                    Pb = init.wel_vle_chk(iWell, nComps+1);
                    kval_0 = init.wel_vle_chk(iWell, 1:nComps);
                    fg_0 = 0.5d0;
                    kflag = 0;  
                    if abs( sum(Zi_wellstream) - 1.0d0 ) < 1.0e-10
                        [fv_ws_SC, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, MolarVolliq_ws_SC, MolarVolvap_ws_SC, ~, ~, ~, ~] ...
                        = Y_FluidFlash(fluid.Psc, Zi_wellstream, fluid.Tsc, fluid, kval_0, fg_0, kflag, Pb);
                    else
                        fv_ws_SC = 0.0d0;
                        MolarVolliq_ws_SC = 0.0d0;
                        MolarVolvap_ws_SC = 0.0d0;
                    end
                                       
                    % Single well production volume rate
                    rst.WOPR(iWell) = - (1 - fv_ws_SC) .* n_HC_moles_well .* MolarVolliq_ws_SC ./ vol_cof;  % Oil rate, STB/day
                    rst.WGPR(iWell) = - fv_ws_SC  .* n_HC_moles_well .* MolarVolvap_ws_SC ./ 1000;          % Gas rate, MSCF/day
                    rst.WWPR(iWell) = - n_wat_moles_well .* fluid.wat_MW ./ fluid.MDen_w_sc ./ vol_cof;     % Water rate, STB/day
                    
                    % Field production volume rate
                    rst.FOPR = rst.FOPR + rst.WOPR(iWell); % Oil rate, STB/day
                    rst.FGPR = rst.FGPR + rst.WGPR(iWell); % Gas rate, MSCF/day
                    rst.FWPR = rst.FWPR + rst.WWPR(iWell); % Water rate, STB/day
                    
                end
                                    
                % Water Injector
                if well.Inje_Wat_Flag(iWell) == 1
                    
                    % Single well injection volume rate
                    rst.WWIR(iWell) = n_wat_moles_well .* fluid.wat_MW ./ fluid.MDen_w_sc ./ vol_cof;                         % Water rate, STB/day
                    
                    % Field injection volume rate
                    rst.FWIR = rst.FWIR + rst.WWIR(iWell); % Water rate, STB/day   
                    
                end
                    
                % Gas Injector
                if well.Inje_Gas_Flag(iWell) == 1
                    
                    % Update Zi_wellstream from initial well schedual
                    Zi_wellstream = well.Gas_Comp(iWell, :);
                    
                    % Surface flash calculation
                    Pb = init.wel_vle_chk(iWell, nComps+1);
                    kval_0 = init.wel_vle_chk(iWell, 1:nComps);
                    fg_0 = 0.5d0;
                    kflag = 0; 
                    
                    [fv_ws_SC, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, MolarVolliq_ws_SC, MolarVolvap_ws_SC, ~, ~, ~, ~] ...
                    = Y_FluidFlash(fluid.Psc, Zi_wellstream, fluid.Tsc, fluid, kval_0, fg_0, kflag, Pb);
                                
                    % Single well injection volume rate
                    rst.WOIR(iWell) = - ( 1 - fv_ws_SC ) .* n_HC_moles_well .* MolarVolliq_ws_SC ./ vol_cof; % Oil rate, STB/day
                    rst.WGIR(iWell) = - fv_ws_SC .* n_HC_moles_well .* MolarVolvap_ws_SC ./ 1000;            % Gas rate, MSCF/day
                    rst.WWIR(iWell) = - n_wat_moles_well .* fluid.wat_MW ./ fluid.MDen_w_sc ./ vol_cof;                        % Water rate, STB/day
                    
                    % Field injection volume rate
                    rst.FOIR = rst.FOIR + rst.WOIR(iWell); % Oil rate, STB/day
                    rst.FGIR = rst.FGIR + rst.WGIR(iWell); % Gas rate, MSCF/day
                    rst.FWIR = rst.FWIR + rst.WWIR(iWell); % Water rate, STB/day
                    
                end
                                
            end  % End well loop
            
            % Cumulative production
            rst.FOPT = FOPT_prev + rst.FOPR * dt_k; % Cumulative oil production, STB
            rst.FGPT = FGPT_prev + rst.FGPR * dt_k; % Cumulative gas production, MSCF
            rst.FWPT = FWPT_prev + rst.FWPR * dt_k; % Cumulative water production, STB
            
            % Cumulative injection 
            rst.FOIT = FOIT_prev + rst.FOIR * dt_k; % Cumulative oil production, STB
            rst.FGIT = FGIT_prev + rst.FGIR * dt_k; % Cumulative gas production, MSCF
            rst.FWIT = FWIT_prev + rst.FWIR * dt_k; % Cumulative water production, STB
            
            % Display results on console
            disp(  strcat('Convergence Reached at:  ', num2str(current_sim_t), ' days of simulation', '   t_step:', num2str(num_time_steps)    ))
            disp(  strcat('CPU time: ', num2str(CPU_time), ' sec')    )
            disp(  strcat('Simulation time: ', num2str(current_sim_t), ' days')    )
            disp(  strcat('dt: ', num2str(dt_k), ' days')    )
            disp(  strcat('Newton iterations: ', num2str(newton_iter))    )
            % Display material balance error
            % Hydrocarbon in the reservoir
            disp(' ')
            disp('Material Balance Error:')
            for i = 1 : nComps
                msg = strcat(fluid.comp_name{i}, ' | L-2 norm :   ', num2str( Balance_hc(i, 1)), ...
                                                 ' , L-inf norm :   ', num2str( Balance_hc(i, 2)), ', unit: lb-mol'); 
                disp(msg);
            end
            % Water in the reservoir
            msg = strcat('H2O | L-2 norm :   ', num2str( Balance_wat(1, 1)), ...
                            ' , L-inf norm :   ', num2str( Balance_wat(1, 2)), ', unit: lb-mol');
            disp(msg);
            % Volume in the reservoir
            msg = strcat('Volume | L-2 norm :   ', num2str( Balance_vol(1, 1)), ...
                            ' , L-inf norm :   ', num2str( Balance_vol(1, 2)), ', unit: cu.ft');
            disp(msg);

            disp(' ')
            disp(  strcat('P_res: ', num2str(rst.P_res_HCPVweigthed), ' psia')    )
            disp(  strcat('So_avg: ', num2str(rst.So_PVweigthed), ' fraction V/V')    )
            disp(  strcat('Sg_avg: ', num2str(rst.Sg_PVweigthed), ' fraction V/V')    )
            disp(  strcat('Sw_avg: ', num2str(rst.Sw_PVweigthed), ' fraction V/V')    )

            disp(  strcat('FOPR: ', num2str( rst.FOPR), ' STB/day')     )
            disp(  strcat('FWPR: ', num2str( rst.FWPR), ' STB/day')     )
            disp(  strcat('FGPR: ', num2str( rst.FGPR), ' MSCF/day')    )
            disp(  strcat('FOIR: ', num2str( rst.FOIR), ' STB/day')     )
            disp(  strcat('FWIR: ', num2str( rst.FWIR), ' STB/day')     )
            disp(  strcat('FGIR: ', num2str( rst.FGIR), ' MSCF/day')    )

            disp(  strcat('FOPT: ', num2str( rst.FOPT), ' STB')     )
            disp(  strcat('FWPT: ', num2str( rst.FWPT), ' STB')     )
            disp(  strcat('FGPT: ', num2str( rst.FGPT), ' MSCF')    )
            disp(  strcat('FOIT: ', num2str( rst.FOIT), ' STB')     )
            disp(  strcat('FWIT: ', num2str( rst.FWIT), ' STB')     )
            disp(  strcat('FGIT: ', num2str( rst.FGIT), ' MSCF')    )

            % Display bottom hole pressure
            for iWell = 1 : nWells
                if well.Prod_Flag(iWell) == 1
                    disp(strcat(well.name{iWell} , ' PROD BHP: ', num2str(prop(next).well_P(iWell)), ' psia'));
                else
                    disp(strcat(well.name{iWell} , ' INJ  BHP: ', num2str(prop(next).well_P(iWell)), ' psia'));
                end                    
            end

            disp('|===========================================================')
            
            % Save dynamic data into a result folder 
            Save_Result(prop, num_time_steps, current_sim_t, rst_fd);  
            % Write output report 
            writeReport(output, rsm_fd, rst, fluid, init, current_sim_t);

            % Time step optimization
            [dt_1, Reduce_dt_counter_1] = Dt_Optimizer(dt_k, fluid, init, well, prop, ...
            sim_ctr, current_sim_t, var_search_dir_0, Reduce_dt_counter);
        
            dt_k = dt_1; 
            Reduce_dt_counter(:) = Reduce_dt_counter(:) + Reduce_dt_counter_1(:);
                                                 
            % Terminate simulation if both oil and gas production is not
            % acceptable, abundan the wells and stop production.
            if  ( rst.FOPR < 1 && rst.FGPR < 1.0E-2  && rst.FWPR < 1 && rst.FOIR < 1 && rst.FGIR < 1.0E-2 && rst.FWIR < 1 ) ...
            || rst.P_res_PVweigthed < 1.001 * min(well.Target_min_BHP)
            % if false
            
                sim_ctr.t_total = current_sim_t;
                
                disp( strcat('...@ Time step: ', num2str(current_sim_t), 'days') );
                
                % Production
                if rst.FOPR < 1
                    disp(  strcat('...OPR below threshold, simulation is ended. Current OPR: ', num2str( rst.FOPR ), ' STB/day')    )
                end
                
                if rst.FGPR < 1.0E-2
                    disp(  strcat('...GPR below threshold, simulation is ended. Current GPR: ', num2str( rst.FGPR ), ' MSCF/day')    )
                end
                
                if rst.FWPR < 1
                    disp(  strcat('...WPR below threshold, simulation is ended. Current WPR: ', num2str( rst.FWPR ), ' STB/day')    )
                end
                
                % Injection
                if abs(rst.FOIR) < 1
                    disp(  strcat('...OIR below threshold, simulation is ended. Current OIR: ', num2str( rst.FOIR ), ' STB/day')    )
                end
                
                if abs(rst.FGIR) < 1.0E-2
                    disp(  strcat('...GIR below threshold, simulation is ended. Current GIR: ', num2str( rst.FGIR ), ' MSCF/day')    )
                end
                
                if abs(rst.FWIR) < 1
                    disp(  strcat('...WIR below threshold, simulation is ended. Current WIR: ', num2str( rst.FWIR ), ' STB/day')    )
                end
                       
                if rst.P_res_PVweigthed < 1.01 * min(well.Target_min_BHP) && ... 
                any(well.Prod_Flag == 1)
                    disp(  strcat('...Average reservoir pressure close to minimum BHP, simulation is ended. Current average reservoir pressure: ', num2str( rst.P_res_PVweigthed ), ' psia')    )
                end
                
                disp(  strcat('...See Schedule_Control routine for additional information')    )
                break;
                
            end
            
        end % End of convergence check 
                        
    end % End of time loop

    average_newtons = num_newton_steps / num_time_steps; 
    msg0 = strcat('                          SUMMARY');
    msg1 = strcat('TOTAL NUMBER OF TIME STEPS: ', num2str(num_time_steps));
    msg2 = strcat('TOTAL NUMBER OF TIME CUTS: ', num2str(Chop_time_counter));
    msg3 = strcat('TOTAL NUMBER OF NEWTON ITERATIONS: ', num2str(num_newton_steps));
    msg4 = strcat('AVERGAE NEWTON ITERATIONS: ', num2str(average_newtons));
    disp(msg0);
    disp(msg1);
    disp(msg2);
    disp(msg3);
    disp(msg4);
    disp('RUN FINISHED')
    
    flag_run_finish_ok = true;
    
    % Plot Jacobian condition number and newton iterations for check
    x1 = 1: size(Newtons, 2); 
    x2 = 1: size(Jacobian_cond, 2);
    ax1 = subplot(2,1,1); % top subplot
    ax2 = subplot(2,1,2); % bottom subplot 
    title('Numerical Performance');
        
    plot(ax1, x1, Newtons, 'b-h', ...
        'LineWidth',1,...
        'MarkerSize',5,...
        'MarkerEdgeColor','m');
    xlabel(ax1, 'Time Step');
    ylabel(ax1, 'Newton Iteration Number');

    plot(ax2, x2,Jacobian_cond, 'k-s', ...
        'LineWidth',1,...
        'MarkerSize',5,...
        'MarkerEdgeColor','r');
    xlabel(ax2, 'Newton Iteration');
    ylabel(ax2, 'Jacobian Condition Number');
    
%     save GKapp_data;

%        disp('y_CO2')
%        prop(2).cell_Zi(:, 1)
% 
%        disp('P')
%        prop(2).cell_P
              
end % End of function