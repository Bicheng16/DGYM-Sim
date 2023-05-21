% Engine for GURU
function GURU(input)
    %% ************************************************************************
    %*************            Part I Environment Setup             ************
    %**************************************************************************

    CPU_time_start = Env_Setup;

    
    %% ************************************************************************
    %*************            Part II Input Scanning               ************
    %**************************************************************************

    % Scan input data file
    [init, fluid, res, rockfluid, sim_ctr, well] = Scan_Input(input);
    
    % Set basic property vectors and check input
    
    [init, fluid, well, rockfluid, flag_run_finish_ok] = Misc_Cals(init, fluid, ...
    res, rockfluid, sim_ctr, well);

    % Terminate code if error in input
    if not(flag_run_finish_ok)
        Run_Finished(CPU_time_start, 0, flag_run_finish_ok);
        return;
    end
    
    %% ************************************************************************
    %*****************   Part III: Initialization    **************************
    %**************************************************************************

    [init, prop, rst, res, fluid, well, well_ctr] = Init_Cals(init, fluid, res, well, rockfluid);
    
    %% ************************************************************************
    %*****************       Part IV: Time Loop      **************************
    %**************************************************************************

    [flag_run_finish_ok, Chop_time_counter] = Time_Run(init, fluid, res, rockfluid, sim_ctr, well, well_ctr, prop, rst, CPU_time_start, input);    

    % Terminate code if error in input
    if not(flag_run_finish_ok)
        Run_Finished(CPU_time_start, Chop_time_counter, flag_run_finish_ok);
        return;
    end    
   
end

