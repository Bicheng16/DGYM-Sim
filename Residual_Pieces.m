function [Node_Residual, Conne_Residual, R_well, W_SS_Residual, perf_P, well_ctr] = ...
          Residual_Pieces(Target_Nodes, Target_Connes, dt_k, init, fluid, res, well, prop, well_ctr)
    % This function is to calculate residual pieces for Jacobian and Residual (RHS).  
    
    slippage_flag = (sum(res.GKap) > 0);         % At least one rock type consider gas slippage 
    fickian_flag  = (sum(res.Fickian) > 0);      % At least one rock type consider Fickian diffusion 
    
    %% Connection-Based Variables 
    % Upstream potential 
    [delta_pot_o, delta_pot_w, delta_pot_g, flag_upstream_o, flag_upstream_w, flag_upstream_g] ...
    = Upstream_Potential(prop, init, Target_Connes);  

    % Convection Flux Mobility
    [a_o, a_w, a_g, a_g_app] = Mobility_Convec_Flux(init, fluid, slippage_flag, prop, Target_Connes, ...
    flag_upstream_o, flag_upstream_w, flag_upstream_g);

    % Diffusion Flux
    [DifFlux_Residual] = DifFlux(init, fluid, fickian_flag, prop, res, Target_Connes);
    
    %% Well control and well residual piece
    [W_SS_Residual, R_well, perf_P, well_ctr] = Residual_Well(init, well, prop, fluid, well_ctr, Target_Nodes);
    
    %% Reservoir residual piece: 
    % Node_Residual of reservoir grid 
    % Conne_Residual including reservoir nodes (reservoir grid to grid)    
    [Node_Residual, Conne_Residual] = Residual_Split(init, fluid, well, prop, Target_Nodes,   ...
    Target_Connes, flag_upstream_o, flag_upstream_g, a_o, a_g, a_w, a_g_app, slippage_flag,   ...
    fickian_flag, DifFlux_Residual, delta_pot_o, delta_pot_g, delta_pot_w, dt_k, W_SS_Residual);

end