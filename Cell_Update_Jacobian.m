% This function is to update cell based properties at newton level
function [dummy, prop] = Cell_Update_Jacobian(Target_Nodes, dummy, init, fluid, res, rockfluid, prop, Perturb_Nw)

    % Define temporary variable
    next = 2;

    %% Boolean for grid properties inheritance
    if nnz(Target_Nodes) == 1
        % Numerical perturbation for reservoir variables
        inherit = true;
        inherit_loc = [ 1 : Target_Nodes-1  Target_Nodes+1 : init.nCells ];
    elseif Target_Nodes == 0
        % Numerical perturbation for well variables
        inherit = true;
        inherit_loc = 1 : init.nCells;
    else
        % RHS calculation
        inherit = false;
        inherit_loc = [];
    end    
    export = (~inherit);  % If no inherit, just export it since it is calculating RHS  
    
    %% Cell-Based Variables 
    % Update newton level properties: prop(next)
    % Calculate fluids properties
    prop = Properties_Fluids_Jacobian(next, prop, fluid, init, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw); 
    
    % Calculate rock properties
    prop = Properties_Rock(next, prop, init, res, rockfluid, dummy, Target_Nodes, inherit, inherit_loc); 

    % Adsorbed gas as secondary variable 
    % [prop] = Gas_Adsorption_Cals(res, init, next, prop);
    prop = Gas_Adsorption_Cals(next, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw);
  
    % Gas slippage apparent permeability multiplier
    prop = Gas_Slippage_Cals(next, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw);

    % General Fickian diffusion
    prop = Gas_GFick_Cals(next, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw);
        
    % Export cell properties if no inheritance is used or if calculating RHS
    dummy = update_inherit_var(export, dummy, next, prop);
    
end