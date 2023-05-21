function [prop] = Gas_GFick_Cals(stat, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw)
    
    % fick_flag: Flag to consider gas Fickian diffusion    
    fick_flag = (sum(res.Fickian) > 0);         % At least one rock type consider gas diffusion 
    
    %% If diffusion is not considered in the whole reservoir, return
    if ~fick_flag
        prop(stat).cell_GFD(:, :, :) = 0.0d0;
        return;
    end
    
    %% If no gas phase exist in the whole reservoir, return
    if sum(prop(stat).cell_fv(:)) == 0
        prop(stat).cell_GFD(:, :, :) = 0.0d0;
        return;
    end    

    %% Inheritance to avoid unnecessary calculation
    if inherit && ~isempty(inherit_loc)
        if(Perturb_Nw)
            % If perturb Nw, slippage is not updated, inherit from dummy as a
            % whole            
            prop(stat).cell_GFD(:, :, :) = dummy.GFD(:, :, :);
            return;
        else
            % Normal inheritance
            prop(stat).cell_GFD(inherit_loc, :, :) = dummy.GFD(inherit_loc, :, :);        
            % Inherit all reservoir properties from dummy when perturbing well pressure 
            if nnz(inherit_loc) == init.nCells
                return
            end                        
        end
    end
    
    %% Update Fickian diffusion coefficient
    nComps = fluid.nComps; 
    Non_Fickian = 0;
    
    % Flag for diffusion coefficient calculation 
    Dc_model = 0;   % 0: from input, fixed for the whole simualtion 
                    % 1: Leahy-Dios and Firoozabadi's model, dynamically
                    % changing
    
    % Cell rock type for gas diffusion calculation
    % Perturb reservoir variables
    if nnz(Target_Nodes) == 1 && Target_Nodes ~= 0
        Target_Nodes_Group = [];
        for iRT = 1 : res.nRockTypes
            dummy_1.iGroup = [];
            Target_Nodes_Group = [Target_Nodes_Group; dummy_1];
        end                
        Target_Nodes_Group(init.cell_type(Target_Nodes)).iGroup = Target_Nodes;          
    else
    % Calculate RHS 
        Target_Nodes_Group = init.cell_group;
    end
    
    % Calcualte gas diffusion
    for iRT = 1 : res.nRockTypes

        % Get cell names belonging to each rock type
        tmp_nodes = Target_Nodes_Group(iRT).iGroup;
        prop(stat).cell_GFD(tmp_nodes, :, :) = 0.0d0;
        
        % No nodes
        if isempty(tmp_nodes)
            continue; 
        end
        
        % No diffusion
        if (res.Fickian(iRT) == Non_Fickian)
            continue;             
        end
        
        % No gas phase 
        if sum(prop(stat).cell_fv(tmp_nodes)) == 0
            continue;
        end
        
        % Initialization and temporary variables 
        nNodes = nnz(tmp_nodes);
        for i = 1 : nNodes 

            nodeID = tmp_nodes(i);
            
            % Diffusivity from input
            if Dc_model == 0
                prop(stat).cell_GFD(nodeID, :, :) = res.GFD(iRT, :, :);
            end
            
            % Diffusivity from thermodynamics
            if Dc_model == 1
                prop(stat).cell_GFD(nodeID, :, :) = 0.0d0;
            end
            
            
        end
                       
    end
        
end