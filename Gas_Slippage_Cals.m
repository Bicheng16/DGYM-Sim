function [prop] = Gas_Slippage_Cals(stat, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw)
    
    % slippage_flag: Flag to consider gas slippage    
    slippage_flag = (sum(res.GKap) > 0);         % At least one rock type consider gas slippage 
    
    %% If slippage is not considered in the whole reservoir, return
    if ~slippage_flag
        prop(stat).cell_GKapp(:, :) = 1.0d0;
        return;
    end
    
    %% If no gas phase exist in the whole reservoir, return
    if sum(prop(stat).cell_fv(:)) == 0
        prop(stat).cell_GKapp(:, :) = 1.0d0;
        return;
    end    

    %% Inheritance to avoid unnecessary calculation
    if inherit && ~isempty(inherit_loc)
        if(Perturb_Nw)
            % If perturb Nw, slippage is not updated, inherit from dummy as a
            % whole            
            prop(stat).cell_GKapp(:, :) = dummy.GKapp(:, :);
            return;
        else
            % Normal inheritance
            prop(stat).cell_GKapp(inherit_loc, :) = dummy.GKapp(inherit_loc, :);        
            % Inherit all reservoir properties from dummy when perturbing well pressure 
            if nnz(inherit_loc) == init.nCells
                return
            end                        
        end
    end
    
    %% Apparent permeability model
    Non_slippage = 0;
    Civan_slippage = 1;

    % lbmol/ft^3 to kmol/m^3: 1.601848973047802d1    
    % Constant: Eq. 4 in SPE 179704
    coef = sqrt(2.0d0) * pi * 1.601848973047802d1 * 6.022d26; 
    dmc2 = coef .* fluid.comp_dmc .^ (2.0d0);
    nComps = fluid.nComps;
         
    % Cell rock type for gas adsorption calculation
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
    
    % Calcualte gas adsorption
    for iRT = 1 : res.nRockTypes

        % Get cell names belonging to each rock type
        tmp_nodes = Target_Nodes_Group(iRT).iGroup;
        prop(stat).cell_GKapp(tmp_nodes, :) = 1.0d0;
        
        % No nodes
        if isempty(tmp_nodes)
            continue; 
        end
        
        % No slippage
        if (res.GKap(iRT) == Non_slippage)
            continue;             
        end
        
        % No gas phase 
        if sum(prop(stat).cell_fv(tmp_nodes)) == 0
            continue;
        end
        
        % Initialization and temporary variables 
        nNodes = nnz(tmp_nodes);
        
        if(res.GKap(iRT) == Civan_slippage) 
            
            lambda = zeros(nNodes, nComps); % Mean free path of gas moleculars
            alphak = zeros(nNodes, nComps);

            dmc2_cell = repmat(dmc2, nNodes, 1);
            MolarDensVap = prop(stat).cell_MolarDensVap(tmp_nodes);
            MolarDensVap_cell = MolarDensVap(:, ones(1, nComps) );  

            lambda(:, :) = 1.0d0 ./ (dmc2_cell .* MolarDensVap_cell);

            % Knudsen number
            lambda = lambda ./ res.dpore(iRT);

            alphak = (128.0d0/15.0d0/pi^2) .* atan(4.0d0 .* lambda .^ 0.4d0);
            
            lambda = (1.0d0 + alphak .* lambda) .* (1.0d0 + 4.0d0 .* lambda ./ (1.0d0 + lambda) );  

            prop(stat).cell_GKapp(tmp_nodes, :) = lambda(:, :); 
            
        else
        
            prop(stat).cell_GKapp(tmp_nodes, :) = 1.0d0;
            
        end
               
    end
        
end