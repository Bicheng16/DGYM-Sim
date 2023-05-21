function [prop] = Gas_Adsorption_Cals(stat, prop, init, res, fluid, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw)
    
    % sorption_flag: Flag to consider adsorption
    ns = 0;
    for iRT = 1 : res.nRockTypes
        if ~ strcmp(res.Sorption{iRT}, 'NA')
            ns = ns + 1;
        end
    end
    
    sorption_flag = (ns > 0);         % At least one rock type consider adsorption 
    
    %% If adsorption is not considered in the whole reservoir
    if ~sorption_flag
        prop(stat).cell_Qai(:, :) = 0.0d0;
        return;
    end
    
    %% If no gas phase exist in the whole reservoir, return
    if sum(prop(stat).cell_fv(:)) == 0
        prop(stat).cell_Qai(:, :) = 0.0d0;
        return;
    end    

    %% Inheritance to avoid unnecessary calculation
    if inherit && ~isempty(inherit_loc)
        if(Perturb_Nw)
            % If perturb Nw, adsorption is not updated, inherit from dummy as a
            % whole            
            prop(stat).cell_Qai(:, :) = dummy.Qai(:, :);
            return;
        else
            % Normal inheritance
            prop(stat).cell_Qai(inherit_loc, :) = dummy.Qai(inherit_loc, :);        
            % Inherit all reservoir properties from dummy when perturbing well pressure 
            if nnz(inherit_loc) == init.nCells
                return
            end                        
        end
    end
    
    %% Extend Langmuir Adsorption Model
    % Define temporary variables 
    nComps = fluid.nComps;     
    MolarVolVap_sc = 379.342;         % Gas molar volume in standard condition, scf/lb-mol. 
         
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
        prop(stat).cell_Qai(tmp_nodes, :) = 0.0d0;
        
        % No nodes
        if isempty(tmp_nodes)
            continue; 
        end
        
        % No desorption
        if strcmp(res.Sorption{iRT}, 'NA')
            continue;             
        end
        
        % No gas phase 
        if sum(prop(stat).cell_fv(tmp_nodes)) == 0
            continue;
        end
        
        % Initialization and temporary variables 
        nNodes = nnz(tmp_nodes);
        Qai = zeros(nNodes, nComps);
        Pi  = zeros(nNodes, nComps);
        Yi  = zeros(nNodes, nComps);
        
        cell_VLi = zeros(nNodes, nComps);
        cell_PLi = zeros(nNodes, nComps);
        cell_RMDen = zeros(nNodes, 1);        

        cell_VLi(:, :) = res.VLi(init.cell_type(tmp_nodes), 1:nComps);
        cell_PLi(:, :) = res.PLi(init.cell_type(tmp_nodes), 1:nComps);
        cell_RMDen(:, 1) = res.RMDen(iRT);
        
        % Extended Langmuir Adsorption Model 
        % Calculate tmp1 = Yi * P / PL_i
        p_rep = repmat(prop(stat).cell_P(tmp_nodes), 1, nComps);
        tmp1 = prop(stat).cell_Yi(tmp_nodes, :) .* p_rep ./ cell_PLi;
        
        % Calculate tmp2 = 1.0 + sum(Yi * P / PL_i, i = 1, ..., nc)
        tmp2 = 1.0 + sum(tmp1, 2);
        
        % Calculate tmp4 = (1.0 - poro) * RMDen / Vgs / [ 1.0 + sum(Yi * P / PL_i, i = 1, ..., nc) ]
        tmp3 = (1.0 - prop(stat).cell_poro(tmp_nodes)) .* cell_RMDen ./ MolarVolVap_sc ./ tmp2;
        tmp4 = repmat(tmp3, 1, nComps);
        
        % Calculate tmp5 = VLi * (Yi * P / PL_i) 
        tmp5 = cell_VLi .* tmp1;   

        % Calculate VLi * { (Yi * P / PL_i) / (1.0 + sum(Yi * P / PL_i, i = 1, ..., nc) }
        prop(stat).cell_Qai(tmp_nodes, :) = tmp5 .* tmp4;  
        
    end
    
    
 end