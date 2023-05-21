function [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Assemble(Jacobian_nnz_est, RHS_0, Node_Residual_0, Conne_Residual_0, ...
          Well_Residual_0, well_ctr, var_search_dir_0, dt_k, dummy, sim_ctr, init, fluid, res, rockfluid, well, prop, Ni_hold,       ...
          Nw_hold, P_hold, T_hold, Pwf_hold)

    % Define temporary variables
    nCells = init.nCells;
    nConnes = init.nConnes;
    nComps = fluid.nComps;
    nWells = well.nWells;
    next = 2;

    % Assign back primary variables
    prop(next).cell_Ni = Ni_hold;
    prop(next).cell_Nw = Nw_hold;
    prop(next).cell_P = P_hold;
    prop(next).cell_T = T_hold;
    prop(next).well_P = Pwf_hold;
    
    % Initialize triplet vector to save sparse Jacobian 
    Triplet_I = zeros(Jacobian_nnz_est, 1);
    Triplet_J = zeros(Jacobian_nnz_est, 1);
    Triplet_Val = zeros(Jacobian_nnz_est, 1); 

    % Lable primary variable from 1 to the end
    iRVar = 0; 

    % Count nonzeros in Jacobian
    nonzeros_JVal = 0;
     
    % Perturb each primary variables in each cell
    Perturb_Pwf = false;
    Perturb_Nw  = false;  % No need flash
    
    Ni_perturb_hold = sim_ctr.rel_perturb * max(Ni_hold(:));
    Nw_perturb_hold = sim_ctr.rel_perturb * max(Nw_hold(:));
    P_perturb_hold  = sim_ctr.rel_perturb * max(P_hold(:));
    
    for iCell = 1 : nCells
        
        % Variable sequence in reservoir variable X_res per reservoir node:
        % [Ni, Nw, P]*nCells, Pwf1, Pwf2, ..., Pwfn 
        
        % Get perturbed node and connection list to calculate residual (RHS) 
        [Target_Nodes, Target_Connes] = GetPool(iCell, nCells, nConnes, init.conne_cell1, init.conne_cell2);        
                        
        %% Perturb Ni
        for iComp = 1 : nComps
            % Variable sequence # 
            iRVar = iRVar + 1;
            Perturb_Nw = false;
            % Variable search direction
            searching_dir = var_search_dir_0(iRVar);
            
            % Ensure perturbation is an exact computer-precision number and improves 
            % computational efficiency:
            Ni_perturb = Ni_perturb_hold;
            temp = prop(next).cell_Ni(iCell, iComp) + Ni_perturb * searching_dir;
            Ni_perturb = (temp - prop(next).cell_Ni(iCell, iComp)) / searching_dir; 
            
            prop(next).cell_Ni(iCell, iComp) = prop(next).cell_Ni(iCell, iComp) + Ni_perturb * searching_dir;
                     
            [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Perturb(iCell, iRVar, Ni_perturb, searching_dir, nonzeros_JVal, ...
             Triplet_I, Triplet_J, Triplet_Val, RHS_0, Node_Residual_0, Conne_Residual_0, Well_Residual_0, well_ctr, dt_k, dummy, init, fluid,   ...
             res, rockfluid, well, prop, Perturb_Pwf, Perturb_Nw, Target_Nodes, Target_Connes);       

            prop(next).cell_Ni = Ni_hold;
            prop(next).cell_Nw = Nw_hold;
            prop(next).cell_P = P_hold;
            prop(next).cell_T = T_hold;
            prop(next).well_P = Pwf_hold;
            
        end   
        
        %% Perturb Nw
        % Variable sequence # 
        iRVar = iRVar + 1; 
        Perturb_Nw = true;
        % Variable search direction
        searching_dir = var_search_dir_0(iRVar);
        
        % Ensure perturbation is an exact computer number and improves 
        % computational efficiency
        Nw_perturb = Nw_perturb_hold; 
        temp = prop(next).cell_Nw(iCell) + Nw_perturb * searching_dir;
        Nw_perturb = (temp - prop(next).cell_Nw(iCell)) / searching_dir; 
                
        prop(next).cell_Nw(iCell) = prop(next).cell_Nw(iCell) + Nw_perturb * searching_dir; 
        
        [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Perturb(iCell, iRVar, Nw_perturb, searching_dir, nonzeros_JVal, ...
         Triplet_I, Triplet_J, Triplet_Val, RHS_0, Node_Residual_0, Conne_Residual_0, Well_Residual_0, well_ctr, dt_k, dummy, init, fluid,   ...
         res, rockfluid, well, prop, Perturb_Pwf, Perturb_Nw, Target_Nodes, Target_Connes);       
        
        prop(next).cell_Ni = Ni_hold;
        prop(next).cell_Nw = Nw_hold;
        prop(next).cell_P = P_hold;
        prop(next).cell_T = T_hold;
        prop(next).well_P = Pwf_hold;
        
        %% Perturb Po
        % Variable sequence # 
        iRVar = iRVar + 1; 
        Perturb_Nw = false;
        % Variable search direction
        searching_dir = var_search_dir_0(iRVar);
        
        % Ensure perturbation is an exact computer number and improves 
        % computational efficiency
        P_perturb = P_perturb_hold;
        temp = prop(next).cell_P(iCell) + P_perturb * searching_dir;
        P_perturb = (temp - prop(next).cell_P(iCell)) / searching_dir;
        
        % Get bubble point pressure, when pressure get very close to Pb,
        % make it cover the Psat
        % Pb = res.Pb(init.cell_type(iCell));
        % if abs(prop(next).cell_P(iCell) - Pb) < 1.0e-2
        %    if prop(next).cell_P(iCell) > Pb && searching_dir == -1
        %         P_perturb = P_perturb + (prop(next).cell_P(iCell) - Pb);
        %    end
        %    if prop(next).cell_P(iCell) < Pb && searching_dir == 1
        %        P_perturb = P_perturb + (Pb - prop(next).cell_P(iCell));
        %    end
        % end
                
        prop(next).cell_P(iCell) = prop(next).cell_P(iCell) + P_perturb * searching_dir; 
        
        [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Perturb(iCell, iRVar, P_perturb, searching_dir, nonzeros_JVal, ...
         Triplet_I, Triplet_J, Triplet_Val, RHS_0, Node_Residual_0, Conne_Residual_0, Well_Residual_0, well_ctr, dt_k, dummy, init, fluid,   ...
         res, rockfluid, well, prop, Perturb_Pwf, Perturb_Nw, Target_Nodes, Target_Connes);       
        
        prop(next).cell_Ni = Ni_hold;
        prop(next).cell_Nw = Nw_hold;
        prop(next).cell_P = P_hold;
        prop(next).cell_T = T_hold;
        prop(next).well_P = Pwf_hold;
        
    end
    
    %% Perturb Pwf
    Perturb_Pwf = true;
    Perturb_Nw  = false;
    for iWell = 1 : nWells
        
        iCell = 0;
        
        % Get perturbed node and connection list to calculate residual (RHS) 
        [Target_Nodes, Target_Connes] = GetPool(iCell, nCells, nConnes, init.conne_cell1, init.conne_cell2);        
        % Variable sequence #
        iRVar = iRVar + 1; 
        % Variable search direction
        searching_dir = var_search_dir_0(iRVar);

        % Ensure perturbation is an exact computer number and improves 
        % computational efficiency
        Pwf_perturb = sim_ctr.rel_perturb * Pwf_hold(iWell);  

        temp = prop(next).well_P(iWell) + Pwf_perturb * searching_dir;
        Pwf_perturb = (temp - prop(next).well_P(iWell)) / searching_dir; 

        prop(next).well_P(iWell) = prop(next).well_P(iWell) + Pwf_perturb * searching_dir; 

        [Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Perturb(iCell, iRVar, Pwf_perturb, searching_dir, nonzeros_JVal, ...
         Triplet_I, Triplet_J, Triplet_Val, RHS_0, Node_Residual_0, Conne_Residual_0, Well_Residual_0, well_ctr, dt_k, dummy, init, fluid,   ...
         res, rockfluid, well, prop, Perturb_Pwf, Perturb_Nw, Target_Nodes, Target_Connes);       

        prop(next).cell_Ni = Ni_hold;
        prop(next).cell_Nw = Nw_hold;
        prop(next).cell_P = P_hold;
        prop(next).cell_T = T_hold;
        prop(next).well_P = Pwf_hold;    
        
    end
    
    % Set back to original values
    prop(next).cell_Ni = Ni_hold;
    prop(next).cell_Nw = Nw_hold;
    prop(next).cell_P = P_hold;
    prop(next).cell_T = T_hold;
    prop(next).well_P = Pwf_hold;

    % Save as sparse for triplet 
    Triplet_I = Triplet_I(1 : nonzeros_JVal);
    Triplet_J = Triplet_J(1 : nonzeros_JVal);
    Triplet_Val = Triplet_Val(1 : nonzeros_JVal);    
        
    Triplet_I = sparse(Triplet_I);
    Triplet_J = sparse(Triplet_J);
    Triplet_Val = sparse(Triplet_Val);    

end