function[Triplet_I, Triplet_J, Triplet_Val, nonzeros_JVal] = Jacobian_Perturb(iCell, iRVar, perturb, searching_dir, nonzeros_JVal, ...
         Triplet_I, Triplet_J, Triplet_Val, RHS_0, Node_Residual_0, Conne_Residual_0, Well_Residual_0, well_ctr, dt_k, dummy, init, fluid,   ...
		 res, rockfluid, well, prop, Perturb_Pwf, Perturb_Nw, Target_Nodes, Target_Connes)
     
    % Define temporary variables
    nComps = fluid.nComps;
    nCells = init.nCells;
    Node1_in_conne = init.conne_cell1;
    Node2_in_conne = init.conne_cell2;

    %% Perturbation Analysis
    % Update cell properties first
    %[~, prop] = Newton_Cell_Update(iCell, dummy, init, fluid, res, rockfluid, prop);
    [~, prop] = Cell_Update_Jacobian(iCell, dummy, init, fluid, res, rockfluid, prop, Perturb_Nw);
                
    % Calculate residual pieces
    [Node_Residual_1, Conne_Residual_1, Well_Residual_1, ~, ~, ~] = ...
     Residual_Pieces(Target_Nodes, Target_Connes, dt_k, init, fluid, res, well, prop, well_ctr);
 
    if Perturb_Pwf
        [RHS_1, ~, ~, ~, ~, ~, ~] = Residual_Assemble( Conne_Residual_1, Node_Residual_1, Well_Residual_1,    ...
        nCells, nComps, init.conne_cell1, init.conne_cell2); 
    
        Jacobian = (RHS_1 - RHS_0) ./ perturb .* searching_dir; 
        
        nonzeros = nnz(Jacobian);
        %ILoc = zeros(nonzeros, 1);
        JLoc = zeros(nonzeros, 1);
        Val = zeros(nonzeros, 1);
        
        [ILoc, ~] = find(Jacobian~=0);
        JLoc(1:nonzeros, 1) = iRVar;
        Val(1:nonzeros, 1) = Jacobian(ILoc, 1); 
        
    else
        %% Reservoir Jacobian Function Evaluation: Rrr = dRr/dXr
        [ILoc_r, JLoc_r, Val_r, nonzeros_r] = ...
        Jacobian_Res_Evaluation( Node_Residual_0, Conne_Residual_0, ...
                             Node_Residual_1, Conne_Residual_1, ...
                             Node1_in_conne,  Node2_in_conne,   ...
                             iRVar, nComps, nCells, perturb,    ...
                             searching_dir, iCell);

        %% Well Jacobian Function Evaluation: Rwr = dRw/dXr
        [ILoc_w, JLoc_w, Val_w, nonzeros_w] = Jacobian_Well_Evaluation(Well_Residual_0, Well_Residual_1, ...
        iRVar, nComps, nCells, perturb, searching_dir, well, iCell);

        %% Combine those two Jacobians
        ILoc = [ILoc_r; ILoc_w];
        JLoc = [JLoc_r; JLoc_w];
        Val = [Val_r; Val_w]; 
        nonzeros = nonzeros_r + nonzeros_w; 
    end
                                          
    %% Transfer to Global Triplet Format
    start = nonzeros_JVal + 1;

    if ~isempty(ILoc) && ~isempty(JLoc) && ~isempty(Val) && nonzeros > 0 
        
        % Check space remaining 
        remain_space = size(Triplet_I, 1) - nonzeros_JVal;   % Avoid no enough space to save the triplet
        if remain_space <= nonzeros
            % Add additional space for safeguard of memory
            Triplet_I = [Triplet_I; zeros([nonzeros-remain_space+1, 1])];
            Triplet_J = [Triplet_J; zeros([nonzeros-remain_space+1, 1])];
            Triplet_Val = [Triplet_Val; zeros([nonzeros-remain_space+1, 1])];
            remain_space
            nonzeros
        end
        
        % Fill in triplet 
        Triplet_I(start : start+nonzeros-1) = ILoc;
        Triplet_J(start : start+nonzeros-1) = JLoc;
        Triplet_Val(start : start+nonzeros-1) = Val; 
        
        % Update remaining space size 
        nonzeros_JVal = nonzeros_JVal + nonzeros; 
        
    end
    
end