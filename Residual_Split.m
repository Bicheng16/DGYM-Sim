function [Node_Residual, Conne_Residual] = Residual_Split(init, fluid, well, prop, ...
    Target_Nodes, Target_Connes, flag_upstream_o, flag_upstream_g, a_o, a_g, a_w,  ...
    a_g_app, slippage_flag, fickian_flag, DifFlux_Residual, delta_pot_o,           ...
    delta_pot_g, delta_pot_w, dt_k, W_SS_Residual)
    
    % Define temporary variables 
    nComps  = fluid.nComps;
    nCells  = init.nCells;
    nConnes = init.nConnes;
    nPerfs  = sum(well.nperf_cell);
    current = 1;
    next    = 2;
        
    % Residual Piece Calculation in Sequence: primary equation: Rc, Rwat, Rvol   
    Node_Residual = zeros(nCells, (nComps+2));
    Conne_Residual = zeros(nConnes, (nComps+1));
        
    % Memory allocation 
    Xi_conne = zeros(nConnes, 1); 
    Yi_conne = zeros(nConnes, 1);   

    Node1_in_conne = init.conne_cell1(Target_Connes);
    Node2_in_conne = init.conne_cell2(Target_Connes);
        
    %% Conne_Residual: include connection residual of [Rc Rw]' * nConnes    
    % Connection residual in Rc 
    for iComp = 1 : nComps
        
        % Upstream hydrocarbon mole fraction 
        Xi_conne(Target_Connes) = prop(next).cell_Xi(Node1_in_conne, iComp) .* flag_upstream_o(Target_Connes) ...
                                + prop(next).cell_Xi(Node2_in_conne, iComp) .* (1-flag_upstream_o(Target_Connes)); 
       
        Yi_conne(Target_Connes) = prop(next).cell_Yi(Node1_in_conne, iComp) .* flag_upstream_g(Target_Connes) ...
                                + prop(next).cell_Yi(Node2_in_conne, iComp) .* (1-flag_upstream_g(Target_Connes)); 

        % Ensure zeros (instead of NaN) when phase volume is zero for 2 neighbor cells
        Xi_conne(isnan(Xi_conne) | isinf(Xi_conne)) = 0;    
        Yi_conne(isnan(Yi_conne) | isinf(Yi_conne)) = 0;
        
        if (slippage_flag)
            Conne_Residual(Target_Connes, iComp) =                                          ...
            + a_o(Target_Connes) .* Xi_conne(Target_Connes) .* delta_pot_o(Target_Connes)   ...
            + a_g_app(Target_Connes, iComp) .* Yi_conne(Target_Connes) .* delta_pot_g(Target_Connes);                   
        else
            Conne_Residual(Target_Connes, iComp) =                                          ...
            + a_o(Target_Connes) .* Xi_conne(Target_Connes) .* delta_pot_o(Target_Connes)   ...
            + a_g(Target_Connes) .* Yi_conne(Target_Connes) .* delta_pot_g(Target_Connes);       
        end  
        
    end
    
    %% Add diffusion flux for hydrocarbon if considered
    if fickian_flag
        Conne_Residual(Target_Connes, 1:nComps) = Conne_Residual(Target_Connes, 1:nComps) + ...
        DifFlux_Residual(Target_Connes, 1:nComps);
    end
    
    % Connection residual in Rw
    Conne_Residual(Target_Connes, nComps+1) = a_w(Target_Connes) .* delta_pot_w(Target_Connes);
    Conne_Residual = Conne_Residual .* dt_k;
    
    %% Node_Residual: include node residual of [Rc Rw]' * nCells 
    %  Node residual in Rc: hydrocarbon component accumulation     
    Vb_rep = init.cell_volume( :, ones(1, nComps) );
           
    Node_Residual(Target_Nodes, 1:nComps) = Vb_rep(Target_Nodes, 1:nComps) ... % Conventional compressed storage term
        .* (prop(next).cell_Ni(Target_Nodes, 1:nComps)- prop(current).cell_Ni(Target_Nodes, 1:nComps) ) ...
        + Vb_rep(Target_Nodes, 1:nComps) ... % Additional adsorption storage term
        .* (prop(next).cell_Qai(Target_Nodes, 1:nComps)- prop(current).cell_Qai(Target_Nodes, 1:nComps) );
                                                          
    % Node residual in Rw: water accumulation
    Node_Residual(Target_Nodes, nComps+1) = init.cell_volume(Target_Nodes) ... 
        .* (prop(next).cell_Nw(Target_Nodes) - prop(current).cell_Nw(Target_Nodes) ); 
                                       
    % Add well residual back to node residual
    % Residual = Accumulation - Source/Sink (injector: +; producer: -) - Influx (P_neighbor - P_center);
    % Formulation (1): doesn't work if a cell is perforated by multi-wells
    % Node_Residual(well.wconne_cell, :) = Node_Residual(well.wconne_cell, :) - W_SS_Residual_1(:, :); 
    % Formulation (2): works fine if a cell is perforated by multi-wells
    % Since perforated cells number is limited, efficiency is ok in forloop
    for iPerf = 1 : nPerfs
        Node_Residual(well.wconne_cell(iPerf), 1:nComps+1) = Node_Residual(well.wconne_cell(iPerf), 1:nComps+1) - W_SS_Residual(iPerf, 1:nComps+1) .* dt_k;     
    end    
    
    % Last Node residual is volume balance 
    % fluid volume - pore volume
    Node_Residual(Target_Nodes, nComps+2) = prop(next).cell_Vfluid(Target_Nodes) .* init.cell_volume(Target_Nodes) - prop(next).cell_PV(Target_Nodes);

    % Sparse Node_Residual, Conne_Residual
    Node_Residual  = sparse(Node_Residual);
    Conne_Residual = sparse(Conne_Residual); 
    
end