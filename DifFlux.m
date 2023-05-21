function [DifFlux_Residual] = DifFlux(init, fluid, fickian_flag, prop, res, Target_Connes)

    % Calculate flux based on Generalized Fick's Law    
    % No diffusion considered, return
    if (~fickian_flag)
        DifFlux_Residual = [];
        return;
    end 
    
    % Define temporary variables 
    nComps = fluid.nComps;
    nConnes = init.nConnes;
    actConnes = nnz(Target_Connes);
    next = 2;    
    
    % Memory allocation  
    delta_yv = zeros(nConnes, nComps-1);
    Lambda_gd = zeros(nConnes, 1);        
    trans_dif   = zeros(nComps-1, nComps-1);
    trans_dif_1 = zeros(nComps-1, nComps-1);
    trans_dif_2 = zeros(nComps-1, nComps-1);
    DifFlux_Residual = zeros(nConnes, nComps);
     
    % Connection
    Node1_in_conne = init.conne_cell1(Target_Connes);
    Node2_in_conne = init.conne_cell2(Target_Connes);
        
    % Gas mole fraction difference 
    delta_yv(Target_Connes, 1 : nComps-1) = ...
    prop(next).cell_Yi(Node2_in_conne, 1 : nComps-1) - prop(next).cell_Yi(Node1_in_conne, 1 : nComps-1);
    
    % Arithmetic mobility: Lambda_gd = rho_g * phi * Sg 
    Lambda_gd(Target_Connes) = (prop(next).cell_MolarDensVap(Node1_in_conne) + prop(next).cell_MolarDensVap(Node2_in_conne)) ...
                            .* (prop(next).cell_poro(Node1_in_conne) + prop(next).cell_poro(Node2_in_conne)) ...
                            .* (prop(next).cell_Sg(Node1_in_conne) + prop(next).cell_Sg(Node2_in_conne)) ./ 8.0d0;
    Lambda_gd(isnan(Lambda_gd) | isinf(Lambda_gd)) = 0.0d0;
    
    % Calculate diffusion transmissibility
    for iConne = 1 : actConnes
        
        ConneID = Target_Connes(iConne);
        
        % Nodes in connection pair  
        Node1 = Node1_in_conne(iConne); 
        Node2 = Node2_in_conne(iConne); 
        
        % If non-diffusion considered in both cell
        if( res.Fickian(init.cell_type(Node1)) == 0 && ...
            res.Fickian(init.cell_type(Node2)) == 0 )    
            continue;
        end

        % Initialize
        trans_dif_1(:, :) = 0.0d0;
        trans_dif_2(:, :) = 0.0d0; 
        trans_dif(:, :)   = 0.0d0;
        
        % Unified consideration for both neighbor / nonneighbor connection 
        % Geometric coef
        trans_dif_1(:, :) = init.conne_geo(ConneID, 1);
        trans_dif_2(:, :) = init.conne_geo(ConneID, 2);
                  
        % Inverse of Cell 1 Half
        trans_dif_1(1:nComps-1, 1:nComps-1) = ...
        1.0d0 ./ ( prop(next).cell_GFD(Node1, 1:nComps-1, 1:nComps-1) .* trans_dif_1(1:nComps-1, 1:nComps-1) );

        % Inverse of Cell 2 Half
        trans_dif_2(1:nComps-1, 1:nComps-1) = ...
        1.0d0 ./ ( prop(next).cell_GFD(Node2, 1:nComps-1, 1:nComps-1) .* trans_dif_2(1:nComps-1, 1:nComps-1) );

        % Eliminate nonsense
        trans_dif_1(isnan(trans_dif_1) | isinf(trans_dif_1)) = 0.0d0;
        trans_dif_2(isnan(trans_dif_2) | isinf(trans_dif_2)) = 0.0d0;

        % Harmonic 
        trans_dif(1:nComps-1, 1:nComps-1) = 1.0d0 ...
        ./ ( trans_dif_1(1:nComps-1, 1:nComps-1) + trans_dif_2(1:nComps-1, 1:nComps-1)); 
       
        % Component 1 to nc-1 
        for iComp = 1 : nComps-1 
            DifFlux_Residual(ConneID, iComp) = sum(trans_dif(iComp, 1:nComps-1) .* delta_yv(ConneID, 1 : nComps-1)); 
        end
        
        % Total flux balance constrain to calculate solute diffusion flux
        DifFlux_Residual(ConneID, nComps) = - sum(DifFlux_Residual(ConneID, 1 : nComps-1)); 
        
        % Ultimately times diffusive mobility
        DifFlux_Residual(ConneID, :) = Lambda_gd(ConneID) .* DifFlux_Residual(ConneID, :);  
       
    end    
    
end