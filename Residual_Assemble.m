    
function [RHS, R_c, R_wat, R_vol, Balance_hc, Balance_wat, Balance_vol] = Residual_Assemble( Conne_Residual, ...
         Node_Residual, R_well, nCells, nComps, Node1_in_conne, Node2_in_conne)
    
    % Assemble residual from node, connection, well residual pieces 
    
	%% Memory preallocation 
    % Reservoir residual 
    R_c     = zeros([nCells, nComps]);
    R_wat   = zeros([nCells, 1]);
    R_vol   = zeros([nCells, 1]);
    
    % Material balance error
    Balance_hc  = zeros([nComps, 2]);  % Hydrocarbon: lb-mol 
    Balance_wat = zeros([1, 2]);       % Water: lb-mol
    Balance_vol = zeros([1, 2]);       % Volume: cu.ft
    
	%% Data scanning 
	%  Reservoir residual 
                            
	R_c(1:nCells, 1:nComps) = Node_Residual(1:nCells, 1:nComps);    % Hydrocarbon accumulation
	R_wat(1:nCells, 1)      = Node_Residual(1:nCells, nComps+1);    % Water accumulation
    R_vol(1:nCells, 1)      = Node_Residual(1:nCells, nComps+2);    % Volume balance
    
	for iCell = 1 : nCells
        %  Hydrocarbon: add flux to accumulation
        R_c(iCell, 1:nComps) = R_c(iCell, 1:nComps) - sum(Conne_Residual((Node1_in_conne == iCell), 1:nComps), 1) ...
                                                    + sum(Conne_Residual((Node2_in_conne == iCell), 1:nComps), 1);
        % Water: add flux to accumulation                              
        R_wat(iCell, 1) = R_wat(iCell, 1) - sum(Conne_Residual((Node1_in_conne == iCell), nComps+1), 1) ...
                                          + sum(Conne_Residual((Node2_in_conne == iCell), nComps+1), 1);
	end

    % Calculate material balance error: column 1 - L-2 norm; column 2 - infinite norm (max absolute) 
    for iComp = 1 : nComps
        Balance_hc(iComp, 1) = norm(R_c(:, iComp), 2);    % Hydrocarbon: lb-mol 
        Balance_hc(iComp, 2) = norm(R_c(:, iComp), inf);
    end
    
    Balance_wat(1, 1) = norm(R_wat(:, 1), 2);             % Water: lb-mol  
    Balance_wat(1, 2) = norm(R_wat(:, 1), inf);

    Balance_vol(1, 1) = norm(R_vol(:, 1), 2);             % Volume: cu.ft  
    Balance_vol(1, 2) = norm(R_vol(:, 1), inf);
    
    %% Compact RHS
    R_c     =	reshape(R_c.',     nComps, nCells);
    R_wat   =	reshape(R_wat,          1, nCells); 
    R_vol   =	reshape(R_vol,          1, nCells);
        
	% Convert to right hand side 
	RHS     =   [R_c; R_wat; R_vol];
	RHS		=	reshape(RHS, nCells*(nComps+2), 1);
    
    % Add well as the last residual
    RHS     =   [RHS; R_well];
    RHS     =   sparse(RHS);
                
end
    
    
    
    
