% This function is to evaluate Well Jacobian and store in triplet format

function [ILoc, JLoc, Val, nonzeros] = Jacobian_Well_Evaluation(Well_Residual_0, Well_Residual_1, ...
          iVar, nComps, nCells, Peturb, searching_dir, well, iCell)

    % Evaluate current nonzero residual pieces for Jacobian and generate
    % triplet format of Jacobian
    % Length of triplet   
    
    % Control parameter: searching_dir is to flag for positive or negative
    % perturbation for different parameters, only 2 values: 1 (positive),
    % -1 (negative). Note: pressure is negatively perturbed.

    ILoc = [];
    JLoc = [];
    Val  = [];  
    nonzeros = 0;
    
    % If the perturbed cell is not perforated, return 
    if isempty( find(well.wconne_cell == iCell, 1) )
        return
    end    
          
    % Calculate reservoir equations (material balance and volume balance) 
    nreqT = nCells * (nComps + 2);
    
    well_Jacobian = (Well_Residual_1  -  Well_Residual_0)  ./ Peturb .* searching_dir; 
    
    % If there is no nonzeros values for well Jacobian, return
    nonzeros = nnz(well_Jacobian);
    if nonzeros == 0
        return
    end
    
    ILoc = zeros([nonzeros, 1]);
    JLoc = zeros([nonzeros, 1]);
    Val  = zeros([nonzeros, 1]);
    
    JLoc(1:nonzeros, 1) = iVar; 
     
    counter = 0;
    for iWell = 1 : well.nWells
        if well_Jacobian(iWell) ~= 0 
            counter = counter + 1;
            ILoc(counter) = nreqT + iWell; 
            Val(counter) = well_Jacobian(iWell); 
        end
    end
    
end
