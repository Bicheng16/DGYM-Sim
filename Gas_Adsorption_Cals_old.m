function [prop] = Gas_Adsorption_Cals(res, init, stat, prop)


    MolarVolVap_sc = 379.342;         % Gas molar volume in standard condition, scf/lb-mol. 
    
    % sorption_flag: Flag to consider adsorption
    ns = 0;
    for iRT = 1 : res.nRockTypes
        if ~ strcmp(res.Sorption{iRT}, 'NA')
            ns = ns + 1;
        end
    end
    
    sorption_flag = (ns > 0);         % At least one rock type consider adsorption 
    
    % If adsorption is not allowed
    if ~sorption_flag
        return;
    end
    
    % If no gas phase exist, return
    if sum(prop(stat).cell_Yi(:)) == 0
        return;
    end
    
    % If adsorption is considered 
    [nCells, nComps] = size(prop(stat).cell_Yi);
    cell_VLi = zeros(nCells, nComps);
    cell_PLi = zeros(nCells, nComps);
    cell_RMDen = zeros(nCells, 1);
    
    cell_VLi(:, :) = res.VLi(init.cell_type(1 : nCells), 1:nComps);
    cell_PLi(:, :) = res.PLi(init.cell_type(1 : nCells), 1:nComps);
    cell_RMDen(:, 1) = res.RMDen(init.cell_type(1 : nCells));
    
    %% Extend Langmuir Model Calculation
    % Calculate tmp1 = Yi * P / PL_i
    p_rep = repmat(prop(stat).cell_P, 1, nComps);
    tmp1 = prop(stat).cell_Yi .* p_rep ./ cell_PLi;
    
    % Calculate tmp2 = 1.0 + sum(Yi * P / PL_i, i = 1, ..., nc)
    tmp2 = 1.0 + sum(tmp1, 2);
    
    % Calculate tmp4 = (1.0 - poro) * RMDen / Vgs / [ 1.0 + sum(Yi * P / PL_i, i = 1, ..., nc) ]
    tmp3 = (1.0 - prop(stat).cell_poro) .* cell_RMDen ./ MolarVolVap_sc ./ tmp2;
    tmp4 = repmat(tmp3, 1, nComps);
    
    % Calculate tmp5 = VLi * (Yi * P / PL_i) 
    tmp5 = cell_VLi .* tmp1;   
    
    % Calculate VLi * { (Yi * P / PL_i) / (1.0 + sum(Yi * P / PL_i, i = 1, ..., nc) }
    prop(stat).cell_Qai(:, :) = tmp5 .* tmp4;  
    
end