% This function is to define different result variables for report 
function [rst] = Result_Var_Cals(stat, rst, fluid, prop, init)

    % Ncells: number of cells in the reservoir 
    % Ncomp: number of hydrocarbon component in the reservoir
    % Nwells: number of wells in the reservoir 

    nComps = fluid.nComps;
    vol_const = 5.61458144; % Convert volume from barrel to cu.ft: 1 barrel = vol_const cu.ft. 5.614581442311112;
    
    %% Calculate reservoir volume variables
    rst.HCPV = sum(prop(stat).cell_HCPV) / vol_const;                     % Hydrocarbon pore volume, RB
    rst.PV = sum(prop(stat).cell_PV) / vol_const;                         % Pore volume, RB
    rst.VWat = sum(prop(stat).cell_PV .* prop(stat).cell_Sw) / vol_const; % Water volume, RB 
    rst.VOil = sum(prop(stat).cell_PV .* prop(stat).cell_So) / vol_const; % Oil volume, RB 
    rst.VGas = sum(prop(stat).cell_PV .* prop(stat).cell_Sg) / vol_const; % Gas volume, RB 

    %% Calculate reservoir average variables    
    rst.P_res_HCPVweigthed = sum(prop(stat).cell_HCPV .* prop(stat).cell_P) / ( rst.HCPV * vol_const);
    rst.P_res_PVweigthed = sum(prop(stat).cell_PV .* prop(stat).cell_P) / ( rst.PV * vol_const); 
    rst.Sw_PVweigthed = rst.VWat / rst.PV; 
    rst.So_PVweigthed = rst.VOil / rst.PV; 
    rst.Sg_PVweigthed = rst.VGas / rst.PV;
        
    %% Store fluid in place 
    % Calc Zi for whole field to flash to SC(separator)
    dummy_cell_vol = repmat(init.cell_volume, 1, nComps);

    % Number of moles of each component in compressed storage at each cell, lb-mole. 
    cell_free_moles = prop(stat).cell_Ni .* dummy_cell_vol;   
    % Total number of moles of each component in the reservoir (compressed storage), lb-mole.
    res_free_moles = sum(cell_free_moles, 1);
    
    % Number of moles of each component in adsorbed storage at each cell, lb-mole.
    cell_ads_moles = prop(stat).cell_Qai .* dummy_cell_vol;
    % Total number of moles of each component in the reservoir (adsorbed storage), lb-mole.
    res_ads_moles = sum(cell_ads_moles, 1); 
    
    % Total amount of hydrocarbon in the reservoir (compressed + adsorbed)
    rst.COMPIP(1, 1:nComps) = res_free_moles(1, 1:nComps) + res_ads_moles(1, 1:nComps);
    
    % Number of moles for free HC in whole field, lb-mole.
    rst.nT_Field_Vector = sum(res_free_moles, 2); %sum( sum( cell_free_moles ) ); 
    % Overall composition from the whole reservoir
    rst.Zi_Field_vector = res_free_moles ./ rst.nT_Field_Vector; 
    
    % Surface flash calculation based on field overall composition   
    Zi_comp = rst.Zi_Field_vector;

    % Fake Pb, assume bubble point pressure is higher than surface condition pressure 
    Pb = fluid.Psc + 1000;
    kval_0 = zeros(1, nComps);
    fg_0 = 1.0d0;
    kflag = 0;
        
    [fv_SC, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, MolarVolliq_SC, MolarVolvap_SC, ~, ~, ~, ~] ...
    = Y_FluidFlash(fluid.Psc, Zi_comp, fluid.Tsc, fluid, kval_0, fg_0, kflag, Pb);
       
    % Oil in place, STB,  
    rst.FOIP = (1 - fv_SC) .* rst.nT_Field_Vector .* MolarVolliq_SC ./ vol_const; 
    
    % Free gas in place, MSCF
    rst.FGIP =  fv_SC  .* rst.nT_Field_Vector .* MolarVolvap_SC ./ 1000; 
    
    % Adsorbed gas in place, MSCF. 379.342 is gas molar volume in standard condition
    rst.FAGIP = sum(res_ads_moles) * 379.342 / 1000;  
    
    % Water in place, STB 
    % lbm of water per cell    
    cell_wat_moles = prop(stat).cell_Nw .* init.cell_volume; 
    % Total amount of water in place, lb-mol
    rst.COMPIP(1, 1+nComps) = sum(cell_wat_moles, 1);
    % Water in place, STB
    rst.FWIP  = sum( cell_wat_moles ) * fluid.wat_MW ./ fluid.MDen_w_sc ./ vol_const; 
     
end
    
