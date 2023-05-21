function [prop] = Init_Fluids(current, prop, fluid, init)

    % Define temporary variables
    nCells = init.nCells;
    nComps = fluid.nComps;
    
    Xc = zeros(nCells, 1);
    Yc = zeros(nCells, 1);
    cell_Bw_inv = zeros(nCells, 1);
    
    % Initialize buffer 
    HC_Vol_fract_oil = zeros(nCells, 1);

    % Step 1: hydrocarbon fluids properties 
    HC_sat_Flag = (prop(current).cell_Sw ~= 1);
    
    % Scanning all cells
    for i = 1 : nCells
        pressure = prop(current).cell_P(i);
        temperature = prop(current).cell_T(i);
        Zi_comp = prop(current).cell_Zi(i, :);
        
        % Get bubble point 
        Pb = init.res_vle_chk(i, nComps+1);
        
        % Get K-factor and fv guess
        if pressure > Pb
            fv_0   = 1.0d-5; 
            kval_0 = init.res_vle_chk(i, 1:nComps);
        elseif pressure <= Pb && isequal(fluid.fluidtype, 'GAS') 
            %fv_0   = prop(1).cell_fv(i);
            %kval_0 = prop(1).cell_Ki(i, :);
            fv_0   = 1.0d-5; 
            kval_0 = init.res_vle_chk(i, 1:nComps);            
        else
            % Table look-up to find K-factor and fv initial guess
            fv_0 = nakeinterp1(fluid.Flash_Table(:, 1), fluid.Flash_Table(:, 2), pressure);
            for iComp = 1 : fluid.nComps
                kval_0(1, iComp) = nakeinterp1(fluid.Flash_Table(:, 1), fluid.Flash_Table(:, iComp+2), pressure); 
            end
        end
        
        % Perform flash calculation
        [fv, xo, yv, kval, visc_o, visc_g, MolarDensHC, MolarDensliq, MolarDensvap, oDen, gDen, MolarVolliq, MolarVolvap, ~, ~, ~, ~] ...
        = Y_FluidFlash(pressure, Zi_comp, temperature, fluid, kval_0, fv_0, 1, Pb);
                
        % Step 1: hydrocarbon properties
        prop(current).cell_Ki(i, :)=                   kval * HC_sat_Flag(i, 1);      % K-factor
        prop(current).cell_fv(i, 1)=                     fv * HC_sat_Flag(i, 1);      % mol fraction
        prop(current).cell_Vsic_o(i, 1)=                 visc_o * HC_sat_Flag(i, 1);  % cP  
        prop(current).cell_Vsic_g(i, 1)=                 visc_g * HC_sat_Flag(i, 1);  % cP        
        prop(current).cell_MolarDensHC(i, 1)=   MolarDensHC * HC_sat_Flag(i, 1);      % lb-mole/ft^3
        prop(current).cell_MolarDensLiq(i, 1)= MolarDensliq * HC_sat_Flag(i, 1);      % lb-mole/ft^3
        prop(current).cell_MolarDensVap(i, 1)= MolarDensvap * HC_sat_Flag(i, 1);      % lb-mole/ft^3
        prop(current).cell_MolarVolLiq(i, 1)  = MolarVolliq * HC_sat_Flag(i, 1);      % ft^3/lb-mole
        prop(current).cell_MolarVolVap(i, 1)  = MolarVolvap * HC_sat_Flag(i, 1);      % ft^3/lb-mole               
        prop(current).cell_MDen_o(i, 1)=                 oDen * HC_sat_Flag(i, 1);    % lbm/ft^3
        prop(current).cell_MDen_g(i, 1)=                 gDen * HC_sat_Flag(i, 1);    % lbm/ft^3        
        prop(current).cell_Xi(i, :)=                     xo .* HC_sat_Flag(i, 1);     % mol fraction
        prop(current).cell_Yi(i, :)=                     yv .* HC_sat_Flag(i, 1);     % mol fraction
  
        % Volumetric fraction of oil assuming on hydrocarbon components present (actual saturation is calculated below)
        HC_Vol_fract_oil(i,1) = (1 - fv) .* MolarVolliq .* MolarDensHC .* HC_sat_Flag(i,1);     
    end    

    prop(current).cell_So = HC_Vol_fract_oil .* (1 - prop(current).cell_Sw);    
    prop(current).cell_Sg = 1 - prop(current).cell_Sw - prop(current).cell_So;
            
    % Water properties nakeinterp1  interp1q spline
    Xc = fluid.Cw .* ( prop(current).cell_P - fluid.Pref_w);
    cell_Bw_inv = 1.0 ./ fluid.Bw_ref .* (1 + Xc .* (1 + 0.5 * Xc) );
    % Water mass density: lb / cu.ft
    prop(current).cell_MDen_w  = fluid.MDen_w_sc .* cell_Bw_inv;
    % Water molar density: lb-mole / cu.ft
    prop(current).cell_MolarDensWat = prop(current).cell_MDen_w ./ fluid.wat_MW; 

    % Water viscosity, cp
    Yc = - fluid.Cvw .* ( prop(current).cell_P - fluid.Pref_w);
    prop(current).cell_Vsic_w = fluid.Vsic_w_ref ./ (1 + Yc .* (1 + 0.5 .* Yc) );     
    
end