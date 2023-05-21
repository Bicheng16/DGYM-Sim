function [ fg, xo, yv, kval, visc_oil, visc_gas, MolarDensHC_Shift, ...
           MolarDensliq_Shift, MolarDensvap_Shift, oDen_Shift,      ...
           gDen_Shift, MolarVolliq_Shift, MolarVolvap_Shift, MWliq, ...
           MWvap, zL_Shift, zV_Shift] = Y_FluidFlash(P, zfrac, T,   ...
           fluid, kval_0, fg_0, kflag, Pb)
           
    % Temporary variables
    FtoR = 459.598;
    Tol = 1.0e-12;
                                                  
    % Input properties from fluid 
    nComps = fluid.nComps;
    comp_MW = fluid.comp_MW;
    comp_Pc = fluid.comp_Pc;
    comp_Tc = fluid.comp_Tc;
    comp_ACF = fluid.comp_ACF;
    comp_Zc = fluid.comp_Zc;
    comp_BIC = fluid.comp_BIC;
    comp_SSHIFT = fluid.comp_SSHIFT;
    LBC_Coeffs = fluid.LBC_Coeffs;
    fluidtype = fluid.fluidtype;
                                             
    %  FF_Variable_Definitions
    T =  T + FtoR;
    comp_Tc = comp_Tc + FtoR;
    illed_zfrac = abs((sum(zfrac)-1)) > Tol; 
    if illed_zfrac
       disp('...WARNING: Fluid overall composition does not at to 1 for at least one cell during the run') 
    end
    
    %% Flash Calculation    
    [fo, fg, xo, yv, kval, oDen, gDen, MWliq, MWvap, MolarVolliq, ...
     MolarVolvap, MolarVol, MolarDensliq, MolarDensvap, zL, zV, Flash_iter, SS_iter, NR_iter, RR_iter] ...
    = Y_VLEflash(nComps, T, P, zfrac, comp_Tc, comp_Pc, comp_ACF, ...
    comp_MW, comp_BIC, kval_0, fg_0, kflag, Pb, fluidtype);
    
    % disp( 'Flash converged !' )
    % disp( strcat( 'Total iteration for flash: ', num2str(Flash_iter) ) )   
    % disp( strcat( 'Iteration for successive method: ', num2str(SS_iter) ) )   
    % disp( strcat( 'Iteration for newton method: ', num2str(NR_iter) ) )
    % disp( strcat( 'Total iteration for Rachford-Rice: ', num2str(RR_iter) ) )       
    
    %% volume translation    
    [ MolarVol_Shift, MolarVolliq_Shift, MolarVolvap_Shift, MolarDensHC_Shift, ...
      MolarDensliq_Shift, MolarDensvap_Shift, zL_Shift, zV_Shift, oDen_Shift, gDen_Shift] ...
    = Y_Volume_Translation(P, T, zfrac, xo, yv, MolarVol, MolarVolliq, ...
      MolarVolvap, MWliq, MWvap, comp_Pc, comp_Tc, comp_SSHIFT);

    [visc_oil] = Y_Viscosity_Calculation(P, T, xo, comp_Zc, comp_Pc, comp_Tc, comp_MW, MWliq, MolarVolliq_Shift, LBC_Coeffs);
    [visc_gas] = Y_Viscosity_Calculation(P, T, yv, comp_Zc, comp_Pc, comp_Tc, comp_MW, MWvap, MolarVolvap_Shift, LBC_Coeffs);

end