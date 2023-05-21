%% Volume Translation (V shift) calculation
%Developed by Ernesto Valbuena at Texas A&M University, College Station. 2013
%Class project for PETE 685 Advanced Phase Behavior II (original code in Excel VBA)

function [ MolarVol_Shift, MolarVolliq_Shift, MolarVolvap_Shift, MolarDensHC_Shift, ...
           MolarDensliq_Shift, MolarDensvap_Shift, zL_Shift, zV_Shift, oDen_Shift, gDen_Shift] ...
         = Y_Volume_Translation(P, T, zfrac, xo, yv, MolarVol, MolarVolliq, ...
           MolarVolvap, MWliq, MWvap, comp_Pc, comp_Tc, comp_SSHIFT)

    Rg = 10.732;
    
    ci = comp_SSHIFT .* (0.0778 .* Rg .* comp_Tc) ./ comp_Pc;  % Using Peng-Robinson EOS
    
    c  = sum(ci .* zfrac);    % Total
    cL = sum(ci .* xo);       % Liquid
    cV = sum(ci .* yv);       % Vapor

    MolarVol_Shift    = MolarVol    - c;
    MolarVolliq_Shift = MolarVolliq - cL;
    MolarVolvap_Shift = MolarVolvap - cV;

    zL_Shift = P .* MolarVolliq_Shift ./ (Rg .* T);
    zV_Shift = P .* MolarVolvap_Shift ./ (Rg .* T);

    MolarDensHC_Shift = 1 ./ MolarVol_Shift;   
    MolarDensHC_Shift(isnan(MolarDensHC_Shift)) = 0;  
    MolarDensHC_Shift(isinf(MolarDensHC_Shift)) = 0;

    if zL_Shift == 0
        
        oDen_Shift = 0;
        gDen_Shift = P .* MWvap ./ (Rg .* T .* zV_Shift);
        MolarDensliq_Shift = 0;
        MolarDensvap_Shift = 1 ./ MolarVolvap_Shift;
        
    elseif zV_Shift == 0 
        
        oDen_Shift = P .* MWliq ./ (Rg .* T .* zL_Shift);
        gDen_Shift = 0;
        MolarDensliq_Shift = 1 ./ MolarVolliq_Shift;
        MolarDensvap_Shift = 0;
        
    else
        
        oDen_Shift = P .* MWliq ./ (Rg .* T .* zL_Shift);
        gDen_Shift = P .* MWvap ./ (Rg .* T .* zV_Shift);
        MolarDensliq_Shift = 1 ./ MolarVolliq_Shift;
        MolarDensvap_Shift = 1 ./ MolarVolvap_Shift;
        
    end
    
end