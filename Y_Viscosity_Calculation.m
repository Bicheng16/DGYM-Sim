
%% Phase (oil/gas) viscosity calculation
% Developed by Ernesto Valbuena at Texas A&M University, College Station. 2013
% Class project for PETE 685 Advanced Phase Behavior II (original code in Excel VBA)
% Based on LBC SPE 915 paper (1964)

function [visc_phase] = Y_Viscosity_Calculation(P, T, xi, comp_Zc, comp_Pc, ...
          comp_Tc, comp_MW, MW_phase, MolarVol_phase, LBC_Coeffs)

    Rg = 10.732;

    Vcriti = comp_Zc .* Rg .* comp_Tc ./ comp_Pc;
    Vm_pc = sum(xi .* Vcriti);  %_pc stands for pseudo critical

    T_Lpc = sum(xi .* comp_Tc);
    p_Lpc = sum(xi .* comp_Pc); 
    
    %Visc_param_i is denoted by greek letter zeta subscript i in the LBC paper
    Visc_param_i = 5.4402 .* comp_Tc .^ (1/6) ./ ( (sqrt(comp_MW)) .* comp_Pc .^ (2/3) );
    
    %Reduced p and T for each component
    Tr_i = T ./ comp_Tc;
    % Pr_i = P ./ comp_Pc;
    
    %Reduced molar density    
    MolDen_red = Vm_pc ./ MolarVol_phase;
    MolDen_red( isnan(MolDen_red) | isinf(MolDen_red) ) = 0; 
        
    % Calc of low pressure viscosity of i'th component and mixture
    Dummy = ( Tr_i <= 1.5 );
    Visc_lowP_i =       Dummy  .* ( 0.00034   .* ( Tr_i .^ 0.94 ) ./ Visc_param_i ) ...
                + ( 1 - Dummy ).* ( 0.0001778 .* ( ( 4.58 .* Tr_i - 1.67 ) .^ 0.625 ) ./ Visc_param_i );
            
    Numer_phase = sum(xi .* Visc_lowP_i .* sqrt(comp_MW));
    Denom_phase = sum(xi .* sqrt(comp_MW));
    
    Visc_lowP = Numer_phase ./ Denom_phase;    
    Visc_lowP( isnan(Visc_lowP) | isinf(Visc_lowP) ) = 0;
 
    %Visc_param_mix Denoted in the LBC paper by greek letter zeta
    Visc_param_phase = 5.4402 .* T_Lpc.^(1/6) ./ ( sqrt(MW_phase) .* p_Lpc.^(2/3) );      
    Visc_param_phase( isnan(Visc_param_phase) | isinf(Visc_param_phase) ) = 0;
        
    % SumLBC_MolDen_red = LBC_Coeffs(1) + LBC_Coeffs(2) .* MolDen_red + LBC_Coeffs(3) .* MolDen_red.^2 + ...
    %                     LBC_Coeffs(4) .* MolDen_red.^3 + LBC_Coeffs(5).*MolDen_red.^4; 
    
    SumLBC_MolDen_red = LBC_Coeffs(1) + MolDen_red .* ( LBC_Coeffs(2) + ...
                        MolDen_red .* ( LBC_Coeffs(3) + MolDen_red .*   ...
                        ( LBC_Coeffs(4) + MolDen_red .* LBC_Coeffs(5) ) ) );  
    
    visc_phase = Visc_lowP + ( ( SumLBC_MolDen_red .^ 4 - 0.0001 ) ./ Visc_param_phase );       
    visc_phase( isnan(visc_phase) | isinf(visc_phase) ) = 0;
    
end





















