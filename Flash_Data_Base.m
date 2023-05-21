function [fluid] = Flash_Data_Base(fluid, res, init)

    % If fluid is dry gas, return and no table provided
    fluidtype = fluid.fluidtype;
    if isequal(fluidtype, 'GAS')
        return;
    end

    % Temporary variables
    nComps = fluid.nComps;
    FtoR = 459.598;
    Tol = 1.0e-12;
    Pb = init.res_vle_chk(1, nComps+1);    
    Pstr = 50;                  % Lower limit of pressure range
    Pend = floor(Pb);           % Upper limit of pressure range: 
                                % the nearest integer <= Pb
    zfrac = res.Zi(1, :);
    T = res.T(1);
                                              
    % Input properties from fluid 
    nComps = fluid.nComps;
    comp_MW = fluid.comp_MW;
    comp_Pc = fluid.comp_Pc;
    comp_Tc = fluid.comp_Tc;
    comp_ACF = fluid.comp_ACF;
    comp_BIC = fluid.comp_BIC;
                                             
    %  FF_Variable_Definitions
    T =  T + FtoR;
    comp_Tc = comp_Tc + FtoR;
    illed_zfrac = abs((sum(zfrac)-1)) > Tol; 
    if illed_zfrac
       disp('...WARNING: Fluid overall composition does not at to 1 for at least one cell during the run') 
    end

    % Assume pressure step size is 5 psia 
    n = (Pend - Pstr) / 5.0d0;   
    n = floor(n);

    Flash_Table = zeros(n, nComps+2);
    delta_P = (Pend - Pstr) / n;
    kflag = 1;
    kval_0 = init.res_vle_chk(1, 1:nComps);
    fg_0 = 1.0e-5;

    for i = 1 : n
        j = n - i + 1;
        P = Pend - (i - 1) * delta_P;

        P = (P < Pb) * P + (P >= Pb) * (Pb - 0.5); 
         
        [~, fg, ~, ~, kval, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] ...         
        = Y_VLEflash(nComps, T, P, zfrac, comp_Tc, comp_Pc, comp_ACF, ...
        comp_MW, comp_BIC, kval_0, fg_0, kflag, Pb, fluidtype);

        Flash_Table(j, 1) = P;
        Flash_Table(j, 2) = fg;
        Flash_Table(j, 3:end) = kval;
        
        fg_0 = fg;
        kval_0 = kval; 
        
    end

    fluid.Flash_Table = Flash_Table;

end