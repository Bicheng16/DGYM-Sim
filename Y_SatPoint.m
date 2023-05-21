function [Psat, kval] = Y_SatPoint(zfrac, T, fluid, Psat_0, kval_0, Flag_Pb)
           
    
    % Temporary variables
    FtoR = 459.598;
    Tol = 1.0e-12;
    Error = 100;
    iter = 0;
    
    Tol_Newton_Switch = 1.0e-4;      
    Iter_Newton_Switch = 3;       
    Flag_Newton_Switch = false;
    Newton_fail = false;
                                                  
    % Input properties from fluid 
    nComps = fluid.nComps;
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

    %% Flag to saturation pressure
    if Flag_Pb
        fg = 0;
    else
        fg = 1;
    end 

    % Quality check: kval_0 and fg_0
    sum_K = sum(kval_0);
    trust_K  = ~( isnan(sum_K) || isinf(sum_K) ) && (sum_K > 1.0d-6); 
    trust_Psat = (Psat_0 >= 0) && (Psat_0 <= 1.0e5);
    
    if (trust_K && trust_Psat)
        % Directly get from input
        % kval_0 = kval_0;
        P = Psat_0;
    else
        % Wilson's correlation to initialize the k-values
        P = 3000.0d0;
        kval_0 = ( comp_Pc ./ P ) .* exp( 5.3727 .* ( 1 + comp_ACF ) .* ( 1 - comp_Tc ./ T ) );        
    end    
    lnKi = log(kval_0);        
    
    %% Sequential scheme first    
    while Error > Tol

        if ( ((Error < Tol_Newton_Switch && Error > Tol) || (iter >= Iter_Newton_Switch)) )  
            Flag_Newton_Switch = true;
            break;
        end
                
        % Phase mole fraction 
        xo_1 = zfrac;
        yv_1 = zfrac .* kval_0;
        
        % Call EOS
        [lnPhiL, ~, ~, dlnPhiL_dP] = Y_CubicDP(nComps, P, T, 1, xo_1, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
        [lnPhiV, ~, ~, dlnPhiV_dP] = Y_CubicDP(nComps, P, T, 0, yv_1, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
        
        % Update kval
        lnKi = lnPhiL - lnPhiV;
        kval_1 = exp(lnKi);
        
        % Get scalar residual and jacobian
        R = sum(zfrac .* kval_1) - 1.0d0;
        dRdP = sum(zfrac .* kval_1 .* (dlnPhiL_dP - dlnPhiV_dP) );
        
        % Update P
        P = P - R / dRdP; 
        
        P_flag = ( ~isnan(P) && ~isnan(P) && P > 0 );  
        if ~P_flag
            Psat = Psat_0; 
            kval = kval_0;
            return;
        end
        
        % 
        kval_0 = kval_1;
        Error = abs(R);
        
        iter = iter + 1; 
               
    end % End sequential loop

    if Flag_Newton_Switch
        
        x_vec = zeros(nComps+1, 1); 
        Error = 100; 
        
        P_0 = 0;
        lnKi_0 = zeros(1, nComps);
        
        P_0 = P;
        lnKi_0(:) = lnKi(:);
                
        %% Fully implicit scheme
        while Error > Tol

            % Continue iteration
            iter = iter + 1;

            % Construct primary variable vector [lnKi, Pb]
            x_vec(1:nComps, 1) = lnKi_0(1, 1:nComps);
            x_vec(nComps+1, 1) = P;

            % Construct RHS and Jacobian
            [RHS, Jacobian] = Newton_SaturationPressure(fg, lnKi_0, zfrac, nComps, T, P, comp_BIC, comp_ACF, comp_Tc, comp_Pc);

            % Check Jacobian
            condJ = cond(Jacobian); 
            Newton_fail = ( condJ < 1.0e-3 || condJ > 1.0e13);                        
            if Newton_fail
                break
            end             
            
            % Solve linear system
            x_vec = x_vec - Jacobian \ RHS;
                                   
            % Construct primary variable vector
            lnKi(1, 1:nComps) = x_vec(1:nComps, 1);
            P = x_vec(nComps+1, 1);
            
            % If nonphysical pressure
            Newton_fail = ( P < 0 || isnan(P) || isinf(P) );
            if Newton_fail
                P = P_0;
                lnKi(:) = lnKi_0(:);
                break
            else
               P_0 = P;
               lnKi_0(:) = lnKi(:);
            end

            % Clean x_vec and Jacobian
            x_vec(:) = 0;

            % Error calculation
            Error = norm(RHS, 2);

        end % Newton loop
        
    end
    
    % If newton failed, go back to sequential again
    if Newton_fail
        Error = 100.0d0;
        while Error > Tol

            % Phase mole fraction 
            xo_1 = zfrac;
            yv_1 = zfrac .* kval_0;

            % Call EOS
            [lnPhiL, ~, ~, dlnPhiL_dP] = Y_CubicDP(nComps, P, T, 1, xo_1, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
            [lnPhiV, ~, ~, dlnPhiV_dP] = Y_CubicDP(nComps, P, T, 0, yv_1, comp_ACF, comp_Tc, comp_Pc, comp_BIC);

            % Update kval
            lnKi = lnPhiL - lnPhiV;
            kval_1 = exp(lnKi);

            % Get scalar residual and jacobian
            R = sum(zfrac .* kval_1) - 1.0d0;
            dRdP = sum(zfrac .* kval_1 .* (dlnPhiL_dP - dlnPhiV_dP) );

            % Update P
            P = P - R / dRdP; 

            % 
            kval_0 = kval_1;
            Error = abs(R);

            iter = iter + 1; 

        end % End sequential loop        
        
    end
    
    % Recalculate K-factor
    kval = exp(lnKi);
    Psat = P;
        
end

function [RHS, Jacobian] = Newton_SaturationPressure(fg, comp_lnKi, zfrac, nComps, T, P, comp_BIC, comp_ACF, comp_Tc, comp_Pc)

    % Define temporary variables
    kval = exp(comp_lnKi);

    % Allocation memory
    RHS = zeros(1+nComps, 1);
    Jacobian = zeros(1+nComps, 1+nComps);

    if fg == 0
        xo = zfrac;
        yv = xo .* kval;  
        % Reduced Rachford-Rice at bubble point
        RR = sum(zfrac .* kval) - 1.0;
    end
    if fg == 1
        yv = zfrac;
        xo = yv ./ kval;
        % Reduced Rachford-Rice at dew-point
        RR = 1.0 - sum(zfrac ./ kval);
    end

    % Calculate fugacities and derivatives
    [lnPhiL, ~, dlnPhiL_dxoj, dlnPhiL_dP] = Y_CubicDP(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
    [lnPhiV, ~, dlnPhiV_dyvj, dlnPhiV_dP] = Y_CubicDP(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
    
    % Thermodynamic equilibrium
    Rif = comp_lnKi + lnPhiV - lnPhiL;
    
    % RHS
    RHS(1:nComps, 1) = Rif(1, 1:nComps);
    RHS(1+nComps, 1) = RR;
    
    %% Construct Jacobian
    
    % Fugacity equilibrium derivatives to lnKi
    temp_1 = 1 + fg .* (kval - 1);
    temp_1 = temp_1 .* temp_1;
    
    dyvj_dlnkj = zfrac .* kval .* (1 - fg) ./ temp_1; 
    dyvj_dlnkj = repmat(dyvj_dlnkj, nComps, 1);
    dlnPhiV_dlnkj = dlnPhiV_dyvj .* dyvj_dlnkj; 

    dxoj_dlnkj = - zfrac .* kval .* fg ./ temp_1;
    dxoj_dlnkj = repmat(dxoj_dlnkj, nComps, 1);
    dlnPhiL_dlnkj = dlnPhiL_dxoj .* dxoj_dlnkj; 

    dRif_dlnkj = eye(nComps) + dlnPhiV_dlnkj - dlnPhiL_dlnkj;
          
    % Fugacity equilibrium derivatives to Pressure
    temp2   = dlnPhiV_dP - dlnPhiL_dP; 
    dRif_dP = temp2;  
        
    % Rachford-Rice derivatives to lnKi
    dRR_dlnkj = kval .* zfrac ./ temp_1; 
    
    % Rachford-Rice derivate to Pressure: zero
    dRR_dP = 0;
            
    % Jacobian
    Jacobian(1:nComps, 1:nComps) = dRif_dlnkj;
    Jacobian(1:nComps, 1+nComps) = dRif_dP(1:nComps); 
    Jacobian(1+nComps, 1:nComps) = dRR_dlnkj(1:nComps);
    Jacobian(1+nComps, 1+nComps) = dRR_dP;    
    
end