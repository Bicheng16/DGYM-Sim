function  [fo, fg, xo, yv, kval, oDen, gDen, MWliq, MWvap, MolarVolliq, ...
     MolarVolvap, MolarVol, MolarDensliq, MolarDensvap, zL, zV, Flash_iter, SS_iter, NR_iter, RR_iter] ...
    = Y_VLEflash(nComps, T, P, zfrac, comp_Tc, comp_Pc, comp_ACF, ...
    comp_MW, comp_BIC, kval_0, fg_0, kflag, Pb, fluidtype)
    
    % Define temporary variables
    Rg = 10.732;
    Flash_iter = 0; % Iteration for Two phase flash (successive + newton)
    SS_iter = 0;    % Iteration for successive method
    NR_iter = 0;    % Iteration for newton method
    RR_iter = 0;    % Iteration for Rachford-Rice Equation solve
    
    % Gas only phase properties
    if ( isequal(fluidtype, 'GAS') ) 

        xo = zeros(size(zfrac));
        yv = zeros(size(zfrac));
        kval = zeros(size(zfrac)); 
        
        fo = 0;
        fg = 1;

        xo(:) = 0;
        yv = zfrac;
        kval(:) = inf;

        MWliq = 0;
        MWvap = sum(yv .* comp_MW);

        zL = 0;        
        [~, zV] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);

        oDen = 0;                           % Oil mass density in cu.ft/lb      
        gDen = P * MWvap / (zV * Rg * T);   % Gas mass density in cu.ft/lb

        MolarVolliq = 0;                    % Oil molar volumes in cu.ft/lb-mole
        MolarVolvap = MWvap / gDen;         % Gas molar volumes in cu.ft/lb-mole        
        MolarVol = MolarVolvap;             % Hydrocarbon molar volumes in cu.ft/lb-mole        

        MolarDensliq = 0;                   % Oil molar densities in lb-mole/cu.ft        
        MolarDensvap = 1 / MolarVolvap;     % Gas molar densities in lb-mole/cu.ft

        SS_iter = 0;
        Flash_iter = 0;

        return;

    end    
       
    % VLE Fluid flash part
    Tol = 1.0e-9;                % Recommend: 1.0e-6; 1.0e-9;
    Tol_Newton_Switch = 1.0e-2;  % Recommend: 1.0e-2
    
    Tol_f = 1.0e-20;              % Single phase judge based on Tol_fg    
    Iter_Newton_Switch = 2;       % Recommend: 2
    Flag_Newton_Switch = false;
    Newton_fail = false;
    Flag_DEM = false; % Flag to use dominant eigenvalue method (DEM), after every 5 SS insert DEM update
    
    iter_max = 100;
    iter = 0;
    Error = 100;
    fg_fake = 500;
        
    % Quality check: kval_0 and fg_0
    sum_K = sum(kval_0);
    trust_K  = ~( isnan(sum_K) || isinf(sum_K) ) && (sum_K > 1.0d-6); 
    trust_fg = (fg_0 >= 0) && (fg_0 <= 1);  
 
    if (trust_K && trust_fg && kflag)
        % Directly get from input
        kval = kval_0;
    else
        % Wilson's correlation to initialize the k-values
        kval = ( comp_Pc ./ P ) .* exp( 5.3727 .* ( 1 + comp_ACF ) .* ( 1 - comp_Tc ./ T ) );
        fg_0 = fg_fake;        
    end

    while ( Error > Tol ) && ( iter <= iter_max ) && (~Flag_Newton_Switch)
                
        if ( ((Error < Tol_Newton_Switch && Error > Tol) || (iter >= Iter_Newton_Switch)) ) ...
            && P < Pb && (~Oil_only) && (~Gas_only)    
            Flag_Newton_Switch = true;
            break;
        end
                
        [fg, xo, yv, Iter_rr] = Y_SS_RachfordRice(zfrac, kval, fg_0);
        RR_iter = RR_iter + Iter_rr; 
        
        fg(fg < 0) = 0;
        
        fg_0 = fg;
        
        fo = 1 - fg;

        % Oil only
        Oil_only  = (abs(fg) < Tol_f);
        
        % Extreme case: P < Pb but calculated fg == 0, switch to newton
        if Oil_only && P < Pb
            Flag_Newton_Switch = true;
            break;            
        end

        if Oil_only  
            
            fo = 1;
            fg = 0;

            xo = zfrac;
            yv(:) = 0;
            kval(:) = 0;

            MWliq = sum(xo .* comp_MW);
            MWvap = 0;

            [~, zL] = Y_Cubic(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);           
            
            zV = 0;

            oDen = P * MWliq / (zL * Rg * T);   % Oil mass density in cu.ft/lb        
            gDen = 0;                           % Gas mass density in cu.ft/lb  

            MolarVolliq = MWliq / oDen;         % Oil molar volumes in cu.ft/lb-mole
            MolarVolvap=0;                      % Gas molar volumes in cu.ft/lb-mole        
            MolarVol = MolarVolliq;             % Hydrocarbon molar volumes in cu.ft/lb-mole        

            MolarDensliq = 1 / MolarVolliq;     % Oil molar densities in lb-mole/cu.ft        
            MolarDensvap = 0 ;                  % Gas molar densities in lb-mole/cu.ft
            
            SS_iter = iter;
            Flash_iter = iter; 

            return;
            
        end
    
        % Gas only
        Gas_only = (abs(fo) < Tol_f);
        if ( Gas_only && P < Pb ) || ( isequal(fluidtype, 'GAS') ) 
            
            fo = 0;
            fg = 1;

            xo(:) = 0;
            yv = zfrac;
            kval(:) = inf;

            MWliq = 0;
            MWvap = sum(yv .* comp_MW);

            zL = 0;        
            [~, zV] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);

            oDen = 0;                           % Oil mass density in cu.ft/lb      
            gDen = P * MWvap / (zV * Rg * T);   % Gas mass density in cu.ft/lb

            MolarVolliq = 0;                    % Oil molar volumes in cu.ft/lb-mole
            MolarVolvap = MWvap / gDen;         % Gas molar volumes in cu.ft/lb-mole        
            MolarVol = MolarVolvap;             % Hydrocarbon molar volumes in cu.ft/lb-mole        

            MolarDensliq = 0;                   % Oil molar densities in lb-mole/cu.ft        
            MolarDensvap = 1 / MolarVolvap;     % Gas molar densities in lb-mole/cu.ft
            
            SS_iter = iter;
            Flash_iter = iter;

            return;
            
        end

        [phiL, ~] = Y_Cubic(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);        
        [phiV, ~] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
        
        % Fugacities
        fiV = phiV .* yv .* P;
        fiL = phiL .* xo .* P;
        
        % Get new K-value
        if Flag_DEM && iter > 1 && mod(iter, 5) == 1
            
            delta = log(fiL ./ fiV);
            b01 = delta * delta_0';
            b11 = norm(delta_0, 2);
            miu1 = - b01 / (b11 * b11);
                                    
            lnKi_1 = log(kval) + delta ./ (1 + miu1); 
            kval_1 = exp(lnKi_1); 

            delta_0 = delta;
                        
        else
            delta_0 = log(fiL ./ fiV);
            kval_1 = kval .*  (fiL ./ fiV);
        end
        
        Error = sum( (fiL ./ fiV - 1.0).^2 );
        % Error = max( abs( (kval_1 ./ kval - 1.0) ) );
        
        kval = kval_1;
        
        iter = iter + 1; 
        
        % When only 1 iteration reach convergence, force it for another
        % loop
        if iter == 1 && Error < Tol
            Error = 100.0;
        end
                
    end
    
    % Count iterations for successive substitution method
    SS_iter = iter;
    
    %% Solve fugacity equilibrium and Rachford Rich together through newton method
    if Flag_Newton_Switch
        
        Error = 100;
        % Define RHS
        x_vec = zeros(nComps+1, 1);

        % Change newton primary variable to lnKi
        lnKi = log(kval); 
                
        while Error > Tol && iter <= iter_max 
                                
            % Construct primary variable vector
            x_vec(1:nComps, 1) = lnKi(1, 1:nComps);
            x_vec(nComps+1, 1) = fg;
                    
            %% Construct RHS and Jacobian 
            [RHS, Jacobian] = Newton_TwoPhaseFlash(fg, lnKi, zfrac, nComps, T, P, comp_BIC, comp_ACF, comp_Tc, comp_Pc);

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
            fg = x_vec(nComps+1, 1);
                       
            % Error calculation
            Error = norm(RHS, 2);
                        
            % Clean x_vec and Jacobian
            x_vec(:) = 0;
                        
            % Continue iteration
            iter = iter + 1;
                    
        end  % Newton loop
    
        % Recalculate K-factor
        kval = exp(lnKi);
        fg = 0 * (fg <= 0) + 1 * (fg >= 1) + fg * (fg > 0 & fg < 1 );
        
    end
    
    % If newton doesn't work go back SS again
    if Newton_fail
        
        Error = 100.0d0;
        
        while ( Error > Tol ) && ( iter <= iter_max )

            [fg, xo, yv, Iter_rr] = Y_SS_RachfordRice(zfrac, kval, fg_0);
            RR_iter = RR_iter + Iter_rr; 

            fg(fg < 0) = 0;

            fg_0 = fg;

            fo = 1 - fg;

            % Oil only
            Oil_only  = (abs(fg) < Tol_f);
            if Oil_only 

                fo = 1;
                fg = 0;

                xo = zfrac;
                yv(:) = 0;
                kval(:) = 0;

                MWliq = sum(xo .* comp_MW);
                MWvap = 0;

                [~, zL] = Y_Cubic(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);           

                zV = 0;

                oDen = P * MWliq / (zL * Rg * T);   % Oil mass density in cu.ft/lb        
                gDen = 0;                           % Gas mass density in cu.ft/lb  

                MolarVolliq = MWliq / oDen;         % Oil molar volumes in cu.ft/lb-mole
                MolarVolvap=0;                      % Gas molar volumes in cu.ft/lb-mole        
                MolarVol = MolarVolliq;             % Hydrocarbon molar volumes in cu.ft/lb-mole        

                MolarDensliq = 1 / MolarVolliq;     % Oil molar densities in lb-mole/cu.ft        
                MolarDensvap = 0 ;                  % Gas molar densities in lb-mole/cu.ft

                SS_iter = iter;
                Flash_iter = iter; 

                return;

            end

            % Gas only
            Gas_only = (abs(fo) < Tol_f);
            if ( Gas_only && P < Pb ) || ( isequal(fluidtype, 'GAS') ) 

                fo = 0;
                fg = 1;

                xo(:) = 0;
                yv = zfrac;
                kval(:) = inf;

                MWliq = 0;
                MWvap = sum(yv .* comp_MW);

                zL = 0;        
                [~, zV] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);

                oDen = 0;                           % Oil mass density in cu.ft/lb      
                gDen = P * MWvap / (zV * Rg * T);   % Gas mass density in cu.ft/lb

                MolarVolliq = 0;                    % Oil molar volumes in cu.ft/lb-mole
                MolarVolvap = MWvap / gDen;         % Gas molar volumes in cu.ft/lb-mole        
                MolarVol = MolarVolvap;             % Hydrocarbon molar volumes in cu.ft/lb-mole        

                MolarDensliq = 0;                   % Oil molar densities in lb-mole/cu.ft        
                MolarDensvap = 1 / MolarVolvap;     % Gas molar densities in lb-mole/cu.ft

                SS_iter = iter;
                Flash_iter = iter;

                return;

            end

            [phiL, ~] = Y_Cubic(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);        
            [phiV, ~] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);

            % Fugacities
            fiV = phiV .* yv .* P;
            fiL = phiL .* xo .* P;

            % Get new K-value
            if Flag_DEM && iter > 1 && mod(iter, 5) == 1

                delta = log(fiL ./ fiV);
                b01 = delta * delta_0';
                b11 = norm(delta_0, 2);
                miu1 = - b01 / (b11 * b11);

                lnKi_1 = log(kval) + delta ./ (1 + miu1); 
                kval_1 = exp(lnKi_1); 

                delta_0 = delta;

            else
                delta_0 = log(fiL ./ fiV);
                kval_1 = kval .*  (fiL ./ fiV);
            end

            %Error = sum( (fiL ./ fiV - 1.0).^2 );
            Error = max( abs( (kval_1 ./ kval - 1.0) ) );

            kval = kval_1;

            iter = iter + 1;        

        end
        % Count Newton iterations
        Flash_iter = iter;
        NR_iter = 0;    
    else    
        % Count Newton iterations
        Flash_iter = iter;
        NR_iter = Flash_iter - SS_iter;
    end
    
    fo = 1.0 - fg;
    [xo, yv] = Y_renorm(fg, kval, zfrac);    

    MWliq = sum(xo .* comp_MW);
    MWvap = sum(yv .* comp_MW);
        
    [~, zL] = Y_Cubic(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);        
    [~, zV] = Y_Cubic(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
        
    oDen = P * MWliq / (zL * Rg * T);
    gDen = P * MWvap / (zV * Rg * T);
    
    MolarVolliq = MWliq / oDen;         % Oil molar volumes in cu.ft/lb-mole    
    MolarVolvap = MWvap / gDen;         % Gas molar volumes in cu.ft/lb-mole
    MolarVol = fo * MolarVolliq + fg * MolarVolvap; % Hydrocarbon molar volumes in cu.ft/lb-mole
        
    MolarDensliq = 1 / MolarVolliq;     % Oil molar densities in lb-mole/cu.ft
    MolarDensvap = 1 / MolarVolvap;     % Gas molar densities in lb-mole/cu.ft
          
end

function [RHS, Jacobian]  = Newton_TwoPhaseFlash(fg, comp_lnKi, zfrac, nComps, T, P, comp_BIC, comp_ACF, comp_Tc, comp_Pc)

    % Define temporary variables
    kval = exp(comp_lnKi);
    
    % Allocation memory
    RHS = zeros(1+nComps, 1);
    Jacobian = zeros(1+nComps, 1+nComps);

    %% Construct RHS
    % Fugacity equilibrium
    [xo, yv] = Y_renorm(fg, kval, zfrac);

    [lnPhiL, ~, dlnPhiL_dxoj] = Y_CubicD(nComps, P, T, 1, xo, comp_ACF, comp_Tc, comp_Pc, comp_BIC);        
    [lnPhiV, ~, dlnPhiV_dyvj] = Y_CubicD(nComps, P, T, 0, yv, comp_ACF, comp_Tc, comp_Pc, comp_BIC);
        
    % Thermodynamic equilibrium
    Rif = comp_lnKi + lnPhiV - lnPhiL; 

    % Rachford-Rice Equation
    temp_0 = zfrac .* (kval - 1) ./ (1 + fg .* (kval - 1) );
    RR = sum(temp_0, 2);
    
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
    
    % Fugacity equilibrium derivatives to fg
    dxoj_dfg = - zfrac .* (kval - 1) ./ temp_1;
    dyvj_dfg = dxoj_dfg .* kval;
    
    dxoj_dfg = repmat(dxoj_dfg, nComps, 1);
    dyvj_dfg = repmat(dyvj_dfg, nComps, 1); 
    
    dlnPhiV_dfg = dlnPhiV_dyvj .* dyvj_dfg; 
    dlnPhiV_dfg = sum(dlnPhiV_dfg, 2);
    
    dlnPhiL_dfg = dlnPhiL_dxoj .* dxoj_dfg;
    dlnPhiL_dfg = sum(dlnPhiL_dfg, 2);
    
    dRif_dfg = dlnPhiV_dfg - dlnPhiL_dfg;  
    
    % Rachford-Rice derivatives to lnKi
    dRR_dlnkj = kval .* zfrac ./ temp_1; 
    
    % Rachford-Rice derivate to fg
    dRR_dfg = - zfrac .* (kval - 1) .* (kval - 1) ./ temp_1;
    dRR_dfg = sum(dRR_dfg, 2);
            
    % Jacobian
    Jacobian(1:nComps, 1:nComps) = dRif_dlnkj;
    Jacobian(1:nComps, 1+nComps) = dRif_dfg(1:nComps); 
    Jacobian(1+nComps, 1:nComps) = dRR_dlnkj(1:nComps);
    Jacobian(1+nComps, 1+nComps) = dRR_dfg;

end
