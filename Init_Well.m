function [well_ctr, prop, well] = Init_Well(well, prop, fluid, init)

    % Define temporary variables
    nWells = well.nWells;                     % Total number of wells
    nPerfs = sum(well.nperf_cell);            % Total number of perforations 
    nComps = fluid.nComps;                    % Total number of hydrocarbon components in reservoir fluids
    current = 1;                              % Status of reservoir cell properties

    % Memory preallocation     
    perf_P = zeros(nPerfs, 1);    
    Hw = zeros(nPerfs, 1);                    % Gravity potential between well reference depth and each perforation
    
    % Well mobility 
    WI_o = zeros(nPerfs, 1);
    WI_w = zeros(nPerfs, 1);
    WI_g = zeros(nPerfs, 1);
        
    % Output different well controls
    well_ctr = alloc_dynamic_well_ctr(nWells, well.wconne_well);
                                                       
    SMolarVolWat_sc = fluid.SMolarVolWat_sc; % Specific molar volume of water at surface condition
                                                      
    % Mobility for producers (Injectors modified individually)
    Oil_Mob = prop(current).cell_kro(well.wconne_cell) ./ prop(current).cell_Vsic_o(well.wconne_cell);         
    Gas_Mob = prop(current).cell_krg(well.wconne_cell) ./ prop(current).cell_Vsic_g(well.wconne_cell);        
    Wat_Mob = prop(current).cell_krw(well.wconne_cell) ./ prop(current).cell_Vsic_w(well.wconne_cell);        
    
    Oil_Mob(isinf(Oil_Mob) | isnan(Oil_Mob)) = 0.0d0;
    Gas_Mob(isinf(Gas_Mob) | isnan(Gas_Mob)) = 0.0d0;
    Wat_Mob(isinf(Wat_Mob) | isnan(Wat_Mob)) = 0.0d0;
             
    for iWell = 1 : nWells
        
        % If the well is inactive, ignore it
        if ~well.Active_Flag(iWell) 
            continue
        end
        
        % Check volume rate constraint  
        OilVolRate_sc = well.TargetOilRate(iWell); 
        OilVolRate_sc(isnan(OilVolRate_sc) | isinf(OilVolRate_sc)) = 0.0d0;
        OilWell_flag = ( OilVolRate_sc > 1.0e-3 && OilVolRate_sc < 1.0e+10 );

        GasVolRate_sc = well.TargetGasRate(iWell);
        GasVolRate_sc(isnan(GasVolRate_sc) | isinf(GasVolRate_sc)) = 0.0d0;
        GasWell_flag = ( GasVolRate_sc > 1.0e-3 && GasVolRate_sc < 1.0e+10 );

        LiquidVolRate_sc = well.TargetLiquidRate(iWell);
        LiquidVolRate_sc(isnan(LiquidVolRate_sc) | isinf(LiquidVolRate_sc)) = 0.0d0;
        LiqWell_flag = ( LiquidVolRate_sc > 1.0e-3 && LiquidVolRate_sc < 1.0e+10 );
        
        WatVolRate_sc = well.TargetWaterRate(iWell);
        WatVolRate_sc(isnan(WatVolRate_sc) | isinf(WatVolRate_sc)) = 0;                      
        WatWell_flag = ( WatVolRate_sc > 1.0e-3 && WatVolRate_sc < 1.0e+10 );
        
        % Check BHP constraint
        effective_BHP = ( well.Target_min_BHP(iWell) >= 14.7 && well.Target_max_BHP(iWell) <= 20000);
        
        active_well_flag = OilWell_flag + GasWell_flag + LiqWell_flag + WatWell_flag + effective_BHP;
        
        if active_well_flag == 0
            continue
        end
                         
        % Number of perforation in iWell
        active_perf = sum(well.wconne_well == iWell);       % Number of perforated cells connected to the iWell 
        
        Perf_Cells = zeros(active_perf, 1);                 % Perforation cell name connected to iWell
        cellP_Hw   = zeros(active_perf, 1);                 % Cell pressure - hydraulic pressure head, psia

        Qo_perf = zeros(nPerfs, 1);                         % Oil molar rate at each perforation 
        Qg_perf = zeros(nPerfs, 1);                         % Gas molar rate at each perforation 
        Qw_perf = zeros(nPerfs, 1);                         % Water molar rate at each perforation
                
        Pwf_effective = zeros(10, 1);
        nnz_p = 0;
        
        % Get perforation cell ID connected to iWell
        Perf_Cells(1 : active_perf) = well.wconne_cell(well.wconne_well == iWell);
        cellP_Hw(:) = prop(current).cell_P(Perf_Cells) - Hw(well.wconne_well == iWell);
                
        %%%%%%%%%%%%%%%%%%%%%  
        %%% Producers %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%
        if well.Prod_Flag(iWell) == 1 

            % Calculate well mobility for the well perforations (upstream: reservoir)
            WI_o(well.wconne_well == iWell, 1) = well.Prod_Flag_Perf(well.wconne_well == iWell, 1)       ...
            .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) .* Oil_Mob(well.wconne_well == iWell, 1) ...
            .* prop(current).cell_MolarDensLiq(Perf_Cells);
        
            WI_w(well.wconne_well == iWell, 1) = well.Prod_Flag_Perf(well.wconne_well == iWell, 1)       ...
            .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) .* Wat_Mob(well.wconne_well == iWell, 1) ...
            .* prop(current).cell_MolarDensWat(Perf_Cells);
        
            WI_g(well.wconne_well == iWell, 1) = well.Prod_Flag_Perf(well.wconne_well == iWell, 1)       ...
            .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) .* Gas_Mob(well.wconne_well == iWell, 1) ...
            .* prop(current).cell_MolarDensVap(Perf_Cells);
        
            % Eliminate values without physical meaning
            WI_o(isnan(WI_o) | isinf(WI_o)) = 0;
            WI_w(isnan(WI_w) | isinf(WI_w)) = 0;
            WI_g(isnan(WI_g) | isinf(WI_g)) = 0;
                       
            %% Check effective well control
            Pwf_o_rate = 0.0d0;
            Pwf_g_rate = 0.0d0;
            Pwf_ow_rate = 0.0d0;
                                   
            HCProd_flag = OilWell_flag ...       % Effective oil rate control
                        + GasWell_flag ...       % Effective gas rate control
                        + LiqWell_flag;          % Effective total liquid rate control
                    
            if HCProd_flag > 0 

                % Calculate well stream composition for the current producer
                % Mobility averaged if Option 1 failed
                dummy_WI_o = WI_o(:, ones(1, nComps)); 
                dummy_WI_g = WI_g(:, ones(1, nComps));                     
                ni = dummy_WI_o(well.wconne_well == iWell, :) .* prop(current).cell_Xi(Perf_Cells, :) ...
                   + dummy_WI_g(well.wconne_well == iWell, :) .* prop(current).cell_Yi(Perf_Cells, :);
                Zi_wellstream = sum(ni, 1) ./ sum(ni(:)); 
                Zi_wellstream(isnan(Zi_wellstream) | isinf(Zi_wellstream)) = 0.0d0;   
                
                % Update producer well stream
                well.Gas_Comp(iWell, :) = Zi_wellstream;
                
                % Surface flash calculation
                Pb = init.wel_vle_chk(iWell, nComps+1);
                kval_0 = init.wel_vle_chk(iWell, 1:nComps);
                fg_0 = 0.5d0;
                kflag = 1;
                
                [fv_ws_SC, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, MolarVolliq_ws_SC, MolarVolvap_ws_SC, ~, ~, ~, ~] ...
                = Y_FluidFlash(fluid.Psc, Zi_wellstream, fluid.Tsc, fluid, kval_0, fg_0, kflag, Pb);
                       
                % Pwf back calculated from surface oil volume rate constrain
                SMolarVolOil_sc = (1.0d0 - fv_ws_SC) * MolarVolliq_ws_SC;  % Specific molar volume of oil at surface condition
                SMolarVolGas_sc = fv_ws_SC * MolarVolvap_ws_SC;          % Specific molar volume of gas at surface condition
                               
                % HC molar rate back calculated from surface oil volume rate constrain
                if OilWell_flag
                    % Bottom hole pressure based on oil rate  
                    Pwf_o_rate = ( SMolarVolOil_sc * sum( WI_o(well.wconne_well == iWell) .* cellP_Hw + ... 
                                                          WI_g(well.wconne_well == iWell) .* ( cellP_Hw + prop(current).cell_pcgo(Perf_Cells) ) ) ...                          
                               -   OilVolRate_sc ) / ( SMolarVolOil_sc * sum( (WI_o + WI_g) .* (well.wconne_well == iWell) ) );
                    Pwf_o_rate(isnan(Pwf_o_rate) | isinf(Pwf_o_rate)) = 0.0d0;      
                    Pwf_o_rate = Pwf_o_rate * (Pwf_o_rate > 14.7d0);  

                    % If Pwf_rate physically make sense, add to Pwf_effective  
                    if Pwf_o_rate > 0.0d0 
                       nnz_p = nnz_p + 1;
                       Pwf_effective(nnz_p, 1) = Pwf_o_rate; 
                    end
                end 

                % HC molar rate back calculated from surface gas volume rate constrain                               
                if GasWell_flag
                    % Bottom hole pressure based on gas rate 
                    Pwf_g_rate = ( SMolarVolGas_sc * sum( WI_o(well.wconne_well == iWell) .* cellP_Hw + ... 
                                               WI_g(well.wconne_well == iWell) .* ( cellP_Hw + prop(current).cell_pcgo(Perf_Cells) ) ) ...                          
                               -   GasVolRate_sc ) / ( SMolarVolGas_sc * sum( (WI_o + WI_g) .* (well.wconne_well == iWell) ) );
                    Pwf_g_rate(isnan(Pwf_g_rate) | isinf(Pwf_g_rate)) = 0.0d0;      
                    Pwf_g_rate = Pwf_g_rate * (Pwf_g_rate > 14.7d0);                                

                    % If Pwf_rate physically make sense, add to Pwf_effective  
                    if Pwf_g_rate > 0.0d0
                       nnz_p = nnz_p + 1; 
                       Pwf_effective(nnz_p, 1) = Pwf_g_rate; 
                    end
                end
                
                % Pwf based on total liquid production rate                                 
                if LiqWell_flag
                    % Bottom hole pressure based on total liquid rate
                    Pwf_ow_rate = ( SMolarVolOil_sc * sum( WI_o(well.wconne_well == iWell) .* cellP_Hw + ... 
                                                WI_g(well.wconne_well == iWell) .* ( cellP_Hw + prop(current).cell_pcgo(Perf_Cells) ) ) ...
                                  + SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) .* ( cellP_Hw - prop(current).cell_pcow(Perf_Cells) ) ) ...
                                  - LiquidVolRate_sc ) / ( SMolarVolOil_sc * sum( (WI_o + WI_g) .* (well.wconne_well == iWell) ) + SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) ) );      

                    Pwf_ow_rate(isnan(Pwf_ow_rate) | isinf(Pwf_ow_rate)) = 0;      
                    Pwf_ow_rate = Pwf_ow_rate * (Pwf_ow_rate > 14.7);                                

                    % If Pwf_rate physically make sense, add to Pwf_effective  
                    if Pwf_ow_rate > 0 
                       nnz_p = nnz_p + 1; 
                       Pwf_effective(nnz_p, 1) = Pwf_ow_rate; 
                    end
                    
                end

            end

            % Water molar rate back calculated from surface water volume rate constrain
            Pwf_w_rate = 0;

            if WatWell_flag > 0                                                                
                % Pwf based on water mass rate 
                Pwf_w_rate = ( SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) .* ( cellP_Hw - prop(current).cell_pcow(Perf_Cells) ) ) ...
                           -   WatVolRate_sc) / ( SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) ) );
                Pwf_w_rate(isnan(Pwf_w_rate) | isinf(Pwf_w_rate)) = 0;      
                Pwf_w_rate = Pwf_w_rate * (Pwf_w_rate > 14.7);                                

                % If Pwf_rate physically make sense, add to Pwf_effective  
                if Pwf_w_rate > 0 
                   nnz_p = nnz_p + 1;
                   Pwf_effective(nnz_p, 1) = Pwf_w_rate; 
                end
                
            end

            % Determine the minimum effective BHP from rate control
            Pwf_rate_min = min(Pwf_effective(1:nnz_p, 1)) * (nnz_p > 0);
            Pwf_rate_min(isempty(Pwf_rate_min)) = 0;

            if Pwf_rate_min >= well.Target_min_BHP(iWell) 
                
                % One of the rate control is effective 
                prop(current).well_P(iWell) = Pwf_rate_min; 
                
                % Oil rate control is effective 
                well_ctr(current).OilRate_Ctr(iWell) = ( Pwf_o_rate == Pwf_rate_min );
                
                % Gas rate control is effective
                well_ctr(current).GasRate_Ctr(iWell) = ( Pwf_g_rate == Pwf_rate_min );

                % Total liquid rate is effective  
                well_ctr(current).LiqRate_Ctr(iWell) = ( Pwf_ow_rate == Pwf_rate_min );

                % Water rate control is effective
                well_ctr(current).WatRate_Ctr(iWell) = ( Pwf_w_rate == Pwf_rate_min );
                
            else
                % BHP control is effective       
                prop(current).well_P(iWell) = well.Target_min_BHP(iWell);
                well_ctr(current).BHP_Ctr(iWell) = true;
            end
            
            % Check the status of each perforation
            perf_P(well.wconne_well == iWell) = prop(current).well_P(iWell) + Hw(well.wconne_well == iWell);        

            % Oil molar rate of each perforation at reservoir condition
            Qo_perf(well.wconne_well == iWell) = WI_o(well.wconne_well == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) - perf_P(well.wconne_well == iWell));
            % Gas molar rate of each perforation at reservoir condition
            Qg_perf(well.wconne_well == iWell) = WI_g(well.wconne_well == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) + prop(current).cell_pcgo(Perf_Cells) - perf_P(well.wconne_well == iWell));
            % Water mass rate of each perforation at reservoir condition
            Qw_perf(well.wconne_well == iWell) = WI_w(well.wconne_well == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) - prop(current).cell_pcow(Perf_Cells) - perf_P(well.wconne_well == iWell));

            % Forbid cross flow if multi-perfs exist for the well
            perf_cross_flow = [find(Qo_perf<0); find(Qg_perf<0); find(Qw_perf<0)];
            well_ctr(current).wconne_well(perf_cross_flow) = 0;
                       
        end % Producer 
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Hydrocarbon Injector %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if well.Inje_Gas_Flag(iWell) == 1 
            
            % Define extra variables for hydrocarbon injector
            MolarDensLiq_inj = zeros(nPerfs, 1);    % Liquid molar density at each perforation to cal. injector mobility
            MolarDensVap_inj = zeros(nPerfs, 1);    % Vapor molar density at each perforation to cal. injector mobility
                                    
            % Get wellstream compositions 
            Zi_wellstream = well.Gas_Comp(iWell, :); 
            
            % Neglect gravity in wellbore first 
            perf_P(well.wconne_well == iWell) = prop(current).well_P(iWell) + Hw(well.wconne_well == iWell);
                        
            % Calculate fluid molar density @ each perforation            
            % Reservoir temperature            
            iperf_T = prop(current).cell_T(1);                                    
            for iPerf = 1 : nPerfs
                if well.wconne_well(iPerf) == iWell
                    iperf_P = perf_P(iPerf);
                    
                    % Calculate fluid density 
                    Pb = init.wel_vle_chk(iWell, nComps+1);
                    kval_0 = init.wel_vle_chk(iWell, 1:nComps);
                    fg_0 = 0.5d0;
                    kflag = 1;                                     
                    [~, ~, ~, ~, ~, ~, ~, MolarDensliq, MolarDensvap, ~, ~, ~, ~, ~, ~, ~, ~] ...
                    = Y_FluidFlash(iperf_P, Zi_wellstream, iperf_T, fluid, kval_0, fg_0, kflag, Pb);  
                
                    % Get perforation phase molar density to calculate well mobility
                    MolarDensLiq_inj(iPerf, 1) = MolarDensliq;   
                    MolarDensVap_inj(iPerf, 1) = MolarDensvap;  
                end
            end
            
            % Calculate total mobility
            Total_Fluid_Mob = (Oil_Mob + Gas_Mob + Wat_Mob);
            % Update injector well mobility
            WI_o(well.wconne_well == iWell, 1) = well.Inje_Gas_Flag_Perf(well.wconne_well == iWell, 1) .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) ...
                                              .* Total_Fluid_Mob(well.wconne_well == iWell, 1) .* MolarDensLiq_inj(well.wconne_well == iWell, 1);    
            WI_g(well.wconne_well == iWell, 1) = well.Inje_Gas_Flag_Perf(well.wconne_well == iWell, 1) .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) ...
                                              .* Total_Fluid_Mob(well.wconne_well == iWell, 1) .* MolarDensVap_inj(well.wconne_well == iWell, 1);    
            WI_o(isnan(WI_o) | isinf(WI_o)) = 0;
            WI_g(isnan(WI_g) | isinf(WI_g)) = 0;            

            %% Check effective well control
            Pwf_o_rate = 0;
            Pwf_g_rate = 0;
            
            HCInje_flag = OilWell_flag + GasWell_flag;                         

            if HCInje_flag > 0
            
                % Surface flash 
                Pb = init.wel_vle_chk(iWell, nComps+1);
                kval_0 = init.wel_vle_chk(iWell, 1:nComps);
                fg_0 = 0.5d0;
                kflag = 0; 
                [fv_ws_SC, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, MolarVolliq_ws_SC, MolarVolvap_ws_SC, ~, ~, ~, ~] ...
                = Y_FluidFlash(fluid.Psc, Zi_wellstream, fluid.Tsc, fluid, kval_0, fg_0, kflag, Pb);
            
                SMolarVolOil_sc = (1.0 - fv_ws_SC) * MolarVolliq_ws_SC;  % Specific molar volume of oil at surface condition
                SMolarVolGas_sc = fv_ws_SC * MolarVolvap_ws_SC;          % Specific molar volume of gas at surface condition 
                        
                % Pwf back calculated from surface oil volume rate constrain                                
                if OilWell_flag  
                    % Bottom hole pressure based on oil rate  
                    Pwf_o_rate = ( SMolarVolOil_sc * sum( WI_o(well.wconne_well == iWell) .* cellP_Hw + ... 
                                                          WI_g(well.wconne_well == iWell) .* ( cellP_Hw + prop(current).cell_pcgo(Perf_Cells) ) ) ...                          
                               +   OilVolRate_sc ) / ( SMolarVolOil_sc * sum( (WI_o + WI_g) .* (well.wconne_well == iWell) ) );
                    Pwf_o_rate(isnan(Pwf_o_rate) | isinf(Pwf_o_rate)) = 0;      
                    Pwf_o_rate = Pwf_o_rate * (Pwf_o_rate > 14.7);  

                    % If Pwf_rate physically make sense, add to Pwf_effective  
                    if Pwf_o_rate > 0 
                       nnz_p = nnz_p + 1;
                       Pwf_effective(nnz_p, 1) = Pwf_o_rate; 
                    end
                end 

                % Pwf back calculated from surface gas volume rate constrain                                
                if GasWell_flag
                    % Bottom hole pressure based on gas rate 
                    Pwf_g_rate = ( SMolarVolGas_sc * sum( WI_o(well.wconne_well == iWell) .* cellP_Hw + ... 
                                                          WI_g(well.wconne_well == iWell) .* ( cellP_Hw + prop(current).cell_pcgo(Perf_Cells) ) ) ...                          
                               +   GasVolRate_sc ) / ( SMolarVolGas_sc * sum( (WI_o + WI_g) .* (well.wconne_well == iWell) ) );
                    Pwf_g_rate(isnan(Pwf_g_rate) | isinf(Pwf_g_rate)) = 0;      
                    Pwf_g_rate = Pwf_g_rate * (Pwf_g_rate > 14.7);                                

                    % If Pwf_rate physically make sense, add to Pwf_effective  
                    if Pwf_g_rate > 0 
                       nnz_p = nnz_p + 1; 
                       Pwf_effective(nnz_p, 1) = Pwf_g_rate; 
                    end
                end  

            end
            
            % Determine the minimum effective BHP from rate control
            Pwf_rate_max = max(Pwf_effective(1:nnz_p, 1)) * (nnz_p > 0);
            Pwf_rate_max(isempty(Pwf_rate_max)) = 1.0E10;
            
            if Pwf_rate_max <= well.Target_max_BHP(iWell) 
                
                % One of the rate control is effective
                prop(current).well_P(iWell) = Pwf_rate_max;
                
                % Oil rate control is effective
                well_ctr(current).OilRate_Ctr(iWell) = ( Pwf_o_rate == Pwf_rate_max );

                % Gas rate control is effective 
                well_ctr(current).GasRate_Ctr(iWell) = ( Pwf_g_rate == Pwf_rate_max );           

            else
                % BHP control is effective
                prop(current).well_P(iWell) = well.Target_max_BHP(iWell);
                well_ctr(current).BHP_Ctr(iWell) = true;
            end 

            % Neglect gravity in wellbore first 
            perf_P(well.wconne_well == iWell) = prop(current).well_P(iWell) + Hw(well.wconne_well == iWell);
                        
            % Oil molar rate of each perforation at reservoir condition
            Qo_perf(well.wconne_well == iWell) = WI_o(well.wconne_well == iWell) ...
            .* ( perf_P(well.wconne_well == iWell) - prop(current).cell_P(Perf_Cells) );
            % Gas molar rate of each perforation at reservoir condition
            Qg_perf(well.wconne_well == iWell) = WI_g(well.wconne_well == iWell) ...
            .* ( perf_P(well.wconne_well == iWell) - ( prop(current).cell_P(Perf_Cells) + prop(current).cell_pcgo(Perf_Cells) ) );

            % Forbid cross flow (injection in producer / production in injector)
            perf_cross_flow = [find(Qo_perf<0); find(Qg_perf<0)];
            well_ctr(current).wconne_well(perf_cross_flow) = 0; 
                                
        end % Hydrocarbon Injector
        
        %% %%%%%%%%%%%%%%%%%%%
        %%% Water Injector %%%
        %%%%%%%%%%%%%%%%%%%%%%  
        if well.Inje_Wat_Flag(iWell) == 1 
            % Define extra variables for hydrocarbon injector
            MolarDensWat_inj = zeros(nPerfs, 1);    % Water molar density at each perforation to cal. injector mobility
            Xc = zeros(active_perf, 1);
            cell_Bw_inv = zeros(active_perf, 1);
                        
            % Neglect gravity in wellbore first 
            perf_P(well.wconne_well == iWell) = prop(current).well_P(iWell) + Hw(well.wconne_well == iWell);            
            
            % Water molar density:  lb-m/cu.ft
            Xc(1 : active_perf) = fluid.Cw .* ( perf_P(well.wconne_well == iWell) - fluid.Pref_w);
            cell_Bw_inv(1 : active_perf) = 1.0 ./ fluid.Bw_ref .* (1 + Xc(1 : active_perf) .* (1 + 0.5 .* Xc(1 : active_perf)) );
            MolarDensWat_inj(well.wconne_well == iWell)  = fluid.MDen_w_sc .* cell_Bw_inv(1 : active_perf) ./ fluid.wat_MW; 
            
            % Calculate total mobility
            Total_Fluid_Mob = (Oil_Mob + Gas_Mob + Wat_Mob);
            % Update injector well mobility
            WI_w(well.wconne_well == iWell, 1) = well.Inje_Wat_Flag_Perf(well.wconne_well == iWell, 1) .* well.Ind_geom_Perf(well.wconne_well == iWell, 1) ...
                                              .* Total_Fluid_Mob(well.wconne_well == iWell, 1) .* MolarDensWat_inj(well.wconne_well == iWell, 1);    

            % Water molar rate back calculated from surface water volume rate constrain
            Pwf_w_rate = 0;
            
            if WatWell_flag > 0                                               
                % Pwf based on water mass rate 
                Pwf_w_rate = ( SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) .* ( cellP_Hw - prop(current).cell_pcow(Perf_Cells) ) ) ...
                           +   WatVolRate_sc) / ( SMolarVolWat_sc * sum( WI_w(well.wconne_well == iWell) ) );
                Pwf_w_rate(isnan(Pwf_w_rate) | isinf(Pwf_w_rate)) = 0;      
                Pwf_w_rate = Pwf_w_rate * (Pwf_w_rate > 14.7);                                                
            end
            
            if Pwf_w_rate <= well.Target_max_BHP(iWell)
                % Water rate control is effective
                prop(current).well_P(iWell) = Pwf_w_rate;
                well_ctr(current).WatRate_Ctr(iWell) = true;
            else
                % BHP control is effective
                prop(current).well_P(iWell) = well.Target_max_BHP(iWell); 
                well_ctr(current).BHP_Ctr(iWell) = true;
            end

            % Neglect gravity in wellbore first 
            perf_P(well.wconne_well == iWell) = prop(current).well_P(iWell) + Hw(well.wconne_well == iWell);
                        
            % Water molar rate of each perforation at reservoir condition
            Qw_perf(well.wconne_well == iWell) = WI_w(well.wconne_well == iWell) ...
            .* ( perf_P(well.wconne_well == iWell) - ( prop(current).cell_P(Perf_Cells) - prop(current).cell_pcow(Perf_Cells) ) );

            % Forbid cross flow (injection in producer / production in injector)
            well_ctr(current).wconne_well(Qw_perf<0) = 0;                 
                        
        end % Water injector
        
    end
    
    % Assign perforation pressure
    prop(current).perf_P = perf_P;
        
end