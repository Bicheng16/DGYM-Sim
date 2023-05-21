function [prop] = Cal_Well_Gradient(well, prop, fluid, well_ctr, init)

    % Define temporary variables 
    nWells = well.nWells;
    nPerfs = sum(well.nperf_cell);
    nComps = fluid.nComps;
    current = 1;

    % Gravity acceleration  
    gc = 32.17 * 0.21584E-3; 
    
    % Output allocation 
    Hw = zeros(nPerfs, 1);
    wconne_well_0 = well.wconne_well;
    % wconne_well_0 = well_ctr(current).wconne_well;

    % Mobility for producers (Injectors modified individually)
    Oil_Mob = prop(current).cell_kro(well.wconne_cell) ./ prop(current).cell_Vsic_o(well.wconne_cell);         
    Gas_Mob = prop(current).cell_krg(well.wconne_cell) ./ prop(current).cell_Vsic_g(well.wconne_cell);        
    Wat_Mob = prop(current).cell_krw(well.wconne_cell) ./ prop(current).cell_Vsic_w(well.wconne_cell);        
    
    Oil_Mob(isinf(Oil_Mob) | isnan(Oil_Mob)) = 0.0d0;
    Gas_Mob(isinf(Gas_Mob) | isnan(Gas_Mob)) = 0.0d0;
    Wat_Mob(isinf(Wat_Mob) | isnan(Wat_Mob)) = 0.0d0;
           
    % Loop well
    for iWell = 1 : nWells
        % Producer only at the begining
        if ~well.Active_Flag(iWell)
            continue
        end
        
        % Number of perforation in iWell
        active_perf     = sum( wconne_well_0 == iWell );
                        
        MDen_mix        = zeros(active_perf, 1);       % Mixture density of wellstream at local perforation     
        Rel_Perf_Depth  = zeros(active_perf, 1);       % Depth between two neighbor Perforations: ft
        Hw_t            = zeros(active_perf, 1);       % Gravity potential between each perforation, psia
        Hwi             = zeros(active_perf, 1);       % Wellbore pressure head between the perforation and the well's bottom hole datum depth
        local_perf_P    = zeros(active_perf, 1);       % Local perforation pressure in each well
        local_perf_DZ   = zeros(active_perf, 1);       % Local perforation depth in each well
        
        % Get perforation cell ID connected to iWell
        Perf_Cells = well.wconne_cell(wconne_well_0 == iWell); 
        
        % Depth from perforated cell center to BHP reference depth
        Rel_Perf_Depth(1:active_perf, 1) = well.Depth_Perf(wconne_well_0 == iWell) - well.Ref_Depth(iWell);
                                             
        % Calculated the depth difference between every two neighbored perforation
        if active_perf > 1 
            Rel_Perf_Depth(2 : active_perf, 1) = Rel_Perf_Depth(2 : active_perf, 1) - Rel_Perf_Depth(1 : active_perf-1, 1);
        end

        local_perf_P(:) = prop(current).perf_P(wconne_well_0 == iWell); 
        local_perf_DZ(:) = well.Depth_Perf(wconne_well_0 == iWell); 
                
        %% %%%%%%%%%%%%%%%%  
        %%% Producers %%%%%
        %%%%%%%%%%%%%%%%%%%
        % 
        if well.Prod_Flag(iWell) == 1 && any(Rel_Perf_Depth) > 0
            
            Qv_o_perf       = zeros(active_perf, 1); 
            Qv_g_perf       = zeros(active_perf, 1); 
            Qv_w_perf       = zeros(active_perf, 1); 

            MDen_o_perf     = zeros(active_perf, 1);
            MDen_g_perf     = zeros(active_perf, 1);
            MDen_w_perf     = zeros(active_perf, 1);
                        
            % Get well stream 
            Zi_wellstream = well.Gas_Comp(iWell, :);
                                              
            % Calculate fluid phase density at each perforation based on wellstream flash and
            % perforation pressure 
            iperf_T = prop(current).cell_T(1);  
            
            % Estimate bubble point pressure of wellstream            
            Pb = init.wel_vle_chk(iWell, nComps+1);
            kval_0 = init.wel_vle_chk(iWell, 1:nComps);
            fg_0 = 0.5d0;  
            kflag = 0;

            % Producer: Upstream: last perf ===> Downstream: first perf
            for iPerf = 1 : active_perf
                                   
                iperf_P = local_perf_P(iPerf);

                % If last timestep the perf is closed, perf_P can be zero
                if iperf_P == 0.0 
                    ups_nnz = nnz( local_perf_P(iPerf+1:end) > 14.7 );
                    % More than 1 perf in upstream with nonzero perf
                    % pressure, interpolation preferred
                    
                    temp = find( local_perf_P > 14.7 );
                    if ups_nnz > 1                        
                        temp = temp(temp > iPerf);

                        % Interpolate iperf_P based on other perf pressures
                        % based on depth from upstream
                        iperf_P = nakeinterp1( local_perf_DZ(temp), ...
                                               local_perf_P(temp), ...
                                               local_perf_DZ(iPerf) );
                    % Only single upstream perf with nnz pressure, inherit
                    elseif ups_nnz == 1
                        temp = temp(temp > iPerf);
                        loc  = min(temp);  
                        iperf_P = local_perf_P(loc);
                    % No upstream perf with nnz pressure, get pressure from
                    % closest downstream perf
                    elseif (ups_nnz == 0 && iPerf > 1) 
                        temp = temp(temp < iPerf);
                        loc  = max(temp);
                        iperf_P = local_perf_P(loc);
                    % Perf is the most downstream, and no upstream perf
                    % with nnz pressure, get pressure from perf cell
                    % pressure
                    else
                        iperf_P = prop(current).cell_P( Perf_Cells(iPerf) );
                    end
                    local_perf_P(iPerf) = iperf_P; 
                end

                % Oil and gas mass densities              
                [~, ~, ~, ~, ~, ~, ~, ~, ~, oDen, gDen, ~, ~, ~, ~, ~, ~] ...
                = Y_FluidFlash(iperf_P, Zi_wellstream, iperf_T, fluid, kval_0, fg_0, kflag, Pb);                

                MDen_o_perf(iPerf, 1) = oDen;   
                MDen_g_perf(iPerf, 1) = gDen; 

                % Water mass density
                Xc = fluid.Cw * ( iperf_P - fluid.Pref_w);
                cell_Bw_inv = 1.0 / fluid.Bw_ref * (1 + Xc * (1 + 0.5 * Xc) );
                MDen_w_perf(iPerf, 1) = fluid.MDen_w_sc * cell_Bw_inv; 
                                                           
            end             
                                               
            % Calculate volumetric rate of each phase
            % Oil volume rate of each perforation at reservoir condition
            Qv_o_perf(wconne_well_0 == iWell) = Oil_Mob(wconne_well_0 == iWell) .* well.Ind_geom_Perf(wconne_well_0 == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) - local_perf_P);
            % Gas volume rate of each perforation at reservoir condition
            Qv_g_perf(wconne_well_0 == iWell) = Gas_Mob(wconne_well_0 == iWell) .* well.Ind_geom_Perf(wconne_well_0 == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) + prop(current).cell_pcgo(Perf_Cells) - local_perf_P);
            % Water volume rate of each perforation at reservoir condition
            Qv_w_perf(wconne_well_0 == iWell) = Wat_Mob(wconne_well_0 == iWell) .* well.Ind_geom_Perf(wconne_well_0 == iWell) ...
            .* (prop(current).cell_P(Perf_Cells) - prop(current).cell_pcow(Perf_Cells) - local_perf_P);
            
            Qv_o_perf(Qv_o_perf < 0.0d0) = 0.0d0;
            Qv_g_perf(Qv_g_perf < 0.0d0) = 0.0d0;
            Qv_w_perf(Qv_w_perf < 0.0d0) = 0.0d0;
                                     
            % Mixture mass density of wellstream at local perforation: lb/ft3
            MDen_mix = ( MDen_o_perf .* Qv_o_perf + MDen_g_perf .* Qv_g_perf + MDen_w_perf .* Qv_w_perf ) ... 
                    ./ ( Qv_o_perf + Qv_g_perf + Qv_w_perf );
                
            MDen_mix(isnan(MDen_mix) | isinf(MDen_mix)) = 0.0;

            % Producer: Upstream: last perf ===> Downstream: first perf
            for iPerf = 1 : active_perf
                                   
                if MDen_mix(iPerf) == 0.0 
                    ups_nnz = nnz( MDen_mix(iPerf+1:end) > 0 );
                    % More than 1 perf in upstream with nonzero perf
                    % mixture density, interpolation preferred
                    
                    temp = find( MDen_mix > 0 );
                    if ups_nnz > 1                        
                        temp = temp(temp > iPerf);

                        % Extrapolate iperf_P based on other perf pressures
                        % based on depth from upstream
                        MDen_mix(iPerf) = nakeinterp1( local_perf_DZ(temp), ...
                                                       MDen_mix(temp), ...
                                                       local_perf_DZ(iPerf) );
                    % Only single upstream perf with nnz density, inherit
                    elseif ups_nnz == 1
                        temp = temp(temp > iPerf);
                        loc  = min(temp);  
                        MDen_mix(iPerf) = MDen_mix(loc);
                    % No upstream perf with nnz density, get density from
                    % closest downstream perf
                    elseif (ups_nnz == 0 && iPerf > 1) 
                        temp = temp(temp < iPerf);
                        loc  = max(temp);
                        MDen_mix(iPerf) = MDen_mix(loc);
                    % Perf is the most downstream, and no upstream perf
                    % with nnz density, get density from perf cell
                    % mixture density
                    else
						iCell = Perf_Cells(iPerf);
                        MDen_mix(iPerf) = ( prop(current).cell_MDen_o(iCell) * prop(current).cell_So(iCell) ... 
						                +   prop(current).cell_MDen_g(iCell) * prop(current).cell_Sg(iCell) ...
                                        +   prop(current).cell_MDen_w(iCell) * prop(current).cell_Sw(iCell) ); 
                                        %/ ( prop(current).cell_So(iCell) + prop(current).cell_Sg(iCell) + prop(current).cell_Sw(iCell) );										
                    end
                    
                end            

            end
            
            % Gravity potential between each two neighbored segment of perforation from heel to toe, psia
            Hw_t(1, 1) = MDen_mix(1, 1);
            Hw_t(2 : active_perf, 1) = 0.5 .* ( MDen_mix(1 : active_perf-1, 1) +  MDen_mix(2 : active_perf, 1) ) .* (active_perf > 1);
            Hw_t(:, 1) = Hw_t(:, 1) .* gc .* Rel_Perf_Depth(:, 1); 
            
            % Gravity potential between each perforated cell and reference depth, psia
            Hwi(1, 1) = Hw_t(1, 1);
            if active_perf > 1
                for iPerf = 2 : active_perf 
                    Hwi(iPerf, 1) = Hwi(iPerf-1, 1) + Hw_t(iPerf, 1);
                end
            end 
            
            % Put it back to global well gradient list
            Hw(wconne_well_0 == iWell) = Hwi; 
        
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %%% Hydrocarbon Injector %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        if well.Inje_Gas_Flag(iWell) == 1 && any(Rel_Perf_Depth) > 0
            
            
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%
        %%% Water Injector %%%
        %%%%%%%%%%%%%%%%%%%%%%  
        if well.Inje_Wat_Flag(iWell) == 1 && any(Rel_Perf_Depth) > 0

            % Mixture density is the water density
            % Injector: Upstream: 1st perf ===> Downstream: last perf
            for iPerf = 1 : active_perf
                                   
                iperf_P = local_perf_P(iPerf);

                % If last timestep the perf is closed, perf_P can be zero
                if iperf_P == 0.0 
                    ups_nnz = nnz( local_perf_P(1:iPerf-1) > 14.7 );
                    % More than 1 perf in upstream with nonzero perf
                    % pressure, interpolation preferred
                    
                    temp = find( local_perf_P > 14.7 );
                    if ups_nnz > 1                        
                        temp = temp(temp < iPerf);

                        % Interpolate iperf_P based on other perf pressures
                        % based on depth from upstream
                        iperf_P = nakeinterp1( local_perf_DZ(temp), ...
                                               local_perf_P(temp), ...
                                               local_perf_DZ(iPerf) );
                    % Only single upstream perf with nnz pressure, inherit
                    elseif ups_nnz == 1
                        temp = temp(temp < iPerf);
                        loc  = max(temp);  
                        iperf_P = local_perf_P(loc);
                    % No upstream perf with nnz pressure, get pressure from
                    % closest downstream perf
                    elseif (ups_nnz == 0 && iPerf < active_perf) 
                        temp = temp(temp > iPerf);
                        loc  = min(temp);
                        iperf_P = local_perf_P(loc);
                    % Perf is the most downstream, and no upstream perf
                    % with nnz pressure, get pressure from perf cell
                    % pressure
                    else
                        iperf_P = prop(current).cell_P( Perf_Cells(iPerf) );
                    end
                    local_perf_P(iPerf) = iperf_P; 
                end

                % Water mass density
                Xc = fluid.Cw * ( iperf_P - fluid.Pref_w);
                cell_Bw_inv = 1.0 / fluid.Bw_ref * (1 + Xc * (1 + 0.5 * Xc) );
                MDen_mix(iPerf, 1) = fluid.MDen_w_sc * cell_Bw_inv; 
                                                           
            end             
                                                            
            % Gravity potential between each two neighbored segment of perforation from toe to heel (injector), psia
            Hw_t(active_perf, 1) = MDen_mix(active_perf, 1);
            Hw_t(active_perf-1 : -1 : 1 , 1) = 0.5 * ( MDen_mix(active_perf : -1 : 2, 1) +  MDen_mix(active_perf-1 : -1 : 1 , 1) ) .* (active_perf > 1);            
            Hw_t(:, 1) = Hw_t(:, 1) .* gc .* Rel_Perf_Depth(:, 1); 
            
            % Gravity potential between each perforated cell and reference depth, psia
            Hwi(1, 1) = Hw_t(1, 1);
            if active_perf > 1
                Hwi(2:active_perf, 1) = Hwi(1 : active_perf - 1, 1) + Hw_t(2 : active_perf, 1);
            end 
            
            % Put it back to global well gradient list
            Hw(wconne_well_0 == iWell) = Hwi;             
                        
        end
                
    end
    
    % Output pressure different between reference depth and each local
    % perforation, psia
    prop(current).perf_Hw(:) = Hw;

end