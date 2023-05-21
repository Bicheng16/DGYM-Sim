% Set basic property vectors and check input
function [init, fluid, well, rockfluid, flag_run_finish_ok] = Misc_Cals(init, fluid, res, rockfluid, sim_ctr, well)


    disp('...Calculating Miscellaneous Data and Checking Input Data')

    % Boolean to terminate simulation 
    flag_run_finish_ok = true;
    
    % Initilize some frequently used variables
    nCells = init.nCells;
    nConnes = init.nConnes;
    nComps = fluid.nComps;
    nWells = well.nWells;
    Rg = 10.732;
    FtoR = 459.598;           % Convert F to R for temperature (+)
    coef1 = 62.4297665d0;     % Convert cu.ft/lbmol to cu.cm/mol (*)
    
    % Calculate non-water moleculars' collision diameter, meter  
    for iComp = 1 : nComps 
        Vc = coef1 * fluid.comp_Zc(iComp) * Rg * (fluid.comp_Tc(iComp) + FtoR) /fluid.comp_Pc(iComp); 
        fluid.comp_dmc(iComp) = 1.0d-10 * 0.809d0 * Vc^(1.0d0/3.0d0);
    end
    
    % Assign values for init struct
    init.cell_Cr(1 : nCells) = res.Cr(init.cell_type(1 : nCells));
    init.cell_actnum = init.cell_poro > 1.0e-6;
    init.cell_P(1 : nCells) = res.P(init.cell_type(1 : nCells));
    init.cell_T(1 : nCells) = res.T(init.cell_type(1 : nCells));
    init.cell_Sw(1 : nCells) = res.Sw(init.cell_type(1 : nCells));
    for i = 1 : nCells
        init.cell_Zi(i, 1 : nComps) = res.Zi(init.cell_type(i), 1 : nComps);
    end
     
    init.conne_dz(1 : nConnes) = init.cell_depth(init.conne_cell2(1 : nConnes)) ...
                               - init.cell_depth(init.conne_cell1(1 : nConnes));
                           
    % Store cell names of each rock type
    cell_list = 1 : init.nCells;
    cell_group = [];
    % Get starting cell number of each rock type
    for iRT = 1 : res.nRockTypes
        tmp1 = cell_list' .* (init.cell_type == iRT); 
        tmp2 = tmp1(tmp1 ~= 0);
        dummy.iGroup = tmp2;
        cell_group = [cell_group; dummy];
    end  
    init.cell_group = cell_group;
    
    % Check key fluid composition
    [~, fluid.key_comp] = max(fluid.comp_MW);
                                                     
    % Modify some data in well struct
    temp_1 = [];
    for iWell = 1 : nWells
        Well_Location = (well.wconne_well == iWell);                % Indices for well locations based on Well_ID
        % Reset flags if the well is inactive 
        if not(well.Active_Flag(iWell))
            well.Prod_Flag(iWell) = 0;
            well.Bond_Flag(iWell) = 0;
            well.Inje_Wat_Flag(iWell) = 0;
            well.Inje_Gas_Flag(iWell) = 0;
            well.Flag_Perf(Well_Location) = 0;
            well.Prod_Flag_Perf(Well_Location) = 0;
            well.Bond_Flag_Perf(Well_Location) = 0;
            well.Inje_Wat_Flag_Perf(Well_Location) = 0;
            well.Inje_Gas_Flag_Perf(Well_Location) = 0;
        end 
        % If gas injector, record the perf cell ID
        if well.Inje_Gas_Flag(iWell)
            temp_2 = well.wconne_cell(well.wconne_well == iWell);
            temp_1 = [temp_1 temp_2];
        end
    end
    well.Depth_Perf(:) = init.cell_depth(well.wconne_cell) .* well.Flag_Perf;     
    well.Ind_geom_Perf = well.Ind_geom_Perf .* well.Flag_Perf;
    well.Inj_Gas_cells = temp_1; 
             
    %% Well check 
    % Check Number_of_Wells and WellID
    m = max(well.wconne_well);
    if m ~= well.nWells
        disp('...ERROR: Number_of_Wells does not match actual number of wells defined in Well_Connect_ID')
        disp('...ERROR: simulation will be terminated')
        flag_run_finish_ok = false;
    end    
    %Check well type flags (only one type can be defined)
    Check_Welltype = well.Prod_Flag + well.Inje_Wat_Flag + well.Inje_Gas_Flag;
    if max(Check_Welltype)>1
        disp('...ERROR: A Well has been set to more than one type of well. Check Well_Prod_Flag, Well_Inje_Wat_Flag, Well_Inje_Gas_Flag, Well_Inje_Oil_Flag')
        disp('...ERROR: simulation will be terminated')
        flag_run_finish_ok = false;    
    end


    if flag_run_finish_ok 
    
        %% Compositional fluid properties check

        %Compositional EOS properties
        % sums compositions row-wise giving sum(Zi) for each cell
        dummy1 = sum(init.cell_Zi, 2);    
        dummy1 = abs(dummy1 - 1.0);
        dummy1 = sum(dummy1 > 1.0e-10);
        if dummy1 ~= 0
            disp('...ERROR: Fluid composition Zi_init does not add to unity for at least one of the cells. Check fluid composition input data')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end  
        
        [~, n] = size(init.cell_Zi);
        if n ~= fluid.nComps
            disp('...ERROR: Composition has not been provided for all components in Zi_init')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end        
        
        [~, n] = size(fluid.comp_MW);
        if n ~= fluid.nComps
            disp('...ERROR: MW has not been provided for all components in MW_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end        
        
        [~, n] = size(fluid.comp_Tc);
        if n ~= fluid.nComps
            disp('...ERROR: Critical temperature has not been provided for all components in Tcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        [~, n] = size(fluid.comp_Pc);
        if n ~= fluid.nComps
            disp('...ERROR: Critical pressure has not been provided for all components in pcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         
        
        [~, n] = size(fluid.comp_Zc);
        if n ~= fluid.nComps
            disp('...ERROR: Critical compressibility factor has not been provided for all components in Zcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         
        
        [~, n] = size(fluid.comp_ACF);
        if n ~= fluid.nComps
            disp('...ERROR: Acentric factor has not been provided for all components in ACF_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        [~, n] = size(fluid.comp_SSHIFT);
        if n ~= fluid.nComps
            disp('...ERROR: Volume shift factor has not been provided for all components in SSHIFT_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         
 
        [~, n] = size(fluid.comp_PARA);
        if n ~= fluid.nComps
            disp('...ERROR: PARACHOR has not been provided for all components in PARA_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         
 
        [m, n] = size(fluid.comp_BIC);
        if (n ~= fluid.nComps) || (m ~= fluid.nComps)
            disp('...ERROR: Binary interaction coefficients has not been provided for all components in BIC_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        [~, n] = size(fluid.LBC_Coeffs);
        if n ~= 5
            disp('...ERROR: All 5 LBC coefficients have not been provided in LBC_Coeffs')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        [~, n] = size(well.Gas_Comp);
        if n ~= fluid.nComps
            disp('...ERROR: Composition of injected gas not defined for all components in Well_Inje_Gas_Comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end        
        
        for i = 1 : nWells
            if well.Inje_Gas_Flag(i) == 1
                dummy1 = sum(well.Gas_Comp(i,:));
                if dummy1 ~= 1
                    disp('...ERROR: Fluid composition for gas injector well (Well_Inje_Gas_Comp) does not add to unity')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end    
            end
        end        

        dummy1 = sum((fluid.comp_SSHIFT ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...WARNING: Volume shift has not been provided for all components using SSHIFT_comp')
            disp('.......     Simulation will continue')
        end        
 
        dummy1 = sum((fluid.comp_MW ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...ERROR: Molecular weight has not been provided for all the components using MW_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end   
        
        dummy1 = sum((fluid.comp_Tc ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...ERROR: Critical temperature has not been provided for all the components using Tcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end 
        
        dummy1 = sum((fluid.comp_Pc ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...ERROR: Critical pressure has not been provided for all the components using Tcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        dummy1 = sum((fluid.comp_Zc ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...ERROR: Critical compressibility factor has not been provided for all the components using Zcrit_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end         

        dummy1 = sum((fluid.comp_ACF ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...ERROR: Acentric factors has not been provided for all the components using ACF_comp')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end          

        dummy1 = sum((fluid.comp_PARA ~= 0));    
        if dummy1 ~= fluid.nComps
            disp('...WARNING: Parachors has not been provided for all the components using PARA_comp')
            disp('.......     Simulation will continue')
        end    
 
        if issymmetric(fluid.comp_BIC) == 0
            disp('...ERROR: Binary interaction coefficients (BIC_comp) values are not symmetrical')
            disp('...ERROR: simulation will be terminated')
            flag_run_finish_ok = false;
        end
        
    
        %% Rock-fluid table check
        
        % Loop through rock type
        for iRT = 1 : res.nRockTypes

            % Check SWFN table
            dummy1 = size(rockfluid(iRT).SWFN(:, 1));
            kr_tab_size = dummy1(1);
            Swco = rockfluid(iRT).SWFN(1, 1);
            Swmax = rockfluid(iRT).SWFN(kr_tab_size, 1);
            if rockfluid(iRT).SWFN(1, 2) ~= 0
                disp('...ERROR: krw at Swco (first value/ second column in SWFN) must be zero')
                disp('...ERROR: simulation will be terminated')
                flag_run_finish_ok = false;
            end  

            for i=1 : kr_tab_size-1

                if rockfluid(iRT).SWFN(i, 2) == 0
                    Swcr = rockfluid(iRT).SWFN(i, 1);
                end        
                dummy1 = rockfluid(iRT).SWFN(i+1, 1) - rockfluid(iRT).SWFN(i, 1);
                dummy2 = rockfluid(iRT).SWFN(i+1, 2) - rockfluid(iRT).SWFN(i, 2);
                dummy3 = rockfluid(iRT).SWFN(i+1, 3) - rockfluid(iRT).SWFN(i, 3);
                if dummy1 < 0
                    disp('...ERROR: Sw must increase or level on SWFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy2 < 0
                    disp('...ERROR: krw must increase or level on SWFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy3 > 0
                    disp('...ERROR: pcow must decrease or level on SWFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end

            end    

            % Check SGFN table
            dummy1 = size(rockfluid(iRT).SGFN(:, 1));
            kr_tab_size = dummy1(1);
            Sgco = rockfluid(iRT).SGFN(1,1);
            Sgmax = rockfluid(iRT).SGFN(kr_tab_size, 1);
            if rockfluid(iRT).SGFN(1, 2) ~= 0
                disp('...ERROR: krg at Sgco (first value/ second column in SGFN) must be zero')
                disp('...ERROR: simulation will be terminated')
                flag_run_finish_ok = false;
            end    

            dummy1 = 1 - Swco;
            if Sgmax ~= dummy1
                disp('...WARNING: Sgmax (in SGFN table) must equal 1-Swco (in SWFN table)')
                disp('......      Simulation will continue. Reviewing kr curves is advised')
            end 

            for i = 1 : kr_tab_size-1

                if rockfluid(iRT).SGFN(i, 2) == 0
                    Sgcr = rockfluid(iRT).SGFN(i, 1);
                end

                dummy1 = rockfluid(iRT).SGFN(i+1, 1) - rockfluid(iRT).SGFN(i, 1);
                dummy2 = rockfluid(iRT).SGFN(i+1, 2) - rockfluid(iRT).SGFN(i, 2);
                dummy3 = rockfluid(iRT).SGFN(i+1, 3) - rockfluid(iRT).SGFN(i, 3);
                if dummy1 < 0
                    disp('...ERROR: Sg must increase or level on SGFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy2 < 0
                    disp('...ERROR: krg must increase or level on SGFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy3 < 0
                    disp('...ERROR: pcgo must increase or level on SGFN')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
            end 

            % Check SOF3 table 
            dummy1 = size(rockfluid(iRT).SOF3(:, 1));
            kr_tab_size = dummy1(1);
            kro_cw = nakeinterp1(rockfluid(iRT).SOF3(:, 1), rockfluid(iRT).SOF3(:, 2), (1 - Swco)); %To be used in STONEII if kro_hat method used in Dynamic_Fluid_Rock_Props_fPressure_3Phase
            Soco = rockfluid(iRT).SOF3(1, 1);
            Somax = rockfluid(iRT).SOF3(kr_tab_size, 1);
            if rockfluid(iRT).SOF3(1, 2) ~= 0 || rockfluid(iRT).SOF3(1, 3) ~= 0 || Soco ~= 0
                disp('...ERROR: Soco, krow, and krog (first values of all three columns in SOF3) must be zero')
                disp('...ERROR: simulation will be terminated')
                flag_run_finish_ok = false;
            end    

            dummy1 = 1 - Swco;
            if Somax ~= dummy1
                disp('...WARNING: Somax (in SOF3 table) must equal 1-Swco (in SWFN table)')
                disp('......      Simulation will continue. Reviewing kr curves is advised')
            end 

            if Somax == 1
                disp('...STONE II method for 3-phase kro_bar will be used. Review PETE 604 Lect1 for additional information')
                STONE2_Type = 1; % 1=uses kro_bar approach.   0=uses kro_hat approach
            else   
                disp('...STONE II method for 3-phase kro_hat will be used. Review PETE 604 Lect1 for additional information')
                STONE2_Type = 0; % 1=uses kro_bar approach.   0=uses kro_hat approach
            end 

            for i=1 : kr_tab_size-1

                if rockfluid(iRT).SOF3(i, 2) == 0 && rockfluid(iRT).SOF3(i, 3) == 0
                    Socr = rockfluid(iRT).SOF3(i, 1);
                end

                if rockfluid(iRT).SOF3(i, 2)==0
                    Sorw = rockfluid(iRT).SOF3(i, 1);
                end

                if rockfluid(iRT).SOF3(i, 3)==0
                    Sorg = rockfluid(iRT).SOF3(i, 1);
                end

                dummy1 = rockfluid(iRT).SOF3(i+1, 1) - rockfluid(iRT).SOF3(i, 1);
                dummy2 = rockfluid(iRT).SOF3(i+1, 2) - rockfluid(iRT).SOF3(i, 2);
                dummy3 = rockfluid(iRT).SOF3(i+1, 3) - rockfluid(iRT).SOF3(i, 3);
                if dummy1 < 0
                    disp('...ERROR: So must increase or level on SOF3')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy2 < 0
                    disp('...ERROR: krow must increase or level on SOF3')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
                if dummy3 < 0
                    disp('...ERROR: krog must increase or level on SOF3')
                    disp('...ERROR: simulation will be terminated')
                    flag_run_finish_ok = false;
                end
            end

            % Assign rockfluid data 
            rockfluid(iRT).method = STONE2_Type;
            rockfluid(iRT).Swco = Swco;
            rockfluid(iRT).kro_cw = kro_cw;
            
        end % End rock type loop
    
    
        %% Check solution method
        if sim_ctr.method == 1
            disp('...Compositional Fully-Implicit Finite Volume Method Selected')
        else
            disp('...ERROR: Invalid Solution Method Selected (Solution_Method in InputData)')
            flag_run_finish_ok = false;
        end
    
    end
    
    
end



















