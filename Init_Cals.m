function [init, prop, rst, res, fluid, well, well_ctr] = Init_Cals(init, fluid, res, well, rockfluid)

    disp('...Initialization');

    % Initilize some frequently used variables
    nCells = init.nCells;
    nConnes = init.nConnes;
    nComps = fluid.nComps;
    nWells = well.nWells;
    nPerfs = sum(well.nperf_cell);
    current = 1; 
     
    % Calculate bubble point in reservoir
    Psat_0 = 2000;
    kval_0 = zeros(1, nComps);
    for iRT = 1 : res.nRockTypes
        cells = init.cell_group(iRT).iGroup;
        nnc   = nnz(cells);
        if isequal(fluid.fluidtype, 'VLE')            
            [Psat, kval] = Y_SatPoint(res.Zi(iRT, :), res.T(iRT), fluid, Psat_0, kval_0, 1);            
            init.res_vle_chk(cells, 1:nComps) = repmat(kval, [nnc, 1]);
            init.res_vle_chk(cells, nComps+1) = Psat;
            init.res_vle_chk(cells, nComps+2) = res.Zi(iRT, fluid.key_comp);
        else 
            init.res_vle_chk(cells, 1:nComps) = 0.0d0;
            init.res_vle_chk(cells, nComps+1) = 1.0d4;
            init.res_vle_chk(cells, nComps+2) = res.Zi(iRT, fluid.key_comp);            
        end           
    end
    
    % Calculate bubble point for well stream composition 
    Psat_0 = 2000;
    for iWell = 1 : nWells
        Gas_Comp = well.Gas_Comp(iWell, :);
        if abs( sum(Gas_Comp) - 1.0 ) > 1.0e-12
            Gas_Comp = res.Zi(1, :);
        end        
        if isequal(fluid.fluidtype, 'VLE') && ~well.Inje_Wat_Flag(iWell)            
            [Psat, kval] = Y_SatPoint(Gas_Comp, res.T(1), fluid, Psat_0, kval_0, 1);            
            init.wel_vle_chk(iWell, 1:nComps) = kval;
            init.wel_vle_chk(iWell, nComps+1) = Psat;
            init.wel_vle_chk(iWell, nComps+2) = Gas_Comp(1, fluid.key_comp);
        else
            init.wel_vle_chk(iWell, 1:nComps) = 0.0d0;
            init.wel_vle_chk(iWell, nComps+1) = 1.0d4;
            init.wel_vle_chk(iWell, nComps+2) = Gas_Comp(1, fluid.key_comp);            
        end           
    end    

    % Prepare data base to correlate K-factor and fv with pressure
    % Pressure range: 100 psia to Pb, step size is around 5 psia
    % Gas injection excluded here
    [fluid] = Flash_Data_Base(fluid, res, init);
    
    % Memory allocation time-step and newton-step (primary & secondary) variable
    prop = alloc_dynamic_properties(nCells, nConnes, nComps, nWells, nPerfs);

    % Assignment for current properties
    prop(current).cell_P(1 : nCells) = init.cell_P(1 : nCells);
    prop(current).cell_T(1 : nCells) = init.cell_T(1 : nCells);
    prop(current).cell_Sw(1 : nCells) = init.cell_Sw(1 : nCells);
    prop(current).cell_Zi(1 : nCells, 1 : nComps) = init.cell_Zi(1 : nCells, 1 : nComps);
    
    prop(current).well_P(well.Prod_Flag == 1) = well.Target_min_BHP(well.Prod_Flag == 1);
    prop(current).well_P(well.Prod_Flag == 0) = well.Target_max_BHP(well.Prod_Flag == 0);
           
    % Initialize fluid system in the whole reservoir
    prop = Init_Fluids(current, prop, fluid, init);
    
    % Initialize rock properties in the whole reservoir
    TargetCells = 1 : nCells;
    dummy = alloc_inherit_var;
    prop = Properties_Rock(current, prop, init, res, rockfluid, dummy, TargetCells, false, []); 
    
    % Initialize gas adsorption in the whole reservoir
    %[prop] = Gas_Adsorption_Cals(res, init, current, prop);
    prop = Gas_Adsorption_Cals(current, prop, init, res, fluid, dummy, TargetCells, false, [], false);

    % Initialize gas slippage apparent permeability multiplier
    prop = Gas_Slippage_Cals(current, prop, init, res, fluid, dummy, TargetCells, false, [], false);
    
    % Initialize generalized Fickian diffusion coefficient 
    prop = Gas_GFick_Cals(current, prop, init, res, fluid, dummy, TargetCells, false, [], false);
    
    % Initialize hydrocarbon and water mass in each cell
    [prop] = Fluid_Mass(prop, current, nComps);

    % Initialize well: get effective well control and well pressure
    % gradient (later)
    [well_ctr, prop, well] = Init_Well(well, prop, fluid, init);
    
    % Calculate well gradient
    [prop] = Cal_Well_Gradient(well, prop, fluid, well_ctr, init);
    
    % Memory allocation for result variable
    rst = alloc_rst_var(nComps, nWells);
    
    % Initialize result variables
    [rst] = Result_Var_Cals(current, rst, fluid, prop, init);
    
    % Assign current time-step variables to initial variables
    init.cell_Zi(1 : nCells, 1 : nComps) = prop(current).cell_Zi(1 : nCells, 1 : nComps);
    init.cell_Ni(1 : nCells, 1 : nComps) = prop(current).cell_Ni(1 : nCells, 1 : nComps);
    init.cell_Nw(1 : nCells) = prop(current).cell_Nw(1 : nCells);
    init.cell_Qai(1 : nCells, 1 : nComps) = prop(current).cell_Qai(1 : nCells, 1 : nComps);

    init.OOIP = rst.FOIP;                        % Initial oil in place, STB 
    init.OGIP = rst.FGIP;                        % Initial gas in place, MSCF
    init.OWIP = rst.FWIP;                        % Initial water in place, STB
    init.OAGIP = rst.FAGIP;                      % Initial adsorbed gas in place, MSCF 
    
end