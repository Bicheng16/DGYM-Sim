function [init, fluid, res, rockfluid, sim_ctr, well] = Scan_Input(input)

    % Initialization
    file_end = false;     % Flag for file end
    nCells = 1;           % Number of cells in the reservoir
    nConnes = 1;          % Number of connections among all cells 
    nComps = 1;           % Number of component in hydrocarbon from reservoir fluid  
    nRockTypes = 1;       % Number of rock type in the reservoir
    nWells = 1;           % Number of wells, default is 1
    
    % Open the file and start reading
    input = strcat(pwd, input);
    fid = fopen(input);
    disp('...Reading Input Data');

    while not(file_end)

        % Read line by line 
        [tline, str1] = scanline(fid);
            
        switch str1
                          
            % Scan grid data                 
            case 'GLOBALCELL'                
                [nCells, cell_name, cell_type, cell_poro, cell_volume, cell_depth] =  Scan_Grid(fid);
                nRockTypes = max(cell_type);
                file_end = false;

            % Scan connection data                 
            case 'CONNECTIONS'  
                [trans_split, nConnes, conne_cell1, conne_cell2, conne_type, conne_ctrans, conne_geo, conne_perm] = Scan_Connection(fid);                                        
                file_end = false;
                
            % Scan fluid PVT data
            case 'PVT' 
                [fluidtype, eos, nComps, comp_name, comp_MW, comp_ACF, comp_Pc, comp_Tc, comp_Zc, ...
                comp_PARA, comp_SSHIFT, comp_BIC, LBC_Coeffs, Psc, Tsc] = Scan_PVT(fid);
                file_end = false;

            % Scan water properties 
            case 'WATER'
                [Pref_w, MDen_w_sc, Cw, Bw_ref, Vsic_w_ref, Cvw] = Scan_Water(fid);
                file_end = false;
                
            % Scan initial data 
            case 'INIT'
                [Cr, P_ref, P, T, Sw, Zi, GKap, dpore, Fickian, GFD, Sorption, RMDen, VLi, PLi] = Scan_InitialData(fid, nRockTypes, nComps);
                file_end = false;
                
            % Scan rock-fluid data 
            case 'ROCKFLUID'
                [rf] = Scan_RockFluidData(fid, nRockTypes);                
                file_end = false;
                
            % Scan numerical control data 
            case 'NUMERICAL'
                [method, solver, t_total, dt_init, dt_max, dt_min, dt_mul, dt_back, ...
                dbhp_thresh, dpp_thresh, tol, iter_max] = Scan_SimCtr(fid);
                file_end = false;
                
            % Scan well data 
            case 'WELL'
                nWells = str2double(tline{2});                
            
                [name, nperf_cell, Prod_Flag, Inje_Wat_Flag, Inje_Gas_Flag,         ... 
                Target_min_BHP, Target_max_BHP, TargetLiquidRate, TargetOilRate,    ...    
                TargetGasRate, TargetWaterRate, Well_Gas_Comp, Ref_Depth,           ...
                Active_Flag, wconne_well, wconne_cell, Flag_Perf, Prod_Flag_Perf,   ...
                Inje_Wat_Flag_Perf, Inje_Gas_Flag_Perf, Ind_geom_Perf]              ...
                = Scan_WellCtr(fid, nWells, nComps);
                        
                file_end = false;
            
            % Scan the end of the file
            case 'ENDF'
                file_end = true;
            
            % Scan any line before end sign that are non-keywords (comment)
            otherwise              
                continue
                
        end
        
    end
    
    % Close file
    fclose(fid);
        
    % Assignment to corresponding struct
    
    % init struct
    init.nCells = nCells;                   % Total number of cells
    init.cell_name = cell_name;             % Cell name
    init.cell_type = cell_type;             % Cell flag for rock type
    init.cell_poro = cell_poro;             % Cell initial porosity 
    init.cell_poro_ref = cell_poro;         % Cell reference porosity
    init.cell_volume = cell_volume;         % Cell volume, cu.ft
    init.cell_depth = cell_depth;           % Cell center depth, ft
    init.cell_Cr = zeros(nCells, 1);        % Rock compressibility for each cell
    init.cell_actnum = zeros(nCells, 1);    % Active cell flag
    init.cell_P = ones(nCells, 1);          % Cell initial pressure, psia
    init.cell_T = ones(nCells, 1);          % Cell initial temperature, F
    init.cell_Sw = ones(nCells, 1);         % Cell initial water saturation
    init.cell_Zi = ones(nCells, nComps);    % Cell initial composition 
    init.cell_Ni = zeros(nCells, nComps);   % Cell initial number of moles of hydrocarbon in unit cell volume, lb-mol/cu.ft
    init.cell_Nw = zeros(nCells, 1);        % Cell initial number of mass of water per unit cell volume, lb/cu.ft
    init.cell_Qai= zeros(nCells, nComps);   % Cell initial number of moles of hydrocarbon adsorbed in unit cell volume, lb-mol/cu.ft
    init.nConnes = nConnes;                 % Total number of connection
    init.trans_split = trans_split;         % Flag for transmissibility format: true(1): split; false(0): bulk;
    init.conne_type = conne_type;           % Connection type
    init.conne_cell1 = conne_cell1;         % Cell 1 name in each connection pair
    init.conne_cell2 = conne_cell2;         % Cell 2 name in each connection pair
    init.conne_ctrans = conne_ctrans;       % Connection transmissibility for Darcy flow, md-ft;
    init.conne_geo = conne_geo;             % Pair of geometric part for transmissiblity, A/L
    init.conne_perm = conne_perm;           % Pair of permeability part for transmissibility, coef * K 
    init.conne_dz = zeros(nConnes, 1);      % Depth difference in each connection pair, ft
    
    % Flash check point
    init.res_vle_chk = zeros(nCells, nComps+2);     % Reservoir VLE check table: [k-factor at Pb, Pb, key component composition]
    init.wel_vle_chk = zeros(nWells, nComps+2);     % Well VLE check table: [k-factor at Pb, Pb, key component composition] 

    % For multiple porosity model
    init.cell_group = [];                   % A struct to classify cells by type
        
    % Initial fluid in place
    init.OOIP = 0.0;                        % Initial oil in place, STB 
    init.OGIP = 0.0;                        % Initial gas in place, MSCF
    init.OAGIP = 0.0;                       % Initial adsorbed gas in place, MSCF
    init.OWIP = 0.0;                        % Initial water in place, STB
    
    % fluid struct
    fluid.fluidtype = fluidtype;            % For reservoir fluid type flag    
    fluid.eos = eos;                        % EOS model: 'PR' - Peng Robinson model
    fluid.nComps = nComps;                  % Number of hydrocarbon components
    fluid.comp_name = comp_name;            % Component name
    fluid.comp_MW = comp_MW;                % Component molecular weight
    fluid.comp_ACF = comp_ACF;              % Component acentric factor
    fluid.comp_Pc  = comp_Pc ;              % Component critical pressure, psia
    fluid.comp_Tc  = comp_Tc ;              % Component critical temperature, F
    fluid.comp_Zc  = comp_Zc ;              % Component critical compressibility
    fluid.comp_PARA = comp_PARA;            % Component parachor factor
    fluid.comp_SSHIFT = comp_SSHIFT;        % Component volume shift coefficients
    fluid.comp_BIC = comp_BIC;              % Component binary interaction coefficients
    fluid.LBC_Coeffs = LBC_Coeffs;          % Lorentz-Bray-Clark viscosity coefficients 
    fluid.Psc = Psc;                        % Pressure under standard condition, psia
    fluid.Tsc  = Tsc ;                      % Temperature under standard condition, F
    fluid.comp_dmc  = zeros(1, nComps);     % Component molecular collision diameter, meter;

    fluid.Pref_w = Pref_w;                  % Water reference pressure, psia
    fluid.MDen_w_sc  = MDen_w_sc ;          % Water density under standard condition, lb/ft^3
    fluid.Cw = Cw;                          % Water compressibility factor, 1/psia
    fluid.Bw_ref = Bw_ref;                  % Water formation volume factor at reference pressure
    fluid.Vsic_w_ref = Vsic_w_ref;          % Water viscosity at reference pressure, cp
    fluid.Cvw = Cvw;                        % Water viscosibility, 1/psia    
    fluid.wat_MW = 18.015d0;                % Water molecular weight, lb-mol/lb;
    fluid.SMolarVolWat_sc = fluid.wat_MW / fluid.MDen_w_sc; % Specific molar volume of water at surface condition

    % Table look-up for K-factor and fv initial guess when key-component
    % doesn't change too much
    fluid.Flash_Table = [];                 % For gas reservoir reservoir, it is empty
    fluid.key_comp = 1;                     % Fluid key component for saturation point check
    
    % res struct 
    res.nRockTypes = nRockTypes;            % Total number of rock types in the reservoir
    res.Cr = Cr;                            % Rock compressibility of each rock type, 1/psia
    res.P_ref = P_ref;                      % Reference pressure of each rock type, psia
    res.P = P;                              % Intial pressure of each rock type, psia
    res.T = T;                              % Initial temperature of each rock type, F
    res.Sw = Sw;                            % Initial water saturation of each rock type
    res.Zi  = Zi;                           % Initial composition of each rock type
    	
    % For gas slippage
    res.GKap = GKap;                        % Gas apparent permeability based on Civan's model
    res.dpore = dpore;                      % Pore diameter, meter
    
    % For general Fickian diffusion
    res.Fickian = Fickian;
    res.GFD = GFD;
    
    % For desorption
    res.Sorption = Sorption;                % Adsorption model in each rock type
    res.RMDen = RMDen;                      % Bulk rock density in each rock type, lb/cu.ft
    res.VLi = VLi;                          % Langmuir volume in each rock type, scf/lb
    res.PLi = PLi;                          % Langmuir pressure in each rock type, psia
                
    % rockfluid struct 
    rockfluid = rf;
    
    % sim_ctr struct
    sim_ctr.method = method;                % Simulation method: 'FIM' - fully implicit method (default); 'IMPEM' - implicit pressure explicit mass (not implimented);
    sim_ctr.solver = solver;                % Matrix solver option: 1 - direct solver; 2 - CG with preconditioner; 3 - GMRES with preconditioner; 4 - BiCGSTAB with preconditioner
    sim_ctr.t_total = t_total;              % Total simulation time, days
    sim_ctr.dt_init = dt_init;              % Initial time step size, days
    sim_ctr.dt_max = dt_max;                % Maximum time step size, days
    sim_ctr.dt_min = dt_min;                % Minimum time step size, days
    sim_ctr.dt_mul = dt_mul;                % Time step size multiplier, default: 1.5
    sim_ctr.iter_max = iter_max;            % Maximum number of newton iteration, default, 10
    sim_ctr.dt_back = dt_back;              % Time step size backward exponent when phase appearace/disappearance, well control switch
    sim_ctr.dbhp_thresh = dbhp_thresh;      % Threshold to reduce time step size for bottom hole pressure drop, psia, default: 100
    sim_ctr.dpp_thresh = dpp_thresh;        % Threshold to reduce time step size for pressure when close to saturation pressure, psia , default: 100     
    sim_ctr.tol = tol;                      % Pressure tolerance, psia
    sim_ctr.tol_fi = 0.001;                 % Fugacity tolerance 
    sim_ctr.rel_perturb = sqrt(eps('double')); % Relative perturbation size
    sim_ctr.dt_array = zeros(1000, 1);      % Array to store 1000 different timestep sizes
    
    % well struct
    well.nWells = nWells;                   % Total number of wells in the reservoir
    well.name = name;                       % Well names (strings)
    well.Active_Flag = Active_Flag;         % Active well flag
    well.Prod_Flag = Prod_Flag;             % Producer flag
    well.Inje_Wat_Flag = Inje_Wat_Flag;     % Water injector flag
    well.Inje_Gas_Flag = Inje_Gas_Flag;     % Gas/Oil injector flag
    well.Target_min_BHP = Target_min_BHP;   % Minimum BHP for producer, psia
    well.Target_max_BHP = Target_max_BHP;   % Maximum BHP for injector, psia
    well.TargetLiquidRate = TargetLiquidRate; % Target liquid rate, STB/Day
    well.TargetOilRate = TargetOilRate;     % Target oil rate, STB/Day
    well.TargetGasRate = TargetGasRate;     % Target gas rate, SCF/Day
    well.TargetWaterRate = TargetWaterRate; % Target water rate, STB/Day
    well.Gas_Comp = Well_Gas_Comp;          % Composition for gas injector
    well.Ref_Depth = Ref_Depth;             % Reference depth for each well, ft
    well.nperf_cell = nperf_cell;           % Number of perforated cells in each well
    well.wconne_well = wconne_well;                   % Well ID connected to perforated cell
    well.wconne_cell = wconne_cell;                   % Cell ID perfoated by wconne_well 
    well.Flag_Perf = Flag_Perf;                       % Flag for perforated cells 
    well.Prod_Flag_Perf = Prod_Flag_Perf;             % Flag for producer perforated cells
    well.Inje_Wat_Flag_Perf = Inje_Wat_Flag_Perf;     % Flag for water injection perforated cells
    well.Inje_Gas_Flag_Perf = Inje_Gas_Flag_Perf;     % Flag for gas injection perforated cells 
    well.Ind_geom_Perf = Ind_geom_Perf;               % Well index for each perforated cells
    well.Depth_Perf = zeros(sum(nperf_cell), 1);      % Perforated cell center depth, ft
    well.Inj_Gas_cells   = [];                        % Gas injection perforation cells
    
    % Extra
    % (1) Unit convert for liquid rate: STB/D to SCF/D
    vol_cof = 5.61458144;
    well.TargetOilRate = well.TargetOilRate .* vol_cof;
    well.TargetLiquidRate = well.TargetLiquidRate .* vol_cof;
    well.TargetWaterRate = well.TargetWaterRate .* vol_cof;
   
end    


%% Function to scan each line in a file and split it by space
function [tline, str1] = scanline(file_id)
   
    % file_id: file ID
    % tline: string array to store the line information
    line = fgets(file_id); 
    tline = strsplit(strtrim(line));
    str1 = tline{1};
    
end


%% Scan unstructured grid data file
function [nCells, cell_name, cell_type, cell_poro, cell_volume, cell_depth] =  Scan_Grid(fid)

    %% Get grid file path
    [~, temppath] = scanline(fid);
    gridpath = strcat(pwd, temppath);                

    %% Load ASCII data
    grid_data = load(gridpath); 
    
    % Get length of arrays
    nCells = size(grid_data,1);               % Number of reservoir cells                                  

    %% Initialization
    cell_name = zeros([nCells, 1]);           % Cell name: record flag for extending the node list in the future
    cell_type = zeros([nCells, 1]);           % Cell type: integer 
    cell_poro = zeros([nCells, 1]);           % Initial porosity of each cell  
    cell_volume = zeros([nCells, 1]);         % Cell bulk volume
    cell_depth = zeros([nCells, 1]);          % Cell center depth
    
    % grid related data
    cell_name(:, 1) = grid_data(:, 1);
    cell_type(:, 1) = grid_data(:, 2);
    cell_poro(:, 1) = grid_data(:, 3);
    cell_volume(:, 1) = grid_data(:, 4);
    cell_depth(:, 1) = grid_data(:, 5);  
    
end


%% Scan unstructured connection data file
function [trans_split, nConnes, conne_cell1, conne_cell2, conne_type, conne_ctrans, conne_geo, conne_perm] = Scan_Connection(fid)

    %% Get connection file path
    [tline, str1] = scanline(fid);
    connectionpath = strcat(pwd, tline{2});                

    %% Trasmissibility format
    % Connection format: T_AB: default, or (T_A, T_B): split format, good for dynamic changing Transmissibility
    trans_split = false;
    if strcmp(str1, 'SPLIT')
        trans_split = true;
    end
    if strcmp(str1, 'NONSPLIT') 
        trans_split = false;
    end
       
    %% Load ASCII data
    conne_data = load(connectionpath); 
    
    % Get length of arrays
    nConnes = size(conne_data, 1);                % Number of cell-cell connection   

    %% Initialization 
    conne_cell1 = zeros([nConnes, 1]);            % Cell 1 in the connection: integer
    conne_cell2 = zeros([nConnes, 1]);            % Cell 2 in the connection: integer
    conne_type = zeros([nConnes, 1]);             % Connection type: Flag for future use
    
    if trans_split
        conne_ctrans = [];                        % Convection (Darcy flow) transmissibility, split format (null)   
        conne_geo = zeros([nConnes, 2]);          % Convection (Darcy flow) geometric part of transmissibility, split format
        conne_perm = zeros([nConnes, 2]);         % Darcy permeability pair (coef * K) for two connected grids, split format   
    else
        conne_ctrans = zeros([nConnes, 1]);       % Convection (Darcy flow) transmissibility, non-split format
        conne_geo = [];
        conne_perm = [];
    end
    
    % Connection related data     
    conne_cell1(:, 1)  = conne_data(:, 1);
    conne_cell2(:, 1) = conne_data(:, 2);
    conne_type(:, 1)   = conne_data(:, 3);
    if trans_split
        conne_geo(:, 1:2) = conne_data(:, 4:5);
        conne_perm(:, 1:2) = conne_data(:, 6:7);
    else
        conne_ctrans(:, 1) = conne_data(:, 4);
    end
    
end


%% Scan PVT parameters
function [fluidtype, eos, nComps, comp_name, comp_MW, comp_ACF, comp_Pc, comp_Tc, comp_Zc, ...
    comp_PARA, comp_SSHIFT, comp_BIC, LBC_Coeffs, Psc, Tsc] = Scan_PVT(fid)

    fluidtype = 'VLE';               % Fluid type: either VLE or GAS, GAS: no need to flash 
    eos = 'PR';                      % EOS model: Default is Peng-Robinson

    % Scan fluid type 
    [tline, str1] = scanline(fid);
    if strcmp(str1, 'FLUID')   
        fluidtype = tline{2};
    end
    
    % Scan EOS model 
    [tline, str1] = scanline(fid);
    if strcmp(str1, 'MODEL')   
        eos = tline{2};
    end
    % Scan number of component
    [tline, str1] = scanline(fid);
    if strcmp(str1, 'NC')
        nComps = str2double(tline{2}); 
    end 

    %% Initialization
    comp_name = cell(1, nComps);     % Component names (CompNames)
    comp_MW = zeros(1, nComps);      % MW of each component (lbm/lbmol)
    comp_ACF = zeros(1, nComps);     % Acentric factor 
    comp_Pc = zeros(1, nComps);      % Pcrit of each component (psia) 
    comp_Tc = zeros(1, nComps);      % Tcrit of each component (°F)
    comp_Zc = zeros(1, nComps);      % Zcrit of each component (dimless, Zc = Vmc*Pc/(R*Tc))
    comp_PARA = zeros(1, nComps);    % Parachor for each component
    comp_SSHIFT = zeros(1, nComps);  % Volume shift (dimless,  s_i = c_i/b_i;       Vm' = Vm - c;      c = sum(c_i * z_i)    )
                                     % Peneloux "A CONSISTENT CORRECTION FOR SRK VOLUMES". Fluid Phase Equilibria, 8 (1982) 7-23
    comp_BIC = zeros(nComps, nComps);% Binary Interaction Coefficients  
    LBC_Coeffs = zeros(1, 5);        % Lorentz-Bray-Clark Viscosity Correlation Coefficients 
    Psc = 0.0;                       % Standard pressure, psia
    Tsc = 0.0;                       % Standard temperature, °F                                   

    module_end = false;
    
    while not(module_end)

        % Read line by line 
        [~, str1] = scanline(fid);
            
        switch str1
            
            % Scan component name
            case 'COMPNAME'
                [tline, ~] = scanline(fid);
                comp_name(1:nComps) = tline(1:nComps);
                module_end = false;
           
            % Scan Lorentz-Bray-Clark Viscosity Correlation Coefficients
            case 'LBC_VISCOEFF'
                [tline, ~] = scanline(fid);
                LBC_Coeffs(1:5) = str2double(tline(1:5));
                module_end = false;                  

            % Scan molecular weight
            case 'MW'
                [tline, ~] = scanline(fid);
                comp_MW(1:nComps) = str2double(tline(1:nComps));
                module_end = false;
                
            % Scan acentric factor 
            case 'AC'
                [tline, ~] = scanline(fid);
                comp_ACF(1:nComps) = str2double(tline(1:nComps));
                module_end = false;
                                
            % Scan critial pressure: psia
            case 'PCRIT'
                [tline, ~] = scanline(fid);
                comp_Pc(1:nComps) = str2double(tline(1:nComps));
                module_end = false;
                
            % Scan critical temperature (°F)
            case 'TCRIT'
                [tline, ~] = scanline(fid);
                comp_Tc(1:nComps) = str2double(tline(1:nComps));
                module_end = false;
                
            % Scan critical compressibility factor 
            case 'ZCRIT'
                [tline, ~] = scanline(fid);
                comp_Zc(1:nComps) = str2double(tline(1:nComps));
                module_end = false;        
                
            % Scan parachor factor for each component 
            case 'PCHOR'
                [tline, ~] = scanline(fid);
                comp_PARA(1:nComps) = str2double(tline(1:nComps));
                module_end = false;
                
            % Scan volume shift factor for each component 
            case 'VSHIF'
                [tline, ~] = scanline(fid);
                comp_SSHIFT(1:nComps) = str2double(tline(1:nComps));
                module_end = false;   
                
            % Scan standard condition, pressure (psia) and temperature (F) 
            case 'STD'
                [tline, ~] = scanline(fid);
                Psc = str2double(tline{1});              % Standard pressure, psia
                Tsc = str2double(tline{2});              % Standard temperature, F                          
                module_end = false;
 
            % Scan binary interaction coefficient 
            case 'BIN'
                for i = 1 : nComps
                    [tline, ~] = scanline(fid);
                    comp_BIC(i, 1:nComps) = str2double(tline(1:nComps)); 
                end
                module_end = false;            
                
            % Scan the end of the module
            case 'ENDPVT'
                module_end = true;

            % Scan any line before end sign that are non-keywords (comment)
            otherwise 
                continue
                
        end
        
    end
    
end


% Scan water properties data 
function [Pref_w, MDen_w_sc, Cw, Bw_ref, Vsic_w_ref, Cvw] = Scan_Water(fid)

    % Initialization
    Pref_w = 14.7;      % Water reference pressure, psia
    MDen_w_sc = 0.0;    % Water density under standard condition, lb/cu.ft
    Cw = 0.0;           % Water compressibility, 1/psi 
    Bw_ref = 1.0;       % Water formation volume factor at reference pressure
    Vsic_w_ref = 0.31;  % Water viscosity at reference pressure, cp
    Cvw = 0.0;          % Water viscosibility, 1/psia

    module_end = false;    
    
    while not(module_end)

        % Read line by line 
        [tline, str1] = scanline(fid);
            
        switch str1
           
            case 'PREFW'
                Pref_w = str2double(tline{2}); 
                module_end = false;

            case 'RHOWSC'
                MDen_w_sc = str2double(tline{2}); 
                module_end = false;
                
            case 'CREFW'
                Cw = str2double(tline{2});                 
                module_end = false;

            case 'FVFREF'
                Bw_ref = str2double(tline{2});                 
                module_end = false;
                
            case 'VREFW'
                Vsic_w_ref = str2double(tline{2});                 
                module_end = false;
                
            case 'CVREFW'
                Cvw = str2double(tline{2});                 
                module_end = false;
                
            case 'ENDWATER'
                module_end = true;
                
            otherwise               
                continue
                
        end
        
    end
    
end


% Scan initial data 
function [Cr, P_ref, P, T, Sw, Zi, GKap, dpore, Fickian, GFD, Sorption, RMDen, VLi, PLi] = Scan_InitialData(fid, nRockTypes, nComps)
    
    % Initialization    
    Cr = zeros(nRockTypes, 1);          % Rock compressibility of different rock types 
    P_ref = zeros(nRockTypes, 1);       % Reference pressure, psia 
    P = zeros(nRockTypes, 1);           % Initial reservoir pressure, psia
    T = zeros(nRockTypes, 1);           % Reservoir temperature, F 
    Sw = zeros(nRockTypes, 1);          % Initial water saturation, fraction
    Zi = zeros(nRockTypes, nComps);     % Initial composition in each rock type
    
    % Control for gas apparent permeability consideration
    GKap = zeros(nRockTypes, 1);        % Gas apparent permeability considering slippage effect
                                        % 0: No slippage
                                        % 1: Civan's model
                                        
    dpore = zeros(nRockTypes, 1);       % Pore diameter, meter
    
    % Control for General Fick Diffusion coefficient: ft^2/Day 
    Fickian = zeros(nRockTypes, 1);      % Flag for diffusion control
                                        % 0: No diffusion 
                                        % 1: General Fickian Diffusion
    GFD = zeros(nRockTypes, nComps-1, nComps-1);
        
    % Control for adsorption/desorption
    Sorption = cell(nRockTypes, 1);     % Adsorption model used in each rock type       
    RMDen = zeros(nRockTypes, 1);       % Rock density for each rock type, lb/scf
    VLi = zeros(nRockTypes, nComps);    % Adsorption (Langmuir) volume for each component, scf/lb
    PLi = zeros(nRockTypes, nComps);    % Adsorption (Langmuir) pressure for each component, psia
        
    % Initialization 
    for iRT = 1 : nRockTypes 
        GKap(iRT) = 0;                  % Default there is no gas slippage 
        Fickian(iRT) = 0;             % Default there is no Fickian diffusion 
        Sorption{iRT, 1} = 'NA';        % Default there is no adsorption
    end
    
    dpore = 1.0d-6;                     % Default pore diameter is 1e-6 meter (1 micron)
       
    % Loop through each rock type
    for iRT = 1 : nRockTypes 
        
        module_end = false;
        [tline, str1] = scanline(fid);
        if strcmp(str1, 'RPT')
            cRT = str2double(tline(2));  % Current rock type
        else
            cRT = iRT;
        end
                
        while not(module_end)

            % Read line by line 
            [tline, str1] = scanline(fid);

            switch str1
                                                
                case 'CMPROCK'
                    Cr(cRT) = str2double(tline(2));                              
                    module_end = false;

                case 'PREF'
                    P_ref(cRT) = str2double(tline(2));                                              
                    module_end = false;

                case 'PRESINIT'
                    P(cRT) = str2double(tline(2));               
                    module_end = false;

                case 'TRES'
                    T(cRT) = str2double(tline(2));
                    module_end = false;

                case 'SWINT'
                    Sw(cRT) = str2double(tline(2));                                                              
                    module_end = false;

                case 'ZIINT'
                    [tline, ~] = scanline(fid);
                    Zi(cRT, 1 : nComps) = str2double(tline(1 : nComps));
                    module_end = false; 

                case 'GKAPP'
                    if(strcmp(tline(2), 'CIVAN'))
                        GKap(cRT) = 1;
                    end
                    module_end = false; 
                    
                case 'PORD'
                    dpore(cRT) = str2double(tline(2)) * 1.0d-9; % To meter                              
                    module_end = false; 
                    
                case 'GFICK'
                    Fickian(iRT) = 1;
                    % Dimension: (nc-1) x (nc-1)
                    for iComp = 1 : nComps - 1 
                        [tline, ~] = scanline(fid);
                        GFD(cRT, iComp, 1 : nComps-1) = str2double(tline(1 : nComps-1));
                    end
                    module_end = false; 
                                        
                case 'SORPTION'
                    Sorption(cRT) = tline(2);
                    module_end = false;
                    
                case 'RMDEN'
                    RMDen(cRT) = str2double(tline(2)); 
                    module_end = false;
                    
                case 'VLI'
                    [tline, ~] = scanline(fid);
                    VLi(cRT, 1 : nComps) = str2double(tline(1 : nComps));
                    module_end = false;
                    
                case 'PLI' 
                    [tline, ~] = scanline(fid);
                    PLi(cRT, 1 : nComps) = str2double(tline(1 : nComps));             
                    module_end = false;
                                  
                case 'ENDRPT'
                    module_end = true;

                case 'ENDINIT'
                    module_end = true;

                otherwise
                    continue

            end

        end
        
    end
    
end


% Scan rock-fluid table data 
function [rf] = Scan_RockFluidData(fid, nRockTypes)

    % Initialization 
    dummy.method = 0;  % Relative permeability model 
    dummy.Swco = 0;    % Connate water saturation 
    dummy.kro_cw = 0;
    dummy.SWFN = [];
    dummy.SGFN = [];
    dummy.SOF3 = [];
    
    rf = [];
    
    % Define multiple rock types
    for iRT = 1 : nRockTypes
        rf = [rf; dummy];
    end

    % Loop through each rock type
    for iRT = 1 : nRockTypes 
        
        [tline, str1] = scanline(fid);
        if strcmp(str1, 'RPT')
            cRT = str2double(tline(2));  % Current rock type
        else
            cRT = iRT;
        end    
            
        module_end = false;
    
        while not(module_end)

            % Read line by line 
            [tline, str1] = scanline(fid);

            switch str1
                
                case 'MODEL'
                    rf(cRT).method = str2double(tline{2}); 
                    
                case 'SWFN'
                    len_SWFN_tab = int8(str2double(tline{2}));
                    rf(cRT).SWFN = zeros(len_SWFN_tab, 3);
                    for i = 1 : len_SWFN_tab
                        [tline, ~] = scanline(fid);
                        rf(cRT).SWFN(i, 1:3) = str2double(tline(1:3));
                    end
                    module_end = false;                

                case 'SGFN'
                    len_SGFN_tab = int8(str2double(tline{2}));
                    rf(cRT).SGFN = zeros(len_SGFN_tab, 3);
                    for i = 1 : len_SGFN_tab
                        [tline, ~] = scanline(fid);
                        rf(cRT).SGFN(i, 1:3) = str2double(tline(1:3)); 
                    end
                    module_end = false;    

                case 'SOF3'
                    len_SOF3_tab = int8(str2double(tline{2}));
                    rf(cRT).SOF3 = zeros(len_SOF3_tab, 3);
                    for i = 1 : len_SOF3_tab
                        [tline, ~] = scanline(fid);
                        rf(cRT).SOF3(i, 1:3) = str2double(tline(1:3));  
                    end
                    module_end = false;                    

                case 'ENDRPT'
                    module_end = true;                    
                    
                case 'ENDROCKFLUID'
                    module_end = true;

                otherwise
                    continue;  
            end
        end % End keyword reading 
        
    end % End rock type loop

end
  

% Scan numerical control data 
function [method, solver, t_total, dt_init, dt_max, dt_min, dt_mul, dt_back, dbhp_thresh, dpp_thresh, tol, iter_max] = Scan_SimCtr(fid)

    % Initialization(default provided)
    method     = 1;       % FIM: 1; IMPEM: 2;
    solver = 1;           % Linear Solver option 
    t_total = 100;        % Time to be simulated, days
    dt_init = 0.5;        % Initial dt for the run, days
    dt_max = 2;           % Maximum allowable delta_t, days
    dt_min = 1E-6;        % Minimum allowable delta_t, days  
    dt_mul = 1.5;         % Multiplier to increase the size of dt if convergence achieved 
    dt_back = 4.0;        % Time step size backward exponent when phase appearace/disappearance, well control switch
    dbhp_thresh = 0.0;    % Threshold to reduce time step size for bottom hole pressure drop, psia
    dpp_thresh = 0.0;     % Threshold to reduce time step size for pressure when close to saturation pressure, psia 
    tol = 0.1;            % Tolerance for residual in newton (mass balance: lb-mol; volume balance: cu.ft)  
    iter_max = 10;        % Max Newton Iterations   
    
    module_end = false;
    
    while not(module_end)
        
        % Read line by line 
        [tline, str1] = scanline(fid);
            
        switch str1
            
            case 'FIM'
                method = 1;
                module_end = false;
                
            case 'IMPEM'
                method = 2;
                module_end = false;

            case 'SOLVER'
                solver = str2double(tline{2});
                module_end = false;
                
            case 'TIME'
                t_total = str2double(tline{2});
                module_end = false;                

            case 'DTINIT'
                dt_init = str2double(tline{2});
                module_end = false;                

            case 'DTMAX'
                dt_max = str2double(tline{2}); 
                module_end = false;
                
            case 'DTMIN'
                dt_min = str2double(tline{2}); 
                module_end = false;                

            case 'DTMULT'
                dt_mul = str2double(tline{2});
                module_end = false;
                
            case 'NETWONMAX'
                iter_max = int8(str2double(tline{2}));
                module_end = false;
                
            case 'DTBACK'
                dt_back = str2double(tline{2});
                module_end = false;
                
            case 'DBHPTHRESH' 
                dbhp_thresh = str2double(tline{2});
                module_end = false;
                
            case 'DPPTHRESH'
                dpp_thresh = str2double(tline{2});
                module_end = false;

            case 'TOL'
                tol = str2double(tline{2});
                module_end = false;
            
            case 'ENDNUMERICAL'
                module_end = true;
                
            otherwise
                continue;
                
        end
        
    end

end


% Scan well control data 
function [name, nperf_cell, Prod_Flag, Inje_Wat_Flag, Inje_Gas_Flag,          ... 
          Target_min_BHP, Target_max_BHP, TargetLiquidRate, TargetOilRate,    ...    
          TargetGasRate, TargetWaterRate, Well_Gas_Comp, Ref_Depth,           ...
          Active_Flag, wconne_well, wconne_cell, Flag_Perf, Prod_Flag_Perf,   ...
          Inje_Wat_Flag_Perf, Inje_Gas_Flag_Perf, Ind_geom_Perf]              ...
          = Scan_WellCtr(fid, nWells, nComps)
      
    % Initialization
    max_nPerf = 1000;                              % Maximum number of perforation in the reservoir
    
    % Size: nWells
    name = cell(nWells, 1);                        % Well name: string 
    nperf_cell = zeros(nWells, 1);                 % Number of perforated cells in each well
    Prod_Flag = zeros(nWells, 1);                  % Flag for producer
    Inje_Wat_Flag = zeros(nWells, 1);              % Flag for water injector 
    Inje_Gas_Flag = zeros(nWells, 1);              % Flag for gas injector
    Target_min_BHP = ones(nWells, 1) .* 14.7;      % Minimum BHP for producer
    Target_max_BHP = ones(nWells, 1) .* 20000;     % Maximum BHP for injector
    TargetLiquidRate = zeros(nWells, 1);           % Surface liquid rate
    TargetOilRate = zeros(nWells, 1);              % Surface oil rate 
    TargetGasRate = zeros(nWells, 1);              % Surface gas rate 
    TargetWaterRate = zeros(nWells, 1);            % Surface water rate 
    Well_Gas_Comp = zeros(nWells, nComps);         % Injector's gas composition
    Ref_Depth = zeros(nWells, 1);                  % Reference depth for BHP in each well
    Active_Flag = zeros(nWells, 1);                % Flag for well status
    
    % Counter for all perforations from different wells 
    it_perf = 0;
    
    for iwell = 1 : nWells
        
        module_end = false;

        while not(module_end)        

            % Read line by line 
            [tline, str1] = scanline(fid);

            switch str1
                
                case 'MAXPERF'
                    max_nPerf = str2double(tline{2});
                    module_end = false;
                    
                    % Allocate memory for perforated cell related properties
                    wconne_well_0 = zeros(max_nPerf, 1);                   % Well ID perforating a specific cell
                    wconne_cell_0 = zeros(max_nPerf, 1);                   % Cell ID corresponding perforation Welll ID
                    Flag_Perf_0 = zeros(max_nPerf, 1);                     % Flag for perforated cells 
                    Prod_Flag_Perf_0 = zeros(max_nPerf, 1);                % Flag for producer perforated cells
                    Inje_Wat_Flag_Perf_0 = zeros(max_nPerf, 1);            % Flag for water injection perforated cells
                    Inje_Gas_Flag_Perf_0 = zeros(max_nPerf, 1);            % Flag for gas injection perforated cells 
                    Ind_geom_Perf_0 = zeros(max_nPerf, 1);                 % Well index for each perforated cells
                      
                case 'PRODUCER'
                    Prod_Flag(iwell) = 1;
                    name{iwell, 1} = tline{2};
                    if strcmp(tline{3}, 'OPEN')
                        Active_Flag(iwell, 1) = 1;
                    end
                    if strcmp(tline{3}, 'SHUTIN')
                        Active_Flag(iwell, 1) = 0;
                    end                    
                    module_end = false;

                case 'INJECTOR'
                    name{iwell, 1} = tline{2};
                    if strcmp(tline{3}, 'OPEN')
                        Active_Flag(iwell, 1) = 1;
                    end
                    if strcmp(tline{3}, 'SHUTIN')
                        Active_Flag(iwell, 1) = 0;
                    end                    
                    module_end = false;
                    
                case 'OILSTB'
                    TargetOilRate(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'GASSCF'
                    TargetGasRate(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'WATSTB'
                    TargetWaterRate(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'LIQSTB'
                    TargetLiquidRate(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'INJCOM'
                    for iComp = 1 : nComps
                        Well_Gas_Comp(iwell, iComp) = str2double(tline{iComp+1});
                    end
                    module_end = false;
                    
                case 'BHPMIN'
                    Target_min_BHP(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'BHPMAX'
                    Target_max_BHP(iwell, 1) = str2double(tline{2});
                    module_end = false;
                                        
                case 'REFDEPTH'
                    Ref_Depth(iwell, 1) = str2double(tline{2});
                    module_end = false;
                    
                case 'PERF'
                    % Number of perforated cell for this well
                    nPerf = int8(str2double(tline{2}));
                    nperf_cell(iwell) = nPerf;
                    
                    for iPerf = 1 : nPerf
                        
                        % Counter perforations globally (in well loop but
                        % not perforation loop)
                        
                        it_perf = it_perf + 1;
                        
                        % Before read 1st perforation, check if perforation 
                        % cell array is empty, allocate memory
                        if isempty(wconne_well_0) && it_perf == 1
                            % Allocate memory for perforated cell related properties
                            wconne_well_0 = zeros(max_nPerf, 1);                   % Well ID perforating a specific cell
                            wconne_cell_0 = zeros(max_nPerf, 1);                   % Cell ID corresponding perforation Welll ID
                            Flag_Perf_0 = zeros(max_nPerf, 1);                     % Flag for perforated cells 
                            Prod_Flag_Perf_0 = zeros(max_nPerf, 1);                % Flag for producer perforated cells
                            Inje_Wat_Flag_Perf_0 = zeros(max_nPerf, 1);            % Flag for water injection perforated cells
                            Inje_Gas_Flag_Perf_0 = zeros(max_nPerf, 1);            % Flag for gas injection perforated cells 
                            Ind_geom_Perf_0 = zeros(max_nPerf, 1);                 % Well index for each perforated cells                            
                        end
                        
                        [tline, ~] = scanline(fid);
                        
                        perfCell = str2double(tline{1});
                        wellIndex = str2double(tline{2});
                        
                        % Assignment
                        wconne_well_0(it_perf, 1) = iwell;
                        wconne_cell_0(it_perf, 1) = perfCell;
                        Ind_geom_Perf_0(it_perf, 1) = wellIndex;                        
                        Flag_Perf_0(it_perf, 1) = 1;
                        
                        if Prod_Flag(iwell) == 1
                            Prod_Flag_Perf_0(it_perf, 1) = 1;
                        end
                                                 
                        % If injector, check if cell is injecting water or
                        % gas Perf_Center_Depth3D
                        if Prod_Flag(iwell) == 0
                            
                            comp_flag = ( sum(Well_Gas_Comp(iwell, :), 2) == 1.0 && ...   % Injecting composition specified
                                          TargetGasRate(iwell, 1) ~= 0           || ...   % Gas or oil rate specified 
                                          TargetOilRate(iwell, 1) ~= 0 );
                            water_flag = (TargetWaterRate(iwell, 1) ~= 0);                % Water rate specified 
                            
                            Inje_Wat_Flag(iwell, 1) = water_flag; 
                            Inje_Gas_Flag(iwell, 1) = comp_flag;  
                            
                            Inje_Wat_Flag_Perf_0(it_perf, 1) = water_flag;
                            Inje_Gas_Flag_Perf_0(it_perf, 1) = comp_flag;
                            
                        end
                                                
                    end
                    
                                        
                case 'ENDWELL'
                    module_end = true;

                otherwise 
                    continue;

            end
        end  % Scanning loop
        
    end % Well loop 
    
    % Shrink array size
    sum_nPerf = sum(nperf_cell);    
    wconne_well = wconne_well_0(1 : sum_nPerf, 1);
    wconne_cell = wconne_cell_0(1 : sum_nPerf, 1);    
    Flag_Perf = Flag_Perf_0(1 : sum_nPerf, 1);                     
    Prod_Flag_Perf = Prod_Flag_Perf_0(1 : sum_nPerf, 1);
    Inje_Wat_Flag_Perf = Inje_Wat_Flag_Perf_0(1 : sum_nPerf, 1);            
    Inje_Gas_Flag_Perf = Inje_Gas_Flag_Perf_0(1 : sum_nPerf, 1);            
    Ind_geom_Perf = Ind_geom_Perf_0(1 : sum_nPerf, 1);                 

end
