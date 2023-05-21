% Function to allocate memory for result variables

function rst = alloc_rst_var(nComps, nWells)

    % Average values of pressure and saturation
    rst.P_res_HCPVweigthed = 0.0;         % Hydrocarbon pore volume weighted pressure, psia
    rst.P_res_PVweigthed = 0.0;           % Pore volume weighted pressure, psia
    rst.So_PVweigthed = 0.0;              % Pore volume weighted oil saturation 
    rst.Sg_PVweigthed = 0.0;              % Pore volume weighted gas saturation 
    rst.Sw_PVweigthed = 0.0;              % Pore volume weighted water saturation 
    
    % Component in place
    rst.COMPIP = zeros(1, nComps+1);      % Component mass in place: 1:nComps: hydrocarbon; nComps+1: water
    rst.PV     = 0.0;                     % Pore volume, RB
    rst.HCPV   = 0.0;                     % Hydrocarbon pore volume, RB 
    rst.VOil   = 0.0;                     % Oil volme, RB 
    rst.VGas   = 0.0;                     % Gas volume, RB
    rst.VWat   = 0.0;                     % Water volume, RB
        
    % Fluid in place @ Standard condition (current)
    rst.FOIP = 0.0;                       % Oil in place, STB
    rst.FGIP = 0.0;                       % Free gas in place, MSCF
    rst.FAGIP = 0.0;                      % Adsorbed gas in place, MSCF
    rst.FWIP = 0.0;                       % Water in place, STB
    
    rst.nT_Field_Vector = 0.0;
    rst.Zi_Field_vector = zeros(1, nComps);
    
    % Time step size 
    rst.dt = 0;

    % Newton iteration
    rst.newtons = 0;

    % Bottom hole pressure 
    rst.well_P = zeros(nWells, 1);         % Bottom hole pressure of each well;

    % Well rate parameters
    % Producer well rate @ standard condition
    rst.WOPR=zeros(nWells, 1);             % Oil producing rate, STB/day
    rst.WGPR=zeros(nWells, 1);             % Gas producing rate, MSCF/day
    rst.WWPR=zeros(nWells, 1);             % Water producing rate, STB/day
    
    % Injector well rate @ standard condition
    rst.WOIR=zeros(nWells, 1);             % Oil producing rate, STB/day      
    rst.WGIR=zeros(nWells, 1);             % Gas producing rate, MSCF/day      
    rst.WWIR=zeros(nWells, 1);             % Water producing rate, STB/day
    
    % Field rate parameters
    % Field total production rate @ standard condition
    rst.FOPR = 0.0;                       % Gas producing rate, MSCF/day
    rst.FGPR = 0.0;                       % Water producing rate, STB/day 
    rst.FWPR = 0.0;                       % Oil producing rate, STB/day 

    % Field total injection rate @ standard condition
    rst.FOIR = 0.0;                       % Gas producing rate, MSCF/day
    rst.FGIR = 0.0;                       % Water producing rate, STB/day 
    rst.FWIR = 0.0;                       % Oil producing rate, STB/day 

    % Field cumulative production amount @ standard condition: ECL Mannual P1526
    rst.FOPT = 0.0;                       % Total oil production, STB
    rst.FGPT = 0.0;                       % Total gas production, MSCF
    rst.FWPT = 0.0;                       % Total water production, STB       
    
    % Field cumulative injection amount @ standard condition: ECL Mannual P1526
    rst.FOIT = 0.0;                       % Total oil injection, STB
    rst.FGIT = 0.0;                       % Total gas injection, MSCF
    rst.FWIT = 0.0;                       % Total water injection, STB   
    
end