% Old time step and newton-level step variable definition
function [prop] = alloc_dynamic_properties(nCells, nConnes, nComps, nWells, nPerfs)

    % nCells: number of cells in the reservoir 
    % nConnes: number of connections in the reservoir
    % nComps: number of hydrocarbon component in the reservoir

    dummy.cell_poro = zeros(nCells, 1);
    dummy.cell_P = zeros(nCells, 1);
    dummy.cell_T = zeros(nCells, 1);
    dummy.cell_Sw = zeros(nCells, 1);
    dummy.cell_So = zeros(nCells, 1);
    dummy.cell_Sg = zeros(nCells, 1);
    dummy.cell_Xi  = zeros(nCells, nComps);
    dummy.cell_Yi  = zeros(nCells, nComps);
    dummy.cell_Zi = zeros(nCells, nComps);
    dummy.cell_Ni = zeros(nCells, nComps);
    dummy.cell_Ki = zeros(nCells, nComps);  % For new PVT SUBOUTINE	
    dummy.cell_No = zeros(nCells, 1);
    dummy.cell_Ng = zeros(nCells, 1);
    dummy.cell_Nw = zeros(nCells, 1); 
    dummy.cell_Vfluid = zeros(nCells, 1);
    dummy.cell_fv = zeros(nCells, 1);
    dummy.cell_MDen_o = zeros(nCells, 1);
    dummy.cell_MDen_g = zeros(nCells, 1);
    dummy.cell_MDen_w = zeros(nCells, 1);
    dummy.cell_MolarDensHC = zeros(nCells, 1);
    dummy.cell_MolarDensLiq = zeros(nCells, 1);
    dummy.cell_MolarDensVap = zeros(nCells, 1);
    % Scale water mass balance equation by water molecular weight
    dummy.cell_MolarDensWat = zeros(nCells, 1);
    dummy.cell_MolarVolLiq = zeros(nCells, 1);
    dummy.cell_MolarVolVap = zeros(nCells, 1);
    dummy.cell_kro  = zeros(nCells, 1);
    dummy.cell_krw  = zeros(nCells, 1);
    dummy.cell_krg  = zeros(nCells, 1);
    dummy.cell_Vsic_o = zeros(nCells, 1);
    dummy.cell_Vsic_g = zeros(nCells, 1);
    dummy.cell_Vsic_w = zeros(nCells, 1);
    dummy.cell_pcow  = zeros(nCells, 1);
    dummy.cell_pcgo  = zeros(nCells, 1);
    dummy.cell_PV = zeros(nCells, 1);
    dummy.cell_HCPV  = zeros(nCells, 1);
    %dummy.conne_ctrans = zeros(nConnes, 1);
    
    % Add for adsorption/desorption, lb-mol/cu.ft 
    dummy.cell_Qai = zeros(nCells, nComps);
    
    % Add for gas phase slippage permeability correction, lb-mol/cu.ft
    dummy.cell_GKapp = zeros(nCells, nComps);
    
    % Add for general Fickian diffusion coefficient for gas phase, ft^2/day
    dummy.cell_GFD = zeros(nCells, nComps-1, nComps-1);
    
    % Add for well variables
    dummy.well_P = zeros(nWells, 1);
    dummy.perf_P = zeros(nPerfs, 1);
    dummy.perf_Hw = zeros(nPerfs, 1);
    
    % Properties storage: (Current; Next)
    prop = [dummy; dummy];
    
end