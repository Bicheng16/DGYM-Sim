% Clear newton-level secondary variables after solve linear system AX = b
function [prop] = zero_dynamic_properties(prop, nCells, nConnes, nComps, nPerfs)

    % Null all secondary variables, reserve values for primary variables.

    % nCells: number of cells in the reservoir 
    % nConnes: number of connections in the reservoir
    % nComps: number of hydrocarbon component in the reservoir

    prop(2).cell_poro = zeros(nCells, 1);
    %prop(2).cell_P = zeros(nCells, 1);
    %prop(2).cell_T = zeros(nCells, 1);
    prop(2).cell_Sw = zeros(nCells, 1);
    prop(2).cell_So = zeros(nCells, 1);
    prop(2).cell_Sg = zeros(nCells, 1);
    prop(2).cell_Xi  = zeros(nCells, nComps);
    prop(2).cell_Yi  = zeros(nCells, nComps);
    prop(2).cell_Zi = zeros(nCells, nComps);
    prop(2).cell_Ki = zeros(nCells, nComps);
    %prop(2).cell_Ni = zeros(nCells, nComps);
    prop(2).cell_No = zeros(nCells, 1);
    prop(2).cell_Ng = zeros(nCells, 1);
    %prop(2).cell_Nw = zeros(nCells, 1);
    prop(2).cell_Vfluid = zeros(nCells, 1);
    prop(2).cell_fv = zeros(nCells, 1);
    prop(2).cell_MDen_o = zeros(nCells, 1);
    prop(2).cell_MDen_g = zeros(nCells, 1);
    prop(2).cell_MDen_w = zeros(nCells, 1);
    prop(2).cell_MolarDensHC = zeros(nCells, 1);
    prop(2).cell_MolarDensLiq = zeros(nCells, 1);
    prop(2).cell_MolarDensVap = zeros(nCells, 1);
    prop(2).cell_MolarDensWat = zeros(nCells, 1);
    prop(2).cell_MolarVolLiq = zeros(nCells, 1);
    prop(2).cell_MolarVolVap = zeros(nCells, 1);
    prop(2).cell_kro  = zeros(nCells, 1);
    prop(2).cell_krw  = zeros(nCells, 1);
    prop(2).cell_krg  = zeros(nCells, 1);
    prop(2).cell_Vsic_o = zeros(nCells, 1);
    prop(2).cell_Vsic_g = zeros(nCells, 1);
    prop(2).cell_Vsic_w = zeros(nCells, 1);
    prop(2).cell_pcow  = zeros(nCells, 1);
    prop(2).cell_pcgo  = zeros(nCells, 1);
    prop(2).cell_PV = zeros(nCells, 1);
    prop(2).cell_HCPV  = zeros(nCells, 1);
    prop(2).conne_ctrans = zeros(nConnes, 1);
    prop(2).cell_Qai= zeros(nCells, nComps);
    prop(2).cell_GKapp = zeros(nCells, nComps);
    prop(2).cell_GFD = zeros(nCells, nComps-1, nComps-1);    
    %prop(2).well_P = zeros(nWells, 1);
    prop(2).perf_P = zeros(nPerfs, 1);
    prop(2).perf_Hw = zeros(nPerfs, 1);
    
        
end