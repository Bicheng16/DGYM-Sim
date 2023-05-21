function [vol_total_error, Balance_vol] = Volumetric_Check(prop, init, stat)

    % This subroutine is to check if totoal fluid volume is equal to pore volume
    
    % Water volume in each cell, ft^3
    vol_water = 0.0 .* prop(stat).cell_P;
    vol_hc    = 0.0 .* vol_water;
    vol_cell  = 0.0 .* vol_water;
    Balance_vol = zeros([1, 2]);
    
    % Total amount of moles of hydrocarbon: moles/pore volume
    N_hc = sum(prop(stat).cell_Ni, 2);
    
    % Calculate water volume, ft^3
    vol_water(:) = prop(stat).cell_Nw(:) .* init.cell_volume(:) ./  prop(stat).cell_MDen_w(:); 
                           
    % Calculate hydrocarbon volume: ft^3
    vol_hc(:)    = N_hc(:) .* init.cell_volume(:) ./  prop(stat).cell_MolarDensHC(:);
            
    % Calculate cell pore volume: ft^3
    vol_cell(:)  = prop(stat).cell_poro(:) .* init.cell_volume(:); 
    
    % Calculate volume error in each cell: ft^3
    vol_error    = vol_cell - (vol_water + vol_hc);
    vol_total_error = sum(vol_cell) - sum(vol_water + vol_hc);
    
    % L-2 norm and infinite norm of volume balance
    Balance_vol(1,1) = norm(vol_error(:, 1), 2);               
    Balance_vol(1,2) = norm(vol_error(:, 1), inf);
    
end