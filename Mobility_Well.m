% Function calculate well mobility and injector total mobility
function [WI_o, WI_w, WI_g] = Mobility_Well(well, prop)
    
    % Define temporary variables
    next = 2;
    
    % Mobility modification for injectors (ECL 2012.2 Manual "Injecting connections" pg 1269)
    Oil_Mob = prop(next).cell_kro(well.wconne_cell) ./ prop(next).cell_Vsic_o(well.wconne_cell);         
    Gas_Mob = prop(next).cell_krg(well.wconne_cell) ./ prop(next).cell_Vsic_g(well.wconne_cell);        
    Wat_Mob = prop(next).cell_krw(well.wconne_cell) ./ prop(next).cell_Vsic_w(well.wconne_cell);        
    
    Oil_Mob(isinf(Oil_Mob) | isnan(Oil_Mob)) = 0;
    Gas_Mob(isinf(Gas_Mob) | isnan(Gas_Mob)) = 0;
    Wat_Mob(isinf(Wat_Mob) | isnan(Wat_Mob)) = 0;
        
    % Calculate well productivity index 
    WI_o = well.Prod_Flag_Perf .* well.Ind_geom_Perf .* Oil_Mob .* prop(next).cell_MolarDensLiq(well.wconne_cell);    
    WI_g = well.Prod_Flag_Perf .* well.Ind_geom_Perf .* Gas_Mob .* prop(next).cell_MolarDensVap(well.wconne_cell);
    WI_w = well.Prod_Flag_Perf .* well.Ind_geom_Perf .* Wat_Mob .* prop(next).cell_MDen_w(well.wconne_cell);
    
    WI_o(isinf(WI_o) | isnan(WI_o)) = 0;
    WI_g(isinf(WI_g) | isnan(WI_g)) = 0;
    WI_w(isinf(WI_w) | isnan(WI_w)) = 0;
         
end