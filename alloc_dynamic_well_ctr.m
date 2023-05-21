function well_ctr = alloc_dynamic_well_ctr(nWells, wconne_well)

    % Allocate memory
    OilRate_Ctr = zeros(nWells, 1);  OilRate_Ctr(:) = false;
    GasRate_Ctr = zeros(nWells, 1);  GasRate_Ctr(:) = false;
    LiqRate_Ctr = zeros(nWells, 1);  LiqRate_Ctr(:) = false;  
    WatRate_Ctr = zeros(nWells, 1);  WatRate_Ctr(:) = false;
    BHP_Ctr     = zeros(nWells, 1);  BHP_Ctr(:) = false;
    
    % Assign to struct
    dummy.OilRate_Ctr = OilRate_Ctr;
    dummy.GasRate_Ctr = GasRate_Ctr;
    dummy.LiqRate_Ctr = LiqRate_Ctr;
    dummy.WatRate_Ctr = WatRate_Ctr;
    dummy.BHP_Ctr = BHP_Ctr;
    dummy.wconne_well = wconne_well;
    
    % Three different levels: Current = 1; Next = 2
    well_ctr = [dummy; dummy; dummy];

end