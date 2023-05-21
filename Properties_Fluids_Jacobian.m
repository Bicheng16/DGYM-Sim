function [prop] = Properties_Fluids_Jacobian(stat, prop, fluid, init, dummy, Target_Nodes, inherit, inherit_loc, Perturb_Nw)

    % Inherit data from previous calculation if you are evaluating
    % single grid blocks     
    if inherit && ~isempty(inherit_loc)
        prop(stat).cell_No(inherit_loc, 1) = dummy.No(inherit_loc, 1); 
        prop(stat).cell_Ng(inherit_loc, 1) = dummy.Ng(inherit_loc, 1); 
        prop(stat).cell_fv(inherit_loc, 1) = dummy.fv(inherit_loc, 1); 
        prop(stat).cell_Vsic_o(inherit_loc, 1) = dummy.Vsic_o(inherit_loc, 1); 
        prop(stat).cell_Vsic_g(inherit_loc, 1) = dummy.Vsic_g(inherit_loc, 1);
        prop(stat).cell_Vsic_w(inherit_loc, 1) = dummy.Vsic_w(inherit_loc, 1);
        prop(stat).cell_MolarDensHC(inherit_loc, 1) = dummy.MolarDensHC(inherit_loc, 1); 
        prop(stat).cell_MolarDensLiq(inherit_loc, 1) = dummy.MolarDensLiq(inherit_loc, 1); 
        prop(stat).cell_MolarDensVap(inherit_loc, 1) = dummy.MolarDensVap(inherit_loc, 1);
        prop(stat).cell_MolarDensWat(inherit_loc, 1) = dummy.MolarDensWat(inherit_loc, 1);
        prop(stat).cell_MolarVolLiq(inherit_loc, 1) = dummy.MolarVolLiq(inherit_loc, 1);
        prop(stat).cell_MolarVolVap(inherit_loc, 1) = dummy.MolarVolVap(inherit_loc, 1); 
        prop(stat).cell_MDen_o(inherit_loc, 1) = dummy.MDen_o(inherit_loc, 1); 
        prop(stat).cell_MDen_g(inherit_loc, 1) = dummy.MDen_g(inherit_loc, 1); 
        prop(stat).cell_MDen_w(inherit_loc, 1) = dummy.MDen_w(inherit_loc, 1);
        prop(stat).cell_Xi(inherit_loc, :) = dummy.Xi(inherit_loc, :); 
        prop(stat).cell_Yi(inherit_loc, :) = dummy.Yi(inherit_loc, :);
        prop(stat).cell_Ki(inherit_loc, :) = dummy.Ki(inherit_loc, :);   % For new PVT algorithm       
        prop(stat).cell_Vfluid(inherit_loc, 1) = dummy.Vfluid(inherit_loc, 1);
        prop(stat).cell_So(inherit_loc, 1) = dummy.So(inherit_loc, 1); 
        prop(stat).cell_Sg(inherit_loc, 1) = dummy.Sg(inherit_loc, 1);
        prop(stat).cell_Sw(inherit_loc, 1) = dummy.Sw(inherit_loc, 1);
    end

    % Define temperary variables        
    N_hc = sum(prop(stat).cell_Ni, 2);    
    prop(stat).cell_Zi = prop(stat).cell_Ni ./ N_hc(:, ones(1, fluid.nComps));
        
    % Inherit all reservoir properties from dummy when perturbing well pressure 
    nCells = init.nCells; % size(prop(stat).cell_P, 1);    
    % If fully inherited from dummy, return.
    if nnz(inherit_loc) == nCells
        return
    end
    
    % Define temperary variables 
    nComps = fluid.nComps;
    Vo = zeros(nCells, 1);
    Vg = zeros(nCells, 1);
    Vw = zeros(nCells, 1);
    Vft = zeros(nCells, 1); 
    Xc = zeros(nCells, 1);
    Yc = zeros(nCells, 1);
    cell_Bw_inv = zeros(nCells, 1);
    nNodes = nnz(Target_Nodes);
    kval_0 = zeros(1, nComps);
    
    if Perturb_Nw
    
        for i = 1 : nNodes
            nodeID = Target_Nodes(i); 
            
            prop(stat).cell_No(nodeID, 1) = dummy.No(nodeID, 1); 
            prop(stat).cell_Ng(nodeID, 1) = dummy.Ng(nodeID, 1);             
            prop(stat).cell_Ki(nodeID, :) = dummy.Ki(nodeID, :);       
            prop(stat).cell_fv(nodeID, 1) = dummy.fv(nodeID, 1); 
            prop(stat).cell_Vsic_o(nodeID, 1) = dummy.Vsic_o(nodeID, 1); 
            prop(stat).cell_Vsic_g(nodeID, 1) = dummy.Vsic_g(nodeID, 1);
            prop(stat).cell_MolarDensHC(nodeID, 1) = dummy.MolarDensHC(nodeID, 1); 
            prop(stat).cell_MolarDensLiq(nodeID, 1) = dummy.MolarDensLiq(nodeID, 1); 
            prop(stat).cell_MolarDensVap(nodeID, 1) = dummy.MolarDensVap(nodeID, 1);
            prop(stat).cell_MolarVolLiq(nodeID, 1) = dummy.MolarVolLiq(nodeID, 1);
            prop(stat).cell_MolarVolVap(nodeID, 1) = dummy.MolarVolVap(nodeID, 1); 
            prop(stat).cell_MDen_o(nodeID, 1) = dummy.MDen_o(nodeID, 1); 
            prop(stat).cell_MDen_g(nodeID, 1) = dummy.MDen_g(nodeID, 1); 
            prop(stat).cell_Xi(nodeID, :) = dummy.Xi(nodeID, :); 
            prop(stat).cell_Yi(nodeID, :) = dummy.Yi(nodeID, :);            
                        
        end
        
    else

        for i = 1 : nNodes 

            nodeID = Target_Nodes(i);
            pressure = prop(stat).cell_P(nodeID);
            Zi_comp  = prop(stat).cell_Zi(nodeID, :);
            temperature = prop(stat).cell_T(nodeID);

            % Get bubble point
            Pb = init.res_vle_chk(nodeID, nComps+1);

            % Get K-factor and fv guess
            if pressure > Pb
                fv_0   = 1.0d-5; 
                kval_0 = init.res_vle_chk(nodeID, 1:nComps);
            elseif pressure <= Pb && isequal(fluid.fluidtype, 'GAS') 
                fv_0   = prop(1).cell_fv(nodeID);
                kval_0 = prop(1).cell_Ki(nodeID, :);
            else
                % Table look-up to find K-factor and fv initial guess
                fv_0 = nakeinterp1(fluid.Flash_Table(:, 1), fluid.Flash_Table(:, 2), pressure);
                for iComp = 1 : nComps
                    kval_0(1, iComp) = nakeinterp1(fluid.Flash_Table(:, 1), fluid.Flash_Table(:, iComp+2), pressure); 
                end
            end

            [fv, xo, yv, kval, visc_o, visc_g, MolarDensHC, MolarDensliq, MolarDensvap, oDen, gDen, MolarVolliq, MolarVolvap, ~, ~, ~, ~] ...
            = Y_FluidFlash(pressure, Zi_comp, temperature, fluid, kval_0, fv_0, 1, Pb);

            % Step 1: calculate amount of moles for oil and gas in unit cell
            % volume, lb-mol/cu.ft
            prop(stat).cell_No(nodeID, 1) = (1 - fv) * N_hc(nodeID, 1);  
            prop(stat).cell_Ng(nodeID, 1) = fv * N_hc(nodeID, 1);  

            % Step 2: hydrocarbon properties

            prop(stat).cell_Ki(nodeID, :)=                   kval;       % ratio
            prop(stat).cell_fv(nodeID, 1)=                     fv;       % mol fraction
            prop(stat).cell_Vsic_o(nodeID, 1)=             visc_o;       % cP  
            prop(stat).cell_Vsic_g(nodeID, 1)=             visc_g;       % cP
            prop(stat).cell_MolarDensHC(nodeID, 1) =  MolarDensHC;       % lb-mole/cu.ft
            prop(stat).cell_MolarDensLiq(nodeID, 1)= MolarDensliq;       % lb-mole/cu.ft
            prop(stat).cell_MolarDensVap(nodeID, 1)= MolarDensvap;       % lb-mole/cu.ft
            prop(stat).cell_MolarVolLiq(nodeID, 1)  = MolarVolliq;
            prop(stat).cell_MolarVolVap(nodeID, 1)  = MolarVolvap;               
            prop(stat).cell_MDen_o(nodeID, 1)=               oDen;       % lbm/cu.ft
            prop(stat).cell_MDen_g(nodeID, 1)=               gDen;       % lbm/cu.ft
            prop(stat).cell_Xi(nodeID, :)=                     xo;       % mol fraction
            prop(stat).cell_Yi(nodeID, :)=                     yv;       % mol fraction

        end
        
    end
            
    % Step 3: water properties        
    % Water mass density:  lb/cu.ft
    Xc(Target_Nodes) = fluid.Cw .* ( prop(stat).cell_P(Target_Nodes) - fluid.Pref_w);
    cell_Bw_inv(Target_Nodes) = 1.0 ./ fluid.Bw_ref .* (1 + Xc(Target_Nodes) .* (1 + 0.5 .* Xc(Target_Nodes)) );
    prop(stat).cell_MDen_w(Target_Nodes)  = fluid.MDen_w_sc .* cell_Bw_inv(Target_Nodes);    
    % Wate molar density: lb-mol/cu.ft
    prop(stat).cell_MolarDensWat(Target_Nodes) = prop(stat).cell_MDen_w(Target_Nodes) ./ fluid.wat_MW;

    % Water viscosity, cp
    Yc(Target_Nodes) = - fluid.Cvw .* ( prop(stat).cell_P(Target_Nodes) - fluid.Pref_w);
    prop(stat).cell_Vsic_w(Target_Nodes) = fluid.Vsic_w_ref ./ (1 + Yc(Target_Nodes) .* (1 + 0.5 .* Yc(Target_Nodes)) );     
    
    % Step 4: calculate oil, gas, water saturation
    Vo(Target_Nodes) = prop(stat).cell_No(Target_Nodes) ./ prop(stat).cell_MolarDensLiq(Target_Nodes);
    Vg(Target_Nodes) = prop(stat).cell_Ng(Target_Nodes) ./ prop(stat).cell_MolarDensVap(Target_Nodes);
    Vw(Target_Nodes) = prop(stat).cell_Nw(Target_Nodes) ./ prop(stat).cell_MolarDensWat(Target_Nodes);
    Vo(isnan(Vo) | isinf(Vo)) = 0.0;
    Vg(isnan(Vg) | isinf(Vg)) = 0.0;
    Vw(isnan(Vw) | isinf(Vw)) = 0.0;
                    
    Vft(Target_Nodes) = Vo(Target_Nodes) + Vg(Target_Nodes) + Vw(Target_Nodes); 
    prop(stat).cell_So(Target_Nodes) = Vo(Target_Nodes) ./ Vft(Target_Nodes); 
    prop(stat).cell_Sg(Target_Nodes) = Vg(Target_Nodes) ./ Vft(Target_Nodes);
    prop(stat).cell_Sw(Target_Nodes) = Vw(Target_Nodes) ./ Vft(Target_Nodes);
    prop(stat).cell_Vfluid(Target_Nodes) = Vft(Target_Nodes); 
    
    % Cut small saturations 
    prop(stat).cell_So(prop(stat).cell_So < 1.0e-10) = 0.0d0;
    prop(stat).cell_Sg(prop(stat).cell_Sg < 1.0e-10) = 0.0d0;
    prop(stat).cell_Sw(prop(stat).cell_Sw < 1.0e-10) = 0.0d0;
                
end