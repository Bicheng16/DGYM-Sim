function [prop] = Fluid_Mass(prop, status, nComps)

    % Calculate water, oil, gas, hydrocarbon component mass in unit cell
    % volume

    % Calculate oil moles in unit cell volume, lb-mol/cu.ft
    prop(status).cell_No = prop(status).cell_poro ...
                        .* prop(status).cell_MolarDensLiq ...
                        .* prop(status).cell_So;

    % Calculate gas moles in unit cell volume, lb-mol/cu.ft
    prop(status).cell_Ng = prop(status).cell_poro ...
                        .* prop(status).cell_MolarDensVap ...
                        .* prop(status).cell_Sg;

    % Calculate water moles in unit cell volume, lb-mol/cu.ft
    prop(status).cell_Nw = prop(status).cell_poro ...
                        .* prop(status).cell_MolarDensWat ...
                        .* prop(status).cell_Sw;

    % Calculate fluid volume / cell volume
    prop(status).cell_Vfluid = prop(status).cell_poro .* ...
                             ( prop(status).cell_So   +  ...
                               prop(status).cell_Sg   +  ...
                               prop(status).cell_Sw );

    dummy_No = repmat(prop(status).cell_No, 1, nComps);
    dummy_Ng = repmat(prop(status).cell_Ng, 1, nComps);
    
    % Calculate moles of hydrocarbon component i in unit cell volume
    prop(status).cell_Ni = dummy_No .* prop(status).cell_Xi ...
                         + dummy_Ng .* prop(status).cell_Yi;
                     
end