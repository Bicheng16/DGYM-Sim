% Calculate mobility for grid-grid connection and well

function [a_o, a_w, a_g, a_g_app] = Mobility_Convec_Flux(init, fluid, slippage_flag, prop, Target_Connes, ...
    flag_upstream_o, flag_upstream_w, flag_upstream_g)

    % Define temperary variables 
    nConnes = init.nConnes;
    nComps = fluid.nComps;
    
    next = 2;

    % Dynamic (fluid-dependant) terms of Convection (Darcy Flow) Transmissibility: for reservoir grid-grid connection
    % Define mobilities cell-cell connection based on connection
    Lambda_o = zeros(nConnes, 1); a_o = zeros(nConnes, 1);  
    Lambda_g = zeros(nConnes, 1); a_g = zeros(nConnes, 1); a_g_app = zeros(nConnes, nComps); 
    Lambda_w = zeros(nConnes, 1); a_w = zeros(nConnes, 1);  
        
    % Upstream Kr*rho/miu (for hydrocarbon, mole fraction also need to be
    % upstreamed, according to Cao Hui(P42, Eq. 2.14)
    Node1_in_conne = init.conne_cell1(Target_Connes);
    Node2_in_conne = init.conne_cell2(Target_Connes);
    
    Lambda_o(Target_Connes) = prop(next).cell_kro(Node1_in_conne) .* prop(next).cell_MolarDensLiq(Node1_in_conne)  ...
                           ./ prop(next).cell_Vsic_o(Node1_in_conne) .* flag_upstream_o(Target_Connes)             ...
                            + prop(next).cell_kro(Node2_in_conne) .* prop(next).cell_MolarDensLiq(Node2_in_conne)  ...
                           ./ prop(next).cell_Vsic_o(Node2_in_conne) .* (1 - flag_upstream_o(Target_Connes));   

    Lambda_g(Target_Connes) = prop(next).cell_krg(Node1_in_conne) .* prop(next).cell_MolarDensVap(Node1_in_conne)  ...
                           ./ prop(next).cell_Vsic_g(Node1_in_conne) .* flag_upstream_g(Target_Connes)             ...
                            + prop(next).cell_krg(Node2_in_conne) .* prop(next).cell_MolarDensVap(Node2_in_conne)  ...
                           ./ prop(next).cell_Vsic_g(Node2_in_conne) .* (1-flag_upstream_g(Target_Connes)); 
                       
    Lambda_w(Target_Connes) = prop(next).cell_krw(Node1_in_conne) .* prop(next).cell_MolarDensWat(Node1_in_conne)  ...
                           ./ prop(next).cell_Vsic_w(Node1_in_conne) .* flag_upstream_w(Target_Connes)             ...
                            + prop(next).cell_krw(Node2_in_conne) .* prop(next).cell_MolarDensWat(Node2_in_conne)  ...
                           ./ prop(next).cell_Vsic_w(Node2_in_conne) .* (1 - flag_upstream_w(Target_Connes)); 


    % Ensure zeros (instead of NaN / inf) when phase volume is zero for 2 neighbor cells
    Lambda_o(isnan(Lambda_o) | isinf(Lambda_o)) = 0;    
    Lambda_g(isnan(Lambda_g) | isinf(Lambda_g)) = 0;
    Lambda_w(isnan(Lambda_w) | isinf(Lambda_w)) = 0;    

    % Full transmissibility (geometric*mobility)
    if(~init.trans_split)
        
        a_o(Target_Connes) = init.conne_ctrans(Target_Connes) .* Lambda_o(Target_Connes); 
        a_g(Target_Connes) = init.conne_ctrans(Target_Connes) .* Lambda_g(Target_Connes);
        a_w(Target_Connes) = init.conne_ctrans(Target_Connes) .* Lambda_w(Target_Connes);
        
    else
        
        % Transmissibility based on intrinsic permeability
        trans = zeros(nConnes, 1); % For water and oil 
        trans(Target_Connes) = 1.0d0 ./ ( init.conne_perm(Target_Connes, 1) .* init.conne_geo(Target_Connes, 1) ) ...
                             + 1.0d0 ./ ( init.conne_perm(Target_Connes, 2) .* init.conne_geo(Target_Connes, 2) );
        trans(Target_Connes) = 1.0d0 ./ trans(Target_Connes);

        % Oil and water mobility without slippage
        a_o(Target_Connes) = trans(Target_Connes) .* Lambda_o(Target_Connes);        
        a_w(Target_Connes) = trans(Target_Connes) .* Lambda_w(Target_Connes);
                
        if(~slippage_flag)
            
            a_g(Target_Connes) = trans(Target_Connes) .* Lambda_g(Target_Connes);
            
        else  
            
            for iComp = 1 : nComps
                
                a_g_app(Target_Connes, iComp) = ...
                  1.0d0 ./ ( init.conne_perm(Target_Connes, 1) .* init.conne_geo(Target_Connes, 1) ...
                          .* prop(next).cell_GKapp(Node1_in_conne, iComp) ) ...
                + 1.0d0 ./ ( init.conne_perm(Target_Connes, 2) .* init.conne_geo(Target_Connes, 2) ... 
                          .* prop(next).cell_GKapp(Node2_in_conne, iComp) );
                      
                a_g_app(Target_Connes, iComp) = Lambda_g(Target_Connes) ./ a_g_app(Target_Connes, iComp);                   
                
            end
            
        end
        
    end
    
end
