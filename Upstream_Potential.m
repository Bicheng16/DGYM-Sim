% This function is to calculate potentials for grid-grid connection 

function [delta_pot_o, delta_pot_w, delta_pot_g, flag_upstream_o, flag_upstream_w, flag_upstream_g] ...
    = Upstream_Potential(prop, init, Target_Connes)

    % Define temperary variables
    nConnes = init.nConnes;
    gc = 32.17 * 0.21584E-3;
    next = 2;
    
    % Memory preallocation 
    rho_o = zeros(nConnes, 1); 
    rho_g = zeros(nConnes, 1);
    rho_w = zeros(nConnes, 1);

    delta_press = zeros(nConnes, 1);
    
    delta_pcgo = zeros(nConnes, 1);
    delta_pcow = zeros(nConnes, 1);

    delta_pot_o = zeros(nConnes, 1);
    delta_pot_g = zeros(nConnes, 1);
    delta_pot_w = zeros(nConnes, 1);

    flag_upstream_o = zeros(nConnes, 1);
    flag_upstream_g = zeros(nConnes, 1);
    flag_upstream_w = zeros(nConnes, 1);

    Node1_in_conne = init.conne_cell1(Target_Connes);
    Node2_in_conne = init.conne_cell2(Target_Connes);    
    
    % Arithmetic average of mass densities to calculate gravity potential in RRConnes_Pool
    % In between upstream weighting and phase pore volume average
    % rho_o(Target_Connes) = 0.5 * (prop(next).cell_MDen_o(Node1_in_conne) + prop(next).cell_MDen_o(Node2_in_conne));
    % rho_g(Target_Connes) = 0.5 * (prop(next).cell_MDen_g(Node1_in_conne) + prop(next).cell_MDen_g(Node2_in_conne));
    % rho_w(Target_Connes) = 0.5 * (prop(next).cell_MDen_w(Node1_in_conne) + prop(next).cell_MDen_w(Node2_in_conne));

    % Saturation average, Equation (7) in SPE 79692
    rho_o(Target_Connes) = ( prop(next).cell_MDen_o(Node1_in_conne) .* prop(next).cell_So(Node1_in_conne)   ...
                        +   prop(next).cell_MDen_o(Node2_in_conne) .* prop(next).cell_So(Node2_in_conne) ) ...
                       ./ ( prop(next).cell_So(Node1_in_conne) + prop(next).cell_So(Node2_in_conne) );
                     
    rho_g(Target_Connes) = ( prop(next).cell_MDen_g(Node1_in_conne) .* prop(next).cell_Sg(Node1_in_conne)   ...
                         +   prop(next).cell_MDen_g(Node2_in_conne) .* prop(next).cell_Sg(Node2_in_conne) ) ...
                        ./ ( prop(next).cell_Sg(Node1_in_conne) + prop(next).cell_Sg(Node2_in_conne) );

    rho_w(Target_Connes) = ( prop(next).cell_MDen_w(Node1_in_conne) .* prop(next).cell_Sw(Node1_in_conne)   ...
                         +   prop(next).cell_MDen_w(Node2_in_conne) .* prop(next).cell_Sw(Node2_in_conne) ) ...
                        ./ ( prop(next).cell_Sw(Node1_in_conne) + prop(next).cell_Sw(Node2_in_conne) );

    % Volume averaged
%     rho_o(Target_Connes) =(  prop(next).cell_MDen_o(Node1_in_conne) .* prop(next).cell_So(Node1_in_conne) ...
%                          .*  prop(next).cell_PV(Node1_in_conne)                                           ...  
%                          +   prop(next).cell_MDen_o(Node2_in_conne) .* prop(next).cell_So(Node2_in_conne) ...
%                          .*  prop(next).cell_PV(Node2_in_conne) )                                         ...
%                          ./( prop(next).cell_So(Node1_in_conne) .*  prop(next).cell_PV(Node1_in_conne)    ...
%                          +   prop(next).cell_So(Node2_in_conne) .*  prop(next).cell_PV(Node2_in_conne) );
% 
%     rho_g(Target_Connes) =(  prop(next).cell_MDen_g(Node1_in_conne) .* prop(next).cell_Sg(Node1_in_conne) ...
%                          .*  prop(next).cell_PV(Node1_in_conne)                                           ...  
%                          +   prop(next).cell_MDen_g(Node2_in_conne) .* prop(next).cell_Sg(Node2_in_conne) ...
%                          .*  prop(next).cell_PV(Node2_in_conne) )                                         ...
%                          ./( prop(next).cell_Sg(Node1_in_conne) .*  prop(next).cell_PV(Node1_in_conne)    ...
%                          +   prop(next).cell_Sg(Node2_in_conne) .*  prop(next).cell_PV(Node2_in_conne) );
% 
%     rho_w(Target_Connes) =(  prop(next).cell_MDen_w(Node1_in_conne) .* prop(next).cell_Sw(Node1_in_conne) ...
%                          .*  prop(next).cell_PV(Node1_in_conne)                                           ...  
%                          +   prop(next).cell_MDen_w(Node2_in_conne) .* prop(next).cell_Sw(Node2_in_conne) ...
%                          .*  prop(next).cell_PV(Node2_in_conne) )                                         ...
%                          ./( prop(next).cell_Sw(Node1_in_conne) .*  prop(next).cell_PV(Node1_in_conne)    ...
%                          +   prop(next).cell_Sw(Node2_in_conne) .*  prop(next).cell_PV(Node2_in_conne) );

    % Check values
    rho_o(isnan(rho_o) | isinf(rho_o)) = 0.0;
    rho_g(isnan(rho_g) | isinf(rho_g)) = 0.0;
    rho_w(isnan(rho_w) | isinf(rho_w)) = 0.0;

    % Pressure difference between cell1 and cell2
    delta_press(Target_Connes) = prop(next).cell_P(Node2_in_conne) - prop(next).cell_P(Node1_in_conne);
    
    % Capillary pressure difference between cell1 and cell2 in RRConnes_Pool
    delta_pcgo(Target_Connes) = prop(next).cell_pcgo(Node2_in_conne) - prop(next).cell_pcgo(Node1_in_conne);
    delta_pcow(Target_Connes) = prop(next).cell_pcow(Node2_in_conne) - prop(next).cell_pcow(Node1_in_conne);
    
    % Potential difference between cell1 and cell2  
    %  Pot_2-Pot_1 = (P_2 - P_1) - (Pcow_2 - Pcow_1) - gc * [(rho_w1 + rho_w2)/2] * (Z_2 - Z_1)
    %              = (P_2 - P_1) - (delta_pcow)      - gc * [rho_w]               *  conne_dz
    %  for oil: capillary presures are zero
    delta_pot_o(Target_Connes) = delta_press(Target_Connes) - gc .* rho_o(Target_Connes) .* init.conne_dz(Target_Connes);
    delta_pot_g(Target_Connes) = delta_press(Target_Connes) + delta_pcgo(Target_Connes) - gc .* rho_g(Target_Connes) .* init.conne_dz(Target_Connes); 
    delta_pot_w(Target_Connes) = delta_press(Target_Connes) - delta_pcow(Target_Connes) - gc .* rho_w(Target_Connes) .* init.conne_dz(Target_Connes);

    % flag upstream indicates the direction of the flow
    flag_upstream_o(Target_Connes) = delta_pot_o(Target_Connes)<=0; % if delta_pot_o<=0, flag_upstream_o=1 (oil flows from cell1 to cell2)
    flag_upstream_g(Target_Connes) = delta_pot_g(Target_Connes)<=0; % if delta_pot_g<=0, flag_upstream_g=1 (gas flows from cell1 to cell2)
    flag_upstream_w(Target_Connes) = delta_pot_w(Target_Connes)<=0; % if delta_pot_w<=0, flag_upstream_w=1 (water flows from cell1 to cell2)

end