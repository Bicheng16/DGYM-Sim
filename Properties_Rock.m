% Function to calculate rock properties 
% Function to calculate rock properties 
function [prop] = Properties_Rock(stat, prop, init, res, rockfluid, dummy, Target_Nodes, inherit, inherit_loc) 

    % Inheritance to avoid unnecessary calculation
    if inherit && ~isempty(inherit_loc)
		prop(stat).cell_poro(inherit_loc) = dummy.poro(inherit_loc); 
        prop(stat).cell_PV(inherit_loc) = dummy.PV(inherit_loc); 
        prop(stat).cell_HCPV(inherit_loc) = dummy.HCPV(inherit_loc);
        prop(stat).cell_krw(inherit_loc) = dummy.krw(inherit_loc);
        prop(stat).cell_krg(inherit_loc) = dummy.krg(inherit_loc);        
		prop(stat).cell_kro(inherit_loc) = dummy.kro(inherit_loc); 		 
		prop(stat).cell_pcow(inherit_loc) = dummy.pcow(inherit_loc); 
		prop(stat).cell_pcgo(inherit_loc) = dummy.pcgo(inherit_loc);		  
    end
    
    % Inherit all reservoir properties from dummy when perturbing well pressure 
    if nnz(inherit_loc) == init.nCells
        return
    end    
    
    % Download variables to local
    P_ref = res.P_ref(init.cell_type(1 : init.nCells));
    
    cell_krow = zeros(init.nCells, 1);
    cell_krog = zeros(init.nCells, 1);
            
    % Calculate porosity, pore volume, and hydrocarbon pore volume
    prop(stat).cell_poro(Target_Nodes) = init.cell_poro_ref(Target_Nodes) ...
    .* exp( init.cell_Cr(Target_Nodes) .* ( prop(stat).cell_P(Target_Nodes) - P_ref(Target_Nodes) ) );                        % Unit of cu.ft

    prop(stat).cell_PV(Target_Nodes) = init.cell_volume(Target_Nodes) .* prop(stat).cell_poro(Target_Nodes);  
        
    prop(stat).cell_HCPV(Target_Nodes) = prop(stat).cell_PV(Target_Nodes) .* ( 1.0 - prop(stat).cell_Sw(Target_Nodes) );      % Unit of cu.ft 

    % Cell rock type for rock fluid calculation
    % Perturb reservoir variables
    if nnz(Target_Nodes) == 1 && Target_Nodes ~= 0
        Target_Nodes_Group = [];
        for iRT = 1 : res.nRockTypes
            dummy_1.iGroup = [];
            Target_Nodes_Group = [Target_Nodes_Group; dummy_1];
        end                
        Target_Nodes_Group(init.cell_type(Target_Nodes)).iGroup = Target_Nodes;          
    else
    % Calculate RHS or perturb Pwf
        Target_Nodes_Group = init.cell_group;
    end
    
    % Calcualte rock-fluid parameters: interp1q, nakeinterp1, spline
    for iRT = 1 : res.nRockTypes
        
        % Get cell names belonging to each rock type
        tmp_nodes = Target_Nodes_Group(iRT).iGroup;
        
        if isempty(tmp_nodes)
           continue; 
        end
        
        % Calculate krw
        prop(stat).cell_krw(tmp_nodes) = nakeinterp1( rockfluid(iRT).SWFN(:, 1), ...
                                                      rockfluid(iRT).SWFN(:, 2), ...
                                                      prop(stat).cell_Sw(tmp_nodes) );
                                             
        % Calculate krg
        prop(stat).cell_krg(tmp_nodes) = nakeinterp1( rockfluid(iRT).SGFN(:, 1), ...
                                                      rockfluid(iRT).SGFN(:, 2), ...
                                                      prop(stat).cell_Sg(tmp_nodes) );
                                             
        % Calculate pcow 
        prop(stat).cell_pcow(tmp_nodes) = nakeinterp1( rockfluid(iRT).SWFN(:, 1), ...
                                                       rockfluid(iRT).SWFN(:, 3), ...
                                                       prop(stat).cell_Sw(tmp_nodes) );
        
        % Calculate pcgo
        prop(stat).cell_pcgo(tmp_nodes) = nakeinterp1( rockfluid(iRT).SGFN(:, 1), ...
                                                       rockfluid(iRT).SGFN(:, 3), ...
                                                       prop(stat).cell_Sg(tmp_nodes) );
                                                   
        % Negative value check 
        prop(stat).cell_krw(prop(stat).cell_krw < 0.0d0) = 0.0d0;
        prop(stat).cell_krg(prop(stat).cell_krg < 0.0d0) = 0.0d0;
        prop(stat).cell_pcow(prop(stat).cell_pcow < 0.0d0) = 0.0d0;
        prop(stat).cell_pcgo(prop(stat).cell_pcgo < 0.0d0) = 0.0d0;
                
        % Calculate krow
        cell_krow(tmp_nodes) = nakeinterp1( rockfluid(iRT).SOF3(:, 1), ...
                                            rockfluid(iRT).SOF3(:, 2), ...
                                            1 - prop(stat).cell_Sw(tmp_nodes) );
    
        cell_krog(tmp_nodes) = nakeinterp1( rockfluid(iRT).SOF3(:, 1), ...
                                            rockfluid(iRT).SOF3(:, 3), ...
                                            1 - prop(stat).cell_Sg(tmp_nodes) - rockfluid(iRT).Swco );
            
        % Different method to calculate three phase relative permeabilities

        % Option 0: rockfluid(iRT).method == 0
        % STONE II kr calcutation using kr_hat (based on tables with conate Sw) kro table Somax < 1 
        % Based on ECL Saturation Functions Tech Description (pg 935) and Aziz and Settari 1979 pg36
        if rockfluid(iRT).method == 0
            prop(stat).cell_kro(tmp_nodes) = ...
            rockfluid(iRT).kro_cw .* ( (cell_krow(tmp_nodes) ./ rockfluid(iRT).kro_cw + prop(stat).cell_krw(tmp_nodes))   ...
                                  .*   (cell_krog(tmp_nodes) ./ rockfluid(iRT).kro_cw + prop(stat).cell_krg(tmp_nodes))   ...
                                     -  prop(stat).cell_krw(tmp_nodes) - prop(stat).cell_krg(tmp_nodes) );
        end

        % Option 1: rockfluid(iRT).method == 1
        % STONE II kr calcutation using kr_bar (based on tables with conate Sw) kro table reaches Somax = 1
        if rockfluid(iRT).method == 1
            prop(stat).cell_kro(tmp_nodes) = ( cell_krow(tmp_nodes) + prop(stat).cell_krw(tmp_nodes) ) ...           
                                          .* ( cell_krog(tmp_nodes) + prop(stat).cell_krg(tmp_nodes) )   ... 
                                           - ( prop(stat).cell_krw(tmp_nodes) + prop(stat).cell_krg(tmp_nodes) );
        end
                                                                     
    end
                    
    % Set as 1 the values of kro from STONE >1
    prop(stat).cell_kro(Target_Nodes) = prop(stat).cell_kro(Target_Nodes) ./ (prop(stat).cell_kro(Target_Nodes) <=1 );

    prop(stat).cell_kro( isnan(prop(stat).cell_kro) | isinf(prop(stat).cell_kro) ) = 1; 

    % Set as zero the values of kro from STONE <0
    prop(stat).cell_kro(Target_Nodes) = prop(stat).cell_kro(Target_Nodes) .* ( prop(stat).cell_kro(Target_Nodes) >= 0 );        

end
