function [init] = Check_PSat(init, well, prop, fluid, res)

    % Define temporary variables
    nCells = init.nCells;
    nWells = well.nWells;
    nComps = fluid.nComps;
    fluidtype = fluid.fluidtype;
    key_comp = fluid.key_comp;
    current = 1;
    Tol = 5.0e-2;
    
    % Gas reservoir, directly return
    if isequal(fluidtype, 'GAS')
        return
    end
    
    N_update = 0;

    % Calculate composition of each cell
    N_hc = sum(prop(current).cell_Ni, 2);    
    prop(current).cell_Zi = prop(current).cell_Ni ./ N_hc(:, ones(1, nComps));
       
    %% Reservoir cell bubble point check
    for iCell = 1 : nCells
        % Check the change of the key composition
        Zi_comp = prop(current).cell_Zi(iCell, :);
        Zi_keycomp = Zi_comp(1, key_comp);
        offset = abs( Zi_keycomp - init.res_vle_chk(iCell, nComps+2) ) / init.res_vle_chk(iCell, nComps+2);  
                   
        % Update bubble point if large offset
        if offset > Tol
            
            % Get initial guess of Pb from neighboring connected cells
            ngh_cell = sort( [ init.conne_cell2(init.conne_cell1 == iCell); ...
                               init.conne_cell1(init.conne_cell2 == iCell) ]); 
            % Exclude gas injection perf cells
            if ~isempty(well.Inj_Gas_cells)
                for i = 1 : nnz(well.Inj_Gas_cells)
                    ngh_cell = ngh_cell( ngh_cell ~= well.Inj_Gas_cells(i) ); 
                end                
            end
              
            if isempty(ngh_cell)
                Psat_0 = init.res_vle_chk(iCell, nComps+1);
                kval_0 = init.res_vle_chk(iCell, 1:nComps);                
            else
                temp = abs(init.res_vle_chk(ngh_cell, nComps+2) - Zi_keycomp); 
                opm_ngh = ngh_cell(temp == min(temp));
                opm_nghc = opm_ngh(1); 
                Psat_0 = init.res_vle_chk(opm_nghc, nComps+1);
                kval_0 = init.res_vle_chk(opm_nghc, 1:nComps);
            end
            
            [Psat, kval] = Y_SatPoint(Zi_comp, prop(current).cell_T(iCell), fluid, Psat_0, kval_0, 1);
            init.res_vle_chk(iCell, 1:nComps) = kval;
            init.res_vle_chk(iCell, nComps+1) = Psat;
            init.res_vle_chk(iCell, nComps+2) = Zi_keycomp; 
            N_update = N_update + 1;
        end
               
    end
    
    %% Producer bubble point check
%     for iWell = 1 : nWells
%         
%         % Only for producer
%         if ~well.Prod_Flag(iWell)
%             continue;
%         end
%         
%         % Check the change of the key composition
%         Zi_comp = well.Gas_Comp(iWell, :);
%         Zi_keycomp = Zi_comp(1, key_comp);
%         offset = abs( Zi_keycomp - init.wel_vle_chk(iWell, nComps+2) ) / init.wel_vle_chk(iWell, nComps+2);
%         
%         % Update bubble point if large offset
%         Psat_0 = init.wel_vle_chk(iWell, nComps+1);
%         if offset > Tol && Psat_0 < 9000.0d0
%                        
%             kval_0 = init.wel_vle_chk(iWell, 1:nComps);
%                         
%             [Psat, kval] = Y_SatPoint(Zi_comp, res.T(1), fluid, Psat_0, kval_0, 1);
%             init.wel_vle_chk(iWell, 1:nComps) = kval;
%             init.wel_vle_chk(iWell, nComps+1) = Psat;
%             init.wel_vle_chk(iWell, nComps+2) = Zi_keycomp;
%             N_update = N_update + 1;
%             
%         end        
%         
%     end
    
    disp( strcat( 'Times of bubble point update: ', num2str(N_update) ) );
    
end