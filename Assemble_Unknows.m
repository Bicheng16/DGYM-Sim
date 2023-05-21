% This function is to assemble unknow vector x_vec, including Ni, Nw, and P

function [x_vec, nVar] = Assemble_Unknows(Ni, Nw, P, Pwf)

    % nVar: number of primary variables in the unknown vector
    % nComps: number of components in hydrocarbon
    % nCells: number of cells in the reservoir
    
    [nCells, nComps] = size(Ni);
    nRVar = nCells * (nComps + 2);
    
    Ni_dummy = reshape(Ni', nComps, nCells);
    Nw_dummy = reshape(Nw, 1, nCells);   
    P_dummy = reshape(P, 1, nCells);
    
    % Reservoir variable
    x_vec = [Ni_dummy; Nw_dummy; P_dummy];
    x_vec = reshape(x_vec, nRVar, 1);  

    % Well variable at the end
    x_vec = [x_vec; Pwf];
    
    nVar = size(x_vec, 1);
end