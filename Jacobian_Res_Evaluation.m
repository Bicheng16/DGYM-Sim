% This function is to evaluate Reservoir Jacobian and store in triplet format

function [ILoc, JLoc, Val, nonzeros] = ...
   Jacobian_Res_Evaluation( node_RESID_0, conne_RESID_0, ...
                            node_RESID_1, conne_RESID_1, ...
                            Node1_in_conne,  Node2_in_conne,   ...
                            iVar, nComps, nCells, Peturb,      ...
                            searching_dir, perturb_node)

    % Evaluate current nonzero residual pieces for Jacobian and generate
    % triplet format of Jacobian
    % Length of triplet   
    
    % Control parameter: searching_dir is to flag for positive or negative
    % perturbation for different parameters, only 2 values: 1 (positive),
    % -1 (negative). Note: pressure is negatively perturbed.
    
    neq = nComps + 2; 
    nonzeros = nnz(node_RESID_1) + nnz(conne_RESID_1);
            
    %% Memory preallocation 
    node_Jacobian = 0.0 * node_RESID_1;
    conne_Jacobian = 0.0 * conne_RESID_1; 
    vector_diagonal = zeros(1, neq);   
    
    ILoc = [];
    JLoc = [];
    Val  = [];
        
    %% Jacobian Piece Evaluation: target only for nnz entries in it 
    
    % Node Jacobian: operate on node_RESID(perturb_node, :) 
    node_Jacobian(perturb_node, :) = (node_RESID_1(perturb_node, :)  -  node_RESID_0(perturb_node, :))  ./ Peturb .* searching_dir;    
    
    % Conne Jacobian: operate on connection including perturb_node
    perturb_conne   = sort( [ find(Node1_in_conne == perturb_node); ...
                              find(Node2_in_conne == perturb_node) ] );
    perturb_conne = perturb_conne';
        
    conne_Jacobian(perturb_conne, :) = (conne_RESID_1(perturb_conne, :)  - conne_RESID_0(perturb_conne, :)) ./ Peturb .* searching_dir;
    
    % Current Node Jacobian block (pieces)
    % Jacobian = JAccumulation - sum(JFlux * flow_dir), 
    % If perturb_node is node1 in connection pair, flow_dir = +1;
    % If perturb_node is node2 in connection pair, flow_dir = -1;
    % well is either part of accumulation (node based part) property   
    % or treated as node, having both accumulation and flux
    
    % In perturbation, perturbed variable belongs to a unique node index
    
    if nnz(node_Jacobian) == 0 && nnz(conne_Jacobian) == 0  
        return
    end

    nnz_offdiag = nnz(conne_Jacobian);   
    
    % Important: diagonal entries for flux Jacobian always lumped 
    % with accumulation Jacobian, extremely dagerous if mistake 
    
    [ID_conne, colJOD] = find(conne_Jacobian~=0);  
    
    % Always has node_Jacobian in block diagonal
    vector_diagonal(1, :) = node_Jacobian(perturb_node, :);  
        
    % Lump flux through connection direction (sign is the key)
    if nnz_offdiag ~= 0 
        for i = 1 : size(perturb_conne,2)
            % Determine perturb_node is node1 or node2 in conne pair  
            % Flux from node1 to node2 is positive, otherwise, negative
            if perturb_node == Node1_in_conne(perturb_conne(i))           
                flow_dir    =   1;  
            end
            if perturb_node == Node2_in_conne(perturb_conne(i))
                flow_dir    =  -1;   
            end
            % Filling vector_diagonal with flux term 
            vector_diagonal(1, 1:nComps+1) = vector_diagonal(1, 1:nComps+1) - conne_Jacobian(perturb_conne(i), 1:nComps+1) * flow_dir;                                               
        end
    end    
            

    % Memory Preallocation
    nnz_diag = nnz(vector_diagonal);
    
    nonzeros = nnz_diag + nnz_offdiag;
    ILoc = zeros([nonzeros, 1]);
    JLoc = zeros([nonzeros, 1]);
    Val  = zeros([nonzeros, 1]);
    
    % Store vector_diagonal in triplet format 
    % Column index decided by which variable you are perturbing 
    % Row index decided by which equation you are evaluating
    [~, col_JD] = find(vector_diagonal~=0);   
        
    if nnz_diag ~= 0
        JLoc(1:nnz_diag) = iVar;
        for i = 1 : nnz_diag
            ILoc(i) = col_JD(i) + (perturb_node-1) * neq;
            Val (i) = vector_diagonal(1, col_JD(i));  
        end        
    end
    
    % Store off_diagonal entries generated during this perturbation
    % Directly from conne_Jacobian, but it has opposite sign compared 
    % to diagonal 
    % (Refer to Cao Hui, disertation (GPRS), page 58)
    
    % Schematic: Flux between node A and B is F(AB)
        
    % % % % % %% % % % % %
    %         %          %
    %   A <<====>>  B    %
    %         %          %
    % % % % % %% % % % % %
    
    % In terms of Jacobian is the following format
     
    %   |                                 |
    %   |  dF(AB)/dX(A)     dF(AB)/dX(B)  |
    %   |                                 |
    %   | -dF(AB)/dX(A)    -dF(AB)/dX(B)  |
    %   |                                 |
    
    start = nnz_diag + 1; 
    if nnz_offdiag ~= 0 
        for i = 1 : nnz_offdiag
            % Locating neighbor node index
            % Flux from node1 to node2 is positive, otherwise, negative
            if perturb_node == Node1_in_conne(ID_conne(i))           
                neigh_node  =  Node2_in_conne(ID_conne(i));
%                 flow_dir    =  1;
                flow_dir    =  -1;
            end
            if perturb_node == Node2_in_conne(ID_conne(i))
                neigh_node  =  Node1_in_conne(ID_conne(i)); 
%                 flow_dir    =  -1; 
                flow_dir    =  1;
            end
            % Filling off-diagonal (upper or lower)
            if  neigh_node <= nCells 
                JLoc(start) = iVar;
                 % Currently considering all nodes has the same number of 
                 % equation, later tracking neq each node is a must.
                ILoc(start) = colJOD(i) + (neigh_node-1) * neq; 
%                 Val (start) = conne_Jacobian(ID_conne(i), colJOD(i)) * flow_dir;
                Val (start) = - conne_Jacobian(ID_conne(i), colJOD(i)) * flow_dir;
                                
                start = start + 1;
            end
            
        end
    end
    
    % Store in sparse since at the end of each triplet vectors might have
    % zeros 
    ILoc = sparse(ILoc); JLoc = sparse(JLoc); Val = sparse(Val); 
    
        
end