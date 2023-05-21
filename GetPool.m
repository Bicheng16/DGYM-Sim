
function [Target_Nodes, Target_Connes] = GetPool(Target_Nodes, nTNodes, nTConne, Node1_in_conne, Node2_in_conne)
    % This function is to get target connection and node for evaluation
    
    % Calculate RHS 
    if size(Target_Nodes, 2) == nTNodes
        
        % Calculate residual pieces for RHS: all nodes and connections
        % involved
        Target_Connes   = 1 : nTConne;  
    end
        
    % Perturb reservoir variables
    if size(Target_Nodes,2) == 1 && Target_Nodes ~= 0
        
        % Calculate residual pieces for Jacobian, perturb by variable, each
        % time one perturb one variable in a single node
        Target_Connes   = sort( [ find(Node1_in_conne==Target_Nodes); ...
                                  find(Node2_in_conne==Target_Nodes) ] );
                            
        Target_Connes   = Target_Connes';
    end
    
    % Perturb well variables
    if Target_Nodes == 0
        Target_Connes   = 1 : nTConne;
        Target_Nodes    = 1 : nTNodes; 
    end

end