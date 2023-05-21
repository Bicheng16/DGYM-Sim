    
% function [x, flag_Chop_dt_k] = LinearSolver(A, b, Sln_Option)
function [x, flag_Chop_dt_k, Jacobian_cond] = LinearSolver(iLoc, jLoc, Value, nonzeros, dim, b, Sln_Option, Jacobian_cond)
    %% Matrix solver library to solver Ax = b
    % iLoc: i index of nonzero entries 
    % jLoc: j index of nonzero entries 
    % Value: value of nonzero entries at location (i, j)
    % dim: dimension of sparse matrix

    %% Triplet configuration to fill in sparse
    A = sparse( iLoc(1:nonzeros),   ...
                jLoc(1:nonzeros),   ...
                Value(1:nonzeros), ...
                dim, dim );

    % Dump this: don't add anything on zero diags
    % Vector with zero-elements in Jac diagonal
    Zero_diag_Elements=(diag(A)==0);
    % Replace zeros in diagonal by ones
    DiagNewElements=(diag(A))+Zero_diag_Elements;
    A(1:(dim+1):dim*dim) = DiagNewElements;
                        
    %Jacobian_cond = [Jacobian_cond cond(full(A))];
    Jacobian_cond = [Jacobian_cond 1];
                            
    flag = 0;
    msgString = '';
    flag_Chop_dt_k = false;

    switch Sln_Option
        
        % Direct solver: back slash, better than inverse 
        case 1
            
            [A, b] = reorderForILU(A, b);
            x = A\b;
        
        % Conjugate gradient with ILU precondition     
        case 2
            
            tol = 1e-6;
            nel = size(A, 1);
            maxit = 1000;
            % PCG with reordering and iLU preconditioner
            [A, b] = reorderForILU(A, b);
            [L,U] = ilu(A,struct('type','ilutp','droptol',1e-4));  % 1e-6
            prec = @(x) U\(L\x);
            [x, flag, relres, iter] = pcg(A, b, [], tol, min(nel, maxit), prec);

        % GMRES with ILU preconditioner 
        case 3
                        
            tol = 1e-6;   % 1e-6;
            nel = size(A, 1);
            maxit = 1000;
            % GMRES with reordering and iLU preconditioner
            [A, b] = reorderForILU(A, b);
            %[L,U] = ilu(A,struct('type','ilutp','droptol',1e-12));
            [L,U] = ilu(A,struct('type','ilutp','droptol',1e-8));
            prec = @(x) U\(L\x);
            [x, flag, relres, iter] = gmres(A, b, [], tol, min(nel, maxit), prec);
       
        case 4
            
            tol = 1e-6;
            nel = size(A, 1);
            maxit = 1000;
            % BiCGSTAB with reordering and iLU preconditioner
            [A, b] = reorderForILU(A, b);
            %[L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
            
            lfil = 3;
            [L, U] = iluk(A, lfil);
            
            prec = @(x) U\(L\x);
            [x, flag, relres, iter]= bicgstab(A, b, tol , min(nel, maxit), prec);
            
        otherwise 
            msg = ['......WARNING: There is no option of this solver implemented in the simulator, please check it'];
            disp(msg); 
            flag_Chop_dt_k=true;
    end
    
    
    %% Check if valid run
    if flag ~= 0
        disp('......WARNING: Matrix solver did not converge. See Solve_Matrix subroutine for additional information')
        flag_Chop_dt_k = true;
        return
    end

end
