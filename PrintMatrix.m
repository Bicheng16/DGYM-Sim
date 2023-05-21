function PrintMatrix(iLoc, jLoc, Value, b)

    format longE;
    nnzeros = nnz(iLoc);
    sizeb = size(b,1);
    
    % Open a file 
    fid = fopen('Jacobian_Output', 'wt');
    
    %
    for i = 1 : nnzeros
        
       fprintf(fid,'%2s %3u %2s %3u %4s %22.15e\n', 'A(', iLoc(i), ', ', jLoc(i), ') = ', Value(i));    
        
    end
    
    fprintf(fid, '\n\n');
    
    for i = 1 : sizeb
        
       fprintf(fid,'%2s %3u %7s %22.15e\n', 'b(', i, ', 1) = ', b(i));    
        
    end
    
end

