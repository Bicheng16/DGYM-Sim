      function [L, U] = iluk(A,lfil) 
%      function [L, U] = iluk(A,lfil) 
%      ILU with level of fill-in of k (ILU(k))
%  Input:
%  A     = sparse matrix 
%  lfil  = integer. Ill-in parameter. entries whose
%          levels-of-fill exceed ffil during the 
%          ILU process are dropped. lfil must be >= 0 
% Output: 
% L, U = factors of ilu(k)-factorization of A - 
%
n = size(A,1) ;
L = speye(n,n); 
U = sparse(n,n);
levs = sparse(n,n); 
%% NOTE: the code uses a shifted definition of lfil : 
lfilp = lfil+1 ; 
% -----------------------------------------------------------------------
%      beginning of main loop.
% ----------------------------------------------------------------------- ;
  for ii = 1:n  %%%% 500 
%
%    unpack row ii of A 
%
     w = A(ii,:);
     jr = find(w) ;
     levii = sparse(1,n) ; 
     levii(jr) = 1;
     jj = 0;
     while (jj < ii)   %% while 
%
% select smallest column index among jw(k), k=jj+1, ..., lenl. 
%
       jj = jj+min(find(w(jj+1:n))) ;
       if (jj >= ii)
          break;
       end 
       jlev = levii(jj) ;
       if (jlev <= lfilp) 
%  
%      combine current row and  row jj.  
%  
       [ir,jr,v] = find(U(jj,:)); 
       fact = w(jj)/v(1); 
       for k=1:length(jr) 
         j = jr(k) ;
         temp = levs(jj,j) + jlev ; 
         if (levii(j) == 0 ) 
             if (temp <= lfilp) 
                 w(j) =  - fact*v(k); 
                 levii(j) = temp;
             end
         else 
             w(j) = w(j) - fact*v(k); 
             levii(j) = min(levii(j),temp); 
         end
       end 
       w(jj) = fact;   
%-----------------------------------------------------------------------
     end %%% if (jlev) 
   end %% while 
   [iw, jw, s] = find(w) ;
   [il,jl,ll] = find(levii); 
   for k=1:length(jw)  
       j = jw(k) ;
       if (j < ii)
          L(ii,j) = s(k)  ; 
       else 
          U(ii,j) = s(k) ;
          levs(ii,j) = ll(k) ; 
       end  %% if
   end   
% ----------------------------------------------------------------------- ;
%      end main loop ;
% ----------------------------------------------------------------------- ;
 end  %% main loop 


