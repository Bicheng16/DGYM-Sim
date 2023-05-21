function [Phi, Z] = Y_Cubic(nComps, P, T, fluid_type, xi, comp_ACF, comp_Tc, comp_Pc, comp_BIC)

    % Define temporary variables 
    Rg = 10.732;
    zz = zeros(1, 3);
    dummy = ones(nComps, nComps); 
        
    % Frequenly used group variables
    RgTc = Rg .* comp_Tc;         % Rg * Tc
    RgT  = Rg .* T;               % Rg * T
    RgTc_Pc = RgTc ./ comp_Pc;    % Rg * Tc / Pc
    T_Tc = T ./ comp_Tc;          % T / Tc

    flag_omega = (comp_ACF > 0.49);
    mi(1, flag_omega == 1) = 0.379642d0 + comp_ACF(flag_omega == 1) ...
                        .* ( 1.48503d0  + comp_ACF(flag_omega == 1) .* ( - 0.164423d0 + 0.016666d0 .* comp_ACF(flag_omega == 1) ) );    
    mi(1, flag_omega == 0) = 0.37464d0  + comp_ACF(flag_omega == 0) .* ( 1.54226 - 0.26992 .* comp_ACF(flag_omega == 0) );  
         
    alphai = (1 + mi .* ( 1 - sqrt(T_Tc) ) );
    alphai = alphai .* alphai;
    
    aci = 0.45724 .* RgTc .* RgTc_Pc ;
    
    ai = aci.* alphai;
    
    bi = 0.0778 .* RgTc_Pc;
    
    %% Applying binary mixing rule
    sqrt_ai = ai .^ 0.5; 
    
    temp1 = xi .* sqrt_ai;
    temp2 = temp1(ones(nComps, 1), :); 
    temp2 = temp2 .* (dummy - comp_BIC);
    
    % Variables
    amix  = temp1 * temp2;
    amix = sum(amix);
    bmix = sum(xi .* bi);
    
    % A and B constant 
    A = (amix * P) ./ (RgT .* RgT);
    B = (bmix * P) ./ RgT;
    
    % Coefficient for PR-Cubic-EOS
    % Constant in PR-EOS
    u = 2;
    w = -1;
    % Cubic: Z^3 + s*Z^2 + q*Z + r = 0;
    s = (u - 1) * B - 1;
    q = A + ( (w - u) * B - u ) * B;
    r = - ( A + (1 + B) * B * w ) * B;
    
    % Solve cubic equation
    [zz, nr] = Y_CubicRoot(s, q, r);
    if (nr <= 0)
        % fluid_type: 0, gas 
        if (fluid_type == 0) 
            Z = zz(1);
        % fluid_type: 1, liquid
        elseif (fluid_type >= 1)
            Z = zz(3);
        end
    else
        Z = zz(1);
    end 
    
    % Get some constant
    c1 = sqrt(u * u - 4 * w);
    c2 = u + c1;
    c3 = u - c1;    

    % Derivatives refer to Cao Hui Dissertation, Appendix C
    temp3 = sum(temp2, 2)';
    damix_dxi = 2 .* sqrt_ai .* temp3; 
    
    % Calculate fugacity coefficient: fi = p * xi * exp(phi)
    bi_b = bi ./ bmix;
    const_ln = log( (2 * Z + B * c2) / (2 * Z + B * c3)); 
    lnPhi = bi_b .* (Z - 1) - log(Z - B) + (A / B / c1) .* (bi_b - damix_dxi ./ amix) .* const_ln;         
    Phi = exp(lnPhi);
    
end