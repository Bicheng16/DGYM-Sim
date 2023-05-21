function [lnPhi, Z, dlnPhi_dxj, dlnPhi_dP] = Y_CubicDP(nComps, P, T, fluid_type, xi, comp_ACF, comp_Tc, comp_Pc, comp_BIC)

    % Define temporary variables 
    Rg = 10.732;
    zz = zeros(1, 3);
        
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
    temp2 = temp2 .* (1.0 - comp_BIC);
    
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
        
    %% Derivatives to xi
    temp3 = sum(temp2, 2)';
    damix_dxi = 2 .* sqrt_ai .* temp3; 
    dbmix_dxi = bi;
    dA_dxi = (A / amix) .* damix_dxi;
    dB_dxi = (B / bmix) .* dbmix_dxi;
    ds_dxi = (u - 1) .* dB_dxi;
    dq_dxi = dA_dxi + (2 * B * (w - u) - u) .* dB_dxi; 
    dr_dxi = - B .* dA_dxi - (A + (2 + 3 * B) * w * B) .* dB_dxi; 
    
    temp4 = - 1.0 / ( ( 3 * Z + 2 * s ) * Z + q );
    dZ_dxi = ( ds_dxi .* Z + dq_dxi ) * Z + dr_dxi; 
    dZ_dxi = temp4 .* dZ_dxi;
    
    % Calculate fugacity coefficient: fi = p * xi * exp(phi)
    bi_b = bi ./ bmix;
    const_ln = log( (2 * Z + B * c2) / (2 * Z + B * c3) ); 
    lnPhi = bi_b .* (Z - 1) - log(Z - B) + (A / B / c1) .* (bi_b - damix_dxi ./ amix) .* const_ln;         
    
    % Derivative of fugacity 
    dlnPhi_dxj = zeros(nComps, nComps);  % Row: i; colum: j;
    
    % Term 1
    dlnPhi_dxj = dlnPhi_dxj + ( bi_b' * dZ_dxi - (bi_b' * dbmix_dxi) .* ( (Z - 1) / bmix ) );
    
    % Term 2
    temp5 = (dZ_dxi - dB_dxi) ./ (Z - B); 
    temp5 = repmat(temp5, nComps, 1);
    dlnPhi_dxj = dlnPhi_dxj - temp5;
    
    % Group term 3:
    temp6 = (bi_b - damix_dxi ./ amix) ./ c1;
    temp6 = temp6';  % [nc, 1]
    
    temp7 = ( dA_dxi ./ B - (A / B / B) .* dB_dxi ) .* const_ln; % [1, nc]
    
    temp8 = (2 .* dZ_dxi + dB_dxi .* c2) ./ (2 * Z + B * c2);
    temp8 = temp8 - (2 .* dZ_dxi + dB_dxi .* c3) ./ (2 * Z + B * c3);
    temp8 = temp8 .* (A / B); % [1, nc]
    
    dlnPhi_dxj = dlnPhi_dxj + temp6 * (temp7 + temp8); 
    
    % Group term 4:
    temp9 = A / (B * c1) * const_ln; % scalar
    temp10 = - (bi_b' * dbmix_dxi) ./ bmix ...
             + (damix_dxi' * damix_dxi) ./ (amix * amix) ...
             - (2 / amix) .* (1 - comp_BIC) .* (sqrt_ai' * sqrt_ai);  % [nc, nc]
    dlnPhi_dxj = dlnPhi_dxj + temp9 .* temp10; 
    
    %% Derivatives to Pressure
    dA_dP = A / P;
    dB_dP = B / P;
    ds_dP = (u - 1) .* dB_dP; 
    dq_dP = dA_dP + (2 * B * (w - u) - u) .* dB_dP;  
    dr_dP = - B .* dA_dP - (A + (2 + 3 * B) * w * B) .* dB_dP; 

    dZ_dP = ( ds_dP .* Z + dq_dP ) * Z + dr_dP; 
    dZ_dP = temp4 .* dZ_dP;
    
    dlnPhi_dP = 0;
    
    % Term 1
    dlnPhi_dP = dlnPhi_dP + bi_b .* dZ_dP - (dZ_dP - dB_dP) / (Z - B);
    
    % Term 2
    temp11 = ( dA_dP / B - dB_dP * A / (B * B) ) * const_ln ...
           + (A / B) * ( (2 * dZ_dP + dB_dP * c2) / (2 * Z + B * c2) ...
                       - (2 * dZ_dP + dB_dP * c3) / (2 * Z + B * c3) ); 
        
    dlnPhi_dP = dlnPhi_dP + (bi_b - damix_dxi ./ amix) .* (temp11 / c1);  
        
end