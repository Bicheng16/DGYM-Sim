function [fg, xo, yv, iter] = Y_SS_RachfordRice(zfrac, kval, fg_0)

    % Temporary variable
    Switch_variable = true;

    fg_fake = 500;
    fg_min = 0;
    fg_max = 1;
    fo_min = 0;
    fo_max = 1;
    
    iter = 0;
    Tol = 1.0e-12;  
    error = 100;
    
    % i
    % Check when fg = 1;
    Flag_deriv = false;
    [~, RHS_2] = RachfordRice_fg(fg_max, zfrac, kval, Flag_deriv);
    if (RHS_2 > 0)
        fg = fg_max;
        % fg = min([fg_0, fg_max]);
        [xo, yv] = Y_renorm(fg, kval, zfrac); 
        return;
    else 
        % Check when fg = 0;
        Flag_deriv = false;
        [~, RHS_1] = RachfordRice_fg(fg_min, zfrac, kval, Flag_deriv);
        if (RHS_1 < 0)
            fg = fg_min;
            % fg = max([fg_0, fg_min]);
            [xo, yv] = Y_renorm(fg, kval, zfrac); 
            return;
        else
            if fg_0 ~= fg_fake
                fg = fg_0;
                fo = 1 - fg;
            else
                % ii
                %% fg pre-analysis
                
                % Lower bound
                tmp = (kval .* zfrac - 1) ./ (kval - 1) .* (kval > 1);
                tmp(tmp < 0 | tmp > 1) = 0;
                tmp_min = max(tmp(tmp > 0));
                tmp_min(isempty(tmp_min)) = 0;
                
                % Upper bound
                tmp = (1 - zfrac) ./ (1 - kval) .* (kval < 1);
                tmp(tmp < 0 | tmp > 1) = 0;
                tmp_max = min(tmp(tmp > 0));
                tmp_max(isempty(tmp_max)) = 0;
                
                % Update fo_max and fo_min
                fg_max = fg_max * (tmp_max == 0) + tmp_max * (tmp_max > 0);
                fg_min = fg_min * (tmp_min == 0) + tmp_min * (tmp_min > 0);
                fg = (fg_max + fg_min) / 2;
                
                %% fo pre-analysis
                fo = 1 - fg;
                
                % Lower bound
                % tmp = (zfrac - kval) ./ (1 - kval) .* (kval < 1);
                % tmp(tmp < 0) = 0;
                % tmp_min = max(tmp(tmp > 0));
                % tmp_min(isempty(tmp_min)) = 0;
                
                % Upper bound
                % tmp = kval .* (1 - zfrac) ./ (kval - 1) .* (kval > 1);
                % tmp(tmp < 0) = 0;
                % tmp_max = min(tmp(tmp > 0));
                % tmp_max(isempty(tmp_max)) = 0;
                
                % Update fo_max and fo_min
                % fo_min = fo_min * (tmp_min == 0) + tmp_min * (tmp_min > 0);                
                % fo_max = fo_max * (tmp_max == 0) + tmp_max * (tmp_max > 0);
                % fo = (fo_max + fo_min) / 2;
                
            end
        end 

        % iii
        Flag_deriv = false;
        [~, RHS] = RachfordRice_fg(fg, zfrac, kval, Flag_deriv);
        if (RHS > 0)
            fg_min = fg;
            fo_max = fo;
            Flag_fo_prim = false;
        else
            fg_max = fg;
            fo_min = fo;
            Flag_fo_prim = true;
        end 
        
        % iv
        if Flag_fo_prim && Switch_variable 
            %
            while (error > Tol)

                iter = iter + 1;

                % Residual and derivative calculation
                [dR, RHS] = RachfordRice_fo(fo, zfrac, kval);

                error = abs(RHS);

                if (RHS > 0)
                    fo_max = fo;
                else
                    fo_min = fo;              
                end 

                fo = fo - RHS / dR;

                if (fo < fo_min || fo > fo_max) 
                    fo = (fo_max + fo_min) / 2;
                end
                
            end
            
            fg = 1 - fo;
            
        else
            %
            while (error > Tol)

                iter = iter + 1;

                % Residual and derivative calculation
                Flag_deriv = true;
                [dR, RHS] = RachfordRice_fg(fg, zfrac, kval, Flag_deriv);

                error = abs(RHS);

                if (RHS > 0)
                    fg_min = fg;
                else
                    fg_max = fg;              
                end 

                fg = fg - RHS / dR;

                if (fg < fg_min || fg > fg_max) 
                    fg = (fg_max + fg_min) / 2;
                end

            end
        end
        
        % Calculate vapor and liquid compositions
        [xo, yv] = Y_renorm(fg, kval, zfrac); 
        
    end
    
end

%% Rachford-Rice Equation and its derivate to phase mole fraction
% Primary variable: fg
function [dRR_dfg, RR] = RachfordRice_fg(fg, zfrac, kval, Flag_deriv)

    % Initialize
    dRR_dfg = 0;    
    RR = 0;
    
    temp_1 = 1 + fg .* (kval - 1);
    
    temp_2 = zfrac .* (kval - 1) ./ temp_1;
    RR = sum(temp_2, 2);
    
    % Return if RHS only
    if Flag_deriv          
        temp_1  = temp_1 .* temp_1;
        dRR_dfg = - zfrac .* (kval - 1) .* (kval - 1) ./ temp_1;
        dRR_dfg = sum(dRR_dfg, 2);
    end
    
end
% Primary variable: fo
function [dRR_dfo, RR] = RachfordRice_fo(fo, zfrac, kval)

    % Initialize
    dRR_dfo = 0;    
    RR = 0;
    
    temp_1 = ( kval - 1 ) ./ ( fo + (1 - fo) .* kval );
    
    temp_2 = zfrac .* temp_1;
    RR     = sum(temp_2, 2);
    
    temp_2  = zfrac .* temp_1 .* temp_1;
    dRR_dfo = sum(temp_2, 2);
    
end
