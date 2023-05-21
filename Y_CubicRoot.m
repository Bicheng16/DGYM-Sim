
function[x, nr] = Y_CubicRoot(a, b, c)

    % Root of cubic equation in the form of: z^3 + a * z^2 + b * z + c = 0
    % Referring to: http://mathworld.wolfram.com/CubicFormula.html 

    te = 0;
    ex = 1/3;
    x(1) = 0;
    x(2) = 0;
    x(3) = 0;
    q = (3 * b - (a ^ 2)) / 9;
    r = (9 * a * b - 27 * c - 2 * a ^ 3) / 54;
    d = q ^ 3 + r ^ 2;
    Tol = 1.0d-8; 

    if (abs(d) <= Tol)
        d = 0;
    end

    if (d == 0)
        nr = 0;

        if (r >= 0)
            s = abs(r) ^ ex;
        else
            s = -abs(r) ^ ex;
        end
        
        z1 = 2 * s - a / 3;
        if (z1 <= Tol)
            z1 = 0;
        end
        
        z2 = -s - a / 3;
        if (z2 <= Tol)
            z2 = 0;
        end
        x(1) = max(z1, z2);
        x(3) = min(z1, z2);
        x(2) = x(3);
        
    elseif (d > 0)
        nr = 1;
        u = r + sqrt(d);
        if (u >= 0)
            s = abs(u) ^ ex;
        else
            s = -abs(u) ^ ex;
        end
        
        if (abs(s) <= Tol)
            s = 0;
        end
        
        u = r -sqrt(d);
        if (u >= 0)
            te = abs(u) ^ ex;
        else
            te = -abs(u) ^ ex;
        end
        
        if (abs(te) <= Tol)
            te = 0;
        end
        
        x(1) = s + te - a / 3;
        x(2) = 0;
        x(3) = 0;
        return
        
    elseif (d < 0)
        nr = -1;
        u = 2 * sqrt(-q);
        vacos = r / sqrt(-q ^ 3);
        if (vacos > 1)
            vacos = 1;
        end

        theta = acos(vacos) / 3;
        if (abs(theta) <= Tol)
            theta = 0;
        end

        z1 = u * cos(theta) - a / 3;
        if (z1 <= Tol)
            z1 = 0;
        end

        z2 = u * cos(theta + 2 * pi / 3) - a / 3;
        if (z2 <= Tol)
            z2 = 0;
        end

        z3 = u * cos(theta + 4 * pi / 3) - a / 3;
        if (z3 <= Tol)
            z3 = 0;
        end

        x(1) = max([z1 max([z2 z3])]);
        x(3) = min([z1 min([z2 z3])]);
        if (x(1) == z1)
            x(2) = max([z2 z3]);
        end
        if (x(1) == z2)
            x(2) = max([z1  z3]) ;
        end
        if (x(1) == z3)
            x(2) = max([z1 z2]) ;
        end
    end

end