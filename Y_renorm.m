function [xo, yv] = Y_renorm(fg, kval, zfrac)

    xo = zfrac ./ (1 + fg .* (kval - 1) );
    yv = kval .* xo;
    
    ul1 = 1 / sum(xo);
    ul2 = 1 / sum(yv);

    xo = ul1 .* xo;
    yv = ul2 .* yv;

end